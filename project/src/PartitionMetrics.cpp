#include "PartitionMetrics.hpp"
#include <mpi.h>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <cmath>

#include <queue>
#include <map>
#include <set>

PartitionMetrics::PartitionMetrics(const std::string outputFolder, const std::string vizFolder) : outputFolder(outputFolder), vizFolder(vizFolder) {
}

void PartitionMetrics::save_partition_info(const Partition& partition, const int actualParts, const CSR& csr) {
    verticesCounter.assign(actualParts, 0);
    vertexImbalance.assign(actualParts, 0.0);
    edgesCounter.assign(actualParts, 0);
    edgeImbalance.assign(actualParts, 0.0);
    boundaryPerPart.assign(actualParts, 0);
    componentsPerPart.assign(actualParts, 0);
    neighbours.clear();

    computeVertexImbalance(partition.partition, actualParts);
    computeEdgeImbalance(partition.partition, actualParts, csr);
    exchange_ghost_parts(partition.partition, actualParts, csr);
    computeBoundaryVertices(partition.partition, actualParts, csr);
    countNeighbours(partition.partition, actualParts, csr);
    computeConnectedComponents(partition.partition, actualParts, csr);
    
    int finalCut = partition.cutEdge;
    if (partition.type == Partition::Type::SCOTCH || partition.type == Partition::Type::ZOLTAN2_SPHYNX || partition.type == Partition::Type::ZOLTAN_PHG) {
        countCutEdges(partition.partition, actualParts, nullptr, finalCut, csr);
    }
    if(csr.rank == 0)
        std::cout << "Edge cut: " << finalCut << std::endl;

    saveMetricsToFile(partition, finalCut, csr.rank, csr.filename);
}

/*определяем к какому процессу относится вершина*/
int PartitionMetrics::findOwnerRank(const int v, const std::vector<int>& vtxdist) {
    int left = 0, right = vtxdist.size() - 1;

    while (left < right - 1) {
        int mid = (left + right) / 2;
        if (v < vtxdist[mid])
            right = mid;
        else
            left = mid;
    }
    return left;
}

/*определяем номер кластера вершины. Если вершина нам не принадлежит, смотрим в ghost_global_ids*/
int PartitionMetrics::getPartitionOfVertex(int v, const CSR& csr, const std::vector<int>& partition) {
    /* 1. Проверяем, не наша ли это вершина */
    if (csr.vtxdist[csr.rank] <= v && v < csr.vtxdist[csr.rank + 1]) {
        return partition[v - csr.vtxdist[csr.rank]];
    }

    /* 2. Бинарный поиск в отсортированном векторе призрачных вершин */
    auto it = std::lower_bound(ghost_global_ids.begin(), ghost_global_ids.end(), v);

    if (it != ghost_global_ids.end() && *it == v) {
        size_t idx = std::distance(ghost_global_ids.begin(), it);
        return ghost_part_ids[idx];
    }

    std::cerr << "Rank " << csr.rank << ": Missing ghost vertex " << v << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 0;
}

void PartitionMetrics::exchange_ghost_parts(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
    if (!graph_structure_initialized) {
        std::vector<std::vector<int>> send_requests(csr.size);
        for (int p = 0; p < csr.size; p++) send_requests[p].reserve(csr.local_rows);

        for (int u = 0; u < (int)partition.size(); u++) {
            for (int e = csr.rows[u]; e < csr.rows[u+1]; e++) {
                int v = csr.cols[e];
                int owner = findOwnerRank(v, csr.vtxdist);
                if (owner != csr.rank) send_requests[owner].push_back(v);
            }
        }

        send_counts.assign(csr.size, 0);
        recv_counts.assign(csr.size, 0);
        for (int p = 0; p < csr.size; p++) {
            std::sort(send_requests[p].begin(), send_requests[p].end());
            send_requests[p].erase(std::unique(send_requests[p].begin(), send_requests[p].end()), send_requests[p].end());
            send_counts[p] = send_requests[p].size();
        }

        MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        send_displs.assign(csr.size, 0);
        recv_displs.assign(csr.size, 0);
        total_send = 0; total_recv = 0;
        for (int p = 0; p < csr.size; p++) {
            send_displs[p] = total_send;
            recv_displs[p] = total_recv;
            total_send += send_counts[p];
            total_recv += recv_counts[p];
        }

        std::vector<int> send_buf(total_send), recv_buf(total_recv);
        for (int p = 0; p < csr.size; p++) {
            std::copy(send_requests[p].begin(), send_requests[p].end(), send_buf.begin() + send_displs[p]);
        }

        MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(), MPI_INT,
                      recv_buf.data(), recv_counts.data(), recv_displs.data(), MPI_INT, MPI_COMM_WORLD);

        export_indices.assign(total_recv, 0);
        for (int i = 0; i < total_recv; i++) {
            export_indices[i] = recv_buf[i] - csr.vtxdist[csr.rank];
        }

        sorting_p.resize(total_send);
        std::iota(sorting_p.begin(), sorting_p.end(), 0);
        std::sort(sorting_p.begin(), sorting_p.end(), [&](int i, int j) {
            return send_buf[i] < send_buf[j];
        });

        ghost_global_ids.resize(total_send);
        for (int i = 0; i < total_send; i++) {
            ghost_global_ids[i] = send_buf[sorting_p[i]];
        }

        graph_structure_initialized = true;
    }

    std::vector<int> send_answers(total_recv);
    for (int i = 0; i < total_recv; i++) {
        send_answers[i] = partition[export_indices[i]];
    }

    std::vector<int> recv_answers(total_send);
    MPI_Alltoallv(send_answers.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                  recv_answers.data(), send_counts.data(), send_displs.data(), MPI_INT, MPI_COMM_WORLD);

    ghost_part_ids.resize(total_send);
    for (int i = 0; i < total_send; i++) {
        ghost_part_ids[i] = recv_answers[sorting_p[i]];
    }
}

void PartitionMetrics::computeVertexImbalance(const std::vector<int>& partition, const int actualParts) {
    std::vector<int> local_counts(actualParts, 0);
    for (int p : partition)
        local_counts[p]++;
    MPI_Allreduce(local_counts.data(), verticesCounter.data(), actualParts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    double total = 0.0;
    for (int v : verticesCounter)
        total += v;

    double avg = total / actualParts;
    for (int i = 0; i < actualParts; i++) {
        vertexImbalance[i] = (verticesCounter[i] / avg - 1.0) * 100.0;
    }
}

void PartitionMetrics::computeEdgeImbalance(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
    std::vector<long long> local_edges(actualParts, 0);
    for (int i = 0; i < csr.local_rows; i++) {
        local_edges[partition[i]] += (csr.rows[i+1] - csr.rows[i]);
    }
    
    MPI_Allreduce(local_edges.data(), edgesCounter.data(), actualParts, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    long long total_e = 0;
    for (long long e : edgesCounter) total_e += e;
    double avg_e = (double)total_e / actualParts;

    for (int i = 0; i < actualParts; i++) {
        edgeImbalance[i] = (avg_e > 0) ? (edgesCounter[i] / avg_e - 1.0) * 100.0 : 0.0;
    }
}


/*вершины, которые имеют соседа в другом кластере*/
void PartitionMetrics::computeBoundaryVertices(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
    std::vector<int> local_boundary(actualParts, 0);

    for (int u = 0; u < partition.size(); u++) {
        int pu = partition[u];
        bool isBoundary = false;

        for (int i = csr.rows[u]; i < csr.rows[u+1]; i++) {
            int v = csr.cols[i];
            int pv = getPartitionOfVertex(v, csr, partition);

            if (pu != pv) {
                isBoundary = true;
                break;
            }
        }

        if (isBoundary)
            local_boundary[pu]++;
    }

    MPI_Allreduce(local_boundary.data(), boundaryPerPart.data(),
                actualParts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

void PartitionMetrics::countCutEdges(const std::vector<int>& partition, const int actualParts, std::vector<int>* cutEdgesPerPart, int& totalCutEdges, const CSR& csr) {
    std::vector<int> local_cut(actualParts, 0);
    int local_total = 0;

    for (int u = 0; u < partition.size(); u++) {
        int pu = partition[u];

        for (int i = csr.rows[u]; i < csr.rows[u+1]; i++) {
            int v = csr.cols[i];
            int pv = getPartitionOfVertex(v, csr, partition);

            if ( pu != pv) {
                local_cut[pu]++;
                local_total++;
            }
        }
    }

    if(cutEdgesPerPart != nullptr)
        MPI_Allreduce(local_cut.data(), cutEdgesPerPart->data(), 
                    actualParts, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    MPI_Allreduce(&local_total, &totalCutEdges,
                1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    totalCutEdges /= 2;
}

bool PartitionMetrics::Edge::operator<(const Edge& other) const {
    return std::tie(p1, p2) < std::tie(other.p1, other.p2);
}

bool PartitionMetrics::Edge::operator==(const Edge& other) const {
    return p1 == other.p1 && p2 == other.p2;
}

/*для каждого кластера считаем число кластеров с которыми тот общается */
void PartitionMetrics::countNeighbours(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
    /* --- 0. Создаем тип данных MPI для структуры Edge --- */
    MPI_Datatype MPI_EDGE_TYPE;
    MPI_Type_contiguous(2, MPI_INT, &MPI_EDGE_TYPE);
    MPI_Type_commit(&MPI_EDGE_TYPE);

    /* 1. Собираем локальные пары (Edge) */
    std::vector<Edge> local_edges;
    for (int u = 0; u < (int)partition.size(); u++) {
        int pu = partition[u];
        for (int i = csr.rows[u]; i < csr.rows[u + 1]; i++) {
            int v = csr.cols[i];
            int pv = getPartitionOfVertex(v, csr, partition);
            if (pu != pv) {
                /* Всегда храним пару в отсортированном виде (min, max), чтобы легче удалять дубликаты */
                local_edges.push_back({std::min(pu, pv), std::max(pu, pv)});
            }
        }
    }

    /* Удаляем локальные дубликаты (внутри текущего MPI-процесса) */
    std::sort(local_edges.begin(), local_edges.end());
    local_edges.erase(std::unique(local_edges.begin(), local_edges.end()), local_edges.end());

    int local_edge_count = (int)local_edges.size();
    std::vector<int> recv_counts(csr.size);
    MPI_Allgather(&local_edge_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> displs(csr.size, 0);
    int total_edges = 0;
    for (int i = 0; i < csr.size; i++) {
        displs[i] = total_edges;
        total_edges += recv_counts[i];
    }

    std::vector<Edge> all_edges(total_edges);

    MPI_Allgatherv(local_edges.data(), local_edge_count, MPI_EDGE_TYPE,
                all_edges.data(), recv_counts.data(), displs.data(), MPI_EDGE_TYPE,
                MPI_COMM_WORLD);


    neighbours.assign(actualParts, std::vector<int>());
    
    /* Сортируем глобальный список ребер и удаляем дубликаты
     (одна и та же граница разделов могла быть обнаружена разными MPI-процессами)*/
    std::sort(all_edges.begin(), all_edges.end());
    all_edges.erase(std::unique(all_edges.begin(), all_edges.end()), all_edges.end());

    /* Заполняем списки смежности для каждого раздела*/
    for (const auto& e : all_edges) {
        neighbours[e.p1].push_back(e.p2);
        neighbours[e.p2].push_back(e.p1);
    }
    
    /* Сортируем списки соседей (опционально, для удобства поиска)*/
    for (auto& list : neighbours) {
        std::sort(list.begin(), list.end());
    }

    MPI_Type_free(&MPI_EDGE_TYPE);
}

void PartitionMetrics::computeConnectedComponents(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
    int n_local = partition.size();
    std::vector<int> local_comp_id(n_local, -1);
    int local_comp_counter = 0;

    /* --- 1. Локальный BFS --- */
    for (int i = 0; i < n_local; ++i) {
        if (local_comp_id[i] != -1) continue;
        int part_id = partition[i];
        std::queue<int> q;
        q.push(i);
        local_comp_id[i] = local_comp_counter;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            for (int e = csr.rows[u]; e < csr.rows[u + 1]; ++e) {
                int v_glob = csr.cols[e];
                if (v_glob >= csr.vtxdist[csr.rank] && v_glob < csr.vtxdist[csr.rank + 1]) {
                    int v_loc = v_glob - csr.vtxdist[csr.rank];
                    if (local_comp_id[v_loc] == -1 && partition[v_loc] == part_id) {
                        local_comp_id[v_loc] = local_comp_counter;
                        q.push(v_loc);
                    }
                }
            }
        }
        local_comp_counter++;
    }

    /* --- 2. ОБМЕН ID КОМПОНЕНТ */
    std::vector<int> send_answers(total_recv);
    for (int i = 0; i < total_recv; i++) {
        send_answers[i] = local_comp_id[export_indices[i]];
    }

    std::vector<int> recv_answers_mpi(total_send);
    MPI_Alltoallv(send_answers.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                  recv_answers_mpi.data(), send_counts.data(), send_displs.data(), MPI_INT, MPI_COMM_WORLD);

    std::vector<int> ghost_local_comp_ids(total_send);
    for (int i = 0; i < total_send; i++) {
        ghost_local_comp_ids[i] = recv_answers_mpi[sorting_p[i]];
    }

    /*  --- 3. Сбор уникальных склеек --- */
    std::set<std::pair<long long, long long>> merge_pairs; 
    std::vector<long long> local_merges; 
    long long offset = 1000000000LL; 

    for (int u = 0; u < n_local; ++u) {
        int pu = partition[u];
        long long u_comp = (long long)csr.rank * offset + local_comp_id[u];

        for (int e = csr.rows[u]; e < csr.rows[u+1]; ++e) {
            int v_glob = csr.cols[e];
            
            if (v_glob < csr.vtxdist[csr.rank] || v_glob >= csr.vtxdist[csr.rank+1]) {
                auto it = std::lower_bound(ghost_global_ids.begin(), ghost_global_ids.end(), v_glob);
                
                if (it != ghost_global_ids.end() && *it == v_glob) {
                    int idx = (int)std::distance(ghost_global_ids.begin(), it);
                    
                    if (ghost_part_ids[idx] == pu) {
                        int v_owner = findOwnerRank(v_glob, csr.vtxdist);
                        long long v_comp = (long long)v_owner * offset + ghost_local_comp_ids[idx];
                        
                        if (u_comp != v_comp && merge_pairs.insert({std::min(u_comp, v_comp), std::max(u_comp, v_comp)}).second) {
                            local_merges.push_back((long long)pu);
                            local_merges.push_back(u_comp);
                            local_merges.push_back(v_comp);
                        }
                    }
                }
            }
        }
    }

    /*  --- 4. Сбор данных на Rank 0 --- */
    int local_merges_count = (int)local_merges.size();
    std::vector<int> all_merge_counts(csr.size);
    MPI_Gather(&local_merges_count, 1, MPI_INT, all_merge_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<long long> global_merges;
    std::vector<int> m_displs(csr.size, 0);
    int total_m_elements = 0;
    if (csr.rank == 0) {
        for (int i = 0; i < csr.size; ++i) {
            m_displs[i] = total_m_elements;
            total_m_elements += all_merge_counts[i];
        }
        global_merges.resize(total_m_elements);
    }

    MPI_Gatherv(local_merges.data(), local_merges_count, MPI_LONG_LONG,
                global_merges.data(), all_merge_counts.data(), m_displs.data(), 
                MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    std::vector<int> local_comps_per_rank_part(actualParts, 0);
    std::vector<std::set<int>> comps_by_part(actualParts);
    for(int i=0; i<n_local; ++i) comps_by_part[partition[i]].insert(local_comp_id[i]);
    for(int p=0; p<actualParts; ++p) local_comps_per_rank_part[p] = comps_by_part[p].size();

    std::vector<int> all_local_comp_counts(actualParts * csr.size);
    MPI_Gather(local_comps_per_rank_part.data(), actualParts, MPI_INT,
               all_local_comp_counts.data(), actualParts, MPI_INT, 0, MPI_COMM_WORLD);

    /*  --- 5. DSU на Rank 0 --- */
    if (csr.rank == 0) {
        componentsPerPart.assign(actualParts, 0);
        for (int p = 0; p < actualParts; ++p) {
            std::map<long long, long long> parent;
            auto find_set = [&](long long v, auto& self) -> long long {
                if (parent.find(v) == parent.end()) return parent[v] = v;
                return (v == parent[v]) ? v : (parent[v] = self(parent[v], self));
            };
            
            int initial_comps = 0;
            for (int r = 0; r < csr.size; ++r) initial_comps += all_local_comp_counts[r * actualParts + p];

            int merges_done = 0;
            for (size_t i = 0; i < global_merges.size(); i += 3) {
                if (global_merges[i] == p) {
                    long long root_u = find_set(global_merges[i+1], find_set);
                    long long root_v = find_set(global_merges[i+2], find_set);
                    if (root_u != root_v) {
                        parent[root_v] = root_u;
                        merges_done++;
                    }
                }
            }
            componentsPerPart[p] = initial_comps - merges_done;
        }
    }
}

std::string PartitionMetrics::getFileName(const std::string& filename) {
    std::filesystem::path p(filename);
    return p.stem().string();
}

void PartitionMetrics::saveMetricsToFile(const Partition& partition, const int cutEdge, const int rank, const std::string filename) {
    if (rank == 0) {
        std::string outfile = outputFolder + "/" + getFileName(filename) + "_" + Partition::getPartitionName(partition.type) 
                                            + "_p" + std::to_string(partition.nparts) + "_" + partition.postfix + "_info";
        std::ofstream f(outfile);
        if (!f.is_open()) {
            std::cout << "ERROR: Failed to open file for writing: " << outfile << std::endl;
            return;
        }

        auto getStats = [&](const auto& vec) {
            if (vec.empty()) return std::make_tuple(0.0, 0.0, 0.0, 0.0);
            
            auto [min_it, max_it] = std::minmax_element(vec.begin(), vec.end());
            double min_val = static_cast<double>(*min_it);
            double max_val = static_cast<double>(*max_it);
            double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
            double avg = sum / vec.size();

            double abs_deviation_sum = 0.0;
            double variance_sum = 0.0;

            for (auto val : vec) {
                double diff = static_cast<double>(val) - avg;
                abs_deviation_sum += std::abs(diff);
                variance_sum += diff * diff;
            }

            double mad = abs_deviation_sum / vec.size();
            double std_dev = std::sqrt(variance_sum / vec.size());

            return std::make_tuple(avg, max_val - min_val, mad, std_dev);
        };

        std::vector<int> neighborCounts(partition.actualParts);
        for(int i=0; i < partition.actualParts; ++i) neighborCounts[i] = neighbours[i].size();

        auto [avgV, diffV, madV, stdV] = getStats(verticesCounter);
        auto [avgE, diffE, madE, stdE] = getStats(edgesCounter);
        auto [avgB, diffB, madB, stdB] = getStats(boundaryPerPart);
        auto [avgN, diffN, madN, stdN] = getStats(neighborCounts);
        auto [avgC, diffC, madC, stdC] = getStats(componentsPerPart);

        long long totalVertices = 0;
        long long totalEdges = 0;

        for (int i = 0; i < partition.actualParts; i++) {
            totalVertices += verticesCounter[i];
            totalEdges += edgesCounter[i];
        }

        f << "======================================================================" << std::endl;
        f << " PARTITION REPORT: " << Partition::getPartitionName(partition.type) << std::endl;
        f << "======================================================================" << std::endl;
        f << "File:             " << filename << std::endl;
        f << "Seed:             " << partition.seed << std::endl;
        f << "Total Edge Vol:   " << totalEdges << std::endl;
        f << "Total Vertex Vol: " << totalVertices << std::endl;
        f << "Target Parts:     " << partition.nparts << std::endl;
        f << "Actual Parts:     " << partition.actualParts << std::endl;
        f << "Target Imbalance: " << (int)(partition.imbalance * 100) << "%" << std::endl;
        f << "----------------------------------------------------------------------" << std::endl;
        f << "Edge Cut:              " << cutEdge << std::endl;
        f << "Time Total (sec):      " << std::fixed << std::setprecision(4) << partition.time_total << std::endl;
        f << "Time Solve Only (sec): " << std::fixed << std::setprecision(4) << partition.time_solve << std::endl;
        f << "Overhead (sec):        " << std::fixed << std::setprecision(4) << (partition.time_total - partition.time_solve) << std::endl;
        f << "----------------------------------------------------------------------" << std::endl;
        f << std::fixed << std::setprecision(2);
        f << "Vertices:  avg: " << std::setw(10) << avgV << " diff: " << std::setw(8) << diffV << " avg_dev: " << std::setw(8) << madV << " std_dev: " << std::setw(8) << stdV << std::endl;
        f << "Edges:     avg: " << std::setw(10) << avgE << " diff: " << std::setw(8) << diffE << " avg_dev: " << std::setw(8) << madE << " std_dev: " << std::setw(8) << stdE << std::endl;
        f << "Boundary:  avg: " << std::setw(10) << avgB << " diff: " << std::setw(8) << diffB << " avg_dev: " << std::setw(8) << madB << " std_dev: " << std::setw(8) << stdB << std::endl;
        f << "Neighbors: avg: " << std::setw(10) << avgN << " diff: " << std::setw(8) << diffN << " avg_dev: " << std::setw(8) << madN << " std_dev: " << std::setw(8) << stdN << std::endl;
        f << "Comps:     avg: " << std::setw(10) << avgC << " diff: " << std::setw(8) << diffC << " avg_dev: " << std::setw(8) << madC << " std_dev: " << std::setw(8) << stdC << std::endl;
        f << "======================================================================" << std::endl;
        f << std::endl;

        f << std::left << std::setw(8)  << "Part" 
        << std::setw(10) << "Vertices" 
        << std::setw(10) << "V-Imb %" 
        << std::setw(12) << "Edges" 
        << std::setw(10) << "E-Imb %" 
        << std::setw(10) << "Boundary" 
        << std::setw(8) << "Neig"
        << std::setw(8) << "Comps" << std::endl;
        f << "----------------------------------------------------------------------" << std::endl;

        f << std::fixed << std::setprecision(2);
        for (int i = 0; i < partition.actualParts; i++) {
            f << std::left << std::setw(8)  << i 
            << std::setw(10) << verticesCounter[i] 
            << std::setw(10) << vertexImbalance[i] 
            << std::setw(12) << edgesCounter[i] 
            << std::setw(10) << edgeImbalance[i] 
            << std::setw(10) << boundaryPerPart[i] 
            << std::setw(8) << neighbours[i].size() 
            << std::setw(8) << componentsPerPart[i] << std::endl;
        }
        f << "======================================================================" << std::endl;

        f.flush();
        f.close();
        std::cout << "Partition report saved to: " << outfile << std::endl;
    }
}

void PartitionMetrics::save_graph_topology(const CSR& csr) {
    std::string filename = vizFolder + "/" + getFileName(csr.filename) + "_topology.csv";
    std::stringstream ss;

    // 1. Каждый процесс готовит свою часть строк в буфер
    for (int i = 0; i < csr.local_rows; i++) {
        int u = csr.start_row + i;
        for (int j = csr.rows[i]; j < csr.rows[i+1]; j++) {
            if (u < csr.cols[j]) {
                ss << u << "," << csr.cols[j] << "\n";
            }
        }
    }

    std::string local_str = ss.str();
    long long local_size = local_str.size();
    long long offset = 0;

    // 2. Вычисляем смещения (префиксная сумма длин буферов)
    MPI_Exscan(&local_size, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (csr.rank == 0) offset = 0;

    // Добавляем длину заголовка к смещению каждого процесса
    std::string header = "source,target\n";
    long long header_len = header.size();
    offset += header_len;

    // 3. Параллельная запись через MPI-IO
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), 
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    MPI_File_set_size(fh, 0); 

    // Только Rank 0 пишет заголовок в самое начало
    if (csr.rank == 0) {
        MPI_File_write_at(fh, 0, header.c_str(), header_len, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    // Все пишут свои куски по вычисленным смещениям
    MPI_File_write_at_all(fh, offset, local_str.c_str(), local_size, MPI_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
}

void PartitionMetrics::save_partition_vec(const Partition& partition, const CSR& csr) {
    std::string filename = vizFolder + "/" + getFileName(csr.filename) + "_" 
                         + Partition::getPartitionName(partition.type) + "_p" 
                         + std::to_string(partition.nparts) + "_" + partition.postfix + "_map.csv";
    
    std::stringstream ss;
    for (int i = 0; i < csr.local_rows; i++) {
        ss << (csr.start_row + i) << "," << partition.partition[i] << "\n";
    }

    std::string local_str = ss.str();
    long long local_size = local_str.size();
    long long offset = 0;

    // Префиксная сумма для определения позиций
    MPI_Exscan(&local_size, &offset, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (csr.rank == 0) offset = 0;

    std::string header = "node_id,partition_id\n";
    long long header_len = header.size();
    offset += header_len;

    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), 
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    MPI_File_set_size(fh, 0); 

    if (csr.rank == 0) {
        MPI_File_write_at(fh, 0, header.c_str(), header_len, MPI_CHAR, MPI_STATUS_IGNORE);
    }

    MPI_File_write_at_all(fh, offset, local_str.c_str(), local_size, MPI_CHAR, MPI_STATUS_IGNORE);

    MPI_File_close(&fh);
}
