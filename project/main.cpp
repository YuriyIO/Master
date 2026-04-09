#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

#include <parmetis.h>
#include <parhip_interface.h>
#include <ptscotch.h>
#include <dkaminpar.h>

/*
./compile.sh
mpirun -np 4 ./main CSR/test 4 Partitions
*/

class CSR {
public:
    const int rank;
    const int size;

    const std::string filename;

    int nrows = 0;
    int nnz   = 0;

    int start_row  = 0;
    int end_row    = 0;
    int local_rows = 0;
    
    std::vector<int> vtxdist;

    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;

    CSR(const int rank, const int size, const std::string& filename) 
        : rank(rank), size(size), filename(filename) 
    {}

    void readGraph() {
        readCSR();
        deleteDiagElems();
        computeVtxDist();
    }

    void printCSR() {
        // Вывод для проверки (после чтения данных)
        for (int u = 0; u < size; u++) {
            if (u == rank) {
                std::cout << "Rank " << rank << " rows [" << start_row 
                        << ", " << end_row << ")\n";
                
                for (int i = 0; i < local_rows; i++) {
                    int global_row = start_row + i;
                    int start = rows[i];
                    int end = rows[i + 1];
                    
                    std::cout << "Row " << global_row << ": ";
                    for (int j = start; j < end; j++) {
                        std::cout << "(" << cols[j] << ", " 
                                << vals[j] << ") ";
                    }
                    std::cout << "\n";
                }
                fflush(stdout);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    void printDiagElems() {
        // Вывод для проверки (после чтения данных)
        for (int u = 0; u < size; u++) {
            if (u == rank) {
                std::cout << "Rank " << rank << " rows [" << start_row 
                        << ", " << end_row << ")\n";
                
                for (int i = 0; i < local_rows; i++) {
                    int global_row = start_row + i;
                    int start = rows[i];
                    int end = rows[i + 1];
                    
                    std::cout << "Row " << global_row << ": ";
                    for (int j = start; j < end; j++) {
                        if(global_row == cols[j])

                            std::cout << "(" << cols[j] << ", " << vals[j] << ") ";
                    }
                    std::cout << "\n";
                }
                fflush(stdout);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    void printVDist() {
        // Вывод для проверки (после чтения данных)
        for (int u = 0; u < size; u++) {
            if (u == rank) {
                std::cout << "Rank: " << rank << std::endl;
                for (int i = 0; i < size+1; i++) {
                    std::cout << vtxdist[i] << " ";
                }
                std::cout << std::endl;
                fflush(stdout);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

private:

    void readCSR() {
        MPI_File fh;
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        
        // Читаем заголовок
        int header[2];
        if (rank == 0) {
            MPI_File_read_at(fh, 0, header, 2, MPI_INT, MPI_STATUS_IGNORE);
        }
        MPI_Bcast(header, 2, MPI_INT, 0, MPI_COMM_WORLD);
        
        nrows = header[0];
        nnz = header[1];
        
        // Вычисляем распределение строк
        int base = nrows / size;
        int rem = nrows % size;
        
        local_rows = base + (rank < rem ? 1 : 0);
        start_row = rank * base + std::min(rank, rem);
        end_row = start_row + local_rows;

        rows.resize(local_rows + 1);
        MPI_Offset rows_offset = 2 * sizeof(int);
        MPI_File_read_at_all(fh,
                            rows_offset + start_row * sizeof(int),
                            rows.data(),
                            local_rows + 1,
                            MPI_INT,
                            MPI_STATUS_IGNORE);
        

        // Убираем смещение
        int global_start_nnz = rows[0];   
        for (int i = 0; i <= local_rows; i++) {
            rows[i] -= global_start_nnz;
        }
        int local_nnz = rows[local_rows];

        // Чтение cols и vals
        cols.resize(local_nnz);
        vals.resize(local_nnz);

        MPI_Offset offset = 2 * sizeof(int);

        MPI_Offset cols_offset = 2 * sizeof(int) + sizeof(int) * (nrows + 1);
        MPI_Offset vals_offset = cols_offset + sizeof(int) * nnz;
        
        MPI_File_read_at_all(fh,
                            cols_offset + global_start_nnz * sizeof(int),
                            cols.data(),
                            local_nnz,
                            MPI_INT,
                            MPI_STATUS_IGNORE);
        
        MPI_File_read_at_all(fh,
                            vals_offset + global_start_nnz * sizeof(double),
                            vals.data(),
                            local_nnz,
                            MPI_DOUBLE,
                            MPI_STATUS_IGNORE);
        
        MPI_File_close(&fh);
    }

    void deleteDiagElems() {
        std::vector<int> new_rows(local_rows + 1, 0);
        std::vector<int> new_cols;
        std::vector<double> new_vals;

        int nnz_count = 0;  // Счётчик новых nnz для локального блока

        for (int i = 0; i < local_rows; i++) {
            int global_row = start_row + i;
            int start = rows[i];
            int end = rows[i + 1];

            for (int j = start; j < end; j++) {
                if (cols[j] != global_row) { // Оставляем только недиагональные элементы
                    new_cols.push_back(cols[j]);
                    new_vals.push_back(vals[j]);
                    nnz_count++;
                }
            }
            new_rows[i + 1] = nnz_count; // Обновляем строку
        }

        cols = std::move(new_cols);
        vals = std::move(new_vals);
        rows = std::move(new_rows);
    }

    void computeVtxDist() {
        std::vector<int> all_local_rows(size);
        MPI_Allgather(&local_rows, 1, MPI_INT,
                    all_local_rows.data(), 1, MPI_INT,
                    MPI_COMM_WORLD);

        vtxdist.resize(size + 1);
        vtxdist[0] = 0;
        for (int i = 0; i < size; i++) {
            vtxdist[i + 1] = vtxdist[i] + all_local_rows[i];
        }
    }

};


class Partition {
public:

    enum class Type {
        PARMETIS,
        SCOTCH,
        KAHIP,
        KAMINPAR,
        ZOLTAN2
    };

    const int nparts;
    int actualParts;
    const double imbalance;   
    std::vector<int> partition;
    int cutEdge;
    int ret;
    Type type;

    Partition(const int nparts, const double imbalance, const int local_rows) : nparts(nparts), imbalance(imbalance) {
        partition.resize(local_rows);
    }

    void run(Type type, const CSR& csr) {
        actualParts = nparts;
        ret = 0;
        this->type = type; 
        cutEdge = -1;
        switch(type) {
            case Type::PARMETIS: run_parmetis(csr); break;
            case Type::SCOTCH:   run_ptscotch(csr);   break;
            case Type::KAHIP:    run_parhip(csr);   break;
            case Type::KAMINPAR:   run_kaminpar(csr);   break;
        }
    }

    void run_parmetis(const CSR& csr) {    
        // ParMETIS параметры
        int wgtflag = 0;  // 0 - нет весов, 1 - веса вершин, 2 - веса рёбер, 3 - оба
        int numflag = 0;  // 0 - C-нумерация (начиная с 0), 1 - Fortran-нумерация (начиная с 1)
        int ncon = 1;     // Количество ограничений (обычно 1)
        int nparts_idx = nparts;  // Количество частей для разбиения
        float ubvec = 1 + imbalance;  // Целевой дисбаланс (5%)
        int options[3] = {1, 0, 0};  // Параметры: 1 - использовать параметры по умолчанию
        
        // // Конвертируем веса рёбер в формат ParMETIS (int)
        // vector<int> adjwgt_idx(adjwgt.size());
        // if (wgtflag == 2 || wgtflag == 3) {
        //     for (size_t i = 0; i < adjwgt.size(); i++) {
        //         // ParMETIS ожидает целые веса, поэтому округляем или масштабируем
        //         adjwgt_idx[i] = static_cast<int>(adjwgt[i] * 1000);  // Масштабируем для точности
        //     }
        // }
        
        // Целевые веса для каждой части (равномерное распределение)
        std::vector<float> tpwgts(nparts, 1.0 / nparts);
        
        if (csr.rank == 0) {
            std::cout << "=== Starting ParMETIS Partitioning ===" << std::endl;
            std::cout << "Number of parts: " << nparts << std::endl;
            std::cout << "Total vertices: " << csr.vtxdist[csr.size] << std::endl;
            std::cout << "Weight flag: " << wgtflag << std::endl;
        }
        
        // Вызов ParMETIS для разбиения графа на nparts частей
        MPI_Comm comm = MPI_COMM_WORLD;

        ret = ParMETIS_V3_PartKway(
            const_cast<int*>(reinterpret_cast<const int*>(csr.vtxdist.data())),
            const_cast<int*>(reinterpret_cast<const int*>(csr.rows.data())),
            const_cast<int*>(reinterpret_cast<const int*>(csr.cols.data())),
            nullptr,  // vwgt (веса вершин) - не используем
            // (wgtflag == 2 || wgtflag == 3) ? adjwgt_idx.data() : nullptr,  // adjwgt
            nullptr,
            &wgtflag,
            &numflag,
            &ncon,
            &nparts_idx,
            tpwgts.data(), 
            &ubvec,
            options,
            &cutEdge,
            reinterpret_cast<int*>(partition.data()),
            &comm
        );
        
        if (ret != METIS_OK) {
            std::string error_msg;
            switch(ret) {
                case METIS_ERROR_MEMORY:
                    error_msg = "Out of memory";
                    break;
                case METIS_ERROR:
                    error_msg = "General error";
                    break;
                default:
                    error_msg = "Unknown error";
            }
            std::cerr << "Process " << csr.rank 
                    << ": ParMETIS error: " << error_msg << std::endl;
        }
        else if (csr.rank == 0) {
            std::cout << "ParMETIS completed successfully!" << std::endl;
            std::cout << "Edge cut: " << cutEdge << std::endl;
        }
    }

    void run_ptscotch(const CSR& csr) {
        SCOTCH_Strat strat;
        SCOTCH_stratInit(&strat);
        SCOTCH_stratDgraphMapBuild(
            &strat,             // указатель на структуру стратегии
            SCOTCH_STRATBALANCE, // стратегия: балансировка нагрузки
            nparts,             // число частей
            1,                  // менять индексы вершин
            imbalance                // imbalance: допустимое отклонение от идеального баланса
        );

        SCOTCH_Dgraph grafdat;
        SCOTCH_dgraphInit(&grafdat, MPI_COMM_WORLD);
        SCOTCH_Num baseval = 0;
        SCOTCH_Num vertlocnbr = csr.local_rows;
        SCOTCH_Num edgelocnbr = csr.cols.size();

        SCOTCH_dgraphBuild(
            &grafdat,
            0,                          // baseval
            vertlocnbr,                    // vertlocnbr
            vertlocnbr,                    // vertlocmax
            (SCOTCH_Num*)csr.rows.data(),               // vertloctab
            NULL,                       // vendloctab
            NULL,                       // veloloctab
            NULL,                       // vlblloctab
            edgelocnbr,                // edgelocnbr
            edgelocnbr,                // edgelocmax
            (SCOTCH_Num*)csr.cols.data(),              // edgeloctab
            NULL,                       // edgegsttab
            NULL                        // edloloctab
        );

        if(SCOTCH_dgraphCheck(&grafdat) != 0)
        {
            if(csr.rank==0)
                std::cout<<"Graph invalid"<<std::endl;
        }

        if (csr.rank == 0) {
            std::cout << "=== Starting SCOTCH Partitioning ===" << std::endl;
            std::cout << "Number of parts: " << nparts << std::endl;
            std::cout << "Total vertices: " << csr.vtxdist[csr.size] << std::endl;
        }

        ret = SCOTCH_dgraphPart(
            &grafdat,
            nparts,
            &strat,
            partition.data()
        );

        if (ret != 0) {
            std::cerr << "Process " << csr.rank << ": SCOTCH returned error code " << ret << std::endl;
        } else if (csr.rank == 0) {
            std::cout << "SCOTCH completed successfully!" << std::endl;
        }
        
        SCOTCH_dgraphExit(&grafdat);
        SCOTCH_stratExit(&strat);
    }

    void run_ptscotch_map(const CSR& csr) {
        MPI_Comm comm = MPI_COMM_WORLD;

        std::vector<SCOTCH_Num> xadj(csr.rows.begin(), csr.rows.end());
        std::vector<SCOTCH_Num> adjncy(csr.cols.begin(), csr.cols.end());

        SCOTCH_Num vertlocnbr = csr.local_rows;
        SCOTCH_Num edgelocnbr = adjncy.size();

        SCOTCH_Dgraph grafdat;
        SCOTCH_dgraphInit(&grafdat, comm);

        SCOTCH_dgraphBuild(
            &grafdat,
            0,
            vertlocnbr,
            vertlocnbr,
            xadj.data(),
            NULL,
            NULL,
            NULL,
            edgelocnbr,
            edgelocnbr,
            adjncy.data(),
            NULL,
            NULL
        );

        if (SCOTCH_dgraphCheck(&grafdat) != 0) {
            if (csr.rank == 0)
                std::cout << "Graph invalid" << std::endl;
        }

        if (csr.rank == 0)
            std::cout << "SCOTCH_archInit" << std::endl;

        // --- architecture ---
        SCOTCH_Arch arch;
        SCOTCH_archInit(&arch);
        SCOTCH_archCmplt(&arch, nparts);

        if (csr.rank == 0)
            std::cout << "SCOTCH_stratInit" << std::endl;

        // --- strategy ---
        SCOTCH_Strat strat;
        SCOTCH_stratInit(&strat);
        SCOTCH_stratDgraphMapBuild(
            &strat,
            SCOTCH_STRATBALANCE,
            nparts,
            1,
            imbalance
        );

        if (csr.rank == 0)
            std::cout << "SCOTCH completed " << std::endl;

        // --- partition ---
        std::vector<SCOTCH_Num> part_scotch(csr.local_rows);

        int ret = SCOTCH_dgraphMap(
            &grafdat,
            &arch,
            &strat,
            part_scotch.data()
        );

        // --- convert to int ---
        for (int i = 0; i < csr.local_rows; i++)
            partition[i] = (int)part_scotch[i];

        if (ret != 0) {
            std::cerr << "Process " << csr.rank
                    << ": SCOTCH returned error code " << ret << std::endl;
        } else if (csr.rank == 0) {
            std::cout << "SCOTCH completed successfully!" << std::endl;
        }

        // --- cleanup ---
        SCOTCH_archExit(&arch);
        SCOTCH_dgraphExit(&grafdat);
        SCOTCH_stratExit(&strat);

    }

    void run_parhip(const CSR& csr) {
        std::vector<idxtype> vtxdist_id(csr.vtxdist.begin(), csr.vtxdist.end());
        std::vector<idxtype> xadj(csr.rows.begin(), csr.rows.end());
        std::vector<idxtype> adjncy(csr.cols.begin(), csr.cols.end());
        std::vector<idxtype> partitionLong(csr.local_rows);

        double imbalance_parhip = imbalance;
        int seed = 42;
        int mode = 2;
        bool suppress_output = true;//печать доп инфы

        MPI_Comm comm = MPI_COMM_WORLD;

        if (csr.rank == 0) {
            std::cout << "=== Starting ParHIP Partitioning ===" << std::endl;
            std::cout << "Number of parts: " << actualParts << std::endl;
            std::cout << "Total vertices: " << csr.vtxdist[csr.size] << std::endl;
        }

        ParHIPPartitionKWay(vtxdist_id.data(), xadj.data(), adjncy.data(),
                                NULL, NULL, &actualParts, &imbalance_parhip,
                                true, 0, FASTSOCIAL,
                                &cutEdge, partitionLong.data(), &comm);

        if (csr.rank == 0) {
            std::cout << "ParHIP completed successfully!" << std::endl;
            std::cout << "Edge cut: " << cutEdge << std::endl;
        }

        partition.assign(partitionLong.begin(), partitionLong.end());
    }

    void run_kaminpar(const CSR& csr) {
        std::vector<kaminpar::dist::GlobalNodeID> vtxdist(csr.vtxdist.begin(), csr.vtxdist.end());
        std::vector<kaminpar::dist::GlobalEdgeID> xadj(csr.rows.begin(), csr.rows.end());
        std::vector<kaminpar::dist::GlobalNodeID> adjncy(csr.cols.begin(), csr.cols.end());
        std::vector<kaminpar::dist::BlockID> partitionLong(csr.local_rows);

        kaminpar::dKaMinPar dist(MPI_COMM_WORLD, 1, kaminpar::dist::create_default_context());

        dist.copy_graph(
            vtxdist,
            xadj,
            adjncy
        );

        dist.set_k(nparts);
        dist.set_uniform_max_block_weights(imbalance);
        dist.set_output_level(kaminpar::OutputLevel::QUIET);

        if (csr.rank == 0) {
            std::cout << "=== Starting KaMinPar Partitioning ===" << std::endl;
            std::cout << "Number of parts: " << nparts << std::endl;
            std::cout << "Total vertices: " << csr.vtxdist[csr.size] << std::endl;
        }

        kaminpar::dist::EdgeWeight cut = dist.compute_partition(partitionLong);

        if (csr.rank == 0) {
            std::cout << "ParHIP completed successfully!" << std::endl;
            std::cout << "Edge cut: " << cut << std::endl;
        }

        cutEdge = cut;
        partition.assign(partitionLong.begin(), partitionLong.end());
    }

    static std::string getPartitionName(Type partitionType) {
        switch (partitionType) {
            case Type::KAHIP:      return "KaHIP";
            case Type::KAMINPAR:   return "KaMinPar";
            case Type::PARMETIS:   return "ParMETIS";
            case Type::SCOTCH:     return "SCOTCH";
            case Type::ZOLTAN2:    return "Zoltan2";
        }
        return "UNKNOWN";
    }

};

class PartitionMetrics {
public:

    const std::string outputFolder;

    std::unordered_map<int,int> ghost_part;

    std::vector<int> verticesCounter;
    std::vector<double> vertexImbalance;

    std::vector<int> boundaryPerPart;

    std::vector<std::unordered_set<int>> neighbours;

    PartitionMetrics(const std::string outputFolder) : outputFolder(outputFolder) {
    }

    void save_partition_info(const Partition partition, const int actualParts, const CSR& csr) {
        verticesCounter = std::vector<int>(actualParts, 0);
        vertexImbalance = std::vector<double>(actualParts, 0.0);
        boundaryPerPart = std::vector<int>(actualParts, 0);
        neighbours = std::vector<std::unordered_set<int>>(actualParts);
        ghost_part.clear(); 

        computeVertexImbalance(actualParts, partition.partition);
        exchange_ghost_parts(csr, partition.partition, actualParts);
        computeBoundaryVertices(actualParts, csr, partition.partition);
        countNeighboursSets(partition.partition, actualParts, csr);
        
        if(partition.type == Partition::Type::SCOTCH) {
            int totalEdgeParts = 0;
            countCutEdges(partition.partition, actualParts, nullptr, totalEdgeParts, csr);
            if(csr.rank == 0)
                std::cout << "Edge cut: " << totalEdgeParts << std::endl;
            saveMetricsToFile(partition.type, partition.nparts, actualParts, totalEdgeParts, csr.rank, csr.filename);
        }
        else
            saveMetricsToFile(partition.type, partition.nparts, actualParts, partition.cutEdge, csr.rank, csr.filename);
    }

private:

    int findOwnerRank(const int v, const std::vector<int>& vtxdist) {
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

    int getPartitionOfVertex(int v, const CSR& csr, const std::vector<int>& partition) {
        if (csr.vtxdist[csr.rank] <= v && v < csr.vtxdist[csr.rank+1]) {//вершина наша
            return partition[v - csr.vtxdist[csr.rank]];
        } 

        auto it = ghost_part.find(v);
        if (it != ghost_part.end())
            return it->second;
        std::cerr << "Missing ghost vertex " << v <<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        
        return 0;
    }

    void computeVertexImbalance(const int actualParts, const std::vector<int>& partition) {
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

    void exchange_ghost_parts(const CSR& csr, const std::vector<int>& partition, const int actualParts) {
        // 1. Собираем запросы
        std::vector<std::vector<int>> send_requests(csr.size);
        for (int p = 0; p < csr.size; p++) {
            send_requests[p].reserve(actualParts);
        }

        for (int u = 0; u < partition.size(); u++) {
            for (int e = csr.rows[u]; e < csr.rows[u+1]; e++) {
                int v = csr.cols[e];

                int owner = findOwnerRank(v, csr.vtxdist);

                if (owner != csr.rank) {
                    send_requests[owner].push_back(v);
                }
            }
        }

        // ❗ можно убрать дубликаты (важно для производительности)
        for (int p = 0; p < csr.size; p++) {
            sort(send_requests[p].begin(), send_requests[p].end());
            send_requests[p].erase(unique(send_requests[p].begin(),
                                        send_requests[p].end()),
                                send_requests[p].end());
        }

        // 2. Сколько отправляем
        std::vector<int> send_counts(csr.size), recv_counts(csr.size);

        for (int p = 0; p < csr.size; p++)
            send_counts[p] = send_requests[p].size();

        // 3. Обмен размерами
        MPI_Alltoall(send_counts.data(), 1, MPI_INT,
                    recv_counts.data(), 1, MPI_INT,
                    MPI_COMM_WORLD);

        // 4. Подготовка буферов
        std::vector<int> send_displs(csr.size, 0), recv_displs(csr.size, 0);

        int total_send = 0, total_recv = 0;

        for (int p = 0; p < csr.size; p++) {
            send_displs[p] = total_send;
            recv_displs[p] = total_recv;
            total_send += send_counts[p];
            total_recv += recv_counts[p];
        }

        std::vector<int> send_buf(total_send), recv_buf(total_recv);

        // 5. Заполняем send_buf
        for (int p = 0; p < csr.size; p++) {
            copy(send_requests[p].begin(),
                send_requests[p].end(),
                send_buf.begin() + send_displs[p]);
        }

        // 6. Обмен вершинами (запрос)
        MPI_Alltoallv(send_buf.data(), send_counts.data(), send_displs.data(), MPI_INT,
                    recv_buf.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                    MPI_COMM_WORLD);

        // 7. Формируем ответы
        std::vector<int> send_answers(total_recv);

        for (int i = 0; i < total_recv; i++) {
            int v = recv_buf[i];

            int local_v = v - csr.vtxdist[csr.rank]; // перевод в локальный индекс
            send_answers[i] = partition[local_v];
        }

        // 8. Буфер для ответов
        std::vector<int> recv_answers(total_send);

        // 9. Обратный обмен (ответы)
        MPI_Alltoallv(send_answers.data(), recv_counts.data(), recv_displs.data(), MPI_INT,
                    recv_answers.data(), send_counts.data(), send_displs.data(), MPI_INT,
                    MPI_COMM_WORLD);

        // 10. Теперь у нас есть part[v] для всех запрошенных вершин

        // Соберём map: global vertex → part
        for (int p = 0; p < csr.size; p++) {
            for (int i = 0; i < send_counts[p]; i++) {
                int idx = send_displs[p] + i;
                int v = send_buf[idx];
                int pv = recv_answers[idx];
                ghost_part[v] = pv;
            }
        }
    }

    //вершины, которые имеют соседа в другом кластере
    void computeBoundaryVertices(const int actualParts, const CSR& csr, const std::vector<int>& partition) {
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

    std::string getFileName(const std::string& filename) {
        std::filesystem::path p(filename);
        return p.stem().string();
    }

    void saveMetricsToFile(const Partition::Type partitionName, const int nparts, const int actualParts, const int cutEdge, const int rank, const std::string filename) {
        if (rank == 0) {
            std::string outfile = outputFolder + "/" + getFileName(filename) + "_" + Partition::getPartitionName(partitionName) + "_" + "parts_" + std::to_string(nparts) + "_info";
            std::ofstream f(outfile);
            
            f << "Parts: " << nparts << std::endl;
            f << "Actual parts: " << actualParts << std::endl;
            f << "Edge cut: " << cutEdge << std::endl;
            f << std::endl;

            for (int i = 0; i < actualParts; i++) 
                f << "Part " << i << ": vertices=" << verticesCounter[i] << ", imbalance=" << vertexImbalance[i] << "%"  << std::endl;
            f << std::endl;

            for (int i = 0; i < actualParts; i++) 
                f << "Part " << i << ": boundary vertices=" << boundaryPerPart[i] << std::endl;
            f << std::endl;

            for (int i = 0; i < actualParts; i++) 
                f << "Part " << i << ": neighbours=" << neighbours[i].size() << std::endl;
            f << std::endl;

            std::cout << "Partition info saved to: " << outfile << std::endl;
        }
    }

    void countCutEdges(const std::vector<int>& partition, const int actualParts, std::vector<int>* cutEdgesPerPart, int& totalCutEdges, const CSR& csr) {
        std::vector<int> local_cut(actualParts, 0);
        int local_total = 0;

        for (int u = 0; u < partition.size(); u++) {
            int pu = partition[u];
            int global_u = csr.vtxdist[csr.rank] + u;

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

    void countNeighboursSets(const std::vector<int>& partition, const int actualParts, const CSR& csr)
    {
        // --- 1. Локальные множества ---
        std::vector<std::unordered_set<int>> local_neighbours(actualParts);

        for (int u = 0; u < partition.size(); u++) {
            int pu = partition[u];

            for (int i = csr.rows[u]; i < csr.rows[u + 1]; i++) {
                int v = csr.cols[i];
                int pv = getPartitionOfVertex(v, csr, partition);

                if (pu != pv) {
                    local_neighbours[pu].insert(pv);
                    local_neighbours[pv].insert(pu);
                }
            }
        }

        // --- 2. Сериализация ---
        // формат: [p, size, list..., p, size, list...]

        std::vector<int> send_buf;

        for (int p = 0; p < actualParts; p++) {
            send_buf.push_back(p);
            send_buf.push_back(local_neighbours[p].size());

            for (int x : local_neighbours[p])
                send_buf.push_back(x);
        }

        int send_size = send_buf.size();

        // --- 3. Обмен размерами ---
        std::vector<int> recv_sizes(csr.size);

        MPI_Allgather(&send_size, 1, MPI_INT,
                    recv_sizes.data(), 1, MPI_INT,
                    MPI_COMM_WORLD);

        // --- 4. Смещения ---
        std::vector<int> displs(csr.size, 0);
        int total_recv = 0;

        for (int i = 0; i < csr.size; i++) {
            displs[i] = total_recv;
            total_recv += recv_sizes[i];
        }

        std::vector<int> recv_buf(total_recv);

        // --- 5. Обмен данными ---
        MPI_Allgatherv(send_buf.data(), send_size, MPI_INT,
                    recv_buf.data(), recv_sizes.data(), displs.data(), MPI_INT,
                    MPI_COMM_WORLD);

        // --- 6. Восстановление глобальных множеств ---
        int idx = 0;

        while (idx < total_recv) {
            int p = recv_buf[idx++];
            int sz = recv_buf[idx++];
            for (int i = 0; i < sz; i++) {
                int neigh_p = recv_buf[idx++];
                neighbours[p].insert(neigh_p);
            }
        }
    }

};

bool fileOpened(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.good()) {
        return false;
    }

    if (std::filesystem::is_directory(filename)) {
        return false;
    }

    return true;
}

bool directory_exists(const std::string& path) {
    std::filesystem::path dir_path(path);
    return std::filesystem::exists(dir_path) && std::filesystem::is_directory(dir_path);
}

bool correctInitParams(int rank, const std::string filename , const int nparts, const std::string outputFolder) {
     if (!fileOpened(filename)) {
        if (rank == 0) {
            std::cerr << "Cannot open file " << filename << std::endl;
        }
        return false;
    }

    if(!directory_exists(outputFolder)) {
        if (rank == 0) {
            std::cerr << "Cannot open directory " << outputFolder <<std:: endl;
        }
        return false;
    }

    if(nparts <= 0) {
        if (rank == 0) {
            std::cerr << "nparts must be positive: " << nparts <<std:: endl;
        }
        return false;
    }

    return true;
}

void start_partitions(const CSR& csr, PartitionMetrics& partitionMetrics, const int nparts) {
    Partition partition(nparts, 0.05, csr.local_rows);

    partition.run(Partition::Type::PARMETIS, csr);
    partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
    if(csr.rank == 0) {
        std::cout << std::endl;
    }

    partition.run(Partition::Type::SCOTCH, csr);
    partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
    if (csr.rank == 0) {
        std::cout << std::endl;
    }

    partition.run(Partition::Type::KAHIP, csr);
    partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
    if(csr.rank == 0) {
        std::cout << std::endl;
    }

    partition.run(Partition::Type::KAMINPAR, csr);
    partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4) {
        if (rank == 0) {
            std::cout << "Usage: mpirun -np N ./program.exe matrixCSR nparts outputFolder\n";
        }
        MPI_Finalize();
        return 0;
    }

    const std::string filename = std::string(argv[1]);
    const int nparts = atoi(argv[2]);
    const std::string outputFolder = argv[3];

    if(!correctInitParams(rank, filename, nparts, outputFolder)) {
        MPI_Finalize();
        return 0;
    }

    CSR csr(rank, size, filename);
    csr.readGraph();

    PartitionMetrics partitionMetrics(outputFolder);

    start_partitions(csr, partitionMetrics, nparts);

    MPI_Finalize();
    return 0;
}
