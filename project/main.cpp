#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include <filesystem>
#include <iomanip>

#include <parmetis.h>
#include <parhip_interface.h>
#include <ptscotch.h>
#include <dkaminpar.h>

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_Core.hpp>
#include <Zoltan2_XpetraCrsGraphAdapter.hpp>
#include <Zoltan2_SphynxProblem.hpp>

#include <zoltan.h>

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

    // Стандартные (без диагонали)
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;

    // С диагональю (для Sphynx)
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;
    std::vector<double> vals_diag;

    CSR(const int rank, const int size, const std::string& filename) 
        : rank(rank), size(size), filename(filename) 
    {}

    void readGraph() {
        readCSR();
        addDiagElems();
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

    void printCSR_diag() {
        // Вывод всей структуры CSR_diag для проверки
        for (int u = 0; u < size; u++) {
            if (u == rank) {
                std::cout << "--- [DIAG VERSION] Rank " << rank << " rows [" << start_row 
                        << ", " << end_row << ") ---\n";
                
                for (int i = 0; i < local_rows; i++) {
                    int global_row = start_row + i;
                    int start = rows_diag[i];     // Используем rows_diag
                    int end = rows_diag[i + 1];   // Используем rows_diag
                    
                    std::cout << "Row " << global_row << ": ";
                    for (int j = start; j < end; j++) {
                        std::cout << "(" << cols_diag[j] << ", "  // Используем cols_diag
                                << vals_diag[j] << ") ";        // Используем vals_diag
                    }
                    std::cout << "\n";
                }
                fflush(stdout);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    void printDiagElems_diag() {
        // Вывод ТОЛЬКО диагональных элементов из структуры CSR_diag
        for (int u = 0; u < size; u++) {
            if (u == rank) {
                std::cout << "--- [ONLY DIAGONALS] Rank " << rank << " rows [" << start_row 
                        << ", " << end_row << ") ---\n";
                
                for (int i = 0; i < local_rows; i++) {
                    int global_row = start_row + i;
                    int start = rows_diag[i];
                    int end = rows_diag[i + 1];
                    
                    std::cout << "Row " << global_row << ": ";
                    for (int j = start; j < end; j++) {
                        // Если индекс столбца совпадает с глобальным индексом строки
                        if (global_row == cols_diag[j]) {
                            std::cout << "(" << cols_diag[j] << ", " << vals_diag[j] << ") ";
                        }
                    }
                    std::cout << "\n";
                }
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

    void addDiagElems() {
        rows_diag.assign(local_rows + 1, 0);
        cols_diag.clear();
        vals_diag.clear();
        cols_diag.reserve(cols.size() + local_rows);
        vals_diag.reserve(vals.size() + local_rows);

        for (int i = 0; i < local_rows; i++) {
            int global_row = start_row + i;
            int start = rows[i];
            int end = rows[i + 1];

            // 1. Ищем, где должна была быть диагональ (точка разрыва)
            int split = start;
            while (split < end && cols[split] < global_row) {
                split++;
            }

            // 2. Проверяем, есть ли она там уже
            if (split < end && cols[split] == global_row) {
                // ДИАГОНАЛЬ ЕСТЬ: Просто копируем всю строку целиком за один раз
                cols_diag.insert(cols_diag.end(), cols.begin() + start, cols.begin() + end);
                vals_diag.insert(vals_diag.end(), vals.begin() + start, vals.begin() + end);
            } else {
                // ДИАГОНАЛИ НЕТ: Копируем блоки "до" и "после", вставляя 1.0 между ними
                
                // Копируем всё ДО диагонали
                cols_diag.insert(cols_diag.end(), cols.begin() + start, cols.begin() + split);
                vals_diag.insert(vals_diag.end(), vals.begin() + start, vals.begin() + split);
                
                // Вставляем саму диагональ
                cols_diag.push_back(global_row);
                vals_diag.push_back(1.0);
                
                // Копируем всё ПОСЛЕ диагонали
                cols_diag.insert(cols_diag.end(), cols.begin() + split, cols.begin() + end);
                vals_diag.insert(vals_diag.end(), vals.begin() + split, vals.begin() + end);
            }
            rows_diag[i + 1] = static_cast<int>(cols_diag.size());
        }
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
        ZOLTAN2_SPHYNX,
        ZOLTAN_PHG
    };

    const int nparts;
    int actualParts;
    const double imbalance;   
    double actualImbalance;   
    std::vector<int> partition;
    int cutEdge;
    int ret;
    Type type;

    int seed;
    bool suppress_output = true;
    std::string ptscotch_strat_name;
    std::string kahip_strat_name;
    std::string kaminpar_strat_name;
    std::string sphynx_problem_type;
    std::string sphynx_preconditioner_type;
    int sphynx_tolerance_pow;

    std::string postfix;

    double time_total;
    double time_solve;

    Partition(const int nparts, const double imbalance, const int local_rows) : nparts(nparts), imbalance(imbalance) {
        partition.resize(local_rows);
        seed = 42; 
    }

    void run(Type type, const CSR& csr) {
        actualParts = nparts;
        actualImbalance = imbalance;
        ret = 0;
        this->type = type; 
        cutEdge = -1;
        postfix = getPartitionPostfix(type);

        if (csr.rank == 0) {
            std::cout << "=== Starting " << getPartitionName(type) << " Partitioning ===" << std::endl;
            std::cout << getPartitionStats(type) << std::endl;
        }

        double loc_total_start, loc_total_end;
        double loc_solve_duration = 0;

        MPI_Barrier(MPI_COMM_WORLD);
        loc_total_start = MPI_Wtime();

        switch(type) {
            case Type::PARMETIS:        run_parmetis(csr, loc_solve_duration);      break;
            case Type::SCOTCH:          run_ptscotch(csr, loc_solve_duration);      break;
            case Type::KAHIP:           run_kahip(csr, loc_solve_duration);        break;
            case Type::KAMINPAR:        run_kaminpar(csr, loc_solve_duration);      break;
            case Type::ZOLTAN2_SPHYNX:  run_zoltan2_sphynx(csr, loc_solve_duration); break;
        }

        loc_total_end = MPI_Wtime();
        double loc_total_duration = loc_total_end - loc_total_start;

        if(type == Type::PARMETIS && ret != METIS_OK) {
            std::string error_msg = (ret == METIS_ERROR_MEMORY) ? "Out of memory" : 
                        (ret == METIS_ERROR) ? "General error" : "Unknown error";
            std::cerr << "ParMETIS error: " << error_msg << "; rank: " << csr.rank << std::endl;
        }
        else if(type == Type::SCOTCH && ret != 0)
            std::cerr << "SCOTCH returned error code " << ret << "; rank: " << csr.rank << std::endl;
        else if(csr.rank == 0)
            std::cout << getPartitionName(type) << " completed successfully!" << std::endl;
        
        MPI_Reduce(&loc_total_duration, &this->time_total, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&loc_solve_duration, &this->time_solve, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    }

    void run_parmetis(const CSR& csr, double& loc_solve_time) {  
        int wgtflag = 0, numflag = 0, ncon = 1, nparts_idx = nparts;
        float ubvec = 1 + imbalance;
        int options[3] = {1, 0, seed};
        std::vector<float> tpwgts(nparts, 1.0 / nparts);
        MPI_Comm comm = MPI_COMM_WORLD;

        double t1 = MPI_Wtime();
        ret = ParMETIS_V3_PartKway(
            const_cast<int*>(csr.vtxdist.data()), const_cast<int*>(csr.rows.data()),
            const_cast<int*>(csr.cols.data()), nullptr, nullptr, &wgtflag, &numflag,
            &ncon, &nparts_idx, tpwgts.data(), &ubvec, options, &cutEdge,
            reinterpret_cast<int*>(partition.data()), &comm
        );
        loc_solve_time = MPI_Wtime() - t1;
    }

    void run_ptscotch(const CSR& csr, double& loc_solve_time) {
        SCOTCH_Strat strat;
        SCOTCH_stratInit(&strat);
        int strat_flag = (ptscotch_strat_name == "SCOTCH_STRATSPEED") ? SCOTCH_STRATSPEED : 
                         (ptscotch_strat_name == "SCOTCH_STRATBALANCE") ? SCOTCH_STRATBALANCE : SCOTCH_STRATQUALITY;

        SCOTCH_stratDgraphMapBuild(&strat, strat_flag, nparts, 1, imbalance);
        SCOTCH_Dgraph grafdat;
        SCOTCH_dgraphInit(&grafdat, MPI_COMM_WORLD);
        SCOTCH_dgraphBuild(&grafdat, 0, csr.local_rows, csr.local_rows, (SCOTCH_Num*)csr.rows.data(),
                           NULL, NULL, NULL, csr.cols.size(), csr.cols.size(), 
                           (SCOTCH_Num*)csr.cols.data(), NULL, NULL);

        double t1 = MPI_Wtime();
        ret = SCOTCH_dgraphPart(&grafdat, nparts, &strat, partition.data());    
        loc_solve_time = MPI_Wtime() - t1;

        SCOTCH_dgraphExit(&grafdat);
        SCOTCH_stratExit(&strat);
    }

    void run_kahip(const CSR& csr, double& loc_solve_time) {
        std::vector<idxtype> vtxdist_id(csr.vtxdist.begin(), csr.vtxdist.end());
        std::vector<idxtype> xadj(csr.rows.begin(), csr.rows.end());
        std::vector<idxtype> adjncy(csr.cols.begin(), csr.cols.end());
        std::vector<idxtype> partitionLong(csr.local_rows);
        MPI_Comm comm = MPI_COMM_WORLD;
        int strat = (kahip_strat_name == "ULTRAFASTMESH") ? 0 : 
                    (kahip_strat_name == "FASTMESH") ? 1 : 2;

        double t1 = MPI_Wtime();
        ParHIPPartitionKWay(vtxdist_id.data(), xadj.data(), adjncy.data(),
                                NULL, NULL, &actualParts, &actualImbalance,
                                suppress_output, seed, strat,
                                &cutEdge, partitionLong.data(), &comm);
        loc_solve_time = MPI_Wtime() - t1;
        partition.assign(partitionLong.begin(), partitionLong.end());
    }

    void run_kaminpar(const CSR& csr, double& loc_solve_time) {
        std::vector<kaminpar::dist::GlobalNodeID> vd(csr.vtxdist.begin(), csr.vtxdist.end());
        std::vector<kaminpar::dist::GlobalEdgeID> xa(csr.rows.begin(), csr.rows.end());
        std::vector<kaminpar::dist::GlobalNodeID> ad(csr.cols.begin(), csr.cols.end());
        std::vector<kaminpar::dist::BlockID> pL(csr.local_rows);

        kaminpar::dist::Context context = (kaminpar_strat_name == "create_default_context") ? 
                                           kaminpar::dist::create_default_context() :
                                          (kaminpar_strat_name == "create_strong_context") ? 
                                           kaminpar::dist::create_strong_context() : 
                                           kaminpar::dist::create_xterapart_context();

        kaminpar::dKaMinPar dist(MPI_COMM_WORLD, 1, context);
        dist.copy_graph(vd, xa, ad);
        dist.set_k(nparts);
        dist.set_uniform_max_block_weights(imbalance);
        if(suppress_output) dist.set_output_level(kaminpar::OutputLevel::QUIET);

        double t1 = MPI_Wtime();
        cutEdge = dist.compute_partition(pL);
        loc_solve_time = MPI_Wtime() - t1;

        partition.assign(pL.begin(), pL.end());
    }

    void run_zoltan2_sphynx(const CSR& csr, double& loc_solve_time) {
        typedef Tpetra::CrsGraph<> graph_t;
        typedef typename graph_t::local_ordinal_type LO;
        typedef typename graph_t::global_ordinal_type GO;
        typedef typename graph_t::node_type Node;
        typedef Tpetra::Map<LO, GO, Node> map_t;
        typedef Zoltan2::XpetraCrsGraphAdapter<graph_t> adapter_t;

        auto tcomm = Tpetra::getDefaultComm();
        GO g_el = csr.vtxdist[csr.size];
        auto map = Teuchos::rcp(new map_t(g_el, csr.local_rows, 0, tcomm));

        Teuchos::Array<size_t> nE(csr.local_rows);
        for (int i = 0; i < csr.local_rows; ++i) nE[i] = csr.rows_diag[i+1] - csr.rows_diag[i];

        auto graph = Teuchos::rcp(new graph_t(map, nE));
        std::vector<GO> buf; 
        for (int i = 0; i < csr.local_rows; ++i) {
            GO g_row = csr.start_row + i;
            size_t r_start = csr.rows_diag[i], r_len = nE[i];
            buf.resize(r_len);
            for (size_t j = 0; j < r_len; ++j) buf[j] = static_cast<GO>(csr.cols_diag[r_start + j]);
            graph->insertGlobalIndices(g_row, Teuchos::ArrayView<const GO>(buf.data(), r_len));
        }
        graph->fillComplete();

        adapter_t adapter(graph, 0);
        Teuchos::ParameterList params;
        params.set("num_global_parts", nparts);

        Teuchos::RCP<Teuchos::ParameterList> sP = Teuchos::rcp(new Teuchos::ParameterList());
        sP->set("sphynx_skip_preprocessing", true);
        sP->set("sphynx_problem_type", sphynx_problem_type);
        sP->set("sphynx_preconditioner_type", sphynx_preconditioner_type);
        sP->set("sphynx_tolerance", 1.0 / std::pow(10.0, sphynx_tolerance_pow));

        Zoltan2::SphynxProblem<adapter_t> problem(&adapter, &params, sP);

        double t1 = MPI_Wtime();
        problem.solve();
        loc_solve_time = MPI_Wtime() - t1;

        auto sol = problem.getSolution();
        const int* pLV = sol.getPartListView();
        for (size_t i = 0; i < (size_t)csr.local_rows; ++i) partition[i] = (int)(pLV[i]);
    }

    static std::string getPartitionName(Type partitionType) {
        switch (partitionType) {
            case Type::KAHIP:      return "KaHIP";
            case Type::KAMINPAR:   return "KaMinPar";
            case Type::PARMETIS:   return "ParMETIS";
            case Type::SCOTCH:     return "SCOTCH";
            case Type::ZOLTAN2_SPHYNX: return "Zoltan2_Sphynx";
            case Type::ZOLTAN_PHG: return "Zoltan_PHG";
        }
        return "UNKNOWN";
    }

    std::string getPartitionStats(Type partitionType) const {
        switch (partitionType) {
            case Type::KAHIP:           return "kahip_strat_name: " + kahip_strat_name;
            case Type::KAMINPAR:        return "kaminpar_strat_name: " + kaminpar_strat_name;
            case Type::PARMETIS:        return "ParMETIS strat: default";
            case Type::SCOTCH:          return "ptscotch_strat_name: " + ptscotch_strat_name;
            case Type::ZOLTAN2_SPHYNX:  return "sphynx_problem_type: " + sphynx_problem_type + "\n" 
                                             + "sphynx_preconditioner_type: " + sphynx_preconditioner_type + "\n" 
                                             + "sphynx_tolerance_pow: " + std::to_string(sphynx_tolerance_pow);
        }
        return "UNKNOWN";
    }

    std::string getPartitionPostfix(Type partitionType) const {
        int imbalancePercent = (int) (imbalance * 100);
        switch (partitionType) {
            case Type::KAHIP:           return kahip_strat_name + "_imbalance_" + std::to_string(imbalancePercent);
            case Type::KAMINPAR:        return kaminpar_strat_name + "_imbalance_" + std::to_string(imbalancePercent);
            case Type::PARMETIS:        return "imbalance_" + std::to_string(imbalancePercent);
            case Type::SCOTCH:          return ptscotch_strat_name + "_imbalance_" + std::to_string(imbalancePercent);
            case Type::ZOLTAN2_SPHYNX:  return sphynx_problem_type + "_" + sphynx_preconditioner_type + "_" + std::to_string(sphynx_tolerance_pow);
        }
        return "UNKNOWN";
    }

};


class PartitionMetrics {
public:

    const std::string outputFolder;

    std::vector<int> ghost_global_ids;
    std::vector<int> ghost_part_ids;

    std::vector<int> verticesCounter;
    std::vector<double> vertexImbalance;

    std::vector<int> boundaryPerPart;

    std::vector<std::vector<int>> neighbours;

    // --- НОВЫЕ МЕТРИКИ ---
    std::vector<long long> edgesPerPart;    // Количество ребер в каждом разделе
    std::vector<double> edgeImbalance;      // Дисбаланс ребер в % (насколько нагружены связи)
    int maxCommDegree;                      // Максимальное количество соседей у одного раздела
    long long totalCommVolume;              // Суммарный объем пересылаемых данных (Total Ghost Nodes)

    PartitionMetrics(const std::string outputFolder) : outputFolder(outputFolder) {
    }

    void save_partition_info(const Partition& partition, const int actualParts, const CSR& csr) {
        verticesCounter.assign(actualParts, 0);
        vertexImbalance.assign(actualParts, 0.0);
        edgesPerPart.assign(actualParts, 0);
        edgeImbalance.assign(actualParts, 0.0);
        boundaryPerPart.assign(actualParts, 0);
        neighbours.clear();
        ghost_global_ids.clear(); 
        ghost_part_ids.clear();

        // 1. Расчет дисбаланса вершин (нагрузка на вычисления CPU)
        computeVertexImbalance(partition.partition, actualParts);
        
        // 2. Расчет дисбаланса ребер (нагрузка на память и локальные операции)
        computeEdgeImbalance(partition.partition, actualParts, csr);
        
        // 3. Обмен ghost-узлами
        exchange_ghost_parts(partition.partition, actualParts, csr);
        
        // 4. Расчет Communication Volume (общая нагрузка на сеть)
        // В MPI-приложениях объем данных равен количеству уникальных "призрачных" вершин
        computeCommunicationVolume(partition.partition, actualParts, csr);

        // 5. Расчет граничных вершин (вершины, инициирующие обмен)
        computeBoundaryVertices(partition.partition, actualParts, csr);
        
        // 6. Расчет графа связности разделов
        countNeighbours(partition.partition, actualParts, csr);
        
        // 7. Максимальная степень связности (латентность сети)
        computeMaxCommDegree();
        
         // 8. Расчет Edge Cut и Locality
        int finalCut = partition.cutEdge;
        if (partition.type == Partition::Type::SCOTCH || partition.type == Partition::Type::ZOLTAN2_SPHYNX) {
            int totalEdgeParts = 0;
            countCutEdges(partition.partition, actualParts, nullptr, totalEdgeParts, csr);
            finalCut = totalEdgeParts;
        }
        if(csr.rank == 0)
            std::cout << "Edge cut: " << finalCut << std::endl;

        // 9. Сохранение итогового отчета
        saveMetricsToFile(partition, finalCut, csr.rank, csr.filename);
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
        // 1. Проверяем, не наша ли это вершина
        if (csr.vtxdist[csr.rank] <= v && v < csr.vtxdist[csr.rank + 1]) {
            return partition[v - csr.vtxdist[csr.rank]];
        }

        // 2. Бинарный поиск в отсортированном векторе призрачных вершин
        auto it = std::lower_bound(ghost_global_ids.begin(), ghost_global_ids.end(), v);

        if (it != ghost_global_ids.end() && *it == v) {
            // Вычисляем индекс найденного элемента
            size_t idx = std::distance(ghost_global_ids.begin(), it);
            return ghost_part_ids[idx];
        }

        std::cerr << "Rank " << csr.rank << ": Missing ghost vertex " << v << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 0;
    }

    void computeVertexImbalance(const std::vector<int>& partition, const int actualParts) {
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

    void exchange_ghost_parts(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
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
        int total_ghosts = send_buf.size(); // это total_send из твоего кода
        ghost_global_ids.resize(total_ghosts);
        ghost_part_ids.resize(total_ghosts);

        // Создаем вектор индексов [0, 1, 2, ..., N-1] для косвенной сортировки
        std::vector<int> p(total_ghosts);
        std::iota(p.begin(), p.end(), 0);

        // Сортируем индексы на основе значений в send_buf
        std::sort(p.begin(), p.end(), [&](int i, int j) {
            return send_buf[i] < send_buf[j];
        });

        // Заполняем итоговые векторы в правильном (отсортированном) порядке
        for (int i = 0; i < total_ghosts; i++) {
            ghost_global_ids[i] = send_buf[p[i]];
            ghost_part_ids[i] = recv_answers[p[i]];
        }
    }

    //вершины, которые имеют соседа в другом кластере
    void computeBoundaryVertices(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
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

    struct Edge {
        int p1, p2;
        // Операторы для правильной сортировки и удаления дубликатов
        bool operator<(const Edge& other) const {
            if (p1 != other.p1) return p1 < other.p1;
            return p2 < other.p2;
        }
        bool operator==(const Edge& other) const {
            return p1 == other.p1 && p2 == other.p2;
        }
    };

    void countNeighbours(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
        // 1. Собираем локальные пары (Edge)
        std::vector<Edge> local_edges;
        for (int u = 0; u < (int)partition.size(); u++) {
            int pu = partition[u];
            for (int i = csr.rows[u]; i < csr.rows[u + 1]; i++) {
                int v = csr.cols[i];
                int pv = getPartitionOfVertex(v, csr, partition);
                if (pu != pv) {
                    local_edges.push_back({std::min(pu, pv), std::max(pu, pv)});
                }
            }
        }

        // Удаляем локальные дубликаты
        std::sort(local_edges.begin(), local_edges.end());
        local_edges.erase(std::unique(local_edges.begin(), local_edges.end()), local_edges.end());

        // 2. Подготовка к MPI_Gatherv
        int local_edge_count = (int)local_edges.size();
        std::vector<int> recv_counts_in_edges;
        if (csr.rank == 0) recv_counts_in_edges.resize(csr.size);

        // Сначала узнаем, сколько Edge-ей у каждого процесса
        MPI_Gather(&local_edge_count, 1, MPI_INT, recv_counts_in_edges.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Массивы, которые нужны только Root-у для приема
        std::vector<int> counts_int;
        std::vector<int> displs_int;
        std::vector<Edge> all_edges;

        if (csr.rank == 0) {
            counts_int.resize(csr.size);
            displs_int.resize(csr.size);
            int total_ints = 0;
            for (int i = 0; i < csr.size; i++) {
                counts_int[i] = recv_counts_in_edges[i] * 2; // Каждое ребро - это 2 int-а
                displs_int[i] = total_ints;
                total_ints += counts_int[i];
            }
            all_edges.resize(total_ints / 2); // Резервируем место под структуры Edge
        }

        // ВНИМАНИЕ: ОДИН ЕДИНСТВЕННЫЙ ВЫЗОВ ДЛЯ ВСЕХ
        // На не-root рангах counts_int.data() и displs_int.data() могут быть nullptr (т.к. векторы пустые)
        // Согласно стандарту MPI, на не-root рангах эти аргументы игнорируются.
        MPI_Gatherv(local_edges.data(), local_edge_count * 2, MPI_INT,
                    all_edges.data(), counts_int.data(), displs_int.data(), MPI_INT,
                    0, MPI_COMM_WORLD);

        // 3. Обработка результатов (только Rank 0)
        if (csr.rank == 0) {
            neighbours.assign(actualParts, std::vector<int>());
            
            // Удаляем дубликаты, пришедшие от разных процессов
            std::sort(all_edges.begin(), all_edges.end());
            all_edges.erase(std::unique(all_edges.begin(), all_edges.end()), all_edges.end());

            for (const auto& e : all_edges) {
                neighbours[e.p1].push_back(e.p2);
                neighbours[e.p2].push_back(e.p1);
            }
            
            // Сортируем списки соседей для каждого раздела
            for (auto& list : neighbours) {
                std::sort(list.begin(), list.end());
            }
        }
    }

    void computeEdgeImbalance(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
        std::vector<long long> local_edges(actualParts, 0);
        for (int i = 0; i < csr.local_rows; i++) {
            local_edges[partition[i]] += (csr.rows[i+1] - csr.rows[i]);
        }
        
        MPI_Allreduce(local_edges.data(), edgesPerPart.data(), actualParts, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        long long total_e = 0;
        for (long long e : edgesPerPart) total_e += e;
        double avg_e = (double)total_e / actualParts;

        for (int i = 0; i < actualParts; i++) {
            edgeImbalance[i] = (avg_e > 0) ? (edgesPerPart[i] / avg_e - 1.0) * 100.0 : 0.0;
        }
    }

    void computeCommunicationVolume(const std::vector<int>& partition, const int actualParts, const CSR& csr) {
        // 1. Для каждого локального раздела (pu) соберем список уникальных глобальных ID вершин (v),
        // которые ему нужны из других разделов.
        std::vector<std::vector<int>> unique_ghosts_per_part(actualParts);

        for (int u = 0; u < csr.local_rows; u++) {
            int pu = partition[u];
            for (int i = csr.rows[u]; i < csr.rows[u + 1]; i++) {
                int v = csr.cols[i];
                int pv = getPartitionOfVertex(v, csr, partition);

                if (pu != pv) {
                    // Разделу pu нужна вершина v из чужого раздела
                    unique_ghosts_per_part[pu].push_back(v);
                }
            }
        }

        // 2. Оставляем только уникальные ID для каждого раздела на текущем ранге
        long long local_total_vol = 0;
        for (int p = 0; p < actualParts; p++) {
            if (unique_ghosts_per_part[p].empty()) continue;

            std::sort(unique_ghosts_per_part[p].begin(), unique_ghosts_per_part[p].end());
            unique_ghosts_per_part[p].erase(
                std::unique(unique_ghosts_per_part[p].begin(), unique_ghosts_per_part[p].end()), 
                unique_ghosts_per_part[p].end()
            );

            local_total_vol += unique_ghosts_per_part[p].size();
        }

        // 3. Суммируем по всем MPI-процессам. 
        // Примечание: Если один раздел (Partition ID) разнесен по разным MPI-рангам, 
        // эта метрика покажет суммарный объем входящих данных для всех частей этого раздела.
        MPI_Allreduce(&local_total_vol, &totalCommVolume, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    }

    void computeMaxCommDegree() {
        int local_max = 0;
        // На Rank 0 переменная neighbours заполнена полностью после countNeighbours
        for (const auto& list : neighbours) {
            local_max = std::max(local_max, (int)list.size());
        }
        maxCommDegree = local_max;
    }

     void saveMetricsToFile(const Partition& partition, const int cutEdge, const int rank, const std::string filename) {
        if (rank == 0) {
            std::string outfile = outputFolder + "/" + getFileName(filename) + "_" + Partition::getPartitionName(partition.type) 
                                                + "_p" + std::to_string(partition.nparts) + "_" + partition.postfix + "_info";
            std::ofstream f(outfile);
            if (!f.is_open()) {
                std::cout << "ERROR: Failed to open file for writing: " << outfile << std::endl;
                return;
            }

            // 1. Вычисляем агрегированную статистику
            int maxBoundary = 0;
            double maxVertexImbalance = 0.0;
            double maxEdgeImbalance = 0.0;
            long long totalVertices = 0;
            long long totalEdges = 0;

            for (int i = 0; i < partition.actualParts; i++) {
                if (std::abs(vertexImbalance[i]) > std::abs(maxVertexImbalance)) 
                    maxVertexImbalance = vertexImbalance[i];
                
                if (std::abs(edgeImbalance[i]) > std::abs(maxEdgeImbalance)) 
                    maxEdgeImbalance = edgeImbalance[i];
                
                if (boundaryPerPart[i] > maxBoundary)
                    maxBoundary = boundaryPerPart[i];
                
                totalVertices += verticesCounter[i];
                totalEdges += edgesPerPart[i];
            }

            f << "======================================================================" << std::endl;
            f << " PARTITION REPORT: " << Partition::getPartitionName(partition.type) << std::endl;
            f << "======================================================================" << std::endl;
            f << "File:             " << filename << std::endl;
            f << partition.getPartitionStats(partition.type) << std::endl;
            f << "Seed:             " << partition.seed << std::endl;
            f << "Total Edge Vol:   " << totalEdges << std::endl;
            f << "Total Vertex Vol: " << totalVertices << std::endl;
            f << "Target Parts:     " << partition.nparts << std::endl;
            f << "Actual Parts:     " << partition.actualParts << std::endl;
            f << "Target Imbalance: " << (partition.imbalance * 100.0) << "%" << std::endl;
            f << "----------------------------------------------------------------------" << std::endl;
            f << "Edge Cut:              " << cutEdge << std::endl;
            f << "Max Neighbour Degree:  " << maxCommDegree << std::endl;
            f << "Max Boundary Degree:   " << maxBoundary << std::endl;
            f << "Total Comm Vol:        " << totalCommVolume << std::endl;
            f << "Max Vertex imbalance:  " << maxVertexImbalance << std::endl;
            f << "Max Edge imbalance:    " << maxEdgeImbalance << std::endl;    
            f << "Time Total (sec):      " << std::fixed << std::setprecision(4) << partition.time_total << std::endl;
            f << "Time Solve Only (sec): " << std::fixed << std::setprecision(4) << partition.time_solve << std::endl;
            f << "Overhead (sec):        " << std::fixed << std::setprecision(4) << (partition.time_total - partition.time_solve) << std::endl;
            f << "======================================================================" << std::endl;
            f << std::endl;

            f << std::left << std::setw(8)  << "Part" 
              << std::setw(12) << "Vertices" 
              << std::setw(12) << "V-Imb %" 
              << std::setw(12) << "Edges" 
              << std::setw(12) << "E-Imb %" 
              << std::setw(12) << "Boundary" 
              << std::setw(10) << "Neighbors" << std::endl;
            f << "----------------------------------------------------------------------" << std::endl;

            f << std::fixed << std::setprecision(2);
            for (int i = 0; i < partition.actualParts; i++) {
                f << std::left << std::setw(8)  << i 
                  << std::setw(12) << verticesCounter[i] 
                  << std::setw(12) << vertexImbalance[i] 
                  << std::setw(12) << edgesPerPart[i] 
                  << std::setw(12) << edgeImbalance[i] 
                  << std::setw(12) << boundaryPerPart[i] 
                  << std::setw(10) << neighbours[i].size() << std::endl;
            }
            f << "======================================================================" << std::endl;

            f.flush();
            f.close();
            std::cout << "Partition report saved to: " << outfile << std::endl;
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


    std::vector<std::string> scotch_strats = {"SCOTCH_STRATQUALITY", "SCOTCH_STRATBALANCE", "SCOTCH_STRATSPEED"};
    for (const auto& strat : scotch_strats) {
        partition.ptscotch_strat_name = strat;
        partition.run(Partition::Type::SCOTCH, csr);
        partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
        if (csr.rank == 0) {
            std::cout << std::endl;
        }
    }

    std::vector<std::string> kahip_strats = {"ULTRAFASTMESH", "FASTMESH"};//{"ULTRAFASTMESH", "FASTMESH", "ECOMESH"}; ECOMESH зависает на маленьких графах
    for (const auto& strat : kahip_strats) {
        partition.kahip_strat_name = strat;
        partition.run(Partition::Type::KAHIP, csr);
        partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
        if (csr.rank == 0) {
            std::cout << std::endl;
        }
    }

    std::vector<std::string> kaminpar_strats = {"create_default_context", "create_strong_context"};// {"create_default_context", "create_strong_context", "create_xterapart_context"}; create_xterapart_context зависает на маленьких графах
    for (const auto& strat : kaminpar_strats) {
        partition.kaminpar_strat_name = strat;
        partition.run(Partition::Type::KAMINPAR, csr);
        partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
        if (csr.rank == 0) {
            std::cout << std::endl;
        }
    }


    std::vector<std::string> sphynx_types = {"combinatorial", "generalized", "normalized"};
    std::vector<std::string> sphynx_preconds = {"muelu", "jacobi", "polynomial"};
    std::vector<int> tolerances = {6}; //{4, 6, 8};

    for (const auto& type : sphynx_types) {
        for (const auto& prec : sphynx_preconds) {
            for (int tol : tolerances) {
                partition.sphynx_problem_type = type;
                partition.sphynx_preconditioner_type = prec;
                partition.sphynx_tolerance_pow = tol;

                partition.run(Partition::Type::ZOLTAN2_SPHYNX, csr);
                partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
                if (csr.rank == 0) {
                    std::cout << std::endl;
                }
            }
        }
    }

}

int main(int argc, char* argv[]) {
    // ScopeGuard инициализирует MPI и Kokkos
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    
    { // Открываем блок, чтобы все объекты удалились до конца main
        int rank, size;
        auto comm = Tpetra::getDefaultComm();
        rank = comm->getRank();
        size = comm->getSize();

        if (argc < 4) {
            if (rank == 0) {
                std::cout << "Usage: mpirun -np N ./program.exe matrixCSR nparts outputFolder\n";
            }
            return 0;
        }

        const std::string filename = std::string(argv[1]);
        const int nparts = atoi(argv[2]);
        const std::string outputFolder = argv[3];

        if(!correctInitParams(rank, filename, nparts, outputFolder)) {
            return 0;
        }

        CSR csr(rank, size, filename);
        csr.readGraph();

        PartitionMetrics partitionMetrics(outputFolder);
        start_partitions(csr, partitionMetrics, nparts);

    } // Здесь вызываются деструкторы csr, partitionMetrics и т.д.

    return 0; 
    // Здесь ScopeGuard автоматически завершит Kokkos, а затем MPI.
}

// int main(int argc, char* argv[]) {
//     MPI_Init(&argc, &argv);
//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     if (argc < 4) {
//         if (rank == 0) {
//             std::cout << "Usage: mpirun -np N ./program.exe matrixCSR nparts outputFolder\n";
//         }
//         MPI_Finalize();
//         return 0;
//     }

//     const std::string filename = std::string(argv[1]);
//     const int nparts = atoi(argv[2]);
//     const std::string outputFolder = argv[3];

//     if(!correctInitParams(rank, filename, nparts, outputFolder)) {
//         MPI_Finalize();
//         return 0;
//     }

//     CSR csr(rank, size, filename);
//     csr.readGraph();

//     PartitionMetrics partitionMetrics(outputFolder);

//     start_partitions(csr, partitionMetrics, nparts);

//     MPI_Finalize();
//     return 0;
// }
