#include "CSR.hpp"
#include "Partition.hpp"
#include "PartitionMetrics.hpp"

#include <Tpetra_Core.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>

bool fileOpened(const std::string& filename) {
    std::ifstream file(filename);
    return file.good() && !std::filesystem::is_directory(filename);
}

bool directory_exists(const std::string& path) {
    std::filesystem::path dir_path(path);
    return std::filesystem::exists(dir_path) && std::filesystem::is_directory(dir_path);
}

bool correctInitParams(const int rank, const std::string filename , const int nparts, const std::string outputFolder, const double imbalance, const int suppress_output_num, const int graph_viz_num, const std::string vizFolder) {
     if (!fileOpened(filename)) {
        if (rank == 0) std::cerr << "Cannot open file " << filename << std::endl;
        return false;
    }
    if(!directory_exists(outputFolder)) {
        if (rank == 0) std::cerr << "Cannot open directory " << outputFolder <<std:: endl;
        return false;
    }
    if(nparts <= 0) {
        if (rank == 0) std::cerr << "nparts must be positive: " << nparts <<std:: endl;
        return false;
    }
    if(imbalance <= 0 || imbalance >= 1) {
        if (rank == 0) std::cerr << "imbalance must be positive and less than 1: " << imbalance <<std:: endl;
        return false;
    }
    if(suppress_output_num != 0 && suppress_output_num != 1) {
        if (rank == 0) std::cerr << "suppress_output must be 0 or 1 but we have: " << suppress_output_num <<std:: endl;
        return false;
    }
    if(graph_viz_num != 0 && graph_viz_num != 1) {
        if (rank == 0) std::cerr << "graph_viz_num must be 0 or 1 but we have: " << graph_viz_num <<std:: endl;
        return false;
    }
    if(!directory_exists(vizFolder)) {
        if (rank == 0) std::cerr << "Cannot open directory " << vizFolder <<std:: endl;
        return false;
    }
    return true;
}

bool partition_exist(const std::string& reportPath, const int rank) {
    int skip = 0;
    if (rank == 0) {
        if (std::filesystem::exists(reportPath)) {
            skip = 1;
        }
    }
    MPI_Bcast(&skip, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return skip == 1;
}

void start_partitions(const CSR& csr, Partition& partition, PartitionMetrics& partitionMetrics, const bool graph_viz) {
    if(graph_viz) {
        if(!partition_exist(partition.getVizTopologyPath(), csr.rank))
            partitionMetrics.save_graph_topology(partition, csr);
        else if (csr.rank == 0)
            std::cout << "[SKIP] " << partition.getVizTopologyPath() << " already exists\n" << std::endl;
    }

    auto run_partition_and_save_stats = [&csr, &partition, &partitionMetrics, graph_viz]() {
        if (!partition_exist(partition.getReportPath(), csr.rank)) {
            partition.run(csr);
            partitionMetrics.save_partition_info(partition, partition.actualParts, csr);
            if (graph_viz) {
                partitionMetrics.save_partition_vec(partition, csr);
            }
            if (csr.rank == 0) {
                std::cout << std::endl;
            }
        }
        else if (csr.rank == 0) {
            std::cout << "[SKIP] " << partition.getReportPath() << " already exists\n" << std::endl;
        }
    };

    /*ParMETIS*/
    partition.type = Partition::Type::PARMETIS;
    run_partition_and_save_stats();

    /*PT-SCOTCH*/
    partition.type = Partition::Type::SCOTCH;
    std::vector<std::string> scotch_strats = {"SCOTCH_STRATQUALITY", "SCOTCH_STRATBALANCE", "SCOTCH_STRATSPEED"};
    for (const auto& strat : scotch_strats) {
        partition.ptscotch_strat_name = strat;
        run_partition_and_save_stats();
    }

    /*ParHIP*/
    partition.type = Partition::Type::KAHIP;
    std::vector<std::string> kahip_strats = {"ULTRAFASTMESH", "FASTMESH"};/*{"ULTRAFASTMESH", "FASTMESH", "ECOMESH"}; ECOMESH зависает на маленьких графах */
    for (const auto& strat : kahip_strats) {
        partition.kahip_strat_name = strat;
        run_partition_and_save_stats();
    }

    /*dKaMinPar*/
    partition.type = Partition::Type::KAMINPAR;
    std::vector<std::string> kaminpar_strats = {"create_default_context", "create_strong_context"};/* {"create_default_context", "create_strong_context", "create_xterapart_context"}; */
    for (const auto& strat : kaminpar_strats) {
        partition.kaminpar_strat_name = strat;
        run_partition_and_save_stats();
    }

    /*Zoltan PHG*/
    partition.type = Partition::Type::ZOLTAN_PHG;
    partition.phg_strat_name = "default";
    run_partition_and_save_stats();


    /*Zoltan2 sphynx*/
    partition.type = Partition::Type::ZOLTAN2_SPHYNX;
    std::vector<std::string> sphynx_types = {"combinatorial", "generalized", "normalized"};/*{"combinatorial", "generalized", "normalized"};*/
    std::vector<std::string> sphynx_preconds = {"muelu"};/*{"muelu", "jacobi", "polynomial"};*/
    std::vector<int> tolerances = {4, 6, 8}; /*{4, 6, 8};*/

    for(const auto& type : sphynx_types) {
        for(const auto& prec : sphynx_preconds) {
            for(int tol : tolerances) {
                partition.sphynx_problem_type = type;
                partition.sphynx_preconditioner_type = prec;
                partition.sphynx_tolerance_pow = tol;
                run_partition_and_save_stats();
            }
        }
    }

    if (csr.rank == 0) {
        std::cout << "Program completed successfully" <<std::endl;
    }

}

int main(int argc, char* argv[]) {
    Tpetra::ScopeGuard tpetraScope(&argc, &argv);
    {
        int rank, size;
        auto comm = Tpetra::getDefaultComm();
        rank = comm->getRank();
        size = comm->getSize();

        if (argc < 8) {
            if (rank == 0) {
                std::cout << "Usage: mpirun -np N ./program.exe matrixCSR nparts outputFolder imbalance(0.05) suppress_output(1) graph_viz(0) vizFolder\n";
            }
            return 0;
        }

        const std::string filename = std::string(argv[1]);
        const int nparts = atoi(argv[2]);
        const std::string outputFolder = argv[3];
        const double imbalance = std::atof(argv[4]);
        const int suppress_output_num = atoi(argv[5]);
        const int graph_viz_num = atoi(argv[6]);
        const std::string vizFolder = argv[7];

        if(!correctInitParams(rank, filename, nparts, outputFolder, imbalance, suppress_output_num, graph_viz_num, vizFolder)) {
            return 0;
        }
        const bool suppress_output = suppress_output_num == 1;
        const bool graph_viz = graph_viz_num == 1;

        CSR csr(rank, size, filename);
        csr.readGraph();

        PartitionMetrics partitionMetrics;

        std::string graphName = std::filesystem::path(filename).stem().string();
        Partition partition(graphName, outputFolder, vizFolder, nparts, imbalance, suppress_output, csr.local_rows);

        start_partitions(csr, partition, partitionMetrics, graph_viz);

    }
    return 0; 
}