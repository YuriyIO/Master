#include "Partition.hpp"

#include <iostream>
#include <cmath>

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
#include <Zoltan2_PartitioningProblem.hpp>


Partition::Partition(const int nparts, const double imbalance, const bool suppress_output, const int local_rows) : 
    nparts(nparts), 
    imbalance(imbalance), 
    suppress_output(suppress_output),
    seed(42) {
    partition.resize(local_rows);
}

void Partition::run(Type type, const CSR& csr) {
    actualParts = nparts;
    actualImbalance = imbalance;
    ret = 0;
    this->type = type; 
    cutEdge = -1;
    postfix = getPartitionPostfix();
    std::fill(partition.begin(), partition.end(), 0);

    if (csr.rank == 0) {
        std::cout << "=== Starting " << getPartitionName() << " Partitioning ===" << std::endl;
        std::cout << getPartitionStats() << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    time_total = MPI_Wtime();

    switch(type) {
        case Type::PARMETIS:        run_parmetis(csr);       break;
        case Type::SCOTCH:          run_ptscotch(csr);       break;
        case Type::KAHIP:           run_kahip(csr);          break;
        case Type::KAMINPAR:        run_kaminpar(csr);       break;
        case Type::ZOLTAN2_SPHYNX:  run_zoltan2_sphynx(csr); break;
        case Type::ZOLTAN_PHG:      run_zoltan_phg(csr);     break;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    time_total = MPI_Wtime() - time_total;

    if(type == Type::PARMETIS && ret != METIS_OK) {
        std::string error_msg = (ret == METIS_ERROR_MEMORY) ? "Out of memory" : 
                    (ret == METIS_ERROR) ? "General error" : "Unknown error";
        std::cerr << "ParMETIS error: " << error_msg << "; rank: " << csr.rank << std::endl;
    }
    else if(type == Type::SCOTCH && ret != 0)
        std::cerr << "SCOTCH returned error code " << ret << "; rank: " << csr.rank << std::endl;
    else if(csr.rank == 0)
        std::cout << getPartitionName() << " completed successfully!" << std::endl;
}

void Partition::run_parmetis(const CSR& csr) {  
    int wgtflag = 0, numflag = 0, ncon = 1, nparts_idx = nparts;
    float ubvec = 1 + imbalance;
    int options[3] = {1, 0, seed};
    std::vector<float> tpwgts(nparts, 1.0 / nparts);
    MPI_Comm comm = MPI_COMM_WORLD;

    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime();
    ret = ParMETIS_V3_PartKway(
        const_cast<int*>(csr.vtxdist.data()), const_cast<int*>(csr.rows.data()),
        const_cast<int*>(csr.cols.data()), nullptr, nullptr, &wgtflag, &numflag,
        &ncon, &nparts_idx, tpwgts.data(), &ubvec, options, &cutEdge,
        reinterpret_cast<int*>(partition.data()), &comm
    );
    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime() - time_solve;
}

void Partition::run_ptscotch(const CSR& csr) {
    SCOTCH_randomSeed(seed); 
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

    MPI_Barrier(MPI_COMM_WORLD);                
    time_solve = MPI_Wtime();
    ret = SCOTCH_dgraphPart(&grafdat, nparts, &strat, partition.data());    
    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime() - time_solve;

    SCOTCH_dgraphExit(&grafdat);
    SCOTCH_stratExit(&strat);
}

void Partition::run_kahip(const CSR& csr) {
    std::vector<idxtype> vtxdist_id(csr.vtxdist.begin(), csr.vtxdist.end());
    std::vector<idxtype> xadj(csr.rows.begin(), csr.rows.end());
    std::vector<idxtype> adjncy(csr.cols.begin(), csr.cols.end());
    std::vector<idxtype> partitionLong(csr.local_rows);
    MPI_Comm comm = MPI_COMM_WORLD;
    int strat = (kahip_strat_name == "ULTRAFASTMESH") ? 0 : 
                (kahip_strat_name == "FASTMESH") ? 1 : 2;

    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime();
    ParHIPPartitionKWay(vtxdist_id.data(), xadj.data(), adjncy.data(),
                            NULL, NULL, &actualParts, &actualImbalance,
                            suppress_output, seed, strat,
                            &cutEdge, partitionLong.data(), &comm);
    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime() - time_solve;

    partition.assign(partitionLong.begin(), partitionLong.end());
}

void Partition::run_kaminpar(const CSR& csr) {
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
    dist.reseed(seed);
    dist.set_uniform_max_block_weights(imbalance);
    if(suppress_output) dist.set_output_level(kaminpar::OutputLevel::QUIET);

    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime();
    cutEdge = dist.compute_partition(pL);
    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime() - time_solve;

    partition.assign(pL.begin(), pL.end());
}

void Partition::run_zoltan2_sphynx(const CSR& csr) {
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

    if (suppress_output) {
        params.set("debug_level", "no_status");
        sP->set("sphynx_verbosity", 0); 
    }

    Zoltan2::SphynxProblem<adapter_t> problem(&adapter, &params, sP);

    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime();
    problem.solve();
    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime() - time_solve;

    auto sol = problem.getSolution();
    const int* pLV = sol.getPartListView();
    for (size_t i = 0; i < (size_t)csr.local_rows; ++i) partition[i] = (int)(pLV[i]);
}

void Partition::run_zoltan_phg(const CSR& csr) {
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
    for (int i = 0; i < csr.local_rows; ++i) nE[i] = csr.rows[i+1] - csr.rows[i];

    auto graph = Teuchos::rcp(new graph_t(map, nE));
    std::vector<GO> buf; 
    for (int i = 0; i < csr.local_rows; ++i) {
        GO g_row = csr.start_row + i;
        size_t r_start = csr.rows[i], r_len = nE[i];
        buf.resize(r_len);
        for (size_t j = 0; j < r_len; ++j) buf[j] = static_cast<GO>(csr.cols[r_start + j]);
        graph->insertGlobalIndices(g_row, Teuchos::ArrayView<const GO>(buf.data(), r_len));
    }
    graph->fillComplete();

    adapter_t adapter(graph, 0);
    
    Teuchos::ParameterList params;
    params.set("algorithm", "phg");
    params.set("num_global_parts", nparts);
    params.set("imbalance_tolerance", 1.0 + imbalance);

    Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters");

    zparams.set("SEED", std::to_string(seed));
    zparams.set("DETERMINISTIC", "TRUE");  
    if(phg_strat_name == "strong") {
        zparams.set("PHG_EDGE_SIZE_THRESHOLD", "1.0");      /*не пропускаем вершины с больщим сичлом связей. порог пропуска по умолчанию 25%*/ 
        /*zparams.set("PHG_COARSENING_METHOD", "ipm");*/
        zparams.set("PHG_COARSENING_METHOD_FAST", "ipm"); 
        zparams.set("PHG_COARSENING_LIMIT", "500");           /*порог остановки сжатия. по умолчанию 100*/
        zparams.set("PHG_REFINEMENT_QUALITY", "5.0");        /* Коэффициент «упорства» при уточнении границ. по умолчанию 1*/
        zparams.set("PHG_REFINEMENT_LOOP_LIMIT", "20");     /*Количество итераций уточнения на каждом уровне. по умолчанию 10*/
    }
    
    if(suppress_output) zparams.set("DEBUG_LEVEL", "0");

    Zoltan2::PartitioningProblem<adapter_t> problem(&adapter, &params);

    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime();
    problem.solve();
    MPI_Barrier(MPI_COMM_WORLD);
    time_solve = MPI_Wtime() - time_solve;
 
    auto sol = problem.getSolution();
    const int* pLV = sol.getPartListView();
    for (size_t i = 0; i < (size_t)csr.local_rows; ++i) partition[i] = (int)(pLV[i]);
}

std::string Partition::getPartitionInfo() const{
    switch (type) {
        case Type::KAHIP:      return "ParHIP " + kahip_strat_name;
        case Type::KAMINPAR:   return "dKaMinPar " + kaminpar_strat_name;
        case Type::PARMETIS:   return "ParMETIS";
        case Type::SCOTCH:     return "PT-SCOTCH " + ptscotch_strat_name;
        case Type::ZOLTAN2_SPHYNX: return "Zoltan2 " + sphynx_problem_type + " " + sphynx_preconditioner_type + " " + std::to_string(sphynx_tolerance_pow);
        case Type::ZOLTAN_PHG: return "Zoltan";
    }
    return "UNKNOWN"; 
}

std::string Partition::getPartitionName() const{
    switch (type) {
        case Type::KAHIP:      return "KaHIP";
        case Type::KAMINPAR:   return "KaMinPar";
        case Type::PARMETIS:   return "ParMETIS";
        case Type::SCOTCH:     return "SCOTCH";
        case Type::ZOLTAN2_SPHYNX: return "Zoltan2";
        case Type::ZOLTAN_PHG: return "Zoltan";
    }
    return "UNKNOWN"; 
}

std::string Partition::getPartitionStats() const {
    switch (type) {
        case Type::KAHIP:           return "kahip_strat_name: " + kahip_strat_name;
        case Type::KAMINPAR:        return "kaminpar_strat_name: " + kaminpar_strat_name;
        case Type::PARMETIS:        return "ParMETIS strat: default";
        case Type::SCOTCH:          return "ptscotch_strat_name: " + ptscotch_strat_name;
        case Type::ZOLTAN2_SPHYNX:  return "sphynx_problem_type: " + sphynx_problem_type + "\n" 
                                        + "sphynx_preconditioner_type: " + sphynx_preconditioner_type + "\n" 
                                        + "sphynx_tolerance_pow: " + std::to_string(sphynx_tolerance_pow);
        case Type::ZOLTAN_PHG:      return "phg_strat_name: " + phg_strat_name;
    }
    return "UNKNOWN";
}

std::string Partition::getPartitionPostfix() const {
    int imbalancePercent = (int) (imbalance * 100);
    switch (type) {
        case Type::KAHIP:           return kahip_strat_name + "_imbalance_" + std::to_string(imbalancePercent);
        case Type::KAMINPAR:        return kaminpar_strat_name + "_imbalance_" + std::to_string(imbalancePercent);
        case Type::PARMETIS:        return "imbalance_" + std::to_string(imbalancePercent);
        case Type::SCOTCH:          return ptscotch_strat_name + "_imbalance_" + std::to_string(imbalancePercent);
        case Type::ZOLTAN2_SPHYNX:  return sphynx_problem_type + "_" + sphynx_preconditioner_type + "_pow_" + std::to_string(sphynx_tolerance_pow);
        case Type::ZOLTAN_PHG:      return "phg_strat_name_" + phg_strat_name + "_imbalance_" + std::to_string(imbalancePercent);
    }
    return "UNKNOWN";
}
