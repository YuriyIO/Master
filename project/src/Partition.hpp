#ifndef PARTITION_HPP
#define PARTITION_HPP

#include "CSR.hpp"
#include <vector>
#include <string>

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

    const std::string graphName;
    const std::string outputFolder;
    const std::string vizFolder;
    const int nparts;
    int actualParts;
    const double imbalance;   
    double actualImbalance;   
    std::vector<int> partition;
    int cutEdge;
    int ret;
    Type type;

    const int seed;
    const bool suppress_output;
    std::string ptscotch_strat_name;
    std::string kahip_strat_name;
    std::string kaminpar_strat_name;
    std::string sphynx_problem_type;
    std::string sphynx_preconditioner_type;
    int sphynx_tolerance_pow;
    std::string phg_strat_name;

    double time_total;
    double time_solve;

    Partition(const std::string graphName, const std::string outputFolder, const std::string vizFolder, const int nparts, const double imbalance, const bool suppress_output,const int local_rows);

    void run(const CSR& csr);

    std::string getPartitionInfo() const;
    std::string getPartitionName() const;
    std::string getPartitionStats() const;
    std::string getPartitionPostfix() const;
    std::string getReportPath() const;
    std::string getVizReportPath() const;
    std::string getVizTopologyPath() const;

private:
    void run_parmetis(const CSR& csr);
    void run_ptscotch(const CSR& csr);
    void run_kahip(const CSR& csr);
    void run_kaminpar(const CSR& csr);
    void run_zoltan2_sphynx(const CSR& csr);
    void run_zoltan_phg(const CSR& csr);
};

#endif