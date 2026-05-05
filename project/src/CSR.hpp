#ifndef CSR_HPP
#define CSR_HPP

#include <mpi.h>
#include <vector>
#include <string>

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

    /* Стандартные (без диагонали) */
    std::vector<int> rows;
    std::vector<int> cols;

    /* С диагональю (для Sphynx) */
    std::vector<int> rows_diag;
    std::vector<int> cols_diag;

    CSR(const int rank, const int size, const std::string& filename);

    void readGraph();
    void printCSR();
    void printDiagElems();
    void printVDist();
    void printCSR_diag();
    void printDiagElems_diag();

private:
    void readCSR();
    void deleteDiagElems();
    void addDiagElems();
    void computeVtxDist();
};

#endif