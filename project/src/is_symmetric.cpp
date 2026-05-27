
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <mpi.h>

class SymmetryChecker {
    int rank, size;
    int nrows, nnz;
    int local_rows, start_row, end_row;
    std::vector<int> rows, cols;
    std::vector<int> vtxdist;
    std::string filename;

public:
    SymmetryChecker(int r, int s, const std::string& fname) 
        : rank(r), size(s), filename(fname) {}

    void readGraph() {
        MPI_File fh;
        if (MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS) {
            if (rank == 0) std::cerr << "Error: Cannot open file " << filename << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        int header[2];
        if (rank == 0) {
            MPI_File_read_at(fh, 0, header, 2, MPI_INT, MPI_STATUS_IGNORE);
        }
        MPI_Bcast(header, 2, MPI_INT, 0, MPI_COMM_WORLD);

        nrows = header[0];
        nnz = header[1];

        int base = nrows / size;
        int rem = nrows % size;
        local_rows = base + (rank < rem ? 1 : 0);
        start_row = rank * base + std::min(rank, rem);
        end_row = start_row + local_rows;

        vtxdist.resize(size + 1);
        std::vector<int> all_local_rows(size);
        MPI_Allgather(&local_rows, 1, MPI_INT, all_local_rows.data(), 1, MPI_INT, MPI_COMM_WORLD);
        vtxdist[0] = 0;
        for (int i = 0; i < size; i++) vtxdist[i + 1] = vtxdist[i] + all_local_rows[i];

        rows.resize(local_rows + 1);
        MPI_Offset rows_offset = 2 * sizeof(int);
        MPI_File_read_at_all(fh, rows_offset + start_row * sizeof(int), rows.data(), local_rows + 1, MPI_INT, MPI_STATUS_IGNORE);

        int global_start_nnz = rows[0];
        for (int i = 0; i <= local_rows; i++) rows[i] -= global_start_nnz;
        int local_nnz = rows[local_rows];

        cols.resize(local_nnz);
        MPI_Offset cols_offset = 2 * sizeof(int) + sizeof(int) * (nrows + 1);
        MPI_File_read_at_all(fh, cols_offset + (MPI_Offset)global_start_nnz * sizeof(int), cols.data(), local_nnz, MPI_INT, MPI_STATUS_IGNORE);

        MPI_File_close(&fh);
    }

    bool checkSymmetry() {
        for (int i = 0; i < local_rows; i++) {
            std::sort(cols.begin() + rows[i], cols.begin() + rows[i + 1]);
        }

        std::vector<std::vector<int>> send_pkts(size);
        for (int i = 0; i < local_rows; i++) {
            int global_i = start_row + i;
            for (int j_idx = rows[i]; j_idx < rows[i + 1]; j_idx++) {
                int global_j = cols[j_idx];

                int target_rank = 0;
                while (target_rank < size - 1 && global_j >= vtxdist[target_rank + 1]) {
                    target_rank++;
                }

                send_pkts[target_rank].push_back(global_j);
                send_pkts[target_rank].push_back(global_i);
            }
        }

        std::vector<int> send_counts(size), send_displs(size, 0);
        std::vector<int> recv_counts(size), recv_displs(size, 0);
        std::vector<int> flat_send;

        for (int i = 0; i < size; i++) {
            send_counts[i] = send_pkts[i].size();
            flat_send.insert(flat_send.end(), send_pkts[i].begin(), send_pkts[i].end());
            if (i > 0) send_displs[i] = send_displs[i - 1] + send_counts[i - 1];
        }

        MPI_Alltoall(send_counts.data(), 1, MPI_INT, recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);

        int total_recv = 0;
        for (int i = 0; i < size; i++) {
            recv_displs[i] = total_recv;
            total_recv += recv_counts[i];
        }

        std::vector<int> flat_recv(total_recv);
        MPI_Alltoallv(flat_send.data(), send_counts.data(), send_displs.data(), MPI_INT,
                      flat_recv.data(), recv_counts.data(), recv_displs.data(), MPI_INT, MPI_COMM_WORLD);

        bool local_symmetric = true;
        for (int i = 0; i < total_recv; i += 2) {
            int r = flat_recv[i];
            int c = flat_recv[i + 1];
            
            int local_r_idx = r - start_row;
            auto row_begin = cols.begin() + rows[local_r_idx];
            auto row_end = cols.begin() + rows[local_r_idx + 1];

            if (!std::binary_search(row_begin, row_end, c)) {
                local_symmetric = false;
                break;
            }
        }

        int global_res;
        int local_res = local_symmetric ? 1 : 0;
        MPI_Allreduce(&local_res, &global_res, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

        return global_res == 1;
    }
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 2) {
        if (rank == 0) std::cout << "Usage: " << argv[0] << " <graph_file.bin>" << std::endl;
        MPI_Finalize();
        return 1;
    }

    std::string filename = argv[1];
    SymmetryChecker checker(rank, size, filename);

    checker.readGraph();
    bool is_sym = checker.checkSymmetry();

    if (rank == 0) {
        if (is_sym) {
            std::cout << "Граф симметричен (Symmetric)" << std::endl;
        } else {
            std::cout << "Граф несимметричен (NOT Symmetric)" << std::endl;
        }
    }

    MPI_Finalize();
    return 0;
}