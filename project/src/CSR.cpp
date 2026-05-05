#include "CSR.hpp"
#include <iostream>
#include <algorithm>

CSR::CSR(const int rank, const int size, const std::string& filename) 
    : rank(rank), size(size), filename(filename) 
{}

void CSR::readGraph() {
    readCSR();
    addDiagElems();
    deleteDiagElems();
    computeVtxDist();
}

void CSR::printCSR() {
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
                    std::cout << "(" << cols[j] << ") ";
                }
                std::cout << "\n";
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void CSR::printDiagElems() {
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

                        std::cout << "(" << cols[j] << ") ";
                }
                std::cout << "\n";
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void CSR::printVDist() {
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

void CSR::printCSR_diag() {
    for (int u = 0; u < size; u++) {
        if (u == rank) {
            std::cout << "--- [DIAG VERSION] Rank " << rank << " rows [" << start_row 
                    << ", " << end_row << ") ---\n";
            
            for (int i = 0; i < local_rows; i++) {
                int global_row = start_row + i;
                int start = rows_diag[i];
                int end = rows_diag[i + 1]; 
                
                std::cout << "Row " << global_row << ": ";
                for (int j = start; j < end; j++) {
                    std::cout << "(" << cols_diag[j] << ") ";
                }
                std::cout << "\n";
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void CSR::printDiagElems_diag() {
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
                    if (global_row == cols_diag[j]) {
                        std::cout << "(" << cols_diag[j] << ") ";
                    }
                }
                std::cout << "\n";
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void CSR::readCSR() {
    MPI_File fh;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    
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

    rows.resize(local_rows + 1);
    MPI_Offset rows_offset = 2 * sizeof(int);
    MPI_File_read_at_all(fh,
                        rows_offset + start_row * sizeof(int),
                        rows.data(),
                        local_rows + 1,
                        MPI_INT,
                        MPI_STATUS_IGNORE);
    

    /* Убираем смещение */
    int global_start_nnz = rows[0];   
    for (int i = 0; i <= local_rows; i++) {
        rows[i] -= global_start_nnz;
    }
    int local_nnz = rows[local_rows];

    cols.resize(local_nnz);

    MPI_Offset cols_offset = 2 * sizeof(int) + sizeof(int) * (nrows + 1);
    
    MPI_File_read_at_all(fh,
                        cols_offset + global_start_nnz * sizeof(int),
                        cols.data(),
                        local_nnz,
                        MPI_INT,
                        MPI_STATUS_IGNORE);
    
    MPI_File_close(&fh);
}

void CSR::deleteDiagElems() {
    std::vector<int> new_rows(local_rows + 1, 0);
    std::vector<int> new_cols;

    new_cols.reserve(cols.size());

    int nnz_count = 0;

    for (int i = 0; i < local_rows; i++) {
        int global_row = start_row + i;
        int start = rows[i];
        int end = rows[i + 1];

        for (int j = start; j < end; j++) {
            if (cols[j] != global_row) { /* Оставляем только недиагональные элементы */
                new_cols.push_back(cols[j]);
                nnz_count++;
            }
        }
        new_rows[i + 1] = nnz_count; /* Обновляем строку */
    }

    cols = std::move(new_cols);
    rows = std::move(new_rows);
}

void CSR::addDiagElems() {
    rows_diag.assign(local_rows + 1, 0);
    cols_diag.clear();
    cols_diag.reserve(cols.size() + local_rows);

    for (int i = 0; i < local_rows; i++) {
        int global_row = start_row + i;
        int start = rows[i];
        int end = rows[i + 1];

        /* 1. Ищем, где должна была быть диагональ (точка разрыва) */
        int split = start;
        while (split < end && cols[split] < global_row) {
            split++;
        }

        /* 2. Проверяем, есть ли она там уже */
        if (split < end && cols[split] == global_row) {
            /* ДИАГОНАЛЬ ЕСТЬ: копируем всю строку целиком за один раз */
            cols_diag.insert(cols_diag.end(), cols.begin() + start, cols.begin() + end);
        } else {
            /* ДИАГОНАЛИ НЕТ: Копируем блоки "до" и "после", вставляя 1.0 между ними */
            
            /* Копируем всё ДО диагонали */
            cols_diag.insert(cols_diag.end(), cols.begin() + start, cols.begin() + split);
            
            /* Вставляем саму диагональ */
            cols_diag.push_back(global_row);
            
            /* Копируем всё ПОСЛЕ диагонали */
            cols_diag.insert(cols_diag.end(), cols.begin() + split, cols.begin() + end);
        }
        rows_diag[i + 1] = static_cast<int>(cols_diag.size());
    }
}

void CSR::computeVtxDist() {
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