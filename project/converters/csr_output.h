#pragma once

struct csr_mmaped_output : public csr_output {
    csr_output::params_t mtx;
    off_t full_len = 0;
    void *ptr = nullptr;
    int *rows = nullptr;
    int *cols = nullptr;
    double *vals = nullptr;
    double *rhs = nullptr;
    double *sol = nullptr;
    virtual void *get_ptr(int _mode = 0) {
        auto mode = (csr_output::mode_t)_mode;
        switch (mode) {
            case csr_output::mode_t::ROWS_PTR: return rows;
            case csr_output::mode_t::COLS_PTR: return cols;
            case csr_output::mode_t::VALS_PTR: return vals;
            case csr_output::mode_t::RHS_PTR: return rhs;
            case csr_output::mode_t::SOL_PTR: return sol;
        }
        return nullptr;
    }
    bool resize_and_zero() {
        off_t bs = 1024 * (off_t)sizeof(unsigned int);
        char *buf = (char *)calloc(bs, 1);
        assert(buf != nullptr);
        for (off_t rest = (off_t)full_len; rest > 0;) {
            size_t chunksize = (size_t)(rest >= bs ? bs : rest);
            size_t n = fwrite(buf, 1, chunksize, fp);
            if (n != chunksize) {
                perror("fwrite");
                free(buf);
                return false;
            }
            rest -= chunksize;
        }
        if (fseek(fp, 0, SEEK_SET) == -1) {
            perror("fseek");
            return false;
        }
        free(buf);
        return true;
    }

    bool map_file() {
        int fd = fileno(fp);
        ptr = mmap(NULL, full_len, PROT_WRITE, MAP_SHARED, fd, 0);
        if (ptr == MAP_FAILED) {
            perror("mmap");
            ptr = nullptr;
            return false;
        }
        return true;
    }

    bool unmap_file() {
        if (ptr && munmap(ptr, full_len)) {
            perror("munmap");
            return false;
        }
        return true;
    }

    virtual bool init_data(void *_p = nullptr) override {
        auto *p = (csr_output::params_t *)_p;
        mtx = *p;

        // Offsets for main csr file blocks
        off_t hdr_len = sizeof(unsigned int) + sizeof(unsigned int);
        off_t matrix_len = (mtx.nrows + 1) * sizeof(int) + mtx.nnz * sizeof(int) + mtx.nnz * sizeof(double);
        off_t rhs_len = mtx.nrows * sizeof(double);
        off_t sol_len = mtx.nrows * sizeof(double);
        full_len = hdr_len + matrix_len;
        if (mtx.with_rhs_vec) {
            full_len += rhs_len;
        }
        if (mtx.with_solution_vec) {
            full_len += sol_len;
        }

        // Create or truncate output file and fill the full_size bytes pf it with zeroes
        if (!resize_and_zero()) {
            fprintf(stderr, "CSR output file: error in reserving space in the output csr matrix file\n");
            return false;
        }

        // Map the zero'ed file into memory space
        if (!map_file()) {
            fprintf(stderr, "CSR output file: error in mapping the output csr matrix file into memory\n");
            return false;
        }

        // Target pointers to main csr file blocks
        rows = (int *)((char *)ptr + hdr_len);
        cols = (int *)(rows + (mtx.nrows + 1));
        vals = (double *)(cols + mtx.nnz);
        rhs = vals + mtx.nnz;
        sol = rhs + mtx.nrows;

        // Fill in the header values
        rows[-2] = mtx.nrows;
        rows[-1] = mtx.nnz;
        return true;
    }

    virtual bool commit() {
        return unmap_file();
    }
};

struct csr_buffered_output : public csr_output {
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> vals;
    std::vector<double> rhs;
    std::vector<double> sol;
    csr_output::params_t mtx;
    virtual void *get_ptr(int _mode = 0) {
        auto mode = (csr_output::mode_t)_mode;
        switch (mode) {
            case csr_output::mode_t::ROWS_PTR: return rows.data();
            case csr_output::mode_t::COLS_PTR: return cols.data();
            case csr_output::mode_t::VALS_PTR: return vals.data();
            case csr_output::mode_t::RHS_PTR: return rhs.data();
            case csr_output::mode_t::SOL_PTR: return sol.data();
        }
        return nullptr;
    }
    virtual bool init_data(void *_p = nullptr) override {
        auto *_mtx = (csr_output::params_t *)_p;
        mtx = *_mtx;
        try {
            rows.resize(mtx.nrows + 1);
            cols.resize(mtx.nnz);
            vals.resize(mtx.nnz);
            if (mtx.with_rhs_vec) {
                rhs.resize(mtx.nrows);
            }
            if (mtx.with_solution_vec) {
                sol.resize(mtx.nrows);
            }
        }
        catch(std::bad_alloc &ex) {
            fprintf(stderr, "Not enough memory when allocating CSR buffers\n");
            return false;
        }
        return true;
    }

    virtual bool commit() {
        fwrite(&mtx.nrows, sizeof(int), 1, fp);
        fwrite(&mtx.nnz, sizeof(int), 1, fp);
        fwrite(rows.data(), sizeof(int), mtx.nrows + 1, fp);
        fwrite(cols.data(), sizeof(int), mtx.nnz, fp);
        fwrite(vals.data(), sizeof(double), mtx.nnz, fp);
        if (mtx.with_rhs_vec) {
            fwrite(rhs.data(), sizeof(double), mtx.nrows, fp);
        }
        if (mtx.with_solution_vec) {
            fwrite(sol.data(), sizeof(double), mtx.nrows, fp);
        }
        fclose(fp);
        return true;
    }
};

