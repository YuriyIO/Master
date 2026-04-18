#pragma once

#include <functional>

struct mtx_input : public input {
    using input::fp;
    std::set<std::string> flags;
    bool checked_getline_header(std::string &line) {
        size_t bsz = 1024;
        char *buf = (char *)calloc(bsz, 1);
        assert(buf != nullptr);
        ssize_t errcode;
        do {
            if ((errcode = getline(&buf, &bsz, fp)) == -1) {
                char msg[101];
                strerror_r(errno, msg, 100);
                fprintf(stderr, "MTX file reading error; error=%s\n", msg);
                return false;
            }
        } while (buf[0] == '%');
        line = buf;
        free(buf);
        return true;
    }
    bool read_header_flags() {
        size_t n = 100;
        char *head = (char *)calloc(n, sizeof(char));
        assert(head != nullptr);
        read_arr<char>(head, 2);
        if (head[0] != '%' || head[1] != '%') {
            fprintf(stderr, "MTX file format error, invalid header\n");
            return false;
        }
        ssize_t errcode;
        if ((errcode = getline(&head, &n, fp)) == -1) {
            char msg[101] = {0,};
            strerror_r(errno, msg, 100);
            fprintf(stderr, "MTX file format error, header read error: error=%s\n", msg);
            return false;
        }
        assert((size_t)errcode == strlen(head));
        head[errcode-1] = 0;
        auto strvec = helpers::str_split(head, ' ');
        for (const auto &s : strvec) {
            flags.insert(s);
        }
        if (flags.count("MatrixMarket")) {
            flags.insert("valid");
        }
        else {
            fprintf(stderr, "MTX file format error, invalid header\n");
            return false;
        }
        return true;
    }
    virtual ~mtx_input() {}
};

struct mtx_matrix_input : public mtx_input {
    using input::fp;
    using mtx_input::flags;
    int nrows = 0;
    int ncols = 0;
    int nnz = 0;
    int nnz_file = 0;
    int ndiag = 0;
    virtual int get_nrows() override { return nrows; }
    virtual int get_ncols() override { return ncols; }
    virtual int get_nnz() override { return nnz; }

    bool traverse_mtx_file(std::function<bool(int, int, double, int)> value_handler) {
        int file_counter = 0;
        size_t bsz = 1024;
        char *buf = (char *)calloc(bsz, 1);
        do {
            int linelen;
            int r, c;
            double v;
            if (!checked_getline_IID(buf, bsz, r, c, v, linelen)) {
                fprintf(stderr, "MTX file getline error, counter=%d\n", file_counter);
                free(buf);
                return false;
            }
			if (!value_handler(r, c, v, linelen)) {
				free(buf);
				return false;
			}
            if (feof(fp)) {
                if (file_counter != nnz_file-1) {
                    fprintf(stderr, "MTX file read error: matrix: unexpected EOF: counter=%d nnz-1=%d\n", file_counter, nnz_file-1);
                    free(buf);
                    return false;
                }
                break;
            }
            file_counter++;
            if (file_counter == nnz_file) {
                break;
            }
        } while(1);
        free(buf);
        return true;
    }

    virtual int get_ndiag() override {
		int nd = 0;
        off_t pos = ftello(fp);
        if (pos == -1) {
            perror("ftello");
            fprintf(stderr, "MTX file read error: matrix: unexpected EOF: cannot ftello input stream\n");
            return -1;
        }

		bool result = traverse_mtx_file(
            [&](int r, int c, double v, int linelen) -> bool { 
                (void)linelen; 
                (void )v;
                if (r == c) 
                    nd++; 
                return true; 
            }
        );

        fseeko(fp, pos, SEEK_SET);
		if (!result)
			nd = -1;
		return nd;
    }
    
    virtual bool read_header(void *_p = nullptr) override {
        (void)_p;
        if (!read_header_flags()) {
            return false;
        }
        if (flags.count("valid") == 0) {
            return false;
        }
        if (flags.count("matrix") == 0) {
            fprintf(stderr, "MTX file: not a matrix file\n");
            return false;
        }
        if (flags.count("real") == 0 && flags.count("integer") == 0 && flags.count("pattern") == 0) {
            fprintf(stderr, "MTX file: 'real', 'integer' or 'pattern' matrix type required\n");
            return false;
        }
        if (flags.count("coordinate") == 0) {
            fprintf(stderr, "MTX file: 'coordinate' MTX type required\n");
            return false;
        }
        std::string line;
        checked_getline_header(line);
        unsigned long long unr, unc, unz;
        int n = sscanf(line.c_str(), "%llu %llu %llu\n", &unr, &unc, &unz);
        if (n != 3) {
            fprintf(stderr, "MTX file format error (matrix dimensions)\n");
            return false;
        }
        nrows = unr;
        ncols = unc;
        nnz_file = unz;
        if (flags.count("symmetric") == 1) {
            ndiag = get_ndiag();
            if (ndiag == -1) {
                return false;
            }
            nnz = nnz_file * 2 - ndiag;
        } else if (flags.count("skew-symmetric") == 1) {
            nnz = nnz_file * 2;
        } else {
            nnz = nnz_file;
        }
        return true;
    }
    virtual ~mtx_matrix_input() {}
};


struct mtx_matrix_direct_input : public mtx_matrix_input {
    using map_item_t = std::pair<int, off_t>;
    std::vector<map_item_t> mtxmap;
    virtual bool init_data(void *params = nullptr) {
        (void)params;
        if (flags.count("symmetric") == 1 || flags.count("skew-symmetric") == 1) {
            fprintf(stderr, "MTX file read error: matrix: 'direct' mode reader cannot handle 'symmetric' MTX-files! Use buffered mode.\n");
            return false;
        }
        if (options & OMIT_SORTING_STAGE)
            return true;
        if (!fill_in_file_map()) {
            return false;
        }
        sort_file_map();
        return true;
    }
    bool fill_in_file_map() {
        mtxmap.reserve(nnz);
        off_t pos = ftello(fp);
        if (pos == -1) {
            perror("ftello");
            fprintf(stderr, "MTX file read error: matrix: unexpected EOF: cannot ftello input stream\n");
            return false;
        }
        bool result = traverse_mtx_file(
            [&](int r, int c, double v, int linelen) -> bool { 
                (void)c;
                (void)v;
                mtxmap.push_back({r, pos});
                pos += linelen;
                return true; 
            }
        );

        return result;
    }

    void sort_file_map() {
        std::stable_sort(mtxmap.begin(), mtxmap.end(), [](map_item_t r, map_item_t l) -> bool { return r.first < l.first; });
    }

    virtual bool copy_to_array(double *) { return true; }
    virtual bool copy_to_csr_arrays(int *rows, int *cols, double *vals) override {
        (void)ncols;
        size_t bsz = 1024;
        char *buf = (char *)calloc(bsz, 1);
        assert(buf != nullptr);
        int curr_row = 0, counter = 0;
        int linelen = 0;
        off_t oldpos = 0;
        bool with_map = (mtxmap.size() > 0);
        bool pattern = (flags.count("pattern") > 0);
        do {
            if (with_map && oldpos + linelen != mtxmap[counter].second) {
                if (fseeko(fp, mtxmap[counter].second, SEEK_SET) == -1) {
                    perror("fseeko");
                    fprintf(stderr, "MTX file read error: matrix: unexpected EOF: cannot fseeko input stream\n");
                }
            }
            if (with_map) {
                oldpos = mtxmap[counter].second;
            }
            int r, c; double v;
            if (pattern) {
                v = 1.0;
                if (!checked_getline_II(buf, bsz, r, c, linelen)) {
                    fprintf(stderr, "MTX file getline error, counter=%d\n", counter);
                    free(buf);
                    return false;
                }
            } else {
                if (!checked_getline_IID(buf, bsz, r, c, v, linelen)) {
                    fprintf(stderr, "MTX file getline error, counter=%d\n", counter);
                    free(buf);
                    return false;
                }
            }
            if (!with_map && feof(fp)) {
                if (counter != nnz_file-1) {
                    fprintf(stderr, "MTX file read error: matrix: unexpected EOF: counter=%d nnz-1=%d\n", counter, nnz_file-1);
                    free(buf);
                    return false;
                }
                assert(counter == nnz_file-1);
                assert(curr_row == nrows);
                rows[nrows] = nnz;
                cols[counter] = c - 1;
                vals[counter] = v;
                break;
            }
            if (curr_row == 0)
                curr_row = r;
            else if (curr_row != r) {
                if (curr_row > r) {
                    fprintf(stderr, "MTX file format error: matrix rows are not in ascending order, counter=%d\n", counter);
                    free(buf);
                    return false;
                }
                curr_row = r;
                rows[r-1] = counter;
            }
            cols[counter] = c - 1;
            vals[counter] = v;
            counter++;
            if (counter == nnz_file) {
                assert(curr_row == nrows);
                rows[nrows] = nnz;
                break;
            }
        } while(1);
        free(buf);
        return true;
    }
    virtual ~mtx_matrix_direct_input() {}
};

struct mtx_matrix_buffered_input : public mtx_matrix_input {
    std::vector<int> r;
    std::vector<int> c;
    std::vector<double> v;
    virtual bool init_data(void *params = nullptr) {
        (void)params;
        try {
            r.resize(nnz);
            c.resize(nnz);
            v.resize(nnz);
        }
        catch(std::bad_alloc &ex) {
            fprintf(stderr, "Not enough memory when allocating MTX buffers\n");
            return false;
        }
		int counter = 0;
        bool result = traverse_mtx_file(
            [&](int r_, int c_, double v_, int linelen) -> bool {
                (void)linelen;
                r[counter] = r_; 
				c[counter] = c_; 
				v[counter] = v_;
				if (flags.count("symmetric") || flags.count("skew-symmetric")) {
					if (r[counter] != c[counter]) {
						++counter;
						r[counter] = c[counter-1];
						c[counter] = r[counter-1];
						v[counter] = v[counter-1];
					}
				}
				++counter;
                return true; 
            }
        );
        for (int i = 0; i < nnz; ++i) {
            --r[i];
            --c[i];
        }
        return result;
    }

    virtual bool copy_to_array(double *) { return true; }
    virtual bool copy_to_csr_arrays(int *rows, int *cols, double *vals) override {
        std::vector<int> row_cntr(nrows, 0);
        for (int i = 0; i < nnz; ++i)
            ++row_cntr[r[i]];
        rows[0] = 0;
        for (int i = 1; i < nrows+1; ++i) {
            rows[i] = rows[i - 1] + row_cntr[i - 1];
            row_cntr[i - 1] = 0;
        }
        for (int i = 0; i < nnz; ++i) {
            size_t idx = rows[r[i]] + (row_cntr[r[i]]++);
            cols[idx] = c[i];
            vals[idx] = v[i];
        }
        return true;
    }
};


struct mtx_vector_input : public mtx_input {
    struct params_t {
        int nrows;
    };
    using input::fp;
    using mtx_input::flags;
    int nrows;
    virtual int get_nrows() override { return nrows; }
    virtual int get_ncols() override { return 0; }
    virtual int get_nnz() override { return 0; }
    virtual int get_ndiag() override { return 0; }
    virtual bool init_data(void *params = nullptr) override { (void)params; return true; }
    virtual bool copy_to_csr_arrays(int *, int *, double *) override { 
        return true; 
    }

    bool read_header(void *_p) {
        auto *p = (params_t *)_p;
        nrows = p->nrows;
        if (!read_header_flags()) {
            return false;
        }
        if (flags.count("valid") == 0) {
            return false;
        }
        if (flags.count("vector") == 0) {
            fprintf(stderr, "MTX file: not a vector file\n");
            return false;
        }
        std::string line;
        checked_getline_header(line);
        unsigned long long nelems; //, x;
        int n = sscanf(line.c_str(), "%llu\n", &nelems);
        if (n != 1) {
            fprintf(stderr, "MTX file format error (vector dimentions)\n");
            return false;
        }
        if ((int)nelems != nrows) {
            fprintf(stderr, "MTX file error: vector size is not equal to matrix rows number\n");
            return false;
        }
        return true;
    }

    virtual bool copy_to_array(double *vec) override {
        int counter = 0;
        size_t bsz = 1024;
        char *buf = (char *)calloc(bsz, 1);
        assert(buf != nullptr);
        do {
            int linelen;
            double v; int n;
            if (!checked_getline_ID(buf, bsz, n, v, linelen)) {
                fprintf(stderr, "MTX file getline error, counter=%d\n", counter);
                free(buf);
                return false;
            }
            if (feof(fp)) {
                if (counter != nrows-1) {
                    fprintf(stderr, "MTX file read error: vector: unexpected EOF: counter=%d nrows-1=%llu\n", counter, (unsigned long long)(nrows-1));
                    return false;
                }
                break;
            }
            vec[counter] = v;
            counter++;
            if (counter == nrows)
                break;
        } while (1);
        free(buf);
        return true;
    }
    virtual ~mtx_vector_input() {}
};

struct mtx_matrixcol_input : public mtx_matrix_buffered_input {
    struct params_t {
        int column;
    };

    int column = 0;

    virtual bool read_header(void *_p = nullptr) override {
        (void)_p;
        if (!read_header_flags()) {
            return false;
        }
        if (flags.count("valid") == 0) {
            return false;
        }
        if (flags.count("matrix") == 0) {
            fprintf(stderr, "MTX file: not a matrix file\n");
            return false;
        }
        if (flags.count("real") == 0) {
            fprintf(stderr, "MTX file: 'real' matrix type required\n");
            return false;
        }
		if (flags.count("general") == 0) {
            fprintf(stderr, "MTX file: 'general' matrix topology required\n");
            return false;
		}
        if (flags.count("coordinate") == 1) {
			std::string line;
			checked_getline_header(line);
			unsigned long long unr, unc, unz;
			int n = sscanf(line.c_str(), "%llu %llu %llu\n", &unr, &unc, &unz);
			if (n != 3) {
				fprintf(stderr, "MTX file format error (matrix dimensions)\n");
				return false;
			}
			nrows = unr;
			ncols = unc;
			nnz_file = unz;
		} else if (flags.count("array") == 1) {
            std::string line;
            checked_getline_header(line);
            unsigned long long unr, unc;
            int n = sscanf(line.c_str(), "%llu %llu\n", &unr, &unc);
            if (n != 2) {
                fprintf(stderr, "MTX file format error (matrix dimensions)\n");
                return false;
            }
            nrows = unr;
            ncols = unc;
            nnz_file = nrows * ncols;
		} else {
			fprintf(stderr, "MTX file: 'coordinate' or 'array' matrix required\n");
		}
        nnz = nnz_file;
        return true;
    }
    
    virtual bool init_data(void *_p = nullptr) {
        auto *p = (params_t *)_p;
        column = p->column;
		if (flags.count("array")) {
            bool found = false;
            try {
                r.resize(nrows);
                c.resize(nrows);
                v.resize(nrows);
            }
            catch(std::bad_alloc &ex) {
                fprintf(stderr, "Not enough memory when allocating MTX buffers\n");
                return false;
            }
			int counter = 0;
			size_t bsz = 1024;
			char *buf = (char *)calloc(bsz, 1);
			assert(buf != nullptr);
			for (int i = 0; i < nnz; i++) {
				int clmn = i / nrows;
				int linelen;
				double val;
				if (!checked_getline_D(buf, bsz, val, linelen)) {
					fprintf(stderr, "MTX file getline error, counter=%d\n", counter);
					free(buf);
					return false;
				}
				if (feof(fp)) {
					if (counter != nnz-1) {
						fprintf(stderr, "MTX file read error: vector: unexpected EOF: counter=%d\n", counter);
						return false;
					}
					break;
				}
				if (clmn == column) {
                    found = true;
					int n = i % nrows;
					r[n] = n;
					c[n] = clmn;
					v[n] = val;
				}
				counter++;
				if (counter == nnz || clmn > column)
					break;
			}
			free(buf);
			nnz = nrows;
			return found;
		} else {
	        return mtx_matrix_buffered_input::init_data();
		}
    }

    virtual bool copy_to_csr_arrays(int *, int *, double *) override {
        assert(0);
    }

    virtual bool copy_to_array(double *vec) override {
        for (int i = 0; i < nnz; ++i) {
            int row = r[i];
            int col = c[i];
            double val = v[i];
            if (column == col) {
                vec[row] = val;
            }
        }
        return true;
    }
    virtual ~mtx_matrixcol_input() {}
};
