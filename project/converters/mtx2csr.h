#pragma once

namespace helpers {

static inline std::vector<std::string> str_split(const std::string &s, char delimiter) {
    std::vector<std::string> result;
    std::string token;
    std::istringstream token_stream(s);
    while (std::getline(token_stream, token, delimiter)) {
        result.push_back(token);
    }
    return result;
}

static FILE *checked_open_file(const std::string &name, const char *mode, const std::string &descr) {
    char buf[101] = {0,};
    FILE *fp = fopen(name.c_str(), mode);
    if (fp == nullptr) {
        strerror_r(errno, buf, 100);
        fprintf(stderr, "Can't open %s file, mode: %s; error=%s\n", descr.c_str(), mode, buf);
        return nullptr;
    }
    return fp;
}

}


#include "fast_double_parser.h"
namespace parsers {

const char *parse_unsigned_int(const char *p, int &r)
{
    int i = 0;
    bool state = false;
    for (; *p != 0; p++) {
        if (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') {
            if (state) {
                break;
            }
        } else {
            state = true;
            if (*p < '0' || *p > '9')
                return nullptr;
            i = i*10 + (*p - '0');
        }
    }
    r = i;
    while (*p && (*p == ' ' || *p == '\t')) p++;
    return p;
}

}

struct sparse_matrix_descriptor {
    int nrows, ncols, nnz;
    bool with_rhs_vec, with_solution_vec;
};

struct input {
    enum { OMIT_SORTING_STAGE = 1, SORT_BY_COLUMN = 2 }; // FIXME SORT_BY_COLUMN to be implemented
    FILE *fp = nullptr;
    char *iobuf = nullptr;
    int options = 0;
    bool open_file(const std::string &name, const std::string &descr) {
        fp = helpers::checked_open_file(name, "r", descr);
        if (fp == nullptr)
            return false;
#if 0
        size_t bufsize = 64;
        iobuf = (char *)calloc(bufsize, 1);
        assert(iobuf != nullptr);
        setvbuf(fp, iobuf, _IOFBF, bufsize); 
#endif        
        return true;
    }

    template <typename T>
        bool read_arr(T *p, size_t n)
        {
            char buf[101] = {0,};
            size_t count = 0;
            while (count != n) {
                size_t res = fread(p, sizeof(T), n, fp);
                if (res != n) {
                    int nerr = ferror(fp);
                    if (nerr) {
                        strerror_r(errno, buf, 100);
                        printf("read_arr: fread failed: error=%s\n", buf);
                        return false;
                    }
                    if (feof(fp))
                        return false;
                }
                count += res;
            }
            return true;
        }

    int parse_IID(const char *p, int &i, int &j, double &d) {
        int n = 0;
        const char *next = p;
        next = parsers::parse_unsigned_int(next, i);
        if (next == nullptr)
            return n;
        ++n;
        next = parsers::parse_unsigned_int(next, j);
        if (next == nullptr)
            return n;
        ++n;
        next = fast_double_parser::parse_number(next, &d);
        if (next == nullptr)
            return n;
        ++n;
        return n;
    }

    int parse_II(const char *p, int &i, int &j) {
        int n = 0;
        const char *next = p;
        next = parsers::parse_unsigned_int(next, i);
        if (next == nullptr)
            return n;
        ++n;
        next = parsers::parse_unsigned_int(next, j);
        if (next == nullptr)
            return n;
        ++n;
        return n;
    }

    int parse_ID(const char *p, int &i, double &d) {
        int n = 0;
        const char *next = p;
        next = parsers::parse_unsigned_int(next, i);
        if (next == nullptr)
            return n;
        ++n;
        next = fast_double_parser::parse_number(next, &d);
        if (next == nullptr)
            return n;
        ++n;
        return n;
    }

    int parse_D(const char *p, double &d) {
        int n = 0;
        const char *next = p;
        next = fast_double_parser::parse_number(next, &d);
        if (next == nullptr)
            return n;
        ++n;
        return n;
    }

    bool checked_getline_IID(char *&buf, size_t &bsz, int &r, int &c, double &v, int &linelen) {
        ssize_t errcode;
        if ((errcode = getline(&buf, &bsz, fp)) == -1) {
            char msg[101];
            strerror_r(errno, msg, 100);
            fprintf(stderr, "Input file reading error; error=%s\n", msg);
            return false;
        }
        linelen = strlen(buf);
        // NOTE: sscanf is 4-10 times slower
        //int retval = sscanf(buf, "%d %d %lf", &r, &c, &v);
        int retval = parse_IID(buf, r, c, v);
        if (retval != 3) {
            fprintf(stderr, "Input file format error (value dimensions or format), read: %d fileds, expected: int, int, double\n", retval);
            fprintf(stderr, "NOTE: failed on a line: ``%s''\n", buf);
            return false;
        }
        return true;
    }

    bool checked_getline_ID(char *&buf, size_t &bsz, int &n, double &v, int &linelen) {
        ssize_t errcode;
        if ((errcode = getline(&buf, &bsz, fp)) == -1) {
            char msg[101];
            strerror_r(errno, msg, 100);
            fprintf(stderr, "Input file reading error; error=%s\n", msg);
            return false;
        }
        linelen = strlen(buf);
        // NOTE: sscanf is 4-10 times slower
        //int retval = sscanf(buf, "%d %lf", &n, &v);
        int retval = parse_ID(buf, n, v);
        if (retval != 2) {
            fprintf(stderr, "Input file format error (value dimensions or format), read: %d fileds, expected: int, double\n", retval);
            fprintf(stderr, "NOTE: failed on a line: ``%s''\n", buf);
            return false;
        }
        return true;
    }

    bool checked_getline_D(char *&buf, size_t &bsz, double &v, int &linelen) {
        ssize_t errcode;
        if ((errcode = getline(&buf, &bsz, fp)) == -1) {
            char msg[101];
            strerror_r(errno, msg, 100);
            fprintf(stderr, "Input file reading error; error=%s\n", msg);
            return false;
        }
        linelen = strlen(buf);
        // NOTE: sscanf is 4-10 times slower
        //int retval = sscanf(buf, "%d %lf", &n, &v);
        int retval = parse_D(buf, v);
        if (retval != 1) {
            fprintf(stderr, "Input file format error (value dimensions or format), read: %d fileds, expected: int, double\n", retval);
            fprintf(stderr, "NOTE: failed on a line: ``%s''\n", buf);
            return false;
        }
        return true;
    }

    bool checked_getline_II(char *&buf, size_t &bsz, int &r, int &c, int &linelen) {
        ssize_t errcode;
        if ((errcode = getline(&buf, &bsz, fp)) == -1) {
            char msg[101];
            strerror_r(errno, msg, 100);
            fprintf(stderr, "Input file reading error; error=%s\n", msg);
            return false;
        }
        linelen = strlen(buf);
        // NOTE: sscanf is 4-10 times slower
        //int retval = sscanf(buf, "%d %d", &r, &c);
        int retval = parse_II(buf, r, c);
        if (retval != 2) {
            fprintf(stderr, "Input file format error (value dimensions or format), read: %d fileds, expected: int, int\n", retval);
            fprintf(stderr, "NOTE: failed on a line: ``%s''\n", buf);
            return false;
        }
        return true;
    }
    virtual int get_nrows() = 0;
    virtual int get_ncols() = 0;
    virtual int get_nnz() = 0;
    virtual int get_ndiag() = 0;
    virtual bool read_header(void *p = nullptr) = 0;
    virtual bool init_data(void *params = nullptr) = 0;
    virtual bool copy_to_array(double *vec) = 0;
    virtual bool copy_to_csr_arrays(int *rows, int *cols, double *vals) = 0;
    virtual ~input() { if (fp) { fclose(fp); fp = nullptr; if (iobuf) free(iobuf); } }
};

struct output {
    FILE *fp = nullptr;
    virtual bool init_data(void *p) = 0;
    virtual void *get_ptr(int mode = 0) = 0;
    virtual bool commit() = 0;
    bool open_file(const std::string &name, const std::string &descr) {
        fp = helpers::checked_open_file(name, "w+", descr);
        if (fp == nullptr)
            return false;
        return true;
    }
};

struct csr_output : public output {
    enum mode_t { ROWS_PTR = 1, COLS_PTR = 2, VALS_PTR = 3, RHS_PTR = 4, SOL_PTR = 5 };
    struct params_t {
        int nrows, ncols, nnz;
        bool with_rhs_vec, with_solution_vec;
    };
    virtual void *get_ptr(int mode = 0) { (void)mode; return nullptr; }
    virtual bool init_data(void *params = nullptr) = 0;
    virtual int *get_rows_ptr() { return (int *)get_ptr(csr_output::mode_t::ROWS_PTR); }
    virtual int *get_cols_ptr() { return (int *)get_ptr(csr_output::mode_t::COLS_PTR); }
    virtual double *get_vals_ptr() { return (double *)get_ptr(csr_output::mode_t::VALS_PTR); }
    virtual double *get_rhs_ptr() { return (double *)get_ptr(csr_output::mode_t::RHS_PTR); }
    virtual double *get_sol_ptr() { return (double *)get_ptr(csr_output::mode_t::SOL_PTR); }
    virtual bool commit() = 0;
};

struct mtx_vec_output {
    double *vec;
    virtual void *get_ptr(int mode = 0) { (void)mode; return vec; }
    virtual bool init_data(void *params = nullptr);
    virtual bool commit();
};

