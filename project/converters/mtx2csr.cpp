#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <assert.h>

#include <string>
#include <vector>
#include <set>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <memory>

#include "mtx2csr.h"
#include "mtx_input.h"
#include "csr_output.h"

#define FATAL(descr) { fprintf(stderr, "FATAL ERROR on stage: %s; source file: %s:%d\n", descr, __FILE__, __LINE__); return 1; }


void usage(char **argv)
{
    auto myname = basename(argv[0]);
    printf("MTX to CRS converter\n"
            "Usage: %s [-mode=[nosort][,direct|buffered]] <input_file> <output_file> [rhs_vec_file [solution_vec_file]]\n", myname);
    printf("\nMode flags:\n  nosort - implies the input MTX file is already sorted by row\n"
            "  direct - conversion method with lower memory consumption\n"
            "  buffered - conversion method with full input and output buffering in RAM\n");
}

bool parse_mode_string(const std::string &str, bool &n, bool &d, bool &b) {
    auto mode = helpers::str_split(str, '=');
    if (mode.size() != 2 || mode[0] != "-mode") {
        return false;
    }
    auto flags = helpers::str_split(mode[1], ',');
    for (const auto &f : flags) {
        if (f == "nosort") {
            n = true;
        } else if (f == "direct") {
            d = true;
        } else if (f == "buffered") {
            b = true;
        } else {
            return false;
        }
    }
    if (b && d)
        return false;
    return true;
}

struct flags_t {
    bool nosort = false;
    bool direct = false;
    bool buffered = false;
    void defaults() { nosort = false; buffered = true; direct = false; }
};

std::shared_ptr<input> get_input(const flags_t &flags) {
    std::shared_ptr<input> retval;
    if (flags.buffered) {
        retval = std::make_shared<mtx_matrix_buffered_input>();
    } else {
        retval = std::make_shared<mtx_matrix_direct_input>();
    }
    if (flags.nosort)
        retval->options |= input::OMIT_SORTING_STAGE;
    return retval;
}

std::shared_ptr<csr_output> get_output(const flags_t &flags) {
    std::shared_ptr<csr_output> retval;
    if (flags.buffered) {
        retval = std::make_shared<csr_buffered_output>();
    } else {
        retval = std::make_shared<csr_mmaped_output>();
    }
    return retval;
}

int main(int argc, char **argv)
{
    if ((argc < 3) || (argc > 6)) {
        usage(argv);
        return 1;
    }
    flags_t flags;
    int args_start = 1;
    if (argv[1][0] == '-') {
        if (!parse_mode_string(argv[1], flags.nosort, flags.direct, flags.buffered)) {
            usage(argv);
            return 1;
        }
        args_start++;
    } else {
        flags.defaults();
    }

    // Open all required files 
    std::shared_ptr<input> mtx = get_input(flags);
    if (std::string(argv[args_start]) == "-") {
        mtx->fp = stdin;
        mtx->options |= input::OMIT_SORTING_STAGE;
    } else {
        if (!mtx->open_file(argv[args_start], "input matrix")) {
            FATAL("opening the input matrix file");
        }
    }
    std::shared_ptr<csr_output> csr = get_output(flags);
    if (!csr->open_file(argv[args_start + 1], "output matrix")) {
        FATAL("opening the output matrix file for writing");
    }
    mtx_vector_input rhs;
    if (argc - args_start >= 3 && !rhs.open_file(argv[args_start + 2], "RHS vector")) {
        FATAL("opening input RHS vector file");
    }
    mtx_vector_input sol;
    if (argc - args_start >= 4 && !sol.open_file(argv[args_start + 3], "solution vector")) {
        FATAL("opening input solution vector file");
    }

    // Read and check headers of all input files
    if (!mtx->read_header()) {
        FATAL("reading and interpreting the input matrix file header");
    }

    mtx_vector_input::params_t pv { mtx->get_nrows() };
    if (rhs.fp != nullptr && (!rhs.read_header(&pv))) {
        FATAL("reading and interpreting the input vector file header");
    }

    if (sol.fp != nullptr && (!sol.read_header(&pv))) {
        FATAL("reading and interpreting the input vector file header");
    }

    csr_output::params_t pm { mtx->get_nrows(), mtx->get_ncols(), mtx->get_nnz(),
                             rhs.fp != nullptr, sol.fp != nullptr };
    if (!csr->init_data(&pm)) {
        FATAL("making output CSR file and mapping it into memory");
    }

    // Read the whole file to save each line position and then sort this map by row number
    if (!mtx->init_data()) {
        FATAL("creating a file map for (presumably) unsorted MTX");
    }

    // Copy matrix data
    if (!mtx->copy_to_csr_arrays(csr->get_rows_ptr(), 
                csr->get_cols_ptr(), 
                csr->get_vals_ptr())) {
        FATAL("copying matrix data from MTX to csr file");
    }

    // Copy RHS vector data
    if (rhs.fp != nullptr && (!rhs.copy_to_array(csr->get_rhs_ptr()))) {
        FATAL("copying vector data from MTX to csr file");
    }

    // Copy solution vector data
    if (sol.fp != nullptr && (!sol.copy_to_array(csr->get_sol_ptr()))) {
        FATAL("copying vector data from MTX to csr file");
    }
    csr->commit();
    return 0;
}
