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
#include "mtx_output.h"

#define FATAL(descr) { fprintf(stderr, "FATAL ERROR on stage: %s; source file: %s:%d\n", descr, __FILE__, __LINE__); return 1; }


void usage(char **argv)
{
    auto myname = basename(argv[0]);
    printf("MTX matrix column to MTX vector converter\n"
            "Usage: %s <input_file> <output_file> [column_number]\n", myname);
}

std::shared_ptr<input> get_input() {
    std::shared_ptr<input> retval;
    retval = std::make_shared<mtx_matrixcol_input>();
    return retval;
}

std::shared_ptr<output> get_output() {
    std::shared_ptr<output> retval;
    retval = std::make_shared<mtx_vector_output>();
    return retval;
}

int main(int argc, char **argv)
{
    if ((argc < 3) || (argc > 4)) {
        usage(argv);
        return 1;
    }
    int args_start = 1;

    // Open all required files 
    std::shared_ptr<input> mtx = get_input();
    if (!mtx->open_file(argv[args_start], "input matrix")) {
        FATAL("opening the input matrix file");
    }
    std::shared_ptr<output> vec = get_output();
    if (!vec->open_file(argv[args_start + 1], "output vector")) {
        FATAL("opening the output matrix file for writing");
    }

    int column = 0;
    if (argc == 4) {
        column = std::stoi(argv[args_start + 2]);
    }

    // Read and check headers of all input files
    if (!mtx->read_header()) {
        FATAL("reading and interpreting the input matrix file header");
    }

    mtx_vector_output::params_t pv { mtx->get_nrows() };
    if (!vec->init_data(&pv)) {
        FATAL("making output CSR file and mapping it into memory");
    }

    // Read the whole file to save each line position and then sort this map by row number
    mtx_matrixcol_input::params_t pm { column };
    if (!mtx->init_data(&pm)) {
        FATAL("init data stage of MTX file processing");
    }

    // Copy matrix data
    if (!mtx->copy_to_array((double *)vec->get_ptr())) {
        FATAL("copying matrix data from MTX matrix to MTX vector file");
    }
    vec->commit();
    return 0;
}
