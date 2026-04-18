#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sys/mman.h>
#include <assert.h>
#include <iostream>




typedef unsigned int uint32_t;


    template <typename T>
bool read_arr(T *p, size_t n, FILE *fp)
{
    size_t res = fread(p, sizeof(T), n, fp);
    if (res != n) {
        int nerr = ferror(fp);
        printf("fread failed: %d errcode\n", nerr);
        return false;
    }
    return true;
}


void usage(const char *exec_name)
{
    char *name = strdup(exec_name);
    char *bname;
    bname = basename(name);
    std::cout << "Usage: " << bname << " [-d significant digits nimber] input_file.csr output_file.amgx" << std::endl;
    free(name);
}


int main(int argc, char **argv)
{
    opterr = 0;
    int c;
    char *file_to_read = NULL;
    char *file_to_write = NULL;
    while ((c = getopt(argc, argv, "")) != -1)
        switch (c)
        {
            case '?':
                if (optopt == 'd')
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else if (isprint(optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
                usage(argv[0]);
                return EXIT_FAILURE;
            default:
                usage(argv[0]);
                return EXIT_FAILURE;
        }
    for (int index = optind; index < argc; index++) {
        if (index - optind == 0) {
            file_to_read = argv[index];
        }
        if (index - optind == 1) {
            file_to_write = argv[index];
        }
        if (index - optind > 1) {
            usage(argv[0]);
            return EXIT_FAILURE;
        }
    }
    if (file_to_read == NULL || file_to_write == NULL) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    std::string filer(file_to_read);
    FILE *fr = fopen(filer.c_str(), "r");
    if (fr == NULL) {
        printf("Can't open input file for reading\n");
        return 1;
    }
    std::string filew(file_to_write);
    FILE *fw = fopen(filew.c_str(), "wb");
    if (fw == NULL) {
        printf("Can't open output file for writing\n");
        return 1;
    }

    unsigned int nrows = 0;
    if (!read_arr<unsigned int>(&nrows, 1, fr)) {
        printf("Input file format error\n");
        return 1;
    }
    unsigned int nonzeros = 0;
    if (!read_arr<unsigned int>(&nonzeros, 1, fr)) {
        printf("Input file format error\n");
        return 1;
    }
    size_t len = (nrows + 1) * sizeof(int) + nonzeros * sizeof(int) + nonzeros * sizeof(double);
    off_t off = ftell(fr);
    fseek(fr, 0, SEEK_END);
    off_t end = ftell(fr);
    assert(end == (off_t)(len + off));
    fseek(fr, 0, SEEK_SET);

    void *ptr = mmap(NULL, end, PROT_READ, MAP_PRIVATE, fileno(fr), 0);
    if (ptr == MAP_FAILED) {
        perror("mmap");
        return 1;
    }
    int *rows = (int *)((char *)ptr + off);
    int *cols = (int *)(rows + (nrows + 1));
    double *vals = (double *)(cols + nonzeros);
    assert(vals + nonzeros - (double *)ptr == (end)/8);

    const char header [] = "%%NVAMGBinary\n";
    const int system_header_size = 9;
    uint32_t matrix_format = 0;
    uint32_t system_flags [] = { (uint32_t)(1), (uint32_t)(0), (uint32_t)(0), matrix_format, (uint32_t)(0), 
        (uint32_t)(1), (uint32_t)(1), (uint32_t)(nrows), (uint32_t)(nonzeros) };    


    fwrite(header, sizeof(char), strlen(header), fw);
    fwrite(system_flags, sizeof(uint32_t), system_header_size, fw);
    fwrite(rows, sizeof(int), nrows + 1, fw);
    fwrite(cols, sizeof(int), nonzeros, fw);
    fwrite(vals, sizeof(double), nonzeros, fw);

    fclose(fr);
    fclose(fw);
    munmap(ptr, end);
    return 0;
}
