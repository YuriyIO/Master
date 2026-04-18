#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <string>
#include <stdlib.h>
#include <sys/mman.h>
#include <assert.h>
#include <iostream>


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
    std::cout << "Usage: " << bname << " [-d significant digits nimber] input_file.csr output_file.mtx" << std::endl;
    free(name);
}


int main(int argc, char **argv)
{
    opterr = 0;
    int c;
    int significant_digits = 20;
    char *file_to_read = NULL;
    char *file_to_write = NULL;
    while ((c = getopt(argc, argv, "d:")) != -1)
        switch (c)
        {
            case 'd':
                significant_digits = atoi(optarg);
                if (significant_digits < 1 || significant_digits > 30) {
                    usage(argv[0]);
                    return EXIT_FAILURE;
                }
                break;
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
    FILE *fp = fopen(filer.c_str(), "r");
    if (fp == NULL) {
        printf("Can't open input file for reading\n");
        return 1;
    }
    std::string filew(file_to_write);
    FILE *fw = fopen(filew.c_str(), "w");
    if (fw == NULL) {
        printf("Can't open output file for writing\n");
        return 1;
    }

    unsigned int nrows = 0;
    if (!read_arr<unsigned int>(&nrows, 1, fp)) {
        printf("Input file format error\n");
        return 1;
    }
    unsigned int nonzeros = 0;
    if (!read_arr<unsigned int>(&nonzeros, 1, fp)) {
        printf("Input file format error\n");
        return 1;
    }
    size_t len = (nrows + 1) * sizeof(int) + nonzeros * sizeof(int) + nonzeros * sizeof(double);
    off_t off = ftell(fp);
    fseek(fp, 0, SEEK_END);
    off_t end = ftell(fp);
    assert(end == (off_t)(len + off));
    fseek(fp, 0, SEEK_SET);

    int fd = fileno(fp);
    void *ptr = mmap(NULL, end, PROT_READ, MAP_PRIVATE, fd, 0);
    if (ptr == MAP_FAILED) {
        perror("mmap");
        return 1;
    }
    int *rows = (int *)((char *)ptr + off);
    int *cols = (int *)(rows + (nrows + 1));
    double *vals = (double *)(cols + nonzeros);
    assert(vals + nonzeros - (double *)ptr == (end)/8);

    fprintf(fw, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fw, " %d %d %d\n", nrows, nrows, nonzeros);
    char format[256];
    sprintf(format, " %%d %%d %%.%de\n", significant_digits);
    for (unsigned int row = 0; row < nrows; row++) {
        int ci_start = rows[row];
        int ci_end = rows[row + 1];
        for (int ci = ci_start; ci < ci_end; ci++) {
            int col = cols[ci];
            double val = vals[ci];

            fprintf(fw, (const char *)format, row+1, col+1, val); 
        }
    }
    fclose(fp);
    fclose(fw);
    munmap(ptr, end);
    return 0;
}
