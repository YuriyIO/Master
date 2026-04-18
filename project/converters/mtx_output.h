#pragma once

struct mtx_vector_output : public output {
    int nrows;
    std::vector<double> vec;
    struct params_t {
        int nrows;
    };
    virtual bool init_data(void *p = nullptr) override {
        params_t *_p = (params_t *)p;
        nrows = _p->nrows;
        vec.resize(nrows, 0);
        return true;
    }

    virtual void *get_ptr(int mode = 0) override {
        (void)mode;
        return &vec[0];
    }

    virtual bool commit() override {
        
        fprintf(fp, "%%%%MatrixMarket vector array real general\n");
        fprintf(fp, "%d\n", nrows);
        int row = 1;
        for (auto val : vec) {
            fprintf(fp, "%d %g\n", row, val);
            row++;
        }
        return true;
    }
};
