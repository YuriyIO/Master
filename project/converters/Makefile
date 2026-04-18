CXXFLAGS += -O3 -Wall -Wextra -std=c++11
#CXXFLAGS += -g -O0 -Wall -Wextra -std=c++11

all: mtx2csr csr2mtx csr2amgx mtxcol2mtxvec

csr2amgx: csr2amgx.cpp
	g++ $(CXXFLAGS) -o $@ $<

mtx2csr: mtx2csr.cpp
	g++ $(CXXFLAGS) -o $@ $<

csr2mtx: csr2mtx.cpp
	g++ $(CXXFLAGS) -o $@ $<

mtxcol2mtxvec: mtxcol2mtxvec.cpp
	g++ $(CXXFLAGS) -o $@ $<

clean:
	rm -f mtx2csr csr2mtx csr2amgx mtxcol2mtxvec

HEADERS=csr_output.h fast_double_parser.h mtx2csr.h mtx_input.h mtx_output.h

mtx2csr csr2mtx csr2amgx mtxcol2mtxvec: $(HEADERS)
