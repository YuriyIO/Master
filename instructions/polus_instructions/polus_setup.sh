scl enable devtoolset-9 bash
module load OpenMPI
export PATH=$HOME/gcc-13.4.0/bin:$PATH
export LD_LIBRARY_PATH=$HOME/gcc-13.4.0/lib64:$LD_LIBRARY_PATH