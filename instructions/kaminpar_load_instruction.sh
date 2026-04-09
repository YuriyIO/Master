# Compiler: C++20-ready GCC or Clang compiler
# Dependencies: CMake, oneAPI TBB, Google Sparsehash (optional), MPI (optional)
# System: Linux (x86, ARM) or macOS (ARM)

# Sparsehash

git clone https://github.com/sparsehash/sparsehash.git

cd sparsehash

./configure --prefix=/mnt/c/Programs/C/VS_Code/local_libs/sparsehash

make -j4

make install

# KaMinPar

git clone https://github.com/KaHIP/KaMinPar.git

cd KaMinPar

cmake -B build --preset=distributed \
  -DCMAKE_PREFIX_PATH="/mnt/c/Programs/C/VS_Code/local_libs/sparsehash" \
  -DCMAKE_INSTALL_PREFIX="/mnt/c/Programs/C/VS_Code/local_libs/KaMinPar"

cmake --build build --parallel 4

cmake --install build