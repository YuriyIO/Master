# Compiler: C++20-ready GCC or Clang compiler
# Dependencies: CMake, oneAPI TBB, Google Sparsehash (optional), MPI (optional)
# System: Linux (x86, ARM) or macOS (ARM)

# Sparsehash

git clone https://github.com/sparsehash/sparsehash.git

cd sparsehash

./configure --prefix=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/sparsehash

make -j8

make install

# KaMinPar

git clone https://github.com/KaHIP/KaMinPar.git

cd KaMinPar

cmake -B build --preset=distributed \
  -DCMAKE_PREFIX_PATH="/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/sparsehash" \
  -DCMAKE_INSTALL_PREFIX="/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/KaMinPar"


cmake -B build --preset=distributed \
  -DCMAKE_C_FLAGS="-mcpu=native" \
  -DCMAKE_CXX_FLAGS="-mcpu=native" \
  -DKAMINPAR_BUILD_TESTS=OFF \
  -DCMAKE_PREFIX_PATH="/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/sparsehash" \
  -DCMAKE_INSTALL_PREFIX="/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/KaMinPar"


find build -type f -exec sed -i 's/-march=native/-mcpu=native/g' {} +

cmake --build build --parallel 8

cmake --install build