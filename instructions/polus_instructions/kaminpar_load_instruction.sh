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

SPARSE_PATH="/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/sparsehash/include"


# Обновляем обертку для C
echo '#!/bin/bash
ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "-march=native" ]]; then
        ARGS+=("-mcpu=native")
    elif [[ "$arg" == "-fcf-protection=full"  "$arg" == "-fcf-protection" ]]; then
        continue
    else
        ARGS+=("$arg")
    fi
done
exec /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/gcc-13.4.0/bin/gcc "-I'$SPARSE_PATH'" "${ARGS[@]}"' > /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/kaminpar_fix/kaminpar_cc

# Обновляем обертку для C++
echo '#!/bin/bash
ARGS=()
for arg in "$@"; do
    if [[ "$arg" == "-march=native" ]]; then
        ARGS+=("-mcpu=native")
    elif [[ "$arg" == "-fcf-protection=full"  "$arg" == "-fcf-protection" ]]; then
        continue
    else
        ARGS+=("$arg")
    fi
done
exec /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/gcc-13.4.0/bin/g++ "-I'$SPARSE_PATH'" "${ARGS[@]}"' > /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/kaminpar_fix/kaminpar_cxx

chmod +x /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/kaminpar_fix/kaminpar_cc
chmod +x /home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/kaminpar_fix/kaminpar_cxx


git clone https://github.com/KaHIP/KaMinPar.git

cd KaMinPar

sed -i 's/#error "Only x64 and ARM are supported"/\/\/#error "PowerPC fallback"/' kaminpar-common/graph_compression/streamvbyte.h

cat << 'EOF' > kaminpar-common/graph_compression/streamvbyte.h
#pragma once
#include <cstdint>
#include <vector>

namespace kaminpar::streamvbyte {
    enum class DifferentialCodingKind { NONE, D1, D2, D4, DM };

    // auto... Args позволяет принимать и типы (NodeID), и значения (false, D1)
    template <typename T, auto... Args>
    struct UniversalStub {
        UniversalStub(auto...) {}
        template <typename L> void decode(L&&) {}
        std::size_t add(auto...) { return 0; }
        void flush() {}
        void finalize() {}
    };

    template <typename T, auto... Args> 
    using StreamVByteDecoder = UniversalStub<T, Args...>;
    
    template <typename T, auto... Args> 
    using StreamVByteEncoder = UniversalStub<T, Args...>;
}

namespace kaminpar {
    template <typename T, auto... Args> 
    using StreamVByteGapAndWeightsDecoder = streamvbyte::UniversalStub<T, Args...>;
    
    template <typename T, auto... Args> 
    using StreamVByteGapDecoder = streamvbyte::UniversalStub<T, Args...>;
    
    template <typename T, auto... Args> 
    using StreamVByteGapAndWeightEncoder = streamvbyte::UniversalStub<T, Args...>;
    
    template <typename T, auto... Args> 
    using StreamVByteGapEncoder = streamvbyte::UniversalStub<T, Args...>;
}
EOF


cmake -B build --preset=distributed \
  -DCMAKE_CXX_COMPILER=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/kaminpar_fix/kaminpar_cxx \
  -DCMAKE_C_COMPILER=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/kaminpar_fix/kaminpar_cc \
  -DCMAKE_PREFIX_PATH="/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/sparsehash" \
  -DCMAKE_INSTALL_PREFIX="/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/KaMinPar"

cmake --build build --parallel 8

cmake --install build

ln -s lib64 lib