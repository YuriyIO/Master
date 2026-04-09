#!/bin/bash

echo "========================================="
echo "Source file: main.cpp"
echo "Output:      main"
echo "========================================="

# Создаем временную директорию
BUILD_DIR="build_temp"
mkdir -p "$BUILD_DIR"

cd "$BUILD_DIR"

# Запускаем CMake и компилируем
cmake .. -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_CXX_FLAGS="-O3"

make

# Проверяем результат
if [ -f "main" ]; then
    cp "main" ..
    cd ..
    echo "✓ Success! Executable: main"
    echo "Run with: mpirun -np N ./main CSR/test 4 Partitions"
    rm -rf "$BUILD_DIR"
else
    cd ..
    echo "✗ Compilation failed!"
    rm -rf "$BUILD_DIR"
    exit 1
fi