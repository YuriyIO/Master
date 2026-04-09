git clone https://github.com/KaHIP/KaHIP

cd KaHIP/

mkdir build

cd build/

# флаги сборки:
# -DCMAKE_BUILD_TYPE=Release
# -DCMAKE_BUILD_TYPE=Debug
# Release → максимальная производительность
# Debug → для отладки (медленнее)

# -DUSE_TCMALLOC=On
# Оптимизации памяти - Используй, если установлен tcmalloc (часть Google tools)

# -D64BITMODE=On
# Работа с большими графами - очень большие графы (много рёбер > 2³¹)

# -DDETERMINISTIC_PARHIP=On
# одинаковый результат при одинаковом seed
# хуже качество разбиения

# можно добавить
# -DCMAKE_INSTALL_PREFIX=/mnt/c/Programs/C/VS_Code/local_libs/KaHIP


cmake ../ -DCMAKE_BUILD_TYPE=Release

make -j4

cd ..