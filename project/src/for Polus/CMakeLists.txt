cmake_minimum_required(VERSION 3.16)
project(MyProgram LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Находим MPI
find_package(MPI REQUIRED)

# --- Настройка путей ---
set(KAMINPAR_PREFIX "/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/KaMinPar")

# Указываем CMake, где искать конфигурационные файлы TBB и KaMinPar
# Используем ваш симлинк 'lib' (или lib64, если симлинк не сработает)
set(TBB_DIR "${KAMINPAR_PREFIX}/lib/cmake/TBB")
set(KaMinPar_DIR "${KAMINPAR_PREFIX}/lib/cmake/KaMinPar")

# Также добавим корень установки в список поиска префиксов
list(APPEND CMAKE_PREFIX_PATH "${KAMINPAR_PREFIX}")
list(APPEND CMAKE_PREFIX_PATH "/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/sparsehash")

# Теперь находим пакеты
find_package(TBB REQUIRED)
find_package(KaMinPar REQUIRED)

# --- СПИСОК ИСХОДНЫХ ФАЙЛОВ ---
set(SOURCE_FILES 
    main.cpp 
    CSR.cpp 
    Partition.cpp 
    PartitionMetrics.cpp
)
set(OUTPUT_NAME "main")

# Пути к библиотекам
set(PREFIX "/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/")

set(PARMETIS_DIR "${PREFIX}/ParMETIS")
set(METIS_DIR "${PREFIX}/METIS")
set(GKLIB_DIR "${PREFIX}/GKlib")
set(KAHIP_DIR "${PREFIX}/KaHIP")
set(SCOTCH_DIR "${PREFIX}/scotch")

# Поиск Trilinos
set(Trilinos_PREFIX "${PREFIX}/Trilinos_install")
list(APPEND CMAKE_PREFIX_PATH ${Trilinos_PREFIX})
set(TRILINOS_SOURCE_DIR "${PREFIX}/Trilinos") 

find_package(Trilinos REQUIRED COMPONENTS 
    Zoltan Zoltan2 Tpetra Teuchos Anasazi Belos Ifpack2 MueLu Amesos2
)

# Создаем исполняемый файл
add_executable(${OUTPUT_NAME} ${SOURCE_FILES})

# Инклюды
target_include_directories(${OUTPUT_NAME} PRIVATE 
    ${Trilinos_INCLUDE_DIRS} 
    ${Trilinos_TPL_INCLUDE_DIRS}
    "${TRILINOS_SOURCE_DIR}/packages/anasazi/src"
    ${PARMETIS_DIR}/include
    ${METIS_DIR}/include
    ${KAHIP_DIR}/include
    ${SCOTCH_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR} # Чтобы файлы видели друг друга
)

# Линковка
target_link_libraries(${OUTPUT_NAME}
    ${Trilinos_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES}
    ${PARMETIS_DIR}/lib/libparmetis.a
    ${METIS_DIR}/lib/libmetis.a
    ${GKLIB_DIR}/lib/libGKlib.a
    KaMinPar::dKaMinPar
    -L${KAHIP_DIR}/lib 
    -L${SCOTCH_DIR}/lib
    -Wl,-rpath,${KAHIP_DIR}/lib
    -lparhip_interface
    -lptscotch
    -lscotch
    -lscotcherr
)


#Сборка библиотек:

# ParMETIS
# #include <parmetis.h>

#   mpicxx -O3 ver1.cpp -o ver1 \
#   -I/mnt/c/Programs/C/VS_Code/local_libs/METIS/include \
#   -I/mnt/c/Programs/C/VS_Code/local_libs/ParMETIS/include \
#   /mnt/c/Programs/C/VS_Code/local_libs/ParMETIS/lib/libparmetis.a \
#   /mnt/c/Programs/C/VS_Code/local_libs/METIS/lib/libmetis.a \
#   /mnt/c/Programs/C/VS_Code/local_libs/GKlib/lib/libGKlib.a

# scotch
# #include <ptscotch.h>

# mpicxx -O3 ver1.cpp -o ver1 \
#     -I/mnt/c/Programs/C/VS_Code/local_libs/scotch/include \
#     -L/mnt/c/Programs/C/VS_Code/local_libs/scotch/lib -lptscotch -lscotch -lscotcherr

#mpirun -np 4 ./ver1 CSR/test 4 Partitions

# KaHIP
# #include <parhip_interface.h>

# mpicxx -O3 ver1.cpp -o ver1 \
#     -I/mnt/c/Programs/C/VS_Code/local_libs/KaHIP/parallel/parallel_src/interface \
#     -L/mnt/c/Programs/C/VS_Code/local_libs/KaHIP/build/parallel/parallel_src \
#     -Wl,-rpath,/mnt/c/Programs/C/VS_Code/local_libs/KaHIP/build/parallel/parallel_src \
#     -lparhip_interface

# KaMinPar
# #include <dkaminpar.h>

# target_link_libraries(${OUTPUT_NAME} KaMinPar::dKaMinPar)
