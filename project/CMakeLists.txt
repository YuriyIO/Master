cmake_minimum_required(VERSION 3.16)
project(MyProgram LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Находим MPI
find_package(MPI REQUIRED)

# Находим KaMinPar
set(KaMinPar_DIR "/mnt/c/Programs/C/VS_Code/local_libs/KaMinPar/lib/cmake/KaMinPar")
find_package(KaMinPar REQUIRED)

# Укажите ваш исходный файл здесь
set(SOURCE_FILE "main.cpp")  # <- Поменяйте на ваше имя файла
set(OUTPUT_NAME "main")   # <- Как назвать executable

# Пути к библиотекам
set(PREFIX "/mnt/c/Programs/C/VS_Code/local_libs")

set(PARMETIS_DIR "${PREFIX}/ParMETIS")
set(METIS_DIR "${PREFIX}/METIS")
set(GKLIB_DIR "${PREFIX}/GKlib")
set(KAHIP_DIR "${PREFIX}/KaHIP")
set(KAHIP_BUILD "${KAHIP_DIR}/build/parallel/parallel_src")
set(SCOTCH_DIR "${PREFIX}/scotch")

# --- 3. Поиск Trilinos ---
set(Trilinos_PREFIX "${PREFIX}/Trilinos_install")
list(APPEND CMAKE_PREFIX_PATH ${Trilinos_PREFIX})
set(TRILINOS_SOURCE_DIR "${PREFIX}/Trilinos") # Путь для Anasazi_src

# 3. Поиск Trilinos
# Указываем необходимые компоненты. Если они найдены, все пути включатся сами. Zoltan - для phg
find_package(Trilinos REQUIRED COMPONENTS 
    Zoltan2 
    Tpetra 
    Teuchos 
    Anasazi 
    Belos 
    Ifpack2 
    MueLu 
    Amesos2
)

# Создаем executable
add_executable(${OUTPUT_NAME} ${SOURCE_FILE})

# Инклюды для трилиноса
target_include_directories(${OUTPUT_NAME} PRIVATE 
    ${Trilinos_INCLUDE_DIRS} 
    ${Trilinos_TPL_INCLUDE_DIRS}
    "${TRILINOS_SOURCE_DIR}/packages/anasazi/src"
)

# Include директории (полностью соответствует вашему bash-скрипту)
include_directories(
    ${PARMETIS_DIR}/include
    ${METIS_DIR}/include
    ${KAHIP_DIR}/parallel/parallel_src/interface
    ${SCOTCH_DIR}/include
)

target_link_libraries(${OUTPUT_NAME}
    ${Trilinos_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES}
    ${PARMETIS_DIR}/lib/libparmetis.a
    ${METIS_DIR}/lib/libmetis.a
    ${GKLIB_DIR}/lib/libGKlib.a
    # Добавляем пути поиска для следующих библиотек
    -L${KAHIP_BUILD}
    -L${SCOTCH_DIR}/lib
    # Добавляем rpath (как в bash-скрипте)
    -Wl,-rpath,${KAHIP_BUILD}
    # Сами библиотеки
    -lparhip_interface
    -lptscotch
    -lscotch
    -lscotcherr
)

target_link_libraries(${OUTPUT_NAME} KaMinPar::dKaMinPar)

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


