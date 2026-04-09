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

# Создаем executable
add_executable(${OUTPUT_NAME} ${SOURCE_FILE})

# Include директории (полностью соответствует вашему bash-скрипту)
include_directories(
    ${PARMETIS_DIR}/include
    ${METIS_DIR}/include
    ${GKLIB_DIR}/include
    ${KAHIP_DIR}/parallel/parallel_src/interface
    ${SCOTCH_DIR}/include
)

target_link_libraries(${OUTPUT_NAME}
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
