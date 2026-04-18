git clone https://github.com/trilinos/Trilinos.git

cd Trilinos

mkdir build

cd build/


  # # --- ВКЛЮЧЕНИЕ ЗОЛТАНА И PHG (v1) --- \
  # -DTrilinos_ENABLE_Zoltan=ON \
  # -DZoltan_ENABLE_PARMETIS=OFF \
  # -DZoltan_ENABLE_SCOTCH=OFF \
  # \
  # # --- ВКЛЮЧЕНИЕ ZOLTAN2 И SPHYNX --- \
  # -DTrilinos_ENABLE_Zoltan2=ON \
  # -DTrilinos_ENABLE_Zoltan2Sphynx=ON \
  # \
  # # --- ЛИНЕЙНАЯ АЛГЕБРА И РЕШАТЕЛИ (Зависимости для Sphynx) --- \
  # -DTrilinos_ENABLE_Tpetra=ON \
  # -DTrilinos_ENABLE_Anasazi=ON \
  # -DTrilinos_ENABLE_Belos=ON \
  # -DTrilinos_ENABLE_Ifpack2=ON \
  # -DTrilinos_ENABLE_MueLu=ON \
  # -DTrilinos_ENABLE_Amesos2=ON \
  # \
  # # --- УПРАВЛЕНИЕ ПАМЯТЬЮ И ПАРАЛЛЕЛИЗМОМ (Kokkos) --- \
  # -DTrilinos_ENABLE_Kokkos=ON \
  # -DKokkos_ENABLE_CUDA=OFF \
  # -DKokkos_ENABLE_OPENMP=ON \
  # -DKokkos_ENABLE_SERIAL=ON \
  # \
  # # --- ТЕСТЫ И ПРИМЕРЫ (Чтобы вы могли запустить Sphynx.exe) --- \
  # -DZoltan2_ENABLE_TESTS=ON \
  # -DZoltan2_ENABLE_EXAMPLES=ON \

cmake \
  -DCMAKE_INSTALL_PREFIX=/mnt/c/Programs/C/VS_Code/local_libs/Trilinos_install \
  -DCMAKE_BUILD_TYPE=RELEASE \
  -DTrilinos_ENABLE_Fortran=OFF \
  -DTPL_ENABLE_MPI=ON \
  -DTrilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF \
  -DTrilinos_SET_GROUP_AND_PERMISSIONS_ON_INSTALL_BASE_DIR:STRING="" \
  -DTrilinos_ENABLE_Zoltan=ON \
  -DZoltan_ENABLE_PARMETIS=OFF \
  -DZoltan_ENABLE_SCOTCH=OFF \
  -DTrilinos_ENABLE_Zoltan2=ON \
  -DTrilinos_ENABLE_Zoltan2Sphynx=ON \
  -DTrilinos_ENABLE_Tpetra=ON \
  -DTrilinos_ENABLE_Anasazi=ON \
  -DTrilinos_ENABLE_Belos=ON \
  -DTrilinos_ENABLE_Ifpack2=ON \
  -DTrilinos_ENABLE_MueLu=ON \
  -DTrilinos_ENABLE_Amesos2=ON \
  -DTrilinos_ENABLE_Kokkos=ON \
  -DKokkos_ENABLE_CUDA=OFF \
  -DKokkos_ENABLE_OPENMP=ON \
  -DKokkos_ENABLE_SERIAL=ON \
  -DZoltan2_ENABLE_TESTS=ON \
  -DZoltan2_ENABLE_EXAMPLES=ON \
  ..

make -j4

make install

# тест:
# cd /mnt/c/Programs/C/VS_Code/local_libs/Trilinos/build/packages/zoltan2/test/sphynx/
# mpirun -np 2 ./Zoltan2_Sphynx.exe
# ожидаем
# UserInputForTests, Matrix type : Laplace3D
# Matrix is distributed.
# UserInputForTests, Create matrix with Laplace3D (and 10 x 10 x 10 mesh)
# UserInputForTests, Implied matrix row coordinates computed
# NumRows     = 1000
# NumNonzeros = 6400
# NumProcs = 2
# NumLocalRows (rank 0) = 500
# Calling solve() 
# precType_ is: jacobi
# Solver type: LOBPCG
# Solver type: LOBPCG
# Done solve() 
# Minimum count:  500 on rank 0
# Maximum count:  500 on rank 0
# Average count:  500
# Total count:    1000 
# Imbalance:     1
# Redistributing matrix...
# Redistributing vectors...
# IntegerVectorTest:  500 == 500 ? OK
# Matvec original...
# Norm of Original matvec prod:       118.217
# Matvec redistributed...
# Norm of Redistributed matvec prod:  118.217
# PASS

# можно добавить в случае ошибки export LD_LIBRARY_PATH=/mnt/c/Programs/C/VS_Code/local_libs/Trilinos_install/lib:$LD_LIBRARY_PATH