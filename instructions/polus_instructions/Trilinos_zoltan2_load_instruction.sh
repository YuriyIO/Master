git clone https://github.com/trilinos/Trilinos.git

cd Trilinos

mkdir build

cd build/

cmake \
  -DCMAKE_INSTALL_PREFIX=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/Trilinos_install \
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

make -j8

make install

cd ..

cd ..

cd Trilinos_install

ln -s lib64 lib

cd ..
