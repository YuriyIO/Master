# GKlib

git clone https://github.com/KarypisLab/GKlib.git

cd GKlib/

make config prefix=/mnt/c/Programs/C/VS_Code/local_libs/GKlib

make install 

# METIS

git clone https://github.com/KarypisLab/METIS.git

cd METIS/

make config cc=gcc prefix=/mnt/c/Programs/C/VS_Code/local_libs/METIS gklib_path=/mnt/c/Programs/C/VS_Code/local_libs/GKlib

make install


# ParMETIS

git clone https://github.com/KarypisLab/ParMETIS.git

cd ParMETIS/

mkdir build && cd build

cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DGKLIB_PATH=/mnt/c/Programs/C/VS_Code/local_libs/GKlib \
  -DMETIS_PATH=/mnt/c/Programs/C/VS_Code/local_libs/METIS \
  -DCMAKE_INSTALL_PREFIX=/mnt/c/Programs/C/VS_Code/local_libs/ParMETIS \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5

make -j4

make install
