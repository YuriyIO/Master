# GKlib

git clone https://github.com/KarypisLab/GKlib.git

cd GKlib/

make config prefix=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/GKlib

make -j8

make install 

ln -s lib64 lib

cd ..

# METIS

git clone https://github.com/KarypisLab/METIS.git

cd METIS/

make config cc=gcc prefix=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/METIS gklib_path=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/GKlib

make -j8

make install

cd ..

# ParMETIS

git clone https://github.com/KarypisLab/ParMETIS.git

cd ParMETIS/

mkdir build && cd build

cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DGKLIB_PATH=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/GKlib \
  -DMETIS_PATH=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/METIS \
  -DCMAKE_INSTALL_PREFIX=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/ParMETIS \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5

make -j8

make install

cd ..
