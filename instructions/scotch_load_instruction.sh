git clone https://gitlab.inria.fr/scotch/scotch.git

cd scotch/ &&  mkdir build && cd build/

cmake .. -DCMAKE_INSTALL_PREFIX=/mnt/c/Programs/C/VS_Code/local_libs/scotch/ 

make -j4

make install
