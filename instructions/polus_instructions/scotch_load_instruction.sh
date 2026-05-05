git clone https://gitlab.inria.fr/scotch/scotch.git

cd scotch/ &&  mkdir build && cd build/

cmake .. -DCMAKE_INSTALL_PREFIX=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/scotch

make -j8

make install

cd ..
