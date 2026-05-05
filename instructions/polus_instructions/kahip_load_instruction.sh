git clone https://github.com/KaHIP/KaHIP

cd KaHIP/

mkdir build

cd build/

cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/home_edu/edu-cmc-sqi23/edu-cmc-sqi23-25/local_libs/KaHIP

make -j8

make install

cd ..
