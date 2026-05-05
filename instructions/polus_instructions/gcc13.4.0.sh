cd my_gcc

wget https://ftp.gnu.org/gnu/gcc/gcc-13.4.0/gcc-13.4.0.tar.gz
tar -xvf gcc-13.4.0.tar.gz
cd gcc-13.4.0

./contrib/download_prerequisites

cd ..
mkdir gcc-build
cd gcc-build


../gcc-13.4.0/configure \
  --prefix=$HOME/gcc-13.4.0 \
  --disable-multilib \
  --enable-languages=c,c++

make -j8

make install

export PATH=$HOME/gcc-13.4.0/bin:$PATH

export LD_LIBRARY_PATH=$HOME/gcc-13.4.0/lib64:$LD_LIBRARY_PATH

gcc --version