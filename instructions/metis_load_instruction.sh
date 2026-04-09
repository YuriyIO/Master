git clone https://github.com/KarypisLab/METIS.git

cd METIS/

# Common configuration options are:
# cc=[compiler]     - The C compiler to use [default is determined by CMake]
# shared=1          - Build a shared library instead of a static one [off by default]
# prefix=[PATH]     - Set the installation prefix [~/local by default]
# gklib_path=[PATH] - Set the prefix path where GKlib has been installed. You can skip
#                     this if GKlib's installation prefix is the same as that of METIS.
# i64=1             - Sets to 64 bits the width of the datatype that will store information
#                     about the vertices and their adjacency lists. 
# r64=1             - Sets to 64 bits the width of the datatype that will store information 
#                     about floating point numbers.

make config cc=gcc prefix=/mnt/c/Programs/C/VS_Code/local_libs/METIS gklib_path=/mnt/c/Programs/C/VS_Code/local_libs/GKlib

make install
