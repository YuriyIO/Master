mpicxx -O3 is_symmetric.cpp -o is_symmetric
./compile.sh
g++ uniteInfo.cpp -O3 -o uniteInfo


mpirun -np 4 ./is_symmetric CSR/test 

mpirun -np 4 ./main CSR/test 4 Partitions 0.05 1 0 Partitions_viz

./uniteInfo Partitions