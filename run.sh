rm -rf out
mkdir out
mpirun -np 16 ./src/Kolmogorov41.out
