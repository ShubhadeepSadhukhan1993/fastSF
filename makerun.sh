cd src
make
cd ..
rm -rf out
mkdir out
mpirun -np 1 ./src/Kolmogorov41.out
