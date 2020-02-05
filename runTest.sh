cd test/test_scalar_2D
mpirun -np 1 ../../src/Kolmogorov41.out
cd ..
cd test_scalar_3D
mpirun -np 1 ../../src/Kolmogorov41.out
cd ..
cd test_velocity_2D
mpirun -np 1 ../../src/Kolmogorov41.out
cd ..
cd test_velocity_2D
mpirun -np 1 ../../src/Kolmogorov41.out
cd ..
python test.py


