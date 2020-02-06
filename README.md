# `Kolmogorov41`

`Kolmogorov41` is an open source hybrid parallel C++ code to compute structure functions for a given velocity or scalar field.

## Getting the Source Code

`Kolmogorov41` is hosted on GitHub. You can download the source code from the following link:

https://github.com/ShubhadeepSadhukhan1993/Kolmogorov41

## Installing `Kolmogorov41`

### Required Libraries

The following libraries are required for installing and running Kolmogorov41:

1. `CMake`
2. `Blitz++`
3. `YAML-cpp`
4. `MPICH`
5. `HDF5`
6. `H5SI`

The instructions to download and install these libraries are provided in the following website:

(http://turbulencehub.org/index.php/codes/tarang/installing-tarang/).

###  Compiling instruction

After downloading `Kolmogorov41`, change into `Kolmogorov41-master/src` directory and run the command `make` in the terminal. An executable named `Kolmogorov41.out` will be created inside the `Kolmogorov41-master/src` folder.

## Testing `Kolmogorov41`
`Kolmogorov41` offers an automated testing process to validate the code. The relevant test scripts can be found in the `tests/` folder of the code. To execute the tesing process, change into `\Kolmogorov41-master` and run the command 

`bash runTest.sh`. 

The code then runs two test cases; these are as follows. 

* In one case, the code will generate a 2D velocity field given by **u** = [*x, y, z*], and compute the structure functions for the given field. For this case, the longitudinal structure functions should equal *l<sup>q</sup>*. 

* In the second case, the code will generate a 2D scalar field given by *T = x + y + z*, and compute the structure functions for the given field. Fir this case, the structure functions should equal *(l<sub>x</sub> + l<sub>z</sub>)<sup>q</sup>*.

For both the cases, `Kolmogorov41` will compare the computed structure functions with the analytical results. If the percentage difference between the two values is less than 10<sup>-10</sup>, the code is deemed to have passed. 

Finally, for visualization purpose, the python script `test/test.py` is invoked. This script generates the plots of the second and third-order longitudinal structure functions versus *l*, and the density plots of the computed second-order scalar structure functions and *(l<sub>x</sub> + l<sub>z</sub>)<sup>2</sup>*. These plots demonstrate that the structure functions are computed accurately. Note that the following python modules are needed to run the test script successfully:

1. `h5py`
2. `numpy`
3. `matplotlib`


## Detailed instruction for running `Kolmogorov41`

This section provides a detailed procedure to execute `Kolmogorov41` for a given velocity or scalar field.

`Kolmogorov41-master` has a folder named `in`. This folder contains the input field files in `hdf5` format, and a parameters file named `para.yaml`. You need to provide the required input parameters in this file. The details of the entries are as follows:


### i) `para.yaml` details

#### `program: grid_switch`

You can enter `true` or `false` 

`true`: Save the structure function output as a function of the difference vector (**l**) in addition to the magnitude of the difference vector (*l*).
 
`false`: Save structure functions as a function of the magnitude of the difference vector (*l*) only.

#### `program: scalar_switch`

You can enter `true` or `false`

`true`: Calculate the structure function for a given scalar field. 

`false`: Calculate the structure function for a given velocity field. 

#### `program: 2D_switch`

You can enter `true` or `false`.

`true`: Calculate the structure function for two dimensional fields. 

`false`: Calculate the structure function for three dimensional fields.

#### `program: Only_logitudinal`

This entry is for structure functions for velocity  fields only. You can enter `true` or `false`.

`true`: Compute only the longitudinal structure function.

`false`: Compute both longitudinal and transverse structure functions.

#### `program: Number_of_OpenMP_processors`

Enter the number of OpenMP processors. Only integer values will be accepted. 

#### `grid: Nx, Ny, Nz`

The number of points along *x*, *y*, and *z* direction respectively of the  grid. Valid for both the vector and scalar fields. 
For two dimensional fields you need to provide `Nx` and `Nz`.

Only integer values will be accepted.

#### `domain_dimension: Lx, Ly, Lz`

Length of the cubical box along *x*, *y*, and *z* direction respectively. 
For two dimensional fields, you need to provide `Lx` and `Lz`. 


#### `structure_function: q1, q2`

The lower and the upper limit of the order of the structure functions to be computed.

#### `test: test_switch`

You can enter `true` or `false`

`true`: For running in test mode. Idealized velocity and scalar fields are generated internally by the code. Computed structure functions are compared with analytical results. The code is PASSED if the percentage difference between the two results is less than `1e-10`.

`false`: The "regular" mode, in which the code reads the fields from the hdf5 files in the `in` folder.

### ii) Files Required:

All the files storing the input fields should be inside the `in` folder.

#### For two dimensional fields

For vector field, two files named as `U.V1r.h5` and `U.V3r.h5` are required. Each file has one dataset.

For scalar field, one file named as `T.Fr.h5` is required. Each file has one dataset.

Size of the array stored in these files should be (`Nx,Nz`). 

*Important:* Dataset name should be the same as the file name. For example, the dataset inside the file `U.V1r.h5` should be named `U.V1r`.

#### For three dimensional fields

For vector field, three files named as `U.V1r.h5`, `U.V2r.h5`, and `U.V3r.h5` are required. Each file has one dataset.

For scalar field, one file named as `T.Fr.h5` is required. Each file has one dataset.

Size of the array stored in these files should be (`Nx, Ny, Nz`). 

*Important:* Dataset name should be the same as the file name. For example, the dataset inside the file `U.V1r.h5` should be named `U.V1r`.


### iii) Running Instructions
Open the terminal change into `Kolmogorov41-master/in` folder. Open `para.yaml` to set all the parameters. Keep all the required files compatible with the parameter file. Now, move out of the `in` folder run the command

`mpirun -np [number of MPI processors] src/Kolmogorov41.out`


### iii) Output Information

#### a) If `grid_switch` is set to `false`:

**Velocity structure functions**:

The logitudinal  structure functions of order `q1` to `q2` are stored in the files `SF.h5` and `SF_perp.h5` respectively as two dimensional arrays. Here, the first index is for different *n*, which ranges from 0 to *N<sub>l</sub>*, where *N<sub>l</sub>* is the number of gridpoints along the diagonal of the domain. The second index is for the order.

**Scalar structure functions**:

The structure functions of order `q1` to `q2` are stored in the files `SF.h5` as two dimensional array. Here, the first index is for different *n*, which ranges from 0 to *N<sub>l</sub>*, where *N<sub>l</sub>* is the number of gridpoints along the diagonal of the domain. The second index is for the order.

#### b) If `grid_switch` is set to `true`

**Velocity structure functions**:

The logitudinal and transverse structure functions of order `q` are stored in the files `SF_Grid_pll`+`q`+`.h5` and `SF_Grid_perp`+`q`+`.h5` respectively as two/three dimensional arrays for two/three dimensional input fields. 

Note: If you only want the logitudinal structure function then it will store the data for positive `lz` only as it saves computation time and computer memory

**Scalar structure functions**:

The structure functions of order `q` are stored in the files `SF_Grid_pll`+`q`+`.h5` as two/three dimensional arrays for two/three dimensional input fields. 

## Documentation

The documentation can be found in `Kolmogorov41-master/docs/index.html`

## License

`Kolmogorov41` is released under the terms of BSD New License.

