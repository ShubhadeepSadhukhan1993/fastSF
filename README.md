# `fastSF`

`fastSF` is an open source hybrid parallel C++ code to compute structure functions for a given velocity or scalar field.

## Getting the Source Code

`fastSF` is hosted on GitHub. You can download the source code from the following link:

https://github.com/ShubhadeepSadhukhan1993/fastSF


## Installing `fastSF`
`fastSF` requires several libraries (listed later in this section) for compiling. These libraries can be installed in any directory as long as the said directory is added to the user $PATH. In this section, we will provide the guidelines to install these libraries in $HOME/local (which is generally recommended).

### Environment Variables
The following environment variables need to be set before compiling `fastSF`. (You may append these lines to `$HOME/.bashrc`):


`export PATH=$HOME/local/bin:$PATH`

`export PKG_CONFIG_DISABLE_UNINSTALLED=true`

`export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig:$PKG_CONFIG_PATH`

`export HDF5_ROOT=$HOME/local`

`export CPATH=$HOME/local/include/:$CPATH`

`export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH`

`export LIBRARY_PATH=$HOME/local/lib:$LIBRARY_PATH`

`export MANPATH=$HOME/local/share/man/:$MANPATH`
`

### Required Libraries

The following libraries are required for installing and running fastSF:
 
* [`Blitz++`](https://github.com/blitzpp/blitz) (Version 1.0.2)-
All array manipulations are performed using the `Blitz++` library. Download `Blitz++` from [here](https://github.com/blitzpp/blitz). After downloading, change to the `blitz-master` directory and enter the following commands

	`CC=gcc CXX=g++ ./configure --prefix=$HOME/local`
	
	`make install`

* [`YAML-cpp`](https://github.com/jbeder/yaml-cpp/releases/tag/release-0.3.0)(Version 0.3.0) - 
	The input parameters are stored in the `para.yaml` file which needs the `YAML-cpp` library to parse. Download `YAML-cpp` from [here](https://github.com/jbeder/yaml-cpp/releases/tag/release-0.3.0). Extract the zip/tar file and change the `yaml-cpp-release-0.3.0` directory. Important: Please ensure that [CMake](https://cmake.org) is installed in your system. Enter the following commands:
	
	`CC=gcc CXX=g++ cmake -DCMAKE_INSTALL_PREFIX=$HOME/local`
	
	`make install`

	
* An `MPI` (Message Passing Interface) Library - 
`fastSF` uses `MPI` for parallelism. The software was tested using [`MPICH`](https://www.mpich.org)(Version 3.3.2), however, any standard `MPI-1` implementation should be sufficient.  Here, we will provide instructions for installing `MPICH`. Download `MPICH` from [here](https://www.mpich.org/downloads/). After extraction, change to `mpich-3.3.2` folder and enter the following:

	`CC=gcc CXX=g++ ./configure --prefix=$HOME/local`
	
	`make install`

* [`HDF5`](http://turbulencehub.org/wp-content/uploads/Download_Files/hdf5-1.8.20.tar.bz2)(Version 1.8.20) -
The output files are written in HDF5 format. Download `HDF5` from [here](http://turbulencehub.org/wp-content/uploads/Download_Files/hdf5-1.8.20.tar.bz2). After extracting the tar file, change to `hdf5-1.8.20` and enter the following:

	`CC=mpicc CXX=mpicxx ./configure --prefix=$HOME/local --enable-parallel --without-zlib`

	`make install`

* [`H5SI`](https://github.com/anandogc/h5si)(Version 1.1.1) - 
This library is used for simplifying the input-output operations of `HDF5`. Download `H5SI` from [here](https://github.com/anandogc/h5si). After downloading the zip file, extract it and change to `h5si-master/trunk`. Important: Please ensure that [CMake](https://cmake.org) is installed in your system. Enter the following:

	`CXX=mpicxx cmake . -DCMAKE_INSTALL_PREFIX=$HOME/local`
	
	`make`
	
	`make install`


IMPORTANT: 

* Note that `fastSF` is not compatible with higher versions `YAML-cpp`(> 0.3.0). 


###  Compiling instruction

After downloading `fastSF`, change into `fastSF/src` directory and run the command `make` in the terminal. An executable named `fastSF.out` will be created inside the `fastSF/src` folder.

## Testing `fastSF`
`fastSF` offers an automated testing process to validate the code. The relevant test scripts can be found in the `tests/` folder of the code. To execute the tesing process, change into `fastSF` and run the command 

`bash runTest.sh`. 

The code then runs four test cases; these are as follows. 

* In the first case, the code will generate a 2D velocity field given by **u** = [*x, z*], and compute the structure functions for the given field. For this case, the longitudinal structure functions should equal *l<sup>q</sup>*. 

* In the second case, the code will generate a 2D scalar field given by *T = x + z*, and compute the structure functions for the given field. For this case, the structure functions should equal *(l<sub>x</sub> + l<sub>z</sub>)<sup>q</sup>*.

* In the third case, the code will generate a 3D velocity field given by **u** = [*x, y, z*], and compute the structure functions for the given field. For this case, the longitudinal structure functions should equal *l<sup>q</sup>*. 

* In the fourth case, the code will generate a 3D scalar field given by *T = x + y + z*, and compute the structure functions for the given field. For this case, the structure functions should equal *(l<sub>x</sub> + l<sub>y</sub> + l<sub>z</sub>)<sup>q</sup>*.

For the above cases, `fastSF` will compare the computed structure functions with the analytical results. If the percentage difference between the two values is less than 10<sup>-10</sup>, the code is deemed to have passed. 

Finally, for visualization purpose, the python script `test/test.py` is invoked. This script generates the plots of the second and third-order longitudinal structure functions versus *l*, and the density plots of the computed second-order scalar structure functions and *(l<sub>x</sub> + l<sub>z</sub>)<sup>2</sup>*. For the 3D scalar field, the density plots of the computed second-order scalar structure functions for *l<sub>y</sub> = 0.5* and *(l<sub>x</sub> + 0.5 + l<sub>z</sub>)<sup>2</sup>* are generated. These plots demonstrate that the structure functions are computed accurately. Note that the following python modules are needed to run the test script successfully:

1. `h5py`
2. `numpy`
3. `matplotlib`


## Detailed instruction for running `fastSF`

This section provides a detailed procedure to execute `fastSF` for a given velocity or scalar field.

### i) Files Required and HDF5 Schema:

All the files storing the input fields needs to be in the `hdf5` format and stored inside the `in` folder. All the input datasets are required to be precisely as per the following schema.

#### For two dimensional fields

For vector field, two datasets files are required:

1. A 2D dataset storing the x-component of the vector / velocity field. This dataset has dimensions (`Nx,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V1r.h5` and the dataset as `U.V1r`.

2. A 2D dataset storing the z-component of the vector / velocity field. This dataset has dimensions (`Nx,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V3r.h5` and the dataset as `U.V3r`. 

For scalar field, one dataset is required:

1. A 2D dataset storing the scalar field. This dataset has dimensions (`Nx,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `T.Fr.h5` and the dataset as `T.Fr`.

#### For three dimensional fields

For vector field, three `hdf5` files are required:

1. A 3D dataset storing the x-component of the vector / velocity field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V1r.h5` and the dataset as `U.V1r`.

2. A 3D dataset storing the y-component of the vector / velocity field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V2r.h5` and the dataset as `U.V2r`.

3. A 3D dataset storing the z-component of the vector / velocity field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `U.V3r.h5` and the dataset as `U.V3r`.

For scalar field, For vector field, one `hdf5` file is required:

1. A 3D dataset storing the scalar field. This dataset has dimensions (`Nx,Ny,Nz`) storing real double-precision floating point values. The names of the dataset and the hdf5 file containing the dataset can be arbitrary, but they need to be specified by the user via command line arguments during the execution of `fastSF`. If the user prefers not to use command line arguments, the hdf5 file by default should be named as `T.Fr.h5` and the dataset as `T.Fr`.


**IMPORTANT:** For vector fields, it must be ensured that the dimensions of `Ux`, `Uy`, and `Uz` are the same, otherwise the code will throw an error.

### ii) `para.yaml` details

The user can specify the relevant parameters via command line or via parameters file. If no command-line options are given, the entries in the parameters file will be taken. First, we will explain how to use the parameters file, and then we will proceed to command-line options.

`fastSF` has a folder named `in`. This folder contains the input field files in `hdf5` format, and a parameters file named `para.yaml`. You need to provide the required input parameters in this file. The details of the entries are as follows:




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

#### `program: Processors_X`

The number of processors in x-direction. Only integer values are accepted. Note that this value should be an integer factor of the total number of processors.

#### `grid: Nx, Ny, Nz` (Only applicable if test switch is set to `true`, in which case the code generates the input fields)

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

### iii) Running Instructions and Command-Line Arguments 
To run `fastSF`, change to `fastSF` directory. Ensure that the input hdf5 files follow the schema described in the previous subsection. If you want all the relevant parameters to be read from "in/para.yaml", you can simply type the following:

`mpirun -np [number of MPI processors] src/fastSF.out`
 
`fastSF` allows the user to pass the input parameters using command line arguments as well. If inputs are provided via the command-line, the corresponding inputs read from the "in/para.yaml" file get overriden. The user can also specify the input and output hdf5 file names via the command-line. The command line arguments are given as follows:

`mpirun -np [number of MPI processors] src/fastSF.out -s [scalar_switch]` 
`-d [2D_switch] -l [Only_longitudinal] -p [Processors_X] -X [Nx] -Y [Ny]`
`-Z [Nz] -x [Lx] -y [Ly] -z [Lz] -1 [q1] -2 [q2] -t [test_switch]`
`-U [Name of the hdf5 file containing the dataset storing Ux]`
`-u [Name of the dataset storing Ux]`
`-V [Name of the hdf5 file containing the dataset storing Uy]`
`-v [Name of the dataset storing Uy]`
`-W [Name of the hdf5 file containing the dataset storing Uz]`
`-w [Name of the dataset storing Uz]`
`-Q [Name of the hdf5 file containing the dataset storing T (scalar)]`
`-q [Name of the dataset storing T (scalar)]`
`-P [Name of the hdf5 file storing the transverse structure functions]`
`-L [Name of the hdf5 file storing the longitudinal structure functions]`
`-h [Help]`

The user need not give all the command line arguments; the arguments that are not provided will be read by the `in/para.yaml` file. For, if the user wants to run `fastSF` with 16 processors with 4 processors in x direction, and wants to compute only the longitudinal structure functions, the following command should be entered:

`mpirun -np 16 ./src/fastSF.out -p 4 -l true`

In this case, the number of processors in the x-direction and the longitudinal structure function switch will be taken via the command line. The rest of the parameters will be taken from the `in/para.yaml` file.

**Note:** `Nx`, `Ny`, and `Nz` should be specified only if test case is "on", in which case the code generates the input fields.

### iv) Output Information

Unless specified otherwise by the user via command-line arguments, the following output files are written by `fastSF`.

**Velocity structure functions**:

The logitudinal and transverse structure functions of order `q` are stored in the files `SF_Grid_pll.h5` and `SF_Grid_perp.h5` respectively. `SF_Grid_pll.h5` and  `SF_Grid_perp.h5` have datasets named `SF_Grid_pll`+`q`  and `SF_Grid_perp`+`q` respectively. These datasets store two/three dimensional arrays for two/three dimensional input fields.

**Scalar structure functions**:

The structure functions of order `q` are stored in the file `SF_Grid_scalar.h5` consisting of the datasets named `SF_Grid_scalar`+`q`. 

## Memory Requirements

The memory requirement per processor for running `fastSF` depends primarily on the resolution of the grid. The memory requirement also depends on the number of orders of the structure functions to be computed, number of processors *P*, and the distribution of processors in *x* and *y* (or *z*) directions. The memory requirement *M* (in bytes) can be estimated as follows:

### Two dimensional scalar field:

*M* = (20 + 4*n*)*N<sub>x</sub>N<sub>z</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 32*P*.

### Three dimensional scalar field:

*M* = (16 + 2*n*)*N<sub>x</sub>N<sub>y</sub>N<sub>z</sub>* + 4*N<sub>x</sub>N<sub>y</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 40*P*.

### Two dimensional vector field:


*M* = (44 + 4*n*)*N<sub>x</sub>N<sub>z</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 32*P*, if only longitudinal structure functions are to be computed:

*M* = (44 + 8*n*)*N<sub>x</sub>N<sub>z</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 40*P*, if both longitudinal and transverse structure functions are to be computed.

### Three dimensional vector field:

*M* = (56 + 2*n*)*N<sub>x</sub>N<sub>y</sub>N<sub>z</sub>* + 4*N<sub>x</sub>N<sub>y</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 40*P*, if only longitudinal structure functions are to be computed.

*M* = (56 + 4*n*)*N<sub>x</sub>N<sub>y</sub>N<sub>z</sub>* + 4*N<sub>x</sub>N<sub>y</sub>* + 8(*N<sub>x</sub>/p<sub>x</sub> + N<sub>z</sub>p<sub>x</sub>/P*) + 48*P*, if both longitudinal and transverse structure functions are to be computed.

In the above expressions, *p<sub>x</sub>* refers to the number of processes in *x* direction and *P* refers to the total number of processors. Note that for large *N<sub>z</sub>*, the first term dominates the remaining terms; thus the memory requirement can be quickly estimated using the first term only. 

## Documentation and Validation

The documentation can be found in `fastSF/docs/index.html`. 

The validation of `fastSF` is reported [here](https://github.com/ShubhadeepSadhukhan1993/fastSF/blob/master/docs/Verification.md).

## Contributions and Bug Reports

We welcome contributions to this project. If you wish to contribute, please create a branch with a [pull request](https://github.com/ShubhadeepSadhukhan1993/fastSF/pulls) and the changes can be discussed there.

If you find a bug in the or errors in the documentation, please open a [new issue](https://github.com/ShubhadeepSadhukhan1993/fastSF/issues/new) in the Github repository and report the bug or the error. Please provide sufficient information for the bug to be reproduced.  

## License

`fastSF` is released under the terms of BSD New License.

