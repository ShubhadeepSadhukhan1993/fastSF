# Hybrid Parallel C++ Code: Kolmogorov41

Kolmogorov41 is an open source hybrid parallel C++ code to compute structure functions for a given velocity or scalar field.

## Getting the Source Code

Kolmogorov41 is hosted on GitHub. You can download the source code from the following link:

https://github.com/ShubhadeepSadhukhan1993/Kolmogorov41

[Click here to download] (https://github.com/ShubhadeepSadhukhan1993/Kolmogorov41)


## Required Libraries

The following libraries are required for installing and running Kolmogorov41:

1. CMake
2. Blitz++
3. YAML-cpp 
4. MPICH
5. HDF5
6. H5SI

The instructions to download and install these libraries are provided in the following website [Click Here] (http://turbulencehub.org/index.php/codes/tarang/installing-tarang/).

##  Installation Procedure
Go into the "src" folder from terminal and run the command "make". An executable named "Kolmogorov41.out" will be created inside the "src" folder.

## Input Procedure
Your input is taken from the "in" folder. A "para.yaml" file inside the "in" folder will contain all the parameters for your code. 



### i) para.yaml Instructions

#### program: grid_switch
true:  It will save the structure function output as a function of the difference vector.
 
false: It will save structure functions as a function of the magnitude of the difference vector.
#### program: scalar_switch
true: Calculate the structure function of a scalar field. 

false: Calculate the structure function of a vector field. 

#### program: 2D_switch
true: Calculate the structure function two dimensional field. 

false: Calculate the structure function three dimensional field.

#### program: Only_logitudinal
It is valid for structure function of vector fields only
true: Calculate the structure functions using the two point difference of the vector field along the difference vector direction

false: Calculate the structure functions using the two point difference of the vector field in both the parallel and perpendicular to the difference vector direction

#### grid: Nx, Ny, Nz 
It is the number of points along x, y, and z direction respectively of the  grid. Valid for both the vector and scalar fields. 
For two dimensional fields you need to provide Nx and Nz. Ny should be set to 1.


#### domain_dimension: Lx, Ly, Lz
Length of the cubical box along x, y, and z direction respectively 


#### structure_function: q1, q2
The lower and the upper limit of the order of the structure functions to be computed

### test: test_switch
true: It will take some input from inside and test the code for some known result implemented inside. (For testing lower grid size is recommended as it takes a lot of time.)

false: It is recommended for calculation from your real data which will be read from the "in" folder

### ii) Files Required:
All the required files should be inside the "in" folder.
#### Two dimension
For vector field, two files named as U.V1r.h5 and U.V3r.h5 are required.

For scalar field, one file named as T.Fr.h5 is required.

Size of the array stored in these files should be (Nx,Nz). Dataset name should be same as the file name.
#### Three dimension
For vector field, three files named as U.V1r.h5, U.V2r.h5, and U.V3r.h5 are required.

For scalar field, one file named as T.Fr.h5 is required.

Size of the array stored in these files should be (Nx,Ny,Nz). Dataset name should be same as the file name.


## Output Information
If grid_switch: false

The logitudinal and transverse structure functions of order q1 to q2 is stored in the files "SF.h5" and "SF_perp.h5" respectively as two dimensional array for both the two and three dimensional input field. Here, first index is for different l and second index is for the order.

If grid_switch: true

The logitudinal and transverse structure functions of order q is stored in the files "SF_Gridq.h5" and "SF_Grid_perpq.h5" respectively as two/three dimensional array for two/three dimensional input fields. 

Note: If you only want the logitudinal structure function then it will store the data for positive lz only as it saves computation time and computer memory


## Running Instructions
Open the terminal where you keep your "in" folder. Open "in/para.yaml" to set all the parameters. Keep all the required files compatible with the parameter file. Now  you run the command

"mpirun -np node relative-path-of-the-executable"




## License

Kolmogorov41 is released under the terms of BSD New License.

