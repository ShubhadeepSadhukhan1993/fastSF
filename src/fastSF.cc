/********************************************************************************************************************************************
 * fastSF
 *
 * Copyright (C) 2020, Mahendra K. Verma
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     1. Redistributions of source code must retain the above copyright
 *        notice, this list of conditions and the following disclaimer.
 *     2. Redistributions in binary form must reproduce the above copyright
 *        notice, this list of conditions and the following disclaimer in the
 *        documentation and/or other materials provided with the distribution.
 *     3. Neither the name of the copyright holder nor the
 *        names of its contributors may be used to endorse or promote products
 *        derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 ********************************************************************************************************************************************
 */

/*! \file fastSF.cc
 *
 *  \brief Code to compute structure functions using velocity or scalar field data.
 *
 *  \author Shubhadeep Sadhukhan, Shashwat Bhattacharya, Mahendra K. Verma
 *  \date Feb 2020
 *  \copyright New BSD License
 *
 ********************************************************************************************************************************************
 */

#include "h5si.h"
#include <yaml-cpp/yaml.h>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include <sstream>
#include <blitz/array.h>
#include <omp.h>
#include <mpi.h>
#include <sys/time.h>
#include <limits.h>
#include <unistd.h>
using namespace std;
using namespace blitz;

//Function declarations
void get_Inputs(int argc, char* argv[]); 
void write_3D(Array<double,3>, string);
void write_4D(Array<double,4>, string);
void read_2D(Array<double,2>, string, string, string);
string int_to_str(int);
bool str_to_bool(string);
void VECTOR_TEST_CASE_3D();
void VECTOR_TEST_CASE_2D();
void SCALAR_TEST_CASE_2D();
void SCALAR_TEST_CASE_3D();
void compute_time_elapsed(timeval, timeval, double&);

void read_3D(Array<double,3>, string, string, string);


void SFunc2D(Array<double,2>, Array<double,2>);


void SFunc_long_2D(Array<double,2>, Array<double,2>);


void SFunc3D(Array<double,3>, Array<double,3>, Array<double,3>);


void SFunc_long_3D(Array<double,3>, Array<double,3>, Array<double,3>);
void Read_Init(Array<double,2>&, Array<double,2>&);
void Read_Init(Array<double,3>&, Array<double,3>&, Array<double,3>&);
void Read_Init(Array<double,2>&);
void Read_Init(Array<double,3>&);



void SF_scalar_3D(Array<double,3>);


void SF_scalar_2D(Array<double,2>);

void Read_fields();
void resize_SFs();
void calc_SFs();
void write_SFs();
void test_cases();
void show_checklist();

bool compare(Array<int,1> , Array<int,1> );
void calculate_grid_spacing();
void resize_input();
void get_input_shape(string , string , string , Array<int,1>&);
void help_command();

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the input 3D scalar field.
 ********************************************************************************************************************************************
 */
Array <double,3> T;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the x-component of the input 3D velocity field.
 ********************************************************************************************************************************************
 */
Array <double,3> V1;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the y-component of the input 3D velocity field.
 ********************************************************************************************************************************************
 */
Array <double,3> V2;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the z-component of the input 3D velocity field.
 ********************************************************************************************************************************************
 */
Array <double,3> V3;

/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the input 2D scalar field.
 ********************************************************************************************************************************************
 */
Array <double,2> T_2D;

/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the x-component of the input 2D velocity field.
 ********************************************************************************************************************************************
 */
Array<double,2> V1_2D;

/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the z-component of the input 2D velocity field.
 ********************************************************************************************************************************************
 */
Array<double,2> V3_2D;


/**
 ********************************************************************************************************************************************
 * \brief   4D array storing the computed longitudinal structure functions as function of the displacement vector.
 *
 *          This array stores the structure functions as function of the displacement vector \f$ \mathbf{l} = (l_x, l_y, l_z )\f$. The fourth
 *          dimension corresponds to the order of the structure functions that are calculated.
 ********************************************************************************************************************************************
 */
Array<double,4> SF_Grid_pll;

/**
 ********************************************************************************************************************************************
 * \brief   4D array storing the computed transverse structure functions as function of the displacement vector.
 *
 *          This array stores the structure functions as function of the displacement vector \f$ \mathbf{l} = (l_x, l_y, l_z )\f$. The fourth
 *          dimension corresponds to the order of the structure functions that are calculated.
 ********************************************************************************************************************************************
 */
Array<double,4> SF_Grid_perp;

/**
 ********************************************************************************************************************************************
 * \brief   4D array storing the computed scalar structure functions as function of the displacement vector.
 *
 *          This array stores the structure functions as function of the displacement vector \f$ \mathbf{l} = (l_x, l_y, l_z )\f$. The fourth
 *          dimension corresponds to the order of the structure functions that are calculated.
 ********************************************************************************************************************************************
 */
Array<double,4> SF_Grid_scalar;


/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the computed longitudinal structure functions as function of the displacement vector.
 *
 *          This array stores the structure functions as function of the displacement vector \f$ \mathbf{l} = (l_x, l_z )\f$. The third
 *          dimension corresponds to the order of the structure functions that are calculated.
 ********************************************************************************************************************************************
 */
Array<double,3> SF_Grid2D_pll;

 /**
 ********************************************************************************************************************************************
 * \brief   3D array storing the computed transverse structure functions as function of the displacement vector.
 *
 *          This array stores the structure functions as function of the displacement vector \f$ \mathbf{l} = (l_x, l_z )\f$. The third
 *          dimension corresponds to the order of the structure functions that are calculated.
 ********************************************************************************************************************************************
 */
Array<double,3> SF_Grid2D_perp;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the computed scalar structure functions as function of the displacement vector.
 *
 *          This array stores the structure functions as function of the displacement vector \f$ \mathbf{l} = (l_x, l_z )\f$. The third
 *          dimension corresponds to the order of the structure functions that are calculated.
 ********************************************************************************************************************************************
 */
Array<double,3> SF_Grid2D_scalar;

/**
 ********************************************************************************************************************************************
 * \brief   This variable decides whether the structure functions are to be calculated using 2D or 3D velocity field data.
 *
 * If "true", then the code will read 2D velocity fields and calculate the corresponding structure functions. Otherwise, it will read 3D velocity
 * fields. Entered by the user
 ********************************************************************************************************************************************
 */
bool two_dimension_switch;

/**
 ********************************************************************************************************************************************
 * \brief    This variable decides whether the scalar or velocity structure functions are to be evaluated.
 *
 * If the value is "true", then scalar structure functions will be evaluated, else the vector (velocity) structure functions will be evaluated.
 ********************************************************************************************************************************************
 */
bool scalar_switch;


/**
 ********************************************************************************************************************************************
 * \brief   This variable decides whether the test cases are required to run.
 *
 * If the value is "true", then the input fields will be generated internally and test cases will be run.
 ********************************************************************************************************************************************
 */
bool test_switch;

/**
 ********************************************************************************************************************************************
 * \brief   Number of gridpoints in the \f$ x \f$ direction.
 ********************************************************************************************************************************************
 */
int Nx;

/**
 ********************************************************************************************************************************************
 * \brief   Number of gridpoints in the \f$ y \f$ direction.
 ********************************************************************************************************************************************
 */
int Ny;

/**
 ********************************************************************************************************************************************
 * \brief   Number of gridpoints in the \f$ z \f$ direction.
 ********************************************************************************************************************************************
 */
int Nz;


/**
 ********************************************************************************************************************************************
 * \brief   The first order of the range of orders of the structure functions to be computed. Entered by the user.
 ********************************************************************************************************************************************
 */
int q1;

/**
 ********************************************************************************************************************************************
 * \brief   The last order of the range of orders of the structure functions to be computed. Entered by the user.
 ********************************************************************************************************************************************
 */
int q2;


/**
 ********************************************************************************************************************************************
 * \brief   This variable decides whether to calculate both transverse and longitudinal structure functions or only the longitudinal structure
 * functions.
 *
 * If the value is false, the code calculates both transverse and longitudinal structure functions. Else, it calculates only the longitudinal
 * structure functions.
 ********************************************************************************************************************************************
 */
bool longitudinal;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the distance between two consecutive gridpoints in the \f$ x \f$ direction.
 *
 ********************************************************************************************************************************************
 */
double dx;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the distance between two consecutive gridpoints in the \f$ y \f$ direction.
 *
 ********************************************************************************************************************************************
 */
double dy;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the distance between two consecutive gridpoints in the \f$ z \f$ direction.
 *
 ********************************************************************************************************************************************
 */
double dz;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the rank of the MPI process.
 *
 ********************************************************************************************************************************************
 */
int rank_mpi;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the length of the domain.
 *
 ********************************************************************************************************************************************
 */
double Lx;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the width of the domain.
 *
 ********************************************************************************************************************************************
 */
double Ly;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the height of the domain.
 *
 ********************************************************************************************************************************************
 */
double Lz;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the number of the MPI process.
 *
 ********************************************************************************************************************************************
 */
int P;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the number of the MPI process along x axis.
 *
 ********************************************************************************************************************************************
 */
int px;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the input file for the x-component of the velocity field.
 *
 ********************************************************************************************************************************************
 */
string UName="U.V1r";

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the input file for the y-component of the velocity field.
 *
 ********************************************************************************************************************************************
 */
string VName="U.V2r";
/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the input file for the z-component of the velocity field.
 *
 ********************************************************************************************************************************************
 */
string WName="U.V3r";

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the input file for scalar field.
 *
 ********************************************************************************************************************************************
 */
 string TName="T.Fr";

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the dataset for the x-component of the velocity field.
 *
 ********************************************************************************************************************************************
 */
string UdName="U.V1r";


/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the dataset for the y-component of the velocity field.
 *
 ********************************************************************************************************************************************
 */
string VdName="U.V2r";


/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the dataset for the z-component of the velocity field.
 *
 ********************************************************************************************************************************************
 */
string WdName="U.V3r";


/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the dataset of the scalar field.
 *
 ********************************************************************************************************************************************
 */
string TdName="T.Fr";

 /**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the output file for parallel structure function of velocity field.
 *
 ********************************************************************************************************************************************
 */
string SF_Grid_pll_name = "SF_Grid_pll";

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the output file for perpendicular structure function of velocity field.
 *
 ********************************************************************************************************************************************
 */
string SF_Grid_perp_name = "SF_Grid_perp";

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the name for the output file for structure function of scalar field.
 *
 ********************************************************************************************************************************************
 */
 string SF_Grid_scalar_name = "SF_Grid_scalar";


/**
 ********************************************************************************************************************************************
 * \brief   The main function of the "fastSF".
 *
 * \param   argc is the number of command-line arguments passed by the user including the name of the program.
 * \param   argv is the array of character pointers listing all the arguments.
 *
 ********************************************************************************************************************************************
 */
int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_mpi);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    
    
    //Initiallizing h5si
    h5::init();
    timeval start_pt, end_pt, start_t, end_t;

    //Record the time of starting of the program
    gettimeofday(&start_t,NULL);

    double elapsedt=0.0;
    double elapsepdt=0.0;
    
    //Get the input parameters
    get_Inputs(argc, argv);
    
    
    
  	



    //Resizing the input fields
    Read_fields();

    if (rank_mpi==0) {
    	cout<<"\nNumber of processors in x direction: "<<px<<endl;
    	if (two_dimension_switch) {
        	cout<<"Number of processors in z direction: "<<P/px<<endl;
    	}
    	else {
        	cout<<"Number of processors in y direction: "<<P/px<<endl;
    	}
  	}  

 	if (px > P) {
        if (rank_mpi==0) {
            cout<<"ERROR! Number of processors in x direction has to be less than or equal to the total number of processors! Aborting.."<<endl;
        }
        h5::finalize();
        MPI_Finalize();
        exit(1);
    }
    if (Nx/2%px != 0) {
        if (rank_mpi==0){
            cout<<"ERROR! Number of processors in x direction should be less or equal to Nx/2 and some power of 2\n Aborting...\n";
        }
        h5::finalize();
        MPI_Finalize();
        exit(1);
    }

    int N2;

    if (two_dimension_switch) {
        N2 = Nz;
    }
    else {
        N2 = Ny;
    }

    if (N2/2%(P/px) != 0) {
        if (rank_mpi==0){
            cout<<"ERROR! Number of processors in y (or z) direction should be less or equal to Ny/2 (or Nz/2) and some power of 2\n Aborting...\n";
        }
        h5::finalize();
        MPI_Finalize();
        exit(1);
    } 


    //Resize the structure function array according to the type of inputs
    resize_SFs();

    
    //Record the time of starting the parallel processing
    gettimeofday(&start_pt,NULL);

    //Calculating the structure functions
    calc_SFs();


    //Record the time of ending of parallel processing
    gettimeofday(&end_pt,NULL);
    
 
    //Write the SF array to disk
    write_SFs();

    if (test_switch){
        test_cases();
    }

    //Record the time when the program ends
    gettimeofday(&end_t,NULL);
    
    compute_time_elapsed(start_t, end_t, elapsedt);
    compute_time_elapsed(start_pt, end_pt, elapsepdt);
    
    
    if (rank_mpi==0) {
        cout<<"\nTime elapsed for the parallel part: "<<elapsepdt<<endl;
        cout<<"\nTotal time elapsed: "<<elapsedt<<endl;
        cout<<"\nProgram ends."<<endl;
   }

    h5::finalize();
    MPI_Finalize();
    return 0;
}

/**
******************************************************************************************************************************
*\brief  Function to calculate dx, dy and dz
*
******************************************************************************************************************************
*/
void calculate_grid_spacing(){
	//Specify the values of dx, dy, and dz
    if (Nx==1){dx=0;}
    else{
      dx=Lx/double(Nx-1);}
    if (Ny==1){dy=0;}
    else{
      dy=Ly/double(Ny-1);}
    if (Nz==1){dz=0;}
    else{
      dz=Lz/double(Nz-1);
  	}
}


/**
*************************************************************************************************************************************
*\brief     Function to covert string to bool
*           
*************************************************************************************************************************************
*/
bool str_to_bool(string s){
	if (s=="true" || s=="1"){
		return true;
	}
	else if (s=="false" || s=="0"){
		return false;
	}
	else{
        if (rank_mpi==0)
            cout<<"Invalid input\n";
		exit(1);
	}
}


/**
********************************************************************************************************************************
*\brief    Function to check whether the file exists or not. Furthermore it saves the shape of the dataset.
*
*\param fold is data path
*\param file is the file name
*\param dset is the dataset name
*\param s is array in which shape will be saved
*
********************************************************************************************************************************
*/
void get_input_shape(string fold, string file, string dset, Array<int,1>& s){
	ifstream file_name(fold+file+".h5");
	int dim;
	s.resize(4);
	if (file_name.is_open()){
    	file_name.close();
    	h5::File f(fold+file+".h5", "r");
    	h5::Dataset data_set;
    	data_set=(f[dset]);
    	dim=data_set.shape().size();
  		if (dim==2){
  			two_dimension_switch=true;
  			Nx=data_set.shape()[0];
  			Ny=1;
  			Nz=data_set.shape()[1];
  		}
  		else{
  			two_dimension_switch=false;
			Nx=data_set.shape()[0];
  			Ny=data_set.shape()[1];
  			Nz=data_set.shape()[2];  			
  		}
  		
  		s(0)=dim;
  		s(1)=Nx;
  		s(2)=Ny;
  		s(3)=Nz;
  		
  		
  	}
  	else{
    	file_name.close();
    	if (rank_mpi==0){
        	cerr<<"\nDesired file does not exist\n\n";
        	show_checklist();
    	}
    	h5::finalize();
        MPI_Finalize();
    	exit(1);
    }
	
}



/**
******************************************************************************************************************************
\brief Resize input arrays
*
******************************************************************************************************************************
*/

void resize_input(){
	if(two_dimension_switch){
        if (scalar_switch) {
            T_2D.resize(Nx, Nz);
        }
        else {
            V1_2D.resize(Nx, Nz);
            V3_2D.resize(Nx, Nz);
        }
        
    }
    else{
        if (scalar_switch) {
            T.resize(Nx, Ny, Nz);
        }
        else {
            V1.resize(Nx,Ny,Nz);
            V2.resize(Nx,Ny,Nz);
            V3.resize(Nx,Ny,Nz);
        }
        
    }
}


/**
********************************************************************************************************************************
*\brief compare two integer arrays
*
*
*\param A is the first array
*\param B is the second array
*
********************************************************************************************************************************
*/
bool compare(Array<int,1> A, Array<int,1> B){
	int N=A.size();
	for (int i=0; i<N; i++){
		if (A(i)!=B(i)){
			return false;
		}
	}
	return true;
}



/**
*************************************************************************************************************************************
*\brief     Function to generate or read the input fields.
*
*
*\param     UName is the name of the hdf5 file and dataset storing the x-component of velocity / vector field. 
*\param     VName is the name of the hdf5 file and dataset storing the y-component of velocity / vector field.
*\param     WName is the name of the hdf5 file and dataset storing the z-component of velocity / vector field.
*\param     TName is the name of the hdf5 file and dataset storing the scalar field.
* 
*************************************************************************************************************************************
*/
void Read_fields() {
    //Defining the input fields
    if (!test_switch){
    	if (rank_mpi==0){
            cout<<"Reading from the hdf5 files\n";
        }
        Array<int,1> s1,s2, s3;
        
        if (two_dimension_switch){
            if (scalar_switch) {
            	get_input_shape("in/", TName, TdName, s1);
                resize_input();
                calculate_grid_spacing();
                read_2D(T_2D,"in/", TName, TdName);
            }
            else {
            	get_input_shape("in/", UName, UdName, s1);
                get_input_shape("in/", WName, WdName, s2);
                
                if (!compare(s1,s2)){
                	if (rank_mpi==0){
            			cerr<<"\nIncompatible dimension data\n\n";
            			show_checklist();
        			}
                	h5::finalize();
        			MPI_Finalize();
        			exit(1);
                }
                
                resize_input();
                calculate_grid_spacing();
                read_2D(V1_2D,"in/", UName, UdName);
                read_2D(V3_2D,"in/", WName, WdName);
            }
        }
        else{
        	
            if (scalar_switch) {
            	get_input_shape("in/", TName, TdName, s1);
            	resize_input();
            	calculate_grid_spacing();
                read_3D(T, "in/", TName, TdName);
            }
            else {
            	
            	get_input_shape("in/", UName, UdName, s1);
            	get_input_shape("in/", VName, VdName, s2);
            	get_input_shape("in/", WName, WdName, s3);
            	
            	if (!compare(s1,s2)){
            		if (rank_mpi==0){
            			cerr<<"\nIncompatible dimension data\n\n";
            			show_checklist();
        			}
        			h5::finalize();
        			MPI_Finalize();
        			exit(1);
                }
                if (!compare(s2,s3)){
                	if (rank_mpi==0){
            			cerr<<"\nIncompatible dimension data\n\n";
            			show_checklist();
        			}
                	h5::finalize();
        			MPI_Finalize();
        			exit(1);
                }
                if (!compare(s3,s1)){
                	if (rank_mpi==0){
            			cerr<<"\nIncompatible dimension data\n\n";
            			show_checklist();
        			}
                	h5::finalize();
        			MPI_Finalize();
        			exit(1);
                }
            	resize_input();
            	calculate_grid_spacing();
                read_3D(V1, "in/", UName, UdName);
                read_3D(V2, "in/", VName, VdName);
                read_3D(V3, "in/", WName, WdName);
            }
        }
    } 
    else {
        if (rank_mpi==0){
            cout<<"\nWARNING: The code is running in TEST mode. It will generate velocity / scalar fields and will take them as inputs.\n";
        }
        resize_input();
        calculate_grid_spacing();
        if (two_dimension_switch) {
            if (scalar_switch) {
                Read_Init(T_2D);
            }
            else {
                Read_Init(V1_2D, V3_2D);
            }
        }
        else {
            if (scalar_switch) {
                Read_Init(T);
            }
            else {
                Read_Init(V1, V2, V3);
            }
        }
    }
}

/**
*************************************************************************************************************************************
*\brief     Function resize the structure function arrays according to the inputs.
*  
*************************************************************************************************************************************
*/
void resize_SFs(){
    if (rank_mpi==0) {
        if (not two_dimension_switch) {
            if (scalar_switch) {
                SF_Grid_scalar.resize(Nx/2, Ny/2, Nz/2, q2-q1+1);
                SF_Grid_scalar = 0; 
            }
            else {
                SF_Grid_pll.resize(Nx/2, Ny/2, Nz/2, q2-q1+1);
                SF_Grid_pll = 0;
                if (not longitudinal) {
                    SF_Grid_perp.resize(Nx/2, Ny/2, Nz/2, q2-q1+1);
                    SF_Grid_perp = 0;
                }
            }
            
        }
        else {
            if (scalar_switch) {
                SF_Grid2D_scalar.resize(Nx/2, Nz/2, q2-q1+1);
                SF_Grid2D_scalar = 0; 
            }
            else {
                SF_Grid2D_pll.resize(Nx/2, Nz/2, q2-q1+1);
                SF_Grid2D_pll = 0; 
                if (not longitudinal) {
                    SF_Grid2D_perp.resize(Nx/2, Nz/2, q2-q1+1);
                    SF_Grid2D_perp = 0;
                }
            }
        }   
    }
}

/**
*************************************************************************************************************************************
*\brief     Function to compute the structure functions based on the inputs provided by the user.
*************************************************************************************************************************************
*/
void calc_SFs() {
    if (two_dimension_switch){
        if (scalar_switch) {
            SF_scalar_2D(T_2D);
        }
        else {
            if (longitudinal) {
                SFunc_long_2D(V1_2D, V3_2D);
            } 
            else {
                SFunc2D(V1_2D, V3_2D);
            }
        }
    }
    
    else {
        if (scalar_switch) {
            SF_scalar_3D(T); 
        }
        else {
            if (longitudinal) {
                SFunc_long_3D(V1, V2, V3);
            }
            else {
                SFunc3D(V1, V2, V3);
            }
        }
    }
}

/**
*************************************************************************************************************************************
*\brief     Function to write the structure function arrays to the disk.
*************************************************************************************************************************************
*/
void write_SFs() {
    if (rank_mpi==0){
        mkdir("out",0777);
       

        if (two_dimension_switch) {
            
            if (scalar_switch){
                write_3D(SF_Grid2D_scalar, SF_Grid_scalar_name);
            }
            else {
                write_3D(SF_Grid2D_pll, SF_Grid_pll_name);    
                if (not longitudinal) {
                    write_3D(SF_Grid2D_perp, SF_Grid_perp_name);
                }
            }
            cout<<"\nWriting completed\n";
        }
        else {
            
            if (scalar_switch){
                write_4D(SF_Grid_scalar, SF_Grid_scalar_name);
            }
            else {
            	write_4D(SF_Grid_pll, SF_Grid_pll_name);    
                if (not longitudinal) {
                    write_4D(SF_Grid_perp, SF_Grid_perp_name);
                }
            }
            cout<<"\nWriting completed\n";
        }
    }
}

/**
*************************************************************************************************************************************
*\brief     Function to test the correctness of the code.
*************************************************************************************************************************************
*/
void test_cases() {
    if(rank_mpi==0){
        cout<<"\nCOMMENCING TESTING OF THE CODE.\n";
        if (scalar_switch){
            if (two_dimension_switch){
                SCALAR_TEST_CASE_2D();
            }
            else{
                SCALAR_TEST_CASE_3D();
            }
        }
        else{
            if (two_dimension_switch){
                VECTOR_TEST_CASE_2D();
            }
            else{
                VECTOR_TEST_CASE_3D();
            }
        }
    }
    
}


/**
*************************************************************************************************************************************
*\brief     Function to obtain the \$ x \$ and \$ y \$ ranks of the processor.
*
*
*\param     rank is the rank of the processor.
*\param     py is the number of processors in \$ y \$ direction.
*\param     rankx is the \$ x \$ rank of the processor.
*\param     ranky is the \$ y \$ rank of the processor.
* 
*************************************************************************************************************************************
*/
void get_rank(int rank, int py, int& rankx, int& ranky){
    ranky=rank%py;
    rankx=(rank-ranky)/py;
}

/**
*************************************************************************************************************************************
*\brief     Function to allocate the index list along a given direction for a particular rank.
*           
*     
*
* \param    index_list stores the list of indices.
* \param    Nx is half of the number of points along the given direction.
* \param    px is the number of processors along the given direction.
* \param    rank is the \$ x \$ or \$ y \$ rank of the processor.
*************************************************************************************************************************************
*/
void compute_index_list(Array<int,1>& index_list, int Nx, int px, int rank){
    int list_size=Nx/px;
    index_list.resize(list_size);
    for (int i=0; i<list_size; i+=2){
        index_list(i)=rank+i*px;
        if (px!=Nx) {
          index_list(i+1)=Nx-1-index_list(i);
        }
    }
}



/**
*************************************************************************************************************************************
*\brief     Function to distribute load equally among all the MPI processors.
*           
*           
* \param    index_list stores the matrix of the indices for different ranks.
* \param    Nx is the number of points along \$ x \$ direction.
* \param    Ny is the number of points along \$ y \$ (or \$ z \$) direction. 
*************************************************************************************************************************************
*/
void compute_index_list(Array<int,3>& index_list, int Nx, int Ny){
    int list_size=(Nx*Ny)/(4*P);
    index_list.resize(list_size,2,P);

    int py=(P/px);
    Array<int,1> x;
    Array<int,1> y;
    
    int rankx,ranky;
    int nx=Nx/(2*px),ny=Ny/(2*py);

    
    for (int rank_id=0; rank_id<P; rank_id++){
        get_rank(rank_id, py, rankx, ranky);
        compute_index_list(x, Nx/2, px, rankx);
        compute_index_list(y, Ny/2, py, ranky);
        for (int i=0; i<nx; i++){
            index_list(Range(ny*i,(i+1)*ny-1),0,rank_id)=x(i);
            index_list(Range(ny*i,(i+1)*ny-1),1,rank_id)=y(Range::all());
        }
    }


}



/**
 ********************************************************************************************************************************************
 * \brief   Test function to validate the calculation of structure functions of 3D velocity field data.
 *
 *          This function validates the calculation of the structure functions computed using 3D velocity field data. The velocity field is
 *          generated as \f$ \mathbf{u} = x \hat{x} + y \hat{y} + z \hat{z} \f$. For such field, the velocity structure functions of order
 *          \f$ q \f$ is given as \f$ S_q^u(l_x, l_y, l_z) = (\sqrt{l_x^2 + l_y^2 + l_z^2})^q \f$. In this function, the theoretical values
 *          obtained from the aforementioned equation are compared with the computed values. If the percentage difference between the two values
 *          is less than \f$1 \times 10^{-10} \f$, the test is passed.
 *
 ********************************************************************************************************************************************
 */
void VECTOR_TEST_CASE_3D()
{	
    double epsilon=1e-10;
    double err1 = 0, err2 = 0;
    double max = 0;
	Array<double,3> test1,test2;

	if (longitudinal==true){
		test1.resize(Nx/2,Ny/2,Nz/2);

		for (int order=0 ; order<=q2-q1; order++){
			string name=int_to_str(order+q1);
			read_3D(test1,"out/",SF_Grid_pll_name,SF_Grid_pll_name+name);
			for (int i=0; i<test1.extent(0); i++){
				double lx=dx*i;
				for (int j=0; j<test1.extent(1); j++){
					double ly=dy*j;
					for (int k=0; k<test1.extent(2); k++){
						double lz=dz*k;
                        if (lx*lx + ly*ly + lz*lz > epsilon) {
                            err1 = abs((test1(i,j,k)-pow(lx*lx+ly*ly+lz*lz,(order+q1)/2.))/pow(lx*lx+ly*ly+lz*lz,(order+q1)/2.));
                        }
                        else {
                            err1 = abs(test1(i,j,k));
                        }

                        if (err1 > max) {
                            max = err1;
                        }


					}
				}
			}
		}

	}
	else{
        cout<<"\nTESTING BOTH TRANSVERSE AND LONGITUDINAL\n";
		test1.resize(Nx/2,Ny/2,Nz/2);
		test2.resize(Nx/2,Ny/2,Nz/2);
		for (int order=0 ; order<=q2-q1; order++){
			string name=int_to_str(order+q1);

			read_3D(test1,"out/",SF_Grid_pll_name,SF_Grid_pll_name+name);
			read_3D(test2,"out/",SF_Grid_perp_name,SF_Grid_perp_name+name);


			for (int i=0; i<test1.extent(0); i++){
				double lx=dx*i;
				for (int j=0; j<test1.extent(1); j++){
					double ly=dy*j;
					for (int k=0; k<test1.extent(2); k++){
						double lz=dz*k;
						if (lx*lx + ly*ly + lz*lz > epsilon) {
                            err1 = abs((test1(i,j,k)-pow(lx*lx+ly*ly+lz*lz,(order+q1)/2.))/pow(lx*lx+ly*ly+lz*lz,(order+q1)/2.));
                        }
                        else {
                            err1 = abs(test1(i,j,k));
                        }
                        err2 = abs(test2(i,j,k));
                        if (err1 > max) {
                            max = err1;
                            
                        }

                        if (err2 > max) {
                            max = err2;
                    
                        }

					}
				}
			}

		}
	}


	if (max > epsilon){
		cout<<"\n\nVECTOR_3D: TEST_FAILED. The structure functions computed numerically using the code do NOT match with the analytically obtained values. \n\n";
	}
	else{
		cout<<"\n\nVECTOR_3D: TEST_PASSED. The structure functions computed numerically using the code match with the analytically obtained values. \n\n";
	}

    cout<<"MAXIMUM PERCENTAGE ERROR: "<<max<<endl<<endl;

}


/**
 ********************************************************************************************************************************************
 * \brief   Test function to validate the calculation of structure functions of 2D velocity field data.
 *
 *          This function validates the calculation of the structure functions computed using 2D velocity field data. The velocity field is
 *          generated as \f$ \mathbf{u} = x \hat{x} + z \hat{z} \f$. For such field, the velocity structure functions of order
 *          \f$ q \f$ is given as \f$ S_q^u(l_x, l_z) = (\sqrt{l_x^2 + l_z^2})^q \f$. In this function, the analytically obtained  values
 *          obtained from the aforementioned equation are compared with the computed values. If the percentage difference between the two values is less
            than \f$1 \times 10^{-10} \f$, the test is passed.
 *
 ********************************************************************************************************************************************
 */
void VECTOR_TEST_CASE_2D()
{	
	double epsilon=1e-10;
	double max=0;
	double err1=0, err2 = 0;
    Array<double,2> test1,test2;
	int count=0;

	if (longitudinal==true){
		test1.resize(Nx/2,Nz/2);

		for (int order=0 ; order<=q2-q1; order++){
			string name=int_to_str(order+q1);
			read_2D(test1,"out/",SF_Grid_pll_name, SF_Grid_pll_name+name);
			for (int i=0; i<test1.extent(0); i++){
				double lx=dx*i;
				for (int k=0; k<test1.extent(1); k++){
					double lz=dz*k;
                    if ((lx*lx + lz*lz)>epsilon) {
                        err1 = abs((test1(i,k)-pow(lx*lx+lz*lz,(order+q1)/2.))/pow(lx*lx+lz*lz,(order+q1)/2.));
                    }
                    else {
                        err1 =  abs(test1(i,k));
                    }

                    if (err1 > max) {
                        max=err1;
                    }

				}
			}
		}

	}
	else{
		test1.resize(Nx/2,Nz/2);
		test2.resize(Nx/2,Nz/2);
		for (int order=0 ; order<=q2-q1; order++){
			string name=int_to_str(order+q1);

			read_2D(test1,"out/",SF_Grid_pll_name,SF_Grid_pll_name+name);
			read_2D(test2,"out/",SF_Grid_perp_name,SF_Grid_perp_name+name);


			for (int i=0; i<test1.extent(0); i++){
				double lx=dx*i;

				for (int k=0; k<test1.extent(1); k++){
					double lz=dz*k;
                    if ((lx*lx + lz*lz)>epsilon) {
                        err1 = abs((test1(i,k)-pow(lx*lx+lz*lz,(order+q1)/2.))/pow(lx*lx+lz*lz,(order+q1)/2.));
                    }
                    else {
                        err1 =  abs(test1(i,k));
                    }

                    err2 = abs(test2(i,k));

                    if (err1 > max) {
                        max = err1;
                    }
                    if (err2 > max) {
                        max = err2;
                    }

				}
			}
		}
	}


	if (max > epsilon){
		cout<<"\n\nVECTOR_2D: TEST_FAILED. The structure functions computed numerically using the code do NOT match with the analytically obtained values. \n\n";
	}
	else{
		cout<<"\n\nVECTOR_2D: TEST_PASSED. The structure functions computed numerically using the code match with the analytically obtained values. \n\n";
	}

    cout<<"MAXIMUM ERROR: "<<max<<endl<<endl;

}


/**
 ********************************************************************************************************************************************
 * \brief   Test function to validate the calculation of structure functions of 2D scalar field data.
 *
 *          This function validates the calculation of the structure functions computed using 2D scalar field data. The scalar field is
 *          generated as \f$ \theta = x + z \f$. For such field, the structure functions of order
 *          \f$ q \f$ is given as \f$ S_q^u(l_x, l_z) = (l_x^2 + l_z^2)^q \f$. In this function, the theoretical values
 *          obtained from the aforementioned equation are compared with the computed values. If the percentage difference between the two values is less
 *          than \f$1 \times 10^{-10} \f$, the test is passed.
 *
 ********************************************************************************************************************************************
 */
void SCALAR_TEST_CASE_2D()
{	double epsilon=1e-10;
	double max=0;
	double err=0;
	Array<double,2> test1;
	int count=0;
	test1.resize(Nx/2,Nz/2);
	for (int order=0 ; order<=q2-q1; order++){
		string name=int_to_str(order+q1);

		read_2D(test1,"out/",SF_Grid_scalar_name, SF_Grid_scalar_name+name);



		for (int i=0; i<test1.extent(0); i++){
			double lx=dx*i;

			for (int k=0; k<test1.extent(1); k++){
				double lz=dz*k;
				if (abs(lx+lz)>epsilon){
					err=abs((test1(i,k)-pow(lx+lz,(order+q1)))/pow(lx+lz,(order+q1)));

				}
				else{
					err=abs(test1(i,k));

				}

                if (err>max) {
                    max = err;
                }
                

			}
            
		}
	}
	if (max > epsilon){
		cout<<"\n\nSCALAR_2D: TEST_FAILED. The structure functions computed numerically using the code do NOT match with the analytically obtained values. \n\n";
	}
	else{
		cout<<"\n\nSCALAR_2D: TEST_PASSED. The structure functions computed numerically using the code match with the analytically obtained values. \n\n";
	}

	cout<<"MAXIMUM ERROR: "<<max<<endl<<endl;
}

/**
 ********************************************************************************************************************************************
 * \brief   Test function to validate the calculation of structure functions of 3D scalar field data.
 *
 *          This function validates the calculation of the structure functions computed using 3D scalar field data. The scalar field is
 *          generated as \f$ \theta = x + y + z \f$. For such field, the structure functions of order
 *          \f$ q \f$ is given as \f$ S_q^u(l_x, l_y, l_z) = (l_x^2 + l_y^2 + l_z^2)^q \f$. In this function, the theoretical values
 *          obtained from the aforementioned equation are compared with the computed values. If the percentage difference between the two values is less
 *          than \f$1 \times 10^{-10} \f$, the test is passed.
 *
 ********************************************************************************************************************************************
 */
void SCALAR_TEST_CASE_3D(){

	double epsilon=1e-10;
	double max=0;
	double err=0;
    Array<double,3> test1;
	int count=0;
	test1.resize(Nx/2,Ny/2,Nz/2);
	for (int order=0 ; order<=q2-q1; order++){
		string name=int_to_str(order+q1);
		read_3D(test1,"out/",SF_Grid_scalar_name, SF_Grid_scalar_name+name);
		for (int i=0; i<test1.extent(0); i++){
			double lx=dx*i;
			for (int j=0; j<test1.extent(1); j++){
				double ly=dy*j;
				for (int k=0; k<test1.extent(2); k++){
					double lz=dz*k;
					if (abs(lx+ly+lz)>epsilon){
						err=abs((test1(i,j,k)-pow(lx+ly+lz,(order+q1)))/pow(lx+ly+lz,(order+q1)));

					}
					else{
						err=abs(test1(i,j,k));

					}

                    if (err>max) {
                        max=err;
                    }
				}
			}
		}
	}
	if (max > epsilon){
		cout<<"\n\nSCALAR_3D: TEST_FAILED. The structure functions computed numerically using the code do NOT match with the analytically obtained values. \n\n";
	}
	else{
		cout<<"\n\nSCALAR_3D: TEST_PASSED. The structure functions computed numerically using the code match with the analytically obtained values. \n\n";

	}

    cout<<"MAXIMUM ERROR: "<<max<<endl<<endl;

}



/**
 ********************************************************************************************************************************************
 * \brief   Function to convert an integer type value to string.
 *
 *
 * \param   number is the integer value of to be converted.
 *
 * \return  The value as a string.
 *
 ********************************************************************************************************************************************
 */
string int_to_str(int number)
{
    stringstream ss;

    ss << number;

    return ss.str();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the time elapsed.
 *
 *          This function computes the time elapsed in seconds between a start point and an end point during the execution of the
 *          program
 *
 * \param start_t is the time corresponding to the start point
 * \param end_t is the time corresponding to the end point
 * \param elapsed stores the time elapsed in seconds
 ********************************************************************************************************************************************
 */
void compute_time_elapsed(timeval start_t, timeval end_t, double& elapsed){
    long elapsed_2 = (end_t.tv_sec-start_t.tv_sec)*1000000u + end_t.tv_usec-start_t.tv_usec;
    elapsed=elapsed_2/1.0e6;
}



/**
 ********************************************************************************************************************************************
 * \brief   Function to write the structure functions as function of \f$ (l_x,l_y,l_z) \f$ as a 3D hdf5 file.
 *
 *          This function reads the structure functions as a 4D array of dimensions \f$ (l_x \times l_y \times l_z \times p)\f$, where \f$ p \f$ is the order of the number of orders of
 *          structure function.
 *          The structure functions of different orders are then stored as separate 3D hdf5 files.
 *
 * \param   A is the 4D array representing the structure functions.
 * \param   file is the name of the hdf5 file and the dataset in which the structure functions are stored.
 * \param   q is the order of the structure function to be stored.
 ********************************************************************************************************************************************
 */
void write_4D(Array<double,4> A, string file) {
  int nx=A(Range::all(),0,0,0).size();
  int ny=A(0,Range::all(),0,0).size();
  int nz=A(0,0,Range::all(),0).size();
  h5::File f("out/"+file+".h5", "w");
  int q = q1;
  Array<double,3> temp(nx,ny,nz);
  while(q<=q2){
      cout<<"Writing "<<q<<" order to file.\n";
      string qstr = int_to_str(q);
      h5::Dataset ds = f.create_dataset(file+qstr, h5::shape(nx,ny,nz), "double");
      temp(Range::all(),Range::all(),Range::all())=A(Range::all(),Range::all(),Range::all(),q-q1);
      ds << temp.data();
      q++;
  }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to write the structure functions as function of lx,lz as a 2D hdf5 file.
 *
 *          This function reads the structure functions as a 3D array of dimensions \f$ (l_x \times l_z \times p) \f$, where \f$ p \f$ is the number of orders of the structure function.
 *          The structure functions of different orders are then stored as separate 2D hdf5 files.
 *
 * \param   A is the 3D array representing the structure functions.
 * \param   file is the name of the hdf5 file and the dataset in which the structure functions are stored.
 * \param   q is the order of the structure function to be stored.
 ********************************************************************************************************************************************
 */
void write_3D(Array<double,3> A, string file) {
  int nx=A(Range::all(),0,0).size();
  int nz=A(0,Range::all(),0).size();
  h5::File f("out/"+file+".h5", "w");
  int q = q1;
  Array<double,2> temp(nx,nz);
  while(q<=q2) {
      cout<<"Writing "<<q<<" order to file.\n";
      string qstr = int_to_str(q);
      h5::Dataset ds = f.create_dataset(file+qstr, h5::shape(nx,nz), "double");
      temp(Range::all(),Range::all())=A(Range::all(),Range::all(),q-q1);
      ds << temp.data();
      q++;
  }
}
/**
 ********************************************************************************************************************************************
 * \brief   Function to show the checklist for proper input files
 ********************************************************************************************************************************************
 */
void show_checklist(){
	
  	cerr<<"Error: Please check the following\n\n";
  	cerr<<"a. 'in' folder contains the input files\n\n";
    cerr<<"b. Unless specified otherwise via command-line arguments, the input files should be as follows:\n";
  	cerr<<"\tCase Vector:\n";
 	cerr<<"\t\tCase 2D: U.V1r.h5, U.V3r.h5\n";
 	cerr<<"\t\tCase 3D: U.V1r.h5, U.V2r.h5, U.V3r.h5\n";
 	cerr<<"\tCase Scalar: \n\t\tT.Fr.h5\n";
    cerr<<"NOTE: The dataset name should be the same as the file name (without the .h5 extension).\n\n";
	cerr<<"c. All the input datasets must have the same dimensions\n";
	cerr<<"Please refer to Readme for details\n\n";
    
    


}

/**
 ********************************************************************************************************************************************
 * \brief   Function to read a 2D field from an hdf5 file.
 *
 *          This function reads an hdf5 file containing a 2D field, which can be the \f$ x \f$ or \f$ z \f$ component of a 2D velocity field. The dimensions of the
 *          2D field is \f$(N_x \times N_z)\f$, where \f$N_x\f$ and \f$N_z\f$ are the number of gridpoints in \f$ x \f$ and \f$ z \f$ directions respectively.
 *          The hdf5 file should have only
 *          one dataset, and the names of the hdf5 file and the dataset must be identical. This function makes use of the H5SI library for reading the
 *          hdf5 file.
 *
 * \param   A is the 2D array to store the field that is read from the file.
 * \param   fold is the name of the folder in which the input files are kept.
 * \param   file is a string storing the name of the file to be read.
 ********************************************************************************************************************************************
 */
void read_2D(Array<double,2> A, string fold, string file, string dset) {
  ifstream file_name(fold+file+".h5");
  h5::File f(fold+file+".h5", "r");
  f[dset] >> A.data();
}




/**
 ********************************************************************************************************************************************
 * \brief   Function to read a 3D field from an hdf5 file.
 *
 *          This function reads an hdf5 file containing a 4D field, which can be the \f$ x \f$, \f$ y \f$, or \f$ z \f$ component of a 3D velocity field. The dimensions of the
 *          3D field is \f$(N_x \times N_y \times N_z)\f$, where \f$N_x\f$, \f$N_y\f$, and \f$N_z\f$ are the number of gridpoints in \f$ x \f$, \f$ y \f$, and \f$ z \f$ directions
 *          respectively. The hdf5 file should
 *          have only one dataset, and the names of the hdf5 file and the dataset must be identical. This function makes use of the H5SI library for
 *          reading the hdf5 file.
 *
 * \param A is the 3D array to store the field that is read from the file.
 * \param fold is the name of the folder in which the input files are kept.
 * \param file is a string storing the name of the file to be read.
 * \param dset is the name of the dataset storing the input field.
 ********************************************************************************************************************************************
 */
void read_3D(Array<double,3> A, string fold, string file, string dset) {
	ifstream file_name(fold+file+".h5");
	h5::File f(fold+file+".h5", "r");
	f[dset] >> A.data();
}


/**
***********************************************************************************************************************************
* \brief Function to outputs the run-time instructions for the user
*
***********************************************************************************************************************************
*/


void help_command(){
	if(rank_mpi==0){
		cout<<"### iii) Running Instructions and Command-Line Arguments \n \
		To run `fastSF`, change to `fastSF` directory. Ensure that the input\n\
		hdf5 files follow the schema described in the previous subsection. \n\
		If you want all the relevant parameters to be read from 'in/para.yaml', \n\
		you can simply type the following: `mpirun -np [number of MPI processors] src/fastSF.out`\n\
		`fastSF` allows the user to pass the input parameters using command line arguments as well.\n \
		If inputs are provided via the command-line, the corresponding inputs read from the 'in/para.yaml'\n\
		file get overriden. The user can also specify the input and output hdf5 file names via the \n\
		command-line. The command line arguments are given as follows:\n\n\n\
		`mpirun -np [number of MPI processors] src/fastSF.out -s [scalar_switch]` \n\
		`-d [2D_switch] -l [Only_longitudinal] -p [Processors_X] -X [Nx] -Y [Ny]`\n\
		`-Z [Nz] -x [Lx] -y [Ly] -z [Lz] -1 [q1] -2 [q2] -t [test_switch]`\n\
		`-U [Name of the hdf5 file storing Ux]`\n\
		`-V [Name of the hdf5 file storing Uy]`\n\
		`-W [Name of the hdf5 file storing Uz]`\n\
		`-Q [Name of the hdf5 file storing T]`\n\
		`-u [Name of the dataset storing Ux]`\n\
		`-v [Name of the dataset storing Uy]`\n\
		`-w [Name of the dataset storing Uz]`\n\
		`-q [Name of the dataset storing T]`\n\
		`-P [Name of the hdf5 file storing the transverse structure functions]`\n\
		`-L [Name of the hdf5 file storing the longitudinal structure functions]\
        `-h [Help]`\n\n\n\
		The user need not give all the command line arguments; the arguments that \n\
		are not provided will be read by the `in/para.yaml` file. For, if the user wants \n\
		to run `fastSF` with 16 processors with 4 processors in x direction, and wants to compute\n\
		only the longitudinal structure functions, the following command should be entered:\n\
		`mpirun -np 16 ./src/fastSF.out -p 4 -l true`\n\
		In this case, the number of processors in the x-direction and the longitudinal\n\
		structure function switch will be taken via the command line. The rest of the parameters \n\
		will be taken from the `in/para.yaml` file.";
	}
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to parse the parameters from the yaml file and from the command-line arguments.
 *
 * \param   argc is the number of command-line arguments passed by the user including the name of the program.
 * \param   argv is the array of character pointers listing all the arguments.
 *
 ********************************************************************************************************************************************
 */
void get_Inputs(int argc, char* argv[]) {
    YAML::Node para;
    ifstream para_yaml,input_field;
    string para_path="in/para.yaml";
    para_yaml.open(para_path.c_str());
  
    if (para_yaml.is_open())
    {
      try
      {
        YAML::Parser parser(para_yaml);
        parser.GetNextDocument(para);
      }
      catch(YAML::ParserException& e)
      {
        cerr << "Global::Parse: Error reading parameter file: \n" << e.what() << endl;
      }
  
    }
    else
    {
      cerr << "Global::Parse: Unable to open '" + para_path + "'." << endl;
      h5::finalize();
      MPI_Finalize();
      exit(1);
    }
    para["program"]["scalar_switch"]>>scalar_switch;
    para["program"]["Only_longitudinal"]>>longitudinal;
    para["program"]["2D_switch"]>>two_dimension_switch;
    para["program"]["Processors_X"]>>px;

    para["test"]["test_switch"]>>test_switch;
    
    if (test_switch){
    	para["grid"]["Nx"]>>Nx;
    	para["grid"]["Ny"]>>Ny;
    	para["grid"]["Nz"]>>Nz;
	}
    para["domain_dimension"]["Lx"]>>Lx;
    para["domain_dimension"]["Ly"]>>Ly;
    para["domain_dimension"]["Lz"]>>Lz;
  
    para["structure_function"]["q1"]>>q1;
    para["structure_function"]["q2"]>>q2;
    
  
    int option;
    while ((option=getopt(argc, argv, "X:Y:Z:1:2:x:y:z:l:d:p:t:s:U:V:W:Q:P:L:M:h:u:v:w:q:"))!=-1){
    	switch(option){
    		case 'h':
    			help_command();
    			exit(1);
    			break;
    		case 'X':
    			Nx=std::stoi(optarg);
    			break;
    		case 'Y':
    			Ny=std::stoi(optarg);
    			break;
    		case 'Z':
    			Nz=std::stoi(optarg);
    			break;
    		case 'x':
    			Lx=std::stod(optarg);
    			break;
    		case 'y':
    			Ly=std::stod(optarg);
    			break;
    		case 'z':
    			Lz=std::stod(optarg);
    			break;
    		case 'p':
    			px=std::stod(optarg);
    			break;
    		case '1':
    			q1=std::stod(optarg);
    			break;
    		case '2':
    			q2=std::stod(optarg);
    			break;
    		case 't':
    			test_switch=str_to_bool(optarg);
    			break;
    		case 's':
    			scalar_switch=str_to_bool(optarg);
    			break;
    		case 'd':
    			two_dimension_switch=str_to_bool(optarg);
    			break;
    		case 'l':
    			longitudinal=str_to_bool(optarg);
    			break;
            case 'U':
                UName = optarg;
                break;
            case 'V':
                VName = optarg;
                break;
            case 'W':
                WName = optarg;
                break;
            case 'Q':
                TName = optarg;
                break;
            case 'u':
                UdName = optarg;
                break;
            case 'v':
                VdName = optarg;
                break;
            case 'w':
                WdName = optarg;
                break;
            case 'q':
                TdName = optarg;
                break;
            case 'P':
                SF_Grid_perp_name = optarg;
                break;
            case 'L':
                SF_Grid_pll_name = optarg;
                break;
            case 'M':
                SF_Grid_scalar_name = optarg;
                break;
            default:
                if (rank_mpi==0){
                    cout<<"\nNo command line options given; reading all the inputs from para.yaml.\n";
                }
    	}
    }

  
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to generate an idealized 3D velocity field.
 *
 *          This function generates the following 3D velocity field.
 *          \f$u_x = x, \quad u_y = y, \quad u_z = z\f$.
 *
 * \param Ux is a 3D array representing the x-component of 3D velocity field.
 * \param Uy is a 3D array representing the y-component of 3D velocity field.
 * \param Uz is a 3D array representing the z-component of 3D velocity field.
 ********************************************************************************************************************************************
 */
void Read_Init(Array<double,3>& Ux, Array<double,3>& Uy, Array<double,3>& Uz){
  if (rank_mpi==0)
  {cout<<"\nGenerating the 3D velocity field: U = [x, y, z] \n";
  }
  for (int i=0; i<Nx; i++){
      for (int j=0; j<Ny; j++){
        for (int k=0; k<Nz; k++){
            Ux(i, j, k) = i*dx;
            Uy(i, j, k) = j*dy;
            Uz(i, j, k) = k*dz;
          }
        }
    }
    if (rank_mpi==0)
    {cout<<"\nField has been generated.\n";
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to generate an idealized 2D velocity field.
 *
 *          This function generates the following 2D velocity field.
 *          \f$u_x = x, \quad u_z = z\f$.
 *
 * \param Ux is a 2D array representing the x-component of 2D velocity field.
 * \param Uz is a 2D array representing the z-component of 2D velocity field.
 ********************************************************************************************************************************************
 */
void Read_Init(Array<double,2>& Ux, Array<double,2>& Uz){
	if (rank_mpi==0){
		cout<<"\nGenerating the 2D velocity field: U = [x, z] \n";
	}
    for (int i=0;i<Nx;i++){
      for (int k=0;k<Nz;k++){
          Ux(i, k) = i*dx;
          Uz(i, k) = k*dz;
       }
  }
  if (rank_mpi==0)
    {cout<<"\nField has been generated.\n";
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to generate an idealized 2D scalar field.
 *
 *          This function generates the following 2D scalar field.
 *          \f$\theta = x + z \f$
 *
 * \param T is a 2D array representing the x-component of 2D velocity field.
 ********************************************************************************************************************************************
 */
void Read_Init(Array<double,2>& T) {
	if (rank_mpi==0){
		cout<<"\nGenerating the scalar field: T = x + z \n";
	}
    for (int i=0;i<Nx;i++){
      for (int k=0;k<Nz;k++){
          T(i, k) = i*dx + k*dz;
       }
  }
  if (rank_mpi==0)
    {cout<<"\nField has been generated.\n";
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to generate an idealized 3D scalar field.
 *
 *          This function generates the following 2D scalar field.
 *          \f$\theta = x + y + z \f$
 *
 * \param T is a 3D array representing the x-component of 2D velocity field.
 ********************************************************************************************************************************************
 */
void Read_Init(Array<double,3>& T) {
	if (rank_mpi==0){
		cout<<"\nGenerating the scalar field: T = x + y + z \n";
	}
    for (int i=0;i<Nx;i++){
      for (int j=0;j<Ny;j++){
          for (int k=0;k<Nz;k++){
              T(i, j, k) = i*dx + j*dy + k*dz;
          }
      }
  }
  if (rank_mpi==0)
    {cout<<"\nField has been generated.\n";
    }
}



/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the longitudinal and transverse structure functions for a 3D velocity field.
 *
 *
 * \param Ux is a 3D array representing the x-component of velocity field
 * \param Uy is a 3D array representing the y-component of velocity field
 * \param Uz is a 3D array representing the z-component of velocity field
 ********************************************************************************************************************************************
 */
void SFunc3D(
        Array<double,3> Ux,
        Array<double,3> Uy,
        Array<double,3> Uz)
{
	if (rank_mpi==0) {
        cout<<"\nComputing longitudinal and transverse S(lx, ly, lz) using 3D velocity field data..\n";
    }
    int c_per_proc = Nx*Ny/(4*P);

    Array<int, 3> index_list;
    compute_index_list(index_list, Nx, Ny);
    Array<double,3> dUx;
    Array<double,3> dUy;
    Array<double,3> dUz;
    Array<double,3> dUpll;
    
    for (int ix=0; ix<c_per_proc; ix++){
        int x=index_list(ix, 0, rank_mpi);
        int y=index_list(ix, 1, rank_mpi);
  		for(int z=0; z<Nz/2; z++){
            
  			dUx.resize(Nx-x,Ny-y,Nz-z);
  			dUy.resize(Nx-x,Ny-y,Nz-z);
  			dUz.resize(Nx-x,Ny-y,Nz-z);
            dUpll.resize(Nx-x,Ny-y,Nz-z);
        		
        	int count=(Nx-x)*(Ny-y)*(Nz-z);
        	double lx=x*dx;
        	double ly=y*dy;
        	double lz=z*dz;

        	double r=sqrt(lx*lx+ly*ly+lz*lz);
        	dUx(Range::all(),Range::all(),Range::all())=Ux(Range(x,Nx-1),Range(y,Ny-1),Range(z,Nz-1))-Ux(Range(0,Nx-x-1),Range(0,Ny-y-1),Range(0,Nz-z-1));
        	dUy(Range::all(),Range::all(),Range::all())=Uy(Range(x,Nx-1),Range(y,Ny-1),Range(z,Nz-1))-Uy(Range(0,Nx-x-1),Range(0,Ny-y-1),Range(0,Nz-z-1));
        	dUz(Range::all(),Range::all(),Range::all())=Uz(Range(x,Nx-1),Range(y,Ny-1),Range(z,Nz-1))-Uz(Range(0,Nx-x-1),Range(0,Ny-y-1),Range(0,Nz-z-1));
        		
        	
            dUpll=(lx*dUx+ly*dUy+lz*dUz)/r;
        
        	dUx=dUx-dUpll*lx/r;
        	dUy=dUy-dUpll*ly/r;
        	dUz=dUz-dUpll*lz/r;

            
        	dUx=pow(dUx*dUx+dUy*dUy+dUz*dUz,0.5);
            
        	for (int p=0; p<=q2-q1; p++){
        		double Spll = sum(pow(dUpll(Range::all(),Range::all(),Range::all()),q1+p))/(count);
        		double Sperp =sum(pow(dUx(Range::all(),Range::all(),Range::all()),q1+p))/(count);
                Array<int, 1> X, Y, Z, p_arr;
                Array<double, 1> Spll_arr, Sperp_arr;
                
                if (rank_mpi==0) {
                    X.resize(P);
                    Y.resize(P);
                    Z.resize(P);
                    p_arr.resize(P);
                    Spll_arr.resize(P);
                    Sperp_arr.resize(P);
                }
        
                MPI_Gather(&x, 1, MPI_INT, X.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&y, 1, MPI_INT, Y.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&z, 1, MPI_INT, Z.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&Spll, 1, MPI_DOUBLE, Spll_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Gather(&p, 1, MPI_INT, p_arr.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&Sperp, 1, MPI_DOUBLE, Sperp_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                if (rank_mpi==0) {
                    for (int i=0; i<P; i++) {
                        SF_Grid_pll(X(i), Y(i), Z(i), p_arr(i)) = Spll_arr(i);
                        SF_Grid_perp(X(i), Y(i), Z(i), p_arr(i)) = Sperp_arr(i);
                    } 
                }

        	}
  		}
  	}
    if (rank_mpi==0) {
        SF_Grid_pll(0,0,0,Range::all())=0;
        SF_Grid_perp(0,0,0,Range::all())=0;
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate only the longitudinal structure functions for a 3D velocity field.
 *
 * \param Ux is a 3D array representing the x-component of velocity field
 * \param Uy is a 3D array representing the y-component of velocity field
 * \param Uz is a 3D array representing the z-component of velocity field
 ********************************************************************************************************************************************
 */
void SFunc_long_3D(
        Array<double,3> Ux,
        Array<double,3> Uy,
        Array<double,3> Uz)
{
if (rank_mpi==0) {
        cout<<"\nComputing longitudinal S(lx, ly, lz) using 3D velocity field data..\n";
    }
    
    int c_per_proc = Nx*Ny/(4*P);

    Array<int, 3> index_list;
    compute_index_list(index_list, Nx, Ny);
    Array<double,3> dUx;
    Array<double,3> dUy;
    Array<double,3> dUz;
    Array<double,3> dUpll;
    
    for (int ix=0; ix<c_per_proc; ix++){
        int x=index_list(ix, 0, rank_mpi);
        int y=index_list(ix, 1, rank_mpi);
  		for(int z=0; z<Nz/2; z++){
  			dUx.resize(Nx-x,Ny-y,Nz-z);
			dUy.resize(Nx-x,Ny-y,Nz-z);
  			dUz.resize(Nx-x,Ny-y,Nz-z);
            dUpll.resize(Nx-x,Ny-y,Nz-z);
        		
    		int count=(Nx-x)*(Ny-y)*(Nz-z);
    		double lx=x*dx;
    		double ly=y*dy;
    		double lz=z*dz;

    		double r=sqrt(lx*lx+ly*ly+lz*lz);
    		dUx(Range::all(),Range::all(),Range::all())=Ux(Range(x,Nx-1),Range(y,Ny-1),Range(z,Nz-1))-Ux(Range(0,Nx-x-1),Range(0,Ny-y-1),Range(0,Nz-z-1));
    		dUy(Range::all(),Range::all(),Range::all())=Uy(Range(x,Nx-1),Range(y,Ny-1),Range(z,Nz-1))-Uy(Range(0,Nx-x-1),Range(0,Ny-y-1),Range(0,Nz-z-1));
    		dUz(Range::all(),Range::all(),Range::all())=Uz(Range(x,Nx-1),Range(y,Ny-1),Range(z,Nz-1))-Uz(Range(0,Nx-x-1),Range(0,Ny-y-1),Range(0,Nz-z-1));

            dUpll=(lx*dUx+ly*dUy+lz*dUz)/r;
        		

    		for (int p=0; p<=q2-q1; p++){
    			double Spll = sum(pow(dUpll(Range::all(),Range::all(),Range::all()),q1+p))/(count);
                Array<int, 1> X, Y, Z, p_arr;
                Array<double, 1> Spll_arr;
                
                if (rank_mpi==0) {
                    X.resize(P);
                    Y.resize(P);
                    Z.resize(P);
                    p_arr.resize(P);
                    Spll_arr.resize(P);
                }
        
                MPI_Gather(&x, 1, MPI_INT, X.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&y, 1, MPI_INT, Y.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&z, 1, MPI_INT, Z.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&Spll, 1, MPI_DOUBLE, Spll_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Gather(&p, 1, MPI_INT, p_arr.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

                if (rank_mpi==0) {
                    for (int i=0; i<P; i++) {
                        SF_Grid_pll(X(i), Y(i), Z(i), p_arr(i)) = Spll_arr(i);
                    }
                }

    		}
  		}
  	}
    if (rank_mpi==0) {
  	    SF_Grid_pll(0,0,0,Range::all())=0;
    }
}




/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate the longitudinal and transverse structure functions for a 2D velocity field.
 *
 *        
 * \param Ux is a 2D array representing the x-component of velocity field
 * \param Uz is a 2D array representing the z-component of velocity field
 ********************************************************************************************************************************************
 */
 void SFunc2D(
         Array<double,2> Ux,
         Array<double,2> Uz)
 {
     if (rank_mpi==0) {
         cout<<"\nComputing longitudinal and transverse S(lx, lz) using 2D velocity field data..\n";
     }

    int p_per_proc = Nx*Nz/(4*P);

    Array<int, 3> index_list;
    compute_index_list(index_list, Nx, Nz);
    Array<double,2> dUz;
    Array<double,2> dUx;
    Array<double,2> dUpll;
    
    for (int ix=0; ix<p_per_proc; ix++){
        int x=index_list(ix, 0, rank_mpi);
        int z=index_list(ix, 1, rank_mpi);
        dUx.resize(Nx-x,Nz-z);
        dUz.resize(Nx-x,Nz-z);
        dUpll.resize(Nx-x,Nz-z);
        int count=(Nx-x)*(Nz-z);
        double lx=x*dx;
        double lz=z*dz;
        double r=sqrt(lx*lx+lz*lz);

        dUx(Range::all(),Range::all())=Ux(Range(x,Nx-1),Range(z,Nz-1))-Ux(Range(0,Nx-x-1),Range(0,Nz-z-1));
        dUz(Range::all(),Range::all())=Uz(Range(x,Nx-1),Range(z,Nz-1))-Uz(Range(0,Nx-x-1),Range(0,Nz-z-1));
        	
        dUpll=(lx*dUx+lz*dUz)/r;
        dUx=dUx-dUpll*lx/r;
        dUz=dUz-dUpll*lz/r;
        dUx=pow(dUx*dUx+dUz*dUz,0.5);

    	for (int p=0; p<=q2-q1; p++){
            double Spll = sum(pow(dUpll(Range::all(),Range::all()),q1+p))/(count);
            double Sperp =sum(pow(dUx(Range::all(),Range::all()),q1+p))/(count);
            Array<int, 1> X, Z, p_arr;
            Array<double, 1> Spll_arr, Sperp_arr;
            
            if (rank_mpi==0) {
                X.resize(P);
                Z.resize(P);
                p_arr.resize(P);
                Spll_arr.resize(P);
                Sperp_arr.resize(P);
            }
        
            MPI_Gather(&x, 1, MPI_INT, X.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&z, 1, MPI_INT, Z.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&Spll, 1, MPI_DOUBLE, Spll_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(&p, 1, MPI_INT, p_arr.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&Sperp, 1, MPI_DOUBLE, Sperp_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            if (rank_mpi==0) {
                for (int i=0; i<P; i++) {
                    SF_Grid2D_pll(X(i), Z(i), p_arr(i)) = Spll_arr(i);
                    SF_Grid2D_perp(X(i), Z(i), p_arr(i)) = Sperp_arr(i);
                } 
            } 
        } 
    }
    if (rank_mpi==0) {
        SF_Grid2D_pll(0,0,Range::all())=0;
        SF_Grid2D_perp(0,0,Range::all())=0;
    }
    
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate only the longitudinal structure functions for a 2D velocity field.
 *
 *         
 * \param Ux is a 2D array representing the x-component of velocity field
 * \param Uz is a 2D array representing the z-component of velocity field
 ********************************************************************************************************************************************
 */
void SFunc_long_2D(
         Array<double,2> Ux,
         Array<double,2> Uz)
 {
     if (rank_mpi==0) {
         cout<<"\nComputing longitudinal S(lx, lz) using 2D velocity field data..\n";
     }

    int p_per_proc = Nx*Nz/(4*P);

    Array<int, 3> index_list;
    compute_index_list(index_list, Nx, Nz);
    Array<double,2> dUz;
    Array<double,2> dUx;
    Array<double,2> dUpll;
    
    for (int ix=0; ix<p_per_proc; ix++){
        int x=index_list(ix, 0, rank_mpi);
        int z=index_list(ix, 1, rank_mpi);
        dUx.resize(Nx-x,Nz-z);
        dUz.resize(Nx-x,Nz-z);
        dUpll.resize(Nx-x,Nz-z);
        int count=(Nx-x)*(Nz-z);
        double lx=x*dx;
        double lz=z*dz;
        double r=sqrt(lx*lx+lz*lz);

        dUx(Range::all(),Range::all())=Ux(Range(x,Nx-1),Range(z,Nz-1))-Ux(Range(0,Nx-x-1),Range(0,Nz-z-1));
        dUz(Range::all(),Range::all())=Uz(Range(x,Nx-1),Range(z,Nz-1))-Uz(Range(0,Nx-x-1),Range(0,Nz-z-1));
            
        dUpll=(lx*dUx+lz*dUz)/r;

        for (int p=0; p<=q2-q1; p++){
            double Spll = sum(pow(dUpll(Range::all(),Range::all()),q1+p))/(count);
            Array<int, 1> X, Z, p_arr;
            Array<double, 1> Spll_arr;
            
            if (rank_mpi==0) {
                X.resize(P);
                Z.resize(P);
                p_arr.resize(P);
                Spll_arr.resize(P);
            }
        
            MPI_Gather(&x, 1, MPI_INT, X.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&z, 1, MPI_INT, Z.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&Spll, 1, MPI_DOUBLE, Spll_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(&p, 1, MPI_INT, p_arr.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (rank_mpi==0) {
                for (int i=0; i<P; i++) {
                    SF_Grid2D_pll(X(i), Z(i), p_arr(i)) = Spll_arr(i);
                } 
            }
        }
    }
    if (rank_mpi==0) {
        SF_Grid2D_pll(0,0,Range::all())=0;
    }
    
}



/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions for a 3D scalar field.
 *
 * \param T is a 3D array representing the scalar field
 ********************************************************************************************************************************************
 */


void SF_scalar_3D(
         Array<double,3> T)
 {
     if (rank_mpi==0) {
         cout<<"\nComputing S(lx, ly, lz) using 3D scalar field data..\n";
     }
    
    int c_per_proc = Nx*Ny/(4*P);

    Array<int, 3> index_list;
    compute_index_list(index_list, Nx, Ny);
    Array<double,3> dT;
    
    for (int ix=0; ix<c_per_proc; ix++){
        int x=index_list(ix, 0, rank_mpi);
        int y=index_list(ix, 1, rank_mpi);
        for(int z=0; z<Nz/2; z++){
            dT.resize(Nx-x,Ny-y,Nz-z);
                
            int count=(Nx-x)*(Ny-y)*(Nz-z);
            double r=sqrt(x*x*dx*dx+y*y*dy*dy+z*z*dz*dz);

            dT(Range::all(),Range::all(),Range::all())=T(Range(x,Nx-1),Range(y,Ny-1),Range(z,Nz-1))-T(Range(0,Nx-x-1),Range(0,Ny-y-1),Range(0,Nz-z-1));
                

            for (int p=0; p<=q2-q1; p++){
                double St = sum(pow(dT(Range::all(),Range::all(),Range::all()),q1+p))/(count);
                Array<int, 1> X, Y, Z, p_arr;
                Array<double, 1> St_arr;
                
                if (rank_mpi==0) {
                    X.resize(P);
                    Y.resize(P);
                    Z.resize(P);
                    p_arr.resize(P);
                    St_arr.resize(P);
                }
        
                MPI_Gather(&x, 1, MPI_INT, X.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&y, 1, MPI_INT, Y.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&z, 1, MPI_INT, Z.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
                MPI_Gather(&St, 1, MPI_DOUBLE, St_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Gather(&p, 1, MPI_INT, p_arr.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

                if (rank_mpi==0) {
                    for (int i=0; i<P; i++) {
                        SF_Grid_scalar(X(i), Y(i), Z(i), p_arr(i)) = St_arr(i);
                    }
                }

            }
                
        }
    }
    if (rank_mpi==0) {
        SF_Grid_scalar(0,0,0,Range::all())=0;
    }


 }

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions for a 2D scalar field.
 *
 *
 * \param T is a 2D array representing the scalar field
 ********************************************************************************************************************************************
 */
void SF_scalar_2D(Array<double,2> T)
 {
     if (rank_mpi==0) {
         cout<<"\nComputing S(lx, lz) using 2D scalar field data..\n";
     }

    int p_per_proc = Nx*Nz/(4*P);

    Array<int, 3> index_list;
    compute_index_list(index_list, Nx, Nz);
    Array<double,2> dT;
    
    for (int ix=0; ix<p_per_proc; ix++){
        int x=index_list(ix, 0, rank_mpi);
        int z=index_list(ix, 1, rank_mpi);
       			
        dT.resize(Nx-x,Nz-z);
        int count=(Nx-x)*(Nz-z);

        dT(Range::all(),Range::all())=T(Range(x,Nx-1),Range(z,Nz-1))-T(Range(0,Nx-x-1),Range(0,Nz-z-1));
        		
        for (int p=0; p<=q2-q1; p++){
            double St = sum(pow(dT(Range::all(),Range::all()),q1+p))/(count);
            Array<int, 1> X, Z, p_arr;
            Array<double, 1> St_arr;
            
            if (rank_mpi==0) {
                X.resize(P);
                Z.resize(P);
                p_arr.resize(P);
                St_arr.resize(P);
            }
        
            MPI_Gather(&x, 1, MPI_INT, X.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&z, 1, MPI_INT, Z.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&St, 1, MPI_DOUBLE, St_arr.data(), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(&p, 1, MPI_INT, p_arr.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (rank_mpi==0) {
                for (int i=0; i<P; i++) {
                    SF_Grid2D_scalar(X(i), Z(i), p_arr(i)) = St_arr(i);
                } 
            }
        }
    }
    if (rank_mpi==0) {
        SF_Grid2D_scalar(0,0,Range::all())=0;
    }
 }

