/********************************************************************************************************************************************
 * Kolmogorov41
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

/*! \file Kolmogorov41.cc
 *
 *  \brief Code to compute structure functions using velocity or scalar field data.
 * 
 *  \author Shashwat Bhattacharya, Shubhadeep Sadhukhan
 *  \date Jan 2020
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
void Read_para(); //Input from the user
void write_2D(Array<double,2>, string);
void write_3D(Array<double,3>, string, int);
void write_4D(Array<double,4>, string, int);
void read_2D(Array<double,2>, string);
string int_to_str(int);

void compute_time_elapsed(timeval, timeval, double&);
double powInt(double, int);
void read_3D(Array<double,3>, string);
void SFunc2D(Array<double,2>, Array<double,2>, Array<double,2>&, Array<double,2>&, Array<double,2>&);
void SFunc2D(Array<double,2>, Array<double,2>, Array<double,2>&, Array<double,2>&, Array<double,2>&, Array<double,3>&, Array<double,3>&, Array<double,3>&);
void SFunc_long_2D(Array<double,2>, Array<double,2>, Array<double,2>&, Array<double,2>&);
void SFunc_long_2D(Array<double,2>, Array<double,2>, Array<double,2>&, Array<double,2>&, Array<double,3>&, Array<double,3>&);
void SFunc3D(Array<double,3>, Array<double,3>, Array<double,3>, Array<double,2>&, Array<double,2>&, Array<double,2>&);
void SFunc3D(Array<double,3>, Array<double,3>, Array<double,3>, Array<double,2>&, Array<double,2>&, Array<double,2>&, Array<double,4>&, Array<double,4>&, Array<double,4>&);
void SFunc_long_3D(Array<double,3>, Array<double,3>, Array<double,3>, Array<double,2>&, Array<double,2>&);
void SFunc_long_3D(Array<double,3>, Array<double,3>, Array<double,3>, Array<double,2>&, Array<double,2>&, Array<double,4>&, Array<double,4>&);
void Read_Init(Array<double,2>&, Array<double,2>&);
void Read_Init(Array<double,3>&, Array<double,3>&, Array<double,3>&);
void Read_Init(Array<double,2>&);
void Read_Init(Array<double,3>&);
double powInt(double, int);
void SF_scalar_3D(Array<double,3>, Array<double,2>&, Array<double,2>&);
void SF_scalar_3D(Array<double,3>, Array<double,2>&, Array<double,2>&, Array<double,4>&, Array<double,4>&);
void SF_scalar_2D(Array<double,2>, Array<double,2>&, Array<double,2>&);
void SF_scalar_2D(Array<double,2>, Array<double,2>&, Array<double,2>&, Array<double,3>&, Array<double,3>&);


complex<double> I=1i;
/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the scalar field.   
 ********************************************************************************************************************************************
 */
Array <double,3> T;
/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the x-component of the 3D velocity field.   
 ********************************************************************************************************************************************
 */
Array <double,3> V1;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the y-component of the 3D velocity field.   
 ********************************************************************************************************************************************
 */
Array <double,3> V2;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the z-component of the 3D velocity field.   
 ********************************************************************************************************************************************
 */
Array <double,3> V3;

/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the 2D scalar field.   
 ********************************************************************************************************************************************
 */
Array <double,2> T_2D;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the x-component of the 3D velocity field.   
 ********************************************************************************************************************************************
 */
Array<double,2> V1_2D;

/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the z-component of the 2D velocity field.   
 ********************************************************************************************************************************************
 */ 
Array<double,2> V3_2D;


/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the computed longitudinal structure functions or scalar structure functions and their corresponding orders.   
 ********************************************************************************************************************************************
 */ 
Array<double,2> SF;

/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the computed transverse structure functions and their corresponding orders.   
 ********************************************************************************************************************************************
 */  
Array<double,2> SF_perp;

/**
 ********************************************************************************************************************************************
 * \brief   2D array storing the values to divide SF and SF_perp for averaging.   
 ********************************************************************************************************************************************
 */  
Array<double,2> counter;

/**
 ********************************************************************************************************************************************
 * \brief   4D array storing the computed longitudinal structure functions as function of (lx, ly, lz, p), where p is the order.   
 ********************************************************************************************************************************************
 */
Array<double,4> SF_Grid_pll;

/**
 ********************************************************************************************************************************************
 * \brief   4D array storing the computed transverse structure functions as function of (lx, ly, lz, p), where p is the order.   
 ********************************************************************************************************************************************
 */
Array<double,4> SF_Grid_perp;

/**
 ********************************************************************************************************************************************
 * \brief   4D array storing the counter array for averaging SF_Grid_pll and SF_Grid_perp
 ********************************************************************************************************************************************
 */
Array<double,4> counter_Grid;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the computed longitudinal structure functions as function of (lx, lz, p), where p is the order.   
 ********************************************************************************************************************************************
 */
Array<double,3> SF_Grid2D_pll;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the computed transverse structure functions as function of (lx, ly, lz, p), where p is the order.   
 ********************************************************************************************************************************************
 */
Array<double,3> SF_Grid2D_perp;

/**
 ********************************************************************************************************************************************
 * \brief   3D array storing the counter array for averaging SF_Grid_pll and SF_Grid_perp
 ********************************************************************************************************************************************
 */
Array<double,3> counter_Grid2D;

/**
 ********************************************************************************************************************************************
 * \brief   This variable decides whether the structure functions are to be calculated using 2D or 3D velocity field data.
 * 
 * If true, then the code will read 2D velocity fields and calculate the corresponding structure functions. Otherwise, it will read 3D velocity
 * fields. Entered by the user  
 ********************************************************************************************************************************************
 */  
bool two_dimension_switch;

/**
 ********************************************************************************************************************************************
 * \brief    This variable decides whether the scalar or velocity structure functions are to be evaluated. 
 *
 * If the value is TRUE, then scalar structure functions will be evaluated, else the vector (velocity) structure functions will be evaluated.
 ********************************************************************************************************************************************
 */
bool scalar_switch;

/**
 ********************************************************************************************************************************************
 * \brief    This variable decides whether the structure functions as function of (lx, ly, lz), or of (lx, lz) in case of 2D, needs to be
 *           calculated in addition to the structure functions as function of (l).
 *
 * If the value is TRUE, then the structure functions as function of (lx, ly, lz) / (lx, lz) will be computed in addtion to the structure 
 * functions as function of (l).
 ********************************************************************************************************************************************
 */
bool grid_switch;

/**
 ********************************************************************************************************************************************
 * \brief   Number of gridpoints in the x direction.    
 ********************************************************************************************************************************************
 */  
int Nx;

/**
 ********************************************************************************************************************************************
 * \brief   Number of gridpoints in the y direction.   
 ********************************************************************************************************************************************
 */   
int Ny;

/**
 ********************************************************************************************************************************************
 * \brief   Number of gridpoints in the z direction.    
 ********************************************************************************************************************************************
 */   
int Nz; 

/**
 ********************************************************************************************************************************************
 * \brief   Number of OpenMP threads.
 ********************************************************************************************************************************************
 */ 
int num; 

/**
 ********************************************************************************************************************************************
 * \brief   Total number of gridpoints in the diagonal direction.    
 ********************************************************************************************************************************************
 */ 
int Nr; 

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

//int q; 
/**
 ********************************************************************************************************************************************
 * \brief   This variable decides whether the code needs to generate velocity field data or read velocity field data from HDF5 files. 
 * 
 * If the value is 0, the code generates velocity field data. Else, it reads from the HDF5 files in the "in" folder.  Entered by the user. 
 ********************************************************************************************************************************************
 */
int field_procedure; 

/**
 ********************************************************************************************************************************************
 * \brief   This variable decides whether to calculate both transverse and longitudinal structure functions or only the longitudinal structure  
 * functions using a less computationally expensive technique. 
 * 
 * If the value is false, the code calculates both transverse and longitudinal structure functions. Else, it calculates only the longitudinal
 * structure functions using less number of iterations.  
 ********************************************************************************************************************************************
 */
bool longitudinal;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the distance between two consecutive gridpoints in the x direction.
 * 
 ********************************************************************************************************************************************
 */
double dx;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the distance between two consecutive gridpoints in the y direction.
 * 
 ********************************************************************************************************************************************
 */
double dy; 

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the distance between two consecutive gridpoints in the z direction.
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
double lx;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the breadth of the domain.
 * 
 ********************************************************************************************************************************************
 */ 
double ly;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the height of the domain.
 * 
 ********************************************************************************************************************************************
 */ 
double lz;

/**
 ********************************************************************************************************************************************
 * \brief   This variable stores the number of the MPI process.
 * 
 ********************************************************************************************************************************************
 */ 
int P;


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$//


/**
 ********************************************************************************************************************************************
 * \brief   The main function of the "Strunc" code.
 *
 *          This function is the main function of the "Strunc" code for computing the velocity and scalar structure functions. The MPI 
 *          decomposition and integration are also carried out in this function.
 *     
 * 
 ********************************************************************************************************************************************
 */
int main(int argc, char *argv[]) {
  num=omp_get_max_threads();//strtol(argv[1],NULL,10);  //omp_get_max_threads();//(shaheen)
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_mpi);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  timeval start_pt, end_pt, start_t, end_t;
  gettimeofday(&start_t,NULL);
  double elapsedt=0.0;
  double elapsepdt=0.0;
  
  Read_para();

  Nr = (int)ceil(sqrt(pow(Nx-1,2)+pow(Ny-1,2)+pow(Nz-1,2)))+1;
  //Resizing the input fields
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

  //Defining the input fields
  if (field_procedure==0){
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
        
  else {
      if (rank_mpi==0){
          cout<<"Reading from the hdf5 files\n";
      }
      if (two_dimension_switch){
          if (scalar_switch) {
              read_2D(T_2D, "T.Fr");
          }
          else {
              read_2D(V1_2D,"U.V1r");
              read_2D(V3_2D,"U.V3r");
          }
      }
      else{
          if (scalar_switch) {
              read_3D(T, "T.Fr");
          }
          else {if (rank_mpi==0) cerr<<"\nFLAG 0\n";
              read_3D(V1, "U.V1r");
              read_3D(V2, "U.V2r");
              read_3D(V3, "U.V3r");
          }
      }
  }
   

  
  SF.resize(Nr,q2-q1+1);
  Array<double, 2> SF_Node(Nr,q2-q1+1);
 
  
  Array<double, 2> SF_Node_perp;
  if (not scalar_switch){
      if (not longitudinal) {
          SF_perp.resize(Nr,q2-q1+1);
          SF_Node_perp.resize(Nr,q2-q1+1);
          SF_perp=0;
          SF_Node_perp=0;
      }
  }

  counter.resize(Nr,q2-q1+1);
  Array<double, 2> counter_Node(Nr,q2-q1+1);

  SF = 0;
  counter = 0;
  SF_Node = 0;
  counter_Node = 0;
  Array<double, 4> SF_Grid_pll_Node;
  Array<double, 4> SF_Grid_perp_Node;
  Array<double, 4> counter_Grid_Node;
  Array<double, 3> SF_Grid2D_pll_Node;
  Array<double, 3> SF_Grid2D_perp_Node;
  Array<double, 3> counter_Grid2D_Node;




  if (not grid_switch) {
      //Initiallizing h5si
      h5::init();

      gettimeofday(&start_pt,NULL);
      
      //Calculating the structure functions
      if (two_dimension_switch){
          
          if (scalar_switch) {
              SF_scalar_2D(T_2D, SF_Node, counter_Node);
          }

          else {
              
              if (longitudinal) {
                  SFunc_long_2D(V1_2D, V3_2D, SF_Node, counter_Node);
              }
              else {
                  SFunc2D(V1_2D, V3_2D, SF_Node, SF_Node_perp, counter_Node);
              }
          }
          
      }
      else{
          
          if (scalar_switch) {
              SF_scalar_3D(T, SF_Node, counter_Node);
          }

          else {
              if (longitudinal) {
        	      SFunc_long_3D(V1, V2, V3, SF_Node, counter_Node);
              }
              else {
        	      SFunc3D(V1, V2, V3, SF_Node, SF_Node_perp, counter_Node);
              }
          }
      }
      gettimeofday(&end_pt,NULL);
  }
  else {

      //Declaring the nodal variables
 
      
      if (not two_dimension_switch) {
          if (scalar_switch) {
              SF_Grid_pll_Node.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
              SF_Grid_pll.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
              counter_Grid_Node.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
              counter_Grid.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
          }
          else {
              if (longitudinal) {
                  SF_Grid_pll_Node.resize(2*Nx-1, 2*Ny-1, Nz, q2-q1+1);
                  SF_Grid_pll.resize(2*Nx-1, 2*Ny-1, Nz, q2-q1+1);
                  counter_Grid_Node.resize(2*Nx-1, 2*Ny-1, Nz, q2-q1+1);
                  counter_Grid.resize(2*Nx-1, 2*Ny-1, Nz, q2-q1+1);
              }
              else {
                  SF_Grid_pll_Node.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
                  SF_Grid_pll.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
                  SF_Grid_perp_Node.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
                  SF_Grid_perp.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
                  counter_Grid_Node.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
                  counter_Grid.resize(2*Nx-1, 2*Ny-1, 2*Nz-1, q2-q1+1);
                  SF_Grid_perp = 0;
                  SF_Grid_perp_Node = 0;
              }
          }

          SF_Grid_pll = 0;
          counter_Grid = 0;
          SF_Grid_pll_Node = 0;
          counter_Grid_Node = 0;
      }
      else {
          if (scalar_switch) {
              SF_Grid2D_pll_Node.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
              SF_Grid2D_pll.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
              counter_Grid2D_Node.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
              counter_Grid2D.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
          }
          else {
              if (longitudinal) {
                  SF_Grid2D_pll_Node.resize(2*Nx-1, Nz, q2-q1+1);
                  SF_Grid2D_pll.resize(2*Nx-1, Nz, q2-q1+1);
                  counter_Grid2D_Node.resize(2*Nx-1, Nz, q2-q1+1);
                  counter_Grid2D.resize(2*Nx-1, Nz, q2-q1+1);
              }
              else {
                  SF_Grid2D_pll_Node.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
                  SF_Grid2D_pll.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
                  SF_Grid2D_perp_Node.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
                  SF_Grid2D_perp.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
                  counter_Grid2D_Node.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
                  counter_Grid2D.resize(2*Nx-1, 2*Nz-1, q2-q1+1);
                  SF_Grid2D_perp = 0;
                  SF_Grid2D_perp_Node = 0;
              }
          }

          SF_Grid2D_pll = 0;
          counter_Grid2D = 0;
          SF_Grid2D_pll_Node = 0;
          counter_Grid2D_Node = 0;
      }


      //Initiallizing h5si
      h5::init();

      gettimeofday(&start_pt,NULL);
      
      //Calculating the structure functions
      if (two_dimension_switch){
          
          if (scalar_switch) {
              SF_scalar_2D(T_2D, SF_Node, counter_Node, SF_Grid2D_pll_Node, counter_Grid2D_Node);
          }

          else {
              if (longitudinal) {
                  SFunc_long_2D(V1_2D, V3_2D, SF_Node, counter_Node, SF_Grid2D_pll_Node, counter_Grid2D_Node);
              }
              else {
                  SFunc2D(V1_2D, V3_2D, SF_Node, SF_Node_perp, counter_Node, SF_Grid2D_pll_Node, SF_Grid2D_perp_Node, counter_Grid2D_Node);
              }
          }
          
      }
      else{
          
          if (scalar_switch) {
              SF_scalar_3D(T, SF_Node, counter_Node, SF_Grid_pll_Node, counter_Grid_Node);
          }

          else {
              if (longitudinal) {
        	      SFunc_long_3D(V1, V2, V3, SF_Node, counter_Node, SF_Grid_pll_Node, counter_Grid_Node);
              }
              else {
        	      SFunc3D(V1, V2, V3, SF_Node, SF_Node_perp, counter_Node, SF_Grid_pll_Node, SF_Grid_perp_Node, counter_Grid_Node);
              }
          }

      }

      gettimeofday(&end_pt,NULL);
  }

  MPI_Reduce(SF_Node.data(), SF.data(), Nr*(q2-q1+1), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
  if (not scalar_switch) {
      if (not longitudinal) {
          MPI_Reduce(SF_Node_perp.data(), SF_perp.data(), Nr*(q2-q1+1), MPI_DOUBLE_PRECISION,MPI_SUM, 0, MPI_COMM_WORLD);
      }
  }
  MPI_Reduce(counter_Node.data(),counter.data(),Nr*(q2-q1+1),MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD);

  
  //Reducing SF Grid
  if (grid_switch) {
      int totPoints;
      if (two_dimension_switch) {
          totPoints = counter_Grid2D_Node.extent(0)*counter_Grid2D_Node.extent(1)*counter_Grid2D_Node.extent(2); 
          MPI_Reduce(SF_Grid2D_pll_Node.data(), SF_Grid2D_pll.data(), totPoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
          MPI_Reduce(counter_Grid2D_Node.data(), counter_Grid2D.data(), totPoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
          if (not scalar_switch) {
              if (not longitudinal) {
                  MPI_Reduce(SF_Grid2D_perp_Node.data(), SF_Grid2D_perp.data(), totPoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
              }
          }
      }

      else {
          totPoints = counter_Grid_Node.extent(0)*counter_Grid_Node.extent(1)*counter_Grid_Node.extent(2)*counter_Grid_Node.extent(3); 
          MPI_Reduce(SF_Grid_pll_Node.data(), SF_Grid_pll.data(), totPoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
          MPI_Reduce(counter_Grid_Node.data(), counter_Grid.data(), totPoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
          if (not scalar_switch) {
              if (not longitudinal) {
                  MPI_Reduce(SF_Grid_perp_Node.data(), SF_Grid_perp.data(), totPoints, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);
              }
          }
      }

  }
   
  gettimeofday(&end_t,NULL);
  compute_time_elapsed(start_t, end_t, elapsedt);
  compute_time_elapsed(start_pt, end_pt, elapsepdt);
  
  if (rank_mpi==0){
      cout<<"\nTime elapsed for parallel part : "<<elapsepdt<<" sec\n";
            
      SF/=counter;
      if (not scalar_switch) {
          if (not longitudinal) {
              SF_perp/=counter;
          }
      }

      if (grid_switch) {
          if (two_dimension_switch) {
              SF_Grid2D_pll /= counter_Grid2D;
              if (not scalar_switch) {
                  if (not longitudinal) {
                      SF_Grid2D_perp /= counter_Grid2D;
                  }
              }
          }
          else {
              SF_Grid_pll /= counter_Grid;
              if (not scalar_switch) {
                  if (not longitudinal) {
                      SF_Grid_perp /= counter_Grid;
                  }
              }
          }
      }
            
      mkdir("out",0777);
            
      cout<<"\nWriting SF as function of l\n";
      write_2D(SF,"SF");
      if (not scalar_switch) {
          if (not longitudinal) {
              write_2D(SF_perp,"SF_perp");
          }
      }
      cout<<"\nWriting completed\n";

      if (grid_switch) {
          int p1 = q1;
          while (p1 <= q2) {
              string name = int_to_str(p1);
              if (two_dimension_switch) {
                  cout<<"\nWriting "<<p1<<" order SF has function of lx and lz\n";
                  write_3D(SF_Grid2D_pll,"SF_Grid_pll"+name, p1);
                  if (not scalar_switch) {
                      if (not longitudinal) {
                          write_3D(SF_Grid2D_perp, "SF_Grid_perp"+name, p1);
                      }
                  }
                  cout<<"\nWriting completed\n";
              }
              else {
                  cout<<"\nWriting "<<p1<<" order SF has function of lx, ly, and ly\n";
                  write_4D(SF_Grid_pll,"SF_Grid_pll"+name, p1);
                  if (not scalar_switch) {
                      if (not longitudinal) {
                          write_4D(SF_Grid_perp, "SF_Grid_perp"+name, p1);
                      }
                  }
                  cout<<"\nWriting completed\n";
              }
              p1++;
          }
      }
                
      cout<<"Total time elapsed : "<<elapsedt<<" sec\n";
            
      cout<<"\nProgram ends\n";
  }

    
  h5::finalize();
  MPI_Finalize();
  return 0;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function which conducts exponentiation with an integer as an exponent. 
 *
 *          This function calculates x raised to the power n, where n is an integer. This function is faster than the standard pow(x,n) function
 *          for n>2. Note that this function cannot accept a non-integer exponent.
 * 
 * \param   x is the base of double-precision floating point datatype.
 * \param   n is the exponent of integer datatype
 * 
 * \return  The value as a string.
 * 
 ********************************************************************************************************************************************
 */
double powInt(double x, int n) {
	double p = 1;
	if (n>0) {
		for (int i=1;i<=n;i++) {
			p = p*x;
		}
	}
	else if (n<0) {
		for (int i=-1;i>=n;i--) {
			p = p/x;
		}
	}
	else {
		p = 1;
	}
	return p;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to convert an integer type value to string.
 *
 *          This function converts an integer type value to a string.
 * 
 * \param   num is the integer value of to be converted.
 * 
 * \return  The value as a string.
 * 
 ********************************************************************************************************************************************
 */
string int_to_str(int num)
{
    stringstream ss;

    ss << num;

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
 * \brief   Function to write the structure functions as function of l as a 2D hdf5 file.
 *
 *          This function reads writes the structure functions as a 2D array of dimensions (N x p). Here, N = \f$\sqrt{Nx^2 + Ny^2 + Nz^2}\f$ is the 
 *          number of equidistant points from 0 to L, where L is the length of the diagonal of the domain in which the structure functions are 
 *          calculated. p is the order of the structure functions. 
 * 
 * \param   A is the 2D array to store the structure functions.
 * \param   file is the name of the hdf5 file and the dataset in which the structure functions are stored.
 ********************************************************************************************************************************************
 */
void write_2D(Array<double,2> A,string file) {
  int np=A(Range::all(),0).size();
  int nr=A(0,Range::all()).size();
  h5::File f("out/"+file+".h5", "w");
  h5::Dataset ds = f.create_dataset(file, h5::shape(np,nr), "double");
  ds << A.data();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to write the structure functions as function of lx,ly,lz as a 3D hdf5 file.
 *
 *          This function reads the structure functions as a 4D array of dimensions (lx x ly x lz x np), where np is the order of the structure function.
 *          The structure functions of different orders are then stored as separate 3D hdf5 files.
 * 
 * \param   A is the 4D array representing the structure functions.
 * \param   file is the name of the hdf5 file and the dataset in which the structure functions are stored.
 * \param   q is the order of the structure function to be stored.
 ********************************************************************************************************************************************
 */
void write_4D(Array<double,4> A, string file,int q) {
  int nx=A(Range::all(),0,0,0).size();
  int ny=A(0,Range::all(),0,0).size();
  int nz=A(0,0,Range::all(),0).size();
  Array<double,3> temp(nx,ny,nz);
  temp(Range::all(),Range::all(),Range::all())=(A(Range::all(),Range::all(),Range::all(),q-q1));
  h5::File f("out/"+file+".h5", "w");
  h5::Dataset ds = f.create_dataset(file, h5::shape(nx,ny,nz), "double");
  ds << temp.data();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to write the structure functions as function of lx,lz as a 2D hdf5 file.
 *
 *          This function reads the structure functions as a 4D array of dimensions (lx x lz x np), where np is the order of the structure function.
 *          The structure functions of different orders are then stored as separate 2D hdf5 files.
 * 
 * \param   A is the 3D array representing the structure functions.
 * \param   file is the name of the hdf5 file and the dataset in which the structure functions are stored.
 * \param   q is the order of the structure function to be stored.
 ********************************************************************************************************************************************
 */
void write_3D(Array<double,3> A, string file,int q) {
  int nx=A(Range::all(),0,0).size();
  int nz=A(0,Range::all(),0).size();
  Array<double,2> temp(nx,nz);
  temp(Range::all(),Range::all())=(A(Range::all(),Range::all(),q-q1));
  h5::File f("out/"+file+".h5", "w");
  h5::Dataset ds = f.create_dataset(file, h5::shape(nx,nz), "double");
  ds << temp.data();
}



/**
 ********************************************************************************************************************************************
 * \brief   Function to read a 2D field from an hdf5 file.
 *
 *          This function reads an hdf5 file containing a 2D field, which can be the x or z component of a 2D velocity field. The dimensions of the
 *          2D field is \f$(N_x \times N_z)\f$, where \f$N_x\f$ and \f$N_z\f$ are the number of gridpoints in x and z directions respectively. 
 *          The hdf5 file should have only 
 *          one dataset, and the names of the hdf5 file and the dataset must be identical. This function make use of the H5SI library for reading the 
 *          hdf5 file.
 * 
 * \param   A is the 2D array to store the field that is read from the file.
 * \param   file is a string storing the name of the file to be read.
 ********************************************************************************************************************************************
 */
void read_2D(Array<double,2> A,string file) {
  h5::File f("in/"+file+".h5", "r");
  f[file] >> A.data();
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to read a 3D field from an hdf5 file.
 *          
 *          This function reads an hdf5 file containing a 4D field, which can be the x y, or z component of a 3D velocity field. The dimensions of the
 *          3D field is \f$(N_x \times N_y \times N_z)\f$, where \f$N_x\f$, \f$N_y\f$, and \f$N_z\f$ are the number of gridpoints in x, y, and z directions 
 *          respectively. The hdf5 file should 
 *          have only one dataset, and the names of the hdf5 file and the dataset must be identical. This function make use of the H5SI library for 
 *          reading the hdf5 file.
 * 
 * \param A is the 3D array to store the field that is read from the file.
 * \param file is a string storing the name of the file to be read.
 ********************************************************************************************************************************************
 */
void read_3D(Array<double,3> A,string file) {
  h5::File f("in/"+file+".h5", "r");
  f[file] >> A.data();
  if (rank_mpi==0) cout << A << endl;
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to open the yaml file and parse the parameters.
 *
 *          The function opens the parameters.yaml file and parses the simulation parameters into its member variables that are publicly
 *          accessible.
 ********************************************************************************************************************************************
 */
void Read_para() {
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
    exit(1);
  }
  para["program"]["field_procedure"]>>field_procedure;
  para["program"]["scalar_switch"]>>scalar_switch;
  para["program"]["Only_longitudinal"]>>longitudinal;
  para["program"]["field_procedure"]>>field_procedure;
  para["program"]["grid_switch"]>>grid_switch;
  para["program"]["2D_switch"]>>two_dimension_switch;
  para["grid"]["Nx"]>>Nx;
  para["grid"]["Ny"]>>Ny;
  para["grid"]["Nz"]>>Nz;
  para["domain_dimension"]["lx"]>>lx;
  para["domain_dimension"]["ly"]>>ly;
  para["domain_dimension"]["lz"]>>lz;
  
  para["structure_function"]["q1"]>>q1;
  para["structure_function"]["q2"]>>q2;
  if (Nx==1){dx=0;}
  else{
    dx=lx/double(Nx-1);}
  if (Ny==1){dy=0;}
  else{
    dy=ly/double(Ny-1);}
  if (Nz==1){dz=0;}
  else{
    dz=lz/double(Nz-1);}

}


/**
 ********************************************************************************************************************************************
 * \brief   Function to assign an exponential function to the 3D velocity field.
 *
 *          This function assigns the following exponential function to the 3D velocity field.
 *          \f$U_x = x, \quad U_y = y, \quad U_z = z\f$.
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
 * \brief   Function to assign an exponential function to the 2D velocity field.
 *
 *          This function assigns the following exponential function to the 2D velocity field.
 *          \f$U_x = x, \quad U_z = z\f$.
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
 * \brief   Function to assign an exponential function to a 2D scalar field.
 *
 *          This function assigns the following exponential function to the scalar field.
 *          \f$T = x + z \f$
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
 * \brief   Function to assign an exponential function to a 3D scalar field.
 *
 *          This function assigns the following exponential function to the scalar field.
 *          \f$T = x + y + z \f$
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
 * \brief   Function to compute the magnitude of a 3D vector A.
 * 
 * \param A is a tiny vector representing the 3D velocity field.
 * \param mag is the variable that stores the magnitude calculated in this function.
 ********************************************************************************************************************************************
 */
void magnitude(TinyVector<double,3> A,double& mag){
    mag=sqrt(A(0)*A(0)+A(1)*A(1)+A(2)*A(2));
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to compute the magnitude of a 2D vector A.
 * 
 * \param A is a tiny vector representing the 2D velocity field.
 * \param mag is the variable that stores the magnitude calculated in this function.
 ********************************************************************************************************************************************
 */
void magnitude(TinyVector<double,2> A, double& mag){
    mag = sqrt(A(0)*A(0)+A(1)*A(1));
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 3D velocity field.
 *
 *          The following function computes both the nodal longitudinal and transverse structure functions of 3D velocity field using six nested for-loops.
 *          The outer three for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner three loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param Ux is a 3D array representing the x-component of velocity field
 * \param Uy is a 3D array representing the y-component of velocity field
 * \param Uz is a 3D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions for a range of orders specified by the user.
 * \param SF_Node_p is a 2D array containing the values of the nodal transverse structure functions for orders q1 to q2.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node and SF_Node_p so as to get the average. 
 ********************************************************************************************************************************************
 */
void SFunc3D(Array<double,3> Ux,
             Array<double,3> Uy,
             Array<double,3> Uz,
             Array<double,2>& SF_Node,
             Array<double,2>& SF_Node_p, 
             Array<double,2>& counter_Node){
    if (rank_mpi==0) {
        cout<<"\nComputing the longitudinal and transverse S(l) using 3D velocity field data..\n";
    }

    double S=0;
    double upll=0;
    double uperp=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);
    
    for (int x1=starPt; x1<endPt; x1++){
# pragma omp parallel for num_threads(num) private(upll,uperp) shared(x1)
        for(int y1=0; y1<Ny; y1++){
            Array<double,2> local_SF(Nr,q2-q1+1),local_SF_p(Nr,q2-q1+1);
            Array<double,2> local_counter(Nr,q2-q1+1);
            local_SF=0;
            local_SF_p=0;
            local_counter=0;
            
            for (int z1=0;z1<Nz;z1++){
                for (int x2=0; x2<Nx; x2++) {
                    for(int y2=0; y2<Ny; y2++){     
                        for(int z2=0;z2<Nz;z2++){
                            int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((y2-y1),2)+pow((z2-z1),2)));
                            int r=(int)Lmag;
                            double rmag;
                            TinyVector<double,3> u1(Ux(x1,y1,z1),Uy(x1,y1,z1),Uz(x1,y1,z1)),u2(Ux(x2,y2,z2),Uy(x2,y2,z2),Uz(x2,y2,z2)),rv(dx*(x2 - x1),dy*(y2 - y1),dz*(z2 - z1)),uperpv;
                            magnitude(rv,rmag);
                            if (r==0){
                                upll=0;
                                uperp=0;
                            }
                            else{
                                upll=dot(u2-u1,rv)/rmag;
                                uperpv=(u2-u1)-rv*(upll)/(rmag);
                                magnitude(uperpv,uperp);
                            }
                            for (int p=0;p<=q2-q1;p++){
                                local_SF(Lmag,p)+= powInt(upll,q1+p); 
                                local_SF_p(Lmag,p)+= powInt(uperp,q1+p);
                                local_counter(Lmag,p)+=1;
                            }
                        }
                    }
                } 
            }
# pragma omp critical
          {
              SF_Node+=local_SF;
              SF_Node_p+=local_SF_p;
              counter_Node+=local_counter;
          } 
        } 
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   A less computationally intensive function to calculate only the longitudinal structure functions using 3D velocity field.
 *
 *          The following function computes the only longitudinal structure function using 3D velocity field data. This function is less computationally 
 *          expensive compared to SFunc3D. The function exploits the fact that \f$ \langle du(l)^p \rangle = \langle du(-l)^p \rangle \f$, where \f$ du(l) \f$ 
 *          is the component parallel to l. Thus, although this function uses six nested loops, the innermost loop starts from z1 instead of 0, where z1 is 
 *          the iteration number of the third for-loop. The rest of the structure is similar to the function SFunc3D. 
 * 
 * \param Ux is a 3D array representing the x-component of velocity field
 * \param Uy is a 3D array representing the y-component of velocity field
 * \param Uz is a 3D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 ********************************************************************************************************************************************
 */ 
void SFunc_long_3D(Array<double,3> Ux,
                   Array<double,3> Uy,
                   Array<double,3> Uz,
                   Array<double,2>& SF_Node,
                   Array<double,2>& counter_Node){
    if (rank_mpi==0) {
        cout<<"\nComputing only longitudinal S(l) using 3D velocity field data..\n";
    }
    double S=0;
    double upll=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);
    
    for (int x1=starPt; x1<endPt; x1++){
# pragma omp parallel for num_threads(num) private(upll) shared(x1)
        for(int y1=0; y1<Ny; y1++){
            Array<double,2> local_SF(Nr,q2-q1+1);
            Array<double,2> local_counter(Nr,q2-q1+1);
            local_SF=0;
            local_counter=0;
            
            for (int z1=0;z1<Nz;z1++){
                for (int x2=0; x2<Nx; x2++) {
                    for(int y2=0; y2<Ny; y2++){     
                        for(int z2=z1;z2<Nz;z2++){
                            int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((y2-y1),2)+pow((z2-z1),2)));
                            int r=(int)Lmag;
                            double rmag;
                            TinyVector<double,3> u1(Ux(x1,y1,z1),Uy(x1,y1,z1),Uz(x1,y1,z1)),u2(Ux(x2,y2,z2),Uy(x2,y2,z2),Uz(x2,y2,z2)),rv(dx*(x2 - x1),dy*(y2 - y1),dz*(z2 - z1));
                            magnitude(rv,rmag);
                            if (r==0){
                                upll=0;
                            }
                            else{
                                upll=dot(u2-u1,rv)/rmag;
                            }
                            for (int p=0;p<=q2-q1;p++){
                                local_SF(Lmag,p)+= powInt(upll,q1+p); 
                                local_counter(Lmag,p)+=1;
                            }
                        }
                    }
                } 
            }
# pragma omp critical
          {
              SF_Node+=local_SF;
              counter_Node+=local_counter;
          } 
        } 
    }
}


/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 2D velocity field.
 *
 *          The following function computes both the nodal longitudinal and transverse structure functions of 2D velocity field using four nested for-loops.
 *          The outer two for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner two loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param Ux is a 2D array representing the x-component of velocity field
 * \param Uz is a 2D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions for a range of orders specified by the user.
 * \param SF_Node_p is a 2D array containing the values of the nodal transverse structure functions for orders q1 to q2.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node and SF_Node_p so as to get the average. 
 ********************************************************************************************************************************************
 */
void SFunc2D(Array<double,2> Ux,
             Array<double,2> Uz, 
             Array<double,2>& SF_Node, 
             Array<double,2>& SF_Node_p, 
             Array<double,2>& counter_Node)
{
	if (rank_mpi==0) {
        cout<<"\nComputing the longitudinal and transverse S(l) using 2D velocity field data..\n";
    }
    double S=0;
    double upll=0;
    double uperp=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);
  
    for (int x1 = starPt; x1 < endPt; x1++){  
# pragma omp parallel for num_threads(num) private(upll,uperp) shared(x1)
        for(int z1=0; z1<Nz; z1++){
            Array<double,2> local_SF(Nr,q2-q1+1),local_SF_p(Nr,q2-q1+1);
            Array<double,2> local_counter(Nr,q2-q1+1);
            local_SF=0;
            local_counter=0;
            for (int x2=0; x2<Nx; x2++) {
                for(int z2=0; z2<Nz; z2++){
                    int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((z2-z1),2)));
                    int r=(int)Lmag;
                    double rmag;
                    TinyVector<double,2> u1(Ux(x1, z1), Uz(x1, z1)), u2(Ux(x2, z2), Uz(x2, z2)), rv(dx*(x2 - x1), dz*(z2 - z1)),uperpv;
                    magnitude(rv,rmag);
                    if (r==0){
                        upll=0;
                        uperp=0;
                    }
                
                    else{
                        upll=dot(u2-u1, rv)/rmag;
                        uperpv=(u2-u1) - rv*(upll)/(rmag);
                        magnitude(uperpv, uperp);
                    }
                    
                    for (int p=0;p<=q2-q1;p++){
                        local_SF(Lmag,p)+= powInt(upll,q1+p);
                        local_SF_p(Lmag,p)+= powInt(uperp,q1+p);
                        local_counter(Lmag,p)+=1;
                    }
                }
            }
# pragma omp critical
            {
                SF_Node+=local_SF;
                SF_Node_p+=local_SF_p;
                counter_Node+=local_counter;
            }   
        } 
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   A less computationally intensive function to calculate only the longitudinal structure functions using 2D velocity field.
 *
 *          The following function computes the only longitudinal structure function using 2D velocity field data. This function is less computationally 
 *          expensive compared to SFunc3D. The function exploits the fact that \f$ \langle du(l)^p \rangle = \langle du(-l)^p \rangle \f$, where \f$ du(l) \f$ 
 *          is the component parallel to l. Thus, although this function uses four nested loops, the innermost loop starts from z1 instead of 0, where z1 is 
 *          the iteration number of the second for-loop. The rest of the structure is similar to the function SFunc2D. 
 * 
 * \param Ux is a 3D array representing the x-component of velocity field
 * \param Uz is a 3D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 ********************************************************************************************************************************************
 */ 
void SFunc_long_2D(
        Array<double,2> Ux, 
        Array<double,2> Uz, 
        Array<double,2>& SF_Node,
        Array<double,2>& counter_Node) 
{
  
    if (rank_mpi==0) {
        cout<<"\nComputing only longitudinal S(l) using 2D velocity field data..\n";
    }
    double S=0;
    double upll=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);

  for (int x1=starPt; x1<endPt; x1++){  
      
     
     # pragma omp parallel for num_threads(num) private(upll)
      for(int z1=0; z1<Nz; z1++){
        Array<double,2> local_SF(Nr,q2-q1+1);
        Array<double,2> local_counter(Nr,q2-q1+1);
        local_SF=0;
        local_counter=0;
        for (int x2=0;x2<Nx;x2++){
          for (int z2=z1; z2<Nz; z2++) {
                  int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((z2-z1),2)));
                  int r=(int)Lmag;
                  double rmag;
                  TinyVector<double,2> u1(Ux(x1,z1),Uz(x1,z1)), 
                                       u2(Ux(x2,z2),Uz(x2,z2)), 
                                       rv(dx*(x2 - x1), dz*(z2 - z1));
                  magnitude(rv,rmag);
                  if (r==0){
                    upll=0;
                  }
                  else{
                      upll=dot(u2-u1,rv)/rmag;
                  }
                  for (int p=0;p<=q2-q1;p++){
                    local_SF(Lmag,p)+= powInt(upll,q1+p);
                    local_counter(Lmag,p)+=1;
                    
                  }
            } 
          }
          # pragma omp critical
          {
            SF_Node += local_SF;
            counter_Node += local_counter;
          }   
      }
  }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 3D scalar field.
 *
 *          The following function computes both the nodal structure functions of 3D scalar field using six nested for-loops.
 *          The outer three for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner three loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param T is a 3D array representing the scalar field 
 * \param SF_Node is a 2D array containing the values of the nodal structure functions for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 ********************************************************************************************************************************************
 */
void SF_scalar_3D(Array<double,3> T, 
                  Array<double,2>& SF_Node, 
                  Array<double,2>& counter_Node){
    if (rank_mpi==0) {
        cout<<"\nComputing S(l) using 3D scalar field data..\n";
    }
    double S=0;
    double dT =0;
 
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);
    
    for (int x1=starPt; x1<endPt; x1++){
# pragma omp parallel for num_threads(num) private(dT) shared(x1)
        for(int y1=0; y1<Ny; y1++){
            Array<double,2> local_SF(Nr,q2-q1+1);
            Array<double,2> local_counter(Nr,q2-q1+1);
            local_SF=0;
            local_counter=0;
            
            for (int z1=0;z1<Nz;z1++){
                for (int x2=0; x2<Nx; x2++) {
                    for(int y2=0; y2<Ny; y2++){     
                        for(int z2=0;z2<Nz;z2++){
                            int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((y2-y1),2)+pow((z2-z1),2)));
                            int r=(int)Lmag;
                            if (r==0){
                                dT=0;
                            }
                            else{
                                dT = T(x2, y2, z2) - T(x1, y1, z1);
                            }
                            for (int p=0;p<=q2-q1;p++){
                                local_SF(Lmag,p)+= powInt(dT,q1+p); 
                                local_counter(Lmag,p)+=1;
                            }
                        }
                    }
                } 
            }
# pragma omp critical
          {
              SF_Node+=local_SF;
              counter_Node+=local_counter;
          } 
        } 
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 2D scalar field.
 *
 *          The following function computes both the nodal structure functions of 2D scalar field using four nested for-loops.
 *          The outer two for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner two loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param T is a 2D array representing the scalar field 
 * \param SF_Node is a 2D array containing the values of the nodal structure functions for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 ********************************************************************************************************************************************
 */
void SF_scalar_2D(Array<double,2> T, 
                  Array<double,2>& SF_Node, 
                  Array<double,2>& counter_Node) {
    if (rank_mpi==0) {
        cout<<"\nComputing S(l) using 2D scalar field data..\n";
    }
    double S=0;
    double dT =0;
 
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);
    
    for (int x1=starPt; x1<endPt; x1++){
	# pragma omp parallel for num_threads(num) private(dT) shared(x1)
        for(int z1=0; z1<Nz; z1++){
            Array<double,2> local_SF(Nr,q2-q1+1);
            Array<double,2> local_counter(Nr,q2-q1+1);
            local_SF=0;
            local_counter=0;
            for (int x2=0; x2<Nx; x2++) {
                for(int z2=0;z2<Nz;z2++){
					int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((z2-z1),2)));
                    int r=(int)Lmag;
                    if (r==0){
						dT=0;
					}
                            
                    else{
						dT = T(x2, z2) - T(x1, z1);
					 }   
                            
                    for (int p=0;p<=q2-q1;p++){
						local_SF(Lmag,p)+= powInt(dT,q1+p);    
                        local_counter(Lmag,p)+=1;
					}
				}
			}                 
# pragma omp critical
          {
              SF_Node+=local_SF;
              counter_Node+=local_counter;
          } 
        } 
    }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 3D velocity field as function of (lx, ly, lz) in addition to the structure functions
 *          as function of (l).
 *
 *          The following function computes both the nodal longitudinal and transverse structure functions of 3D velocity field using six nested for-loops.
 *          The outer three for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner three loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param Ux is a 3D array representing the x-component of velocity field
 * \param Uy is a 3D array representing the y-component of velocity field
 * \param Uz is a 3D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions as function of (l) for a range of orders specified by the user.
 * \param SF_p_Node is a 2D array containing the values of the nodal transverse structure functions for as function of (l) for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node and SF_Node_p so as to get the average. 
 * \param SF_Grid_pll_Node is a 4D array containing the values of the nodal longitudinal structure functions as function of (lx,ly,lz) for a range of orders specified by the user.
 * \param SF_Grid_perp_Node is a 4D array containing the values of the nodal transverse structure functions for as function of (lx,ly,lz) for a range of orders specified by the user.
 * \param counter_Grid_Node is a 4D array containing the numbers for dividing the values of SF_Grid_pll_Node and SF_Grid_perp_Node_p so as to get the average. 
 ********************************************************************************************************************************************
 */
void SFunc3D(
        Array<double,3> Ux, 
        Array<double,3> Uy, 
        Array<double,3> Uz, 
        Array<double,2>& SF_Node,
        Array<double,2>& SF_p_Node,
        Array<double,2>& counter_Node, 
        Array<double,4>& SF_Grid_pll_Node, 
        Array<double,4>& SF_Grid_perp_Node, 
        Array<double,4>& counter_Grid_Node
        )
{
    if (rank_mpi==0) {
        cout<<"\nComputing longitudinal and transverse S(l) and S(lx, ly, lz) using 3D velocity field data..\n";
    }
    double S=0;
    double upll=0;
    double uperp=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);

  for (int x1=starPt; x1<endPt; x1++){  
      
     
     # pragma omp parallel for num_threads(num) private(upll,uperp)
      for(int y1=0; y1<Ny; y1++){
        Array<double,2> local_SF(Nr,q2-q1+1);
        Array<double,2> local_SF_p(Nr,q2-q1+1);
        Array<double,2> local_counter(Nr,q2-q1+1);
        Array<double,4> local_SF_Grid_pll(2*Nx-1,2*Ny-1,2*Nz-1,q2-q1+1);
        Array<double,4> local_SF_Grid_perp(2*Nx-1,2*Ny-1,2*Nz-1,q2-q1+1);
        Array<double,4> local_counter_Grid(2*Nx-1,2*Ny-1,2*Nz-1,q2-q1+1);
        local_SF=0;
        local_SF_p=0;
        local_counter=0;
        local_counter_Grid=0;
        local_SF_Grid_pll=0;
        local_SF_Grid_perp=0;
        for (int z1=0;z1<Nz;z1++){
          for (int x2=0; x2<Nx; x2++) {
              for(int y2=0; y2<Ny; y2++){     
                for(int z2=0;z2<Nz;z2++){
                  int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((y2-y1),2)+pow((z2-z1),2)));
                  int r=(int)Lmag;
                  double rmag;
                  TinyVector<double,3> u1(Ux(x1,y1,z1),Uy(x1,y1,z1),Uz(x1,y1,z1)),
                                       u2(Ux(x2,y2,z2),Uy(x2,y2,z2),Uz(x2,y2,z2)),
                                       rv(dx*(x2 - x1),dy*(y2 - y1),dz*(z2 - z1)),
                                       uperpv;
                  magnitude(rv,rmag);
                  if (r==0){
                    upll=0;
                    uperp=0;
                  }
                  else{
                      upll=dot(u2-u1,rv)/rmag;
                      uperpv=(u2-u1)-rv*(upll)/(rmag);
                      magnitude(uperpv,uperp);
                  }
                  for (int p=0;p<=q2-q1;p++){
                    local_SF(Lmag,p) += powInt(upll,q1+p);
                    local_SF_p(Lmag,p) += powInt(uperp,q1+p);
                    local_counter(Lmag,p) += 1;
                    local_SF_Grid_pll(Nx+(x2-x1-1), Ny+(y2-y1-1), Nz+(z2-z1),p) += pow(upll,q1+p);
                    local_SF_Grid_perp(Nx+(x2-x1-1), Ny+(y2-y1-1), Nz+(z2-z1),p) += pow(uperp,q1+p);
                    local_counter_Grid(Nx+(x2-x1-1), Ny+(y2-y1-1), Nz+(z2-z1),p) += 1;
                    
                  }
                }
              }
            } 
          }
          # pragma omp critical
          {
            SF_Node += local_SF;
            SF_p_Node += local_SF_p;
            counter_Node += local_counter;
            SF_Grid_pll_Node += local_SF_Grid_pll;
            SF_Grid_perp_Node += local_SF_Grid_perp;
            counter_Grid_Node += local_counter_Grid;
          
          }   
      }
  }
}

/**
 ********************************************************************************************************************************************
 * \brief   A less computationally intensive function to calculate only the longitudinal structure functions as functions of (lx,ly,lz) in
 *          addition to function of (l) using 3D velocity field.
 *
 *          The following function computes the only longitudinal structure function using 3D velocity field data. This function is less computationally 
 *          expensive compared to SFunc3D. The function exploits the fact that \f$ \langle du(l)^p \rangle = \langle du(-l)^p \rangle \f$, where \f$ du(l) \f$ 
 *          is the component parallel to l. Thus, although this function uses six nested loops, the innermost loop starts from z1 instead of 0, where z1 is 
 *          the iteration number of the third for-loop. The rest of the structure is similar to the function SFunc3D. 
 * 
 * \param Ux is a 3D array representing the x-component of velocity field
 * \param Uy is a 3D array representing the y-component of velocity field
 * \param Uz is a 3D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 * \param SF_Grid_pll_Node is a 4D array containing the values of the nodal longitudinal structure functions as function of (lx,ly,lz) for a range of orders specified by the user.
 * \param counter_Grid_Node is a 4D array containing the numbers for dividing the values of SF_Grid_pll_Node and SF_Grid_perp_Node_p so as to get the average. 
 ********************************************************************************************************************************************
 */
void SFunc_long_3D(
        Array<double,3> Ux, 
        Array<double,3> Uy, 
        Array<double,3> Uz, 
        Array<double,2>& SF_Node,
        Array<double,2>& counter_Node, 
        Array<double,4>& SF_Grid_pll_Node, 
        Array<double,4>& counter_Grid_Node
        )
{
  
    if (rank_mpi==0) {
        cout<<"\nComputing only longitudinal S(l) and S(lx, ly, lz) using 3D velocity field data..\n";
    }
    double S=0;
    double upll=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);

  for (int x1=starPt; x1<endPt; x1++){  
      
     
     # pragma omp parallel for num_threads(num) private(upll)
      for(int y1=0; y1<Ny; y1++){
        Array<double,2> local_SF(Nr,q2-q1+1);
        Array<double,2> local_counter(Nr,q2-q1+1);
        Array<double,4> local_SF_Grid_pll(2*Nx-1,2*Ny-1,Nz,q2-q1+1);
        Array<double,4> local_counter_Grid(2*Nx-1,2*Ny-1,Nz,q2-q1+1);
        local_SF=0;
        local_counter=0;
        local_counter_Grid=0;
        local_SF_Grid_pll=0;
        for (int z1=0;z1<Nz;z1++){
          for (int x2=0; x2<Nx; x2++) {
              for(int y2=0; y2<Ny; y2++){     
                for(int z2=z1;z2<Nz;z2++){
                  int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((y2-y1),2)+pow((z2-z1),2)));
                  int r=(int)Lmag;
                  double rmag;
                  TinyVector<double,3> u1(Ux(x1,y1,z1),Uy(x1,y1,z1),Uz(x1,y1,z1)),
                                       u2(Ux(x2,y2,z2),Uy(x2,y2,z2),Uz(x2,y2,z2)),
                                       rv(dx*(x2 - x1),dy*(y2 - y1),dz*(z2 - z1));
                  magnitude(rv,rmag);
                  if (r==0){
                    upll=0;
                  }
                  else{
                      upll=dot(u2-u1,rv)/rmag;
                  }
                  for (int p=0;p<=q2-q1;p++){
                    local_SF(Lmag,p) += powInt(upll,q1+p);
                    local_counter(Lmag,p) += 1;
                    local_SF_Grid_pll(Nx+(x2-x1-1), Ny+(y2-y1-1), (z2-z1),p) += pow(upll,q1+p);
                    local_counter_Grid(Nx+(x2-x1-1), Ny+(y2-y1-1), (z2-z1),p) += 1;
                    
                  }
                }
              }
            } 
          }
          # pragma omp critical
          {
            SF_Node += local_SF;
            counter_Node += local_counter;
            SF_Grid_pll_Node += local_SF_Grid_pll;
            counter_Grid_Node += local_counter_Grid;
          
          }   
      }
  }
}




/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 2D velocity field as function of (lx, lz) in addition to the structure functions
 *          as function of (l).
 *
 *          The following function computes both the nodal longitudinal and transverse structure functions of 2D velocity field using six nested for-loops.
 *          The outer three for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner three loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param Ux is a 2D array representing the x-component of velocity field
 * \param Uz is a 2D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions as function of (l) for a range of orders specified by the user.
 * \param SF_p_Node is a 2D array containing the values of the nodal transverse structure functions for as function of (l) for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node and SF_Node_p so as to get the average. 
 * \param SF_Grid2D_pll_Node is a 3D array containing the values of the nodal longitudinal structure functions as function of (lx,lz) for a range of orders specified by the user.
 * \param SF_Grid2D_perp_Node is a 3D array containing the values of the nodal transverse structure functions for as function of (lx,lz) for a range of orders specified by the user.
 * \param counter_Grid2D_Node is a 3D array containing the numbers for dividing the values of SF_Grid_pll_Node and SF_Grid_perp_Node_p so as to get the average. 
 ********************************************************************************************************************************************
 */
void SFunc2D(
        Array<double,2> Ux, 
        Array<double,2> Uz, 
        Array<double,2>& SF_Node,
        Array<double,2>& SF_p_Node,
        Array<double,2>& counter_Node, 
        Array<double,3>& SF_Grid2D_pll_Node, 
        Array<double,3>& SF_Grid2D_perp_Node, 
        Array<double,3>& counter_Grid2D_Node)
{
    if (rank_mpi==0) {
        cout<<"\nComputing longitudinal and transverse S(l) and S(lx, lz) using 2D velocity field data..\n"; 
    }
    double S=0;
    double upll=0;
    double uperp=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);


  for (int x1=starPt; x1<endPt; x1++){  
      
     
     # pragma omp parallel for num_threads(num) private(upll,uperp)
      for(int z1=0; z1<Nz; z1++){
        Array<double,2> local_SF(Nr,q2-q1+1);
        Array<double,2> local_SF_p(Nr,q2-q1+1);
        Array<double,2> local_counter(Nr,q2-q1+1);
        Array<double,3> local_SF_Grid_pll(2*Nx-1,2*Nz-1,q2-q1+1); 
        Array<double,3> local_SF_Grid_perp(2*Nx-1,2*Nz-1,q2-q1+1);
        Array<double,3> local_counter_Grid(2*Nx-1,2*Nz-1,q2-q1+1);
        local_SF=0;
        local_SF_p=0;
        local_counter=0;
        local_counter_Grid=0;
        local_SF_Grid_pll=0;
        local_SF_Grid_perp=0;
        for (int x2=0;x2<Nx;x2++){
          for (int z2=0; z2<Nz; z2++) {
                  int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((z2-z1),2)));
                  int r=(int)Lmag;
                  double rmag;
                  TinyVector<double,2> u1(Ux(x1,z1),Uz(x1,z1)), u2(Ux(x2,z2),Uz(x2,z2)), rv(dx*(x2 - x1), dz*(z2 - z1)),uperpv;
                  magnitude(rv,rmag);
                  if (r==0){
                    upll=0;
                    uperp=0;
                  }
                  else{
                      upll=dot(u2-u1,rv)/rmag;
                      uperpv=(u2-u1)-rv*(upll)/(rmag);
                      magnitude(uperpv,uperp);
                  }

                  for (int p=0;p<=q2-q1;p++){
                    local_SF(Lmag,p)+= powInt(upll,q1+p);
                    local_SF_p(Lmag,p)+= powInt(uperp,q1+p);
                    local_counter(Lmag,p)+=1;
                    local_SF_Grid_pll(Nx+(x2-x1-1), Nz+(z2-z1),p) += powInt(upll,q1+p);
                    local_SF_Grid_perp(Nx+(x2-x1-1), Nz+(z2-z1),p) += powInt(uperp,q1+p);
                    local_counter_Grid(Nx+(x2-x1-1), Nz+(z2-z1),p) += 1;
                    
                  }
            } 
          }
          # pragma omp critical
          {
            SF_Node += local_SF;
            SF_p_Node += local_SF_p;
            counter_Node += local_counter;
            SF_Grid2D_pll_Node += local_SF_Grid_pll;
            SF_Grid2D_perp_Node += local_SF_Grid_perp;
            counter_Grid2D_Node += local_counter_Grid;
          
          } 

      }
  }

}


/**
 ********************************************************************************************************************************************
 * \brief   A less computationally intensive function to calculate only the longitudinal structure functions as functions of (lx,lz) in
 *          addition to function of (l) using 2D velocity field.
 *
 *          The following function computes the only longitudinal structure function using 2D velocity field data. This function is less computationally 
 *          expensive compared to SFunc2D. The function exploits the fact that \f$ \langle du(l)^p \rangle = \langle du(-l)^p \rangle \f$, where \f$ du(l) \f$ 
 *          is the component parallel to l. Thus, although this function uses six nested loops, the innermost loop starts from z1 instead of 0, where z1 is 
 *          the iteration number of the third for-loop. The rest of the structure is similar to the function SFunc3D. 
 * 
 * \param Ux is a 2D array representing the x-component of velocity field
 * \param Uz is a 2D array representing the z-component of velocity field 
 * \param SF_Node is a 2D array containing the values of the nodal longitudinal structure functions for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 * \param SF_Grid2D_pll_Node is a 3D array containing the values of the nodal longitudinal structure functions as function of (lx,lz) for a range of orders specified by the user.
 * \param counter_Grid2D_Node is a 3D array containing the numbers for dividing the values of SF_Grid_pll_Node so as to get the average. 
 ********************************************************************************************************************************************
 */
void SFunc_long_2D(
        Array<double,2> Ux, 
        Array<double,2> Uz, 
        Array<double,2>& SF_Node,
        Array<double,2>& counter_Node, 
        Array<double,3>& SF_Grid2D_pll_Node, 
        Array<double,3>& counter_Grid2D_Node)
{
    if (rank_mpi==0) {
        cout<<"\nComputing only longitudinal S(l) and S(lx, lz) using 2D velocity field data..\n";
    }
    double S=0;
    double upll=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);

  for (int x1=starPt; x1<endPt; x1++){  
      
     
     # pragma omp parallel for num_threads(num) private(upll)
      for(int z1=0; z1<Nz; z1++){
        Array<double,2> local_SF(Nr,q2-q1+1);
        Array<double,2> local_counter(Nr,q2-q1+1);
        Array<double,3> local_SF_Grid_pll(2*Nx-1,Nz,q2-q1+1); 
        Array<double,3> local_counter_Grid(2*Nx-1,Nz,q2-q1+1);
        local_SF=0;
        local_counter=0;
        local_counter_Grid=0;
        local_SF_Grid_pll=0;
        for (int x2=0;x2<Nx;x2++){
          for (int z2=z1; z2<Nz; z2++) {
                  int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((z2-z1),2)));
                  int r=(int)Lmag;
                  double rmag;
                  TinyVector<double,2> u1(Ux(x1,z1),Uz(x1,z1)), 
                                       u2(Ux(x2,z2),Uz(x2,z2)), 
                                       rv(dx*(x2 - x1), dz*(z2 - z1));
                  magnitude(rv,rmag);
                  if (r==0){
                    upll=0;
                  }
                  else{
                      upll=dot(u2-u1,rv)/rmag;
                  }
                  for (int p=0;p<=q2-q1;p++){
                    local_SF(Lmag,p)+= powInt(upll,q1+p);
                    local_counter(Lmag,p)+=1;
                    local_SF_Grid_pll(Nx+(x2-x1-1), (z2-z1),p) += powInt(upll,q1+p);
                    local_counter_Grid(Nx+(x2-x1-1),(z2-z1),p) += 1;
                    
                  }
            } 
          }
          # pragma omp critical
          {
            SF_Node += local_SF;
            counter_Node += local_counter;
            SF_Grid2D_pll_Node += local_SF_Grid_pll;
            counter_Grid2D_Node += local_counter_Grid;
          
          }   
      }
  }
}





//SCALAR STRUCTURE FUNCTION 3D GRID
/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 3D scalar field as function of (lx, ly, lz) in addition to the structure functions
 *          as function of (l).
 *
 *          The following function computes both the nodal structure functions of 3D scalar field using six nested for-loops.
 *          The outer three for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner three loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param T is a 3D array representing the scalar field 
 * \param SF_Node is a 2D array containing the values of the nodal structure functions as function of (l) for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 * \param SF_Grid_Node is a 4D array containing the values of the nodal structure functions as function of (lx,ly,lz) for a range of orders specified by the user.
 * \param counter_Grid_Node is a 4D array containing the numbers for dividing the values of SF_Grid_pll_Node so as to get the average. 
 ********************************************************************************************************************************************
 */
void SF_scalar_3D(Array<double,3> T, 
                  Array<double,2>& SF_Node,
                  Array<double,2>& counter_Node,
                  Array<double,4>& SF_Grid_Node, 
                  Array<double,4>& counter_Grid_Node
              ) {
    if (rank_mpi==0) {
        cout<<"\nComputing S(l) and S(lx, ly, lz) using 3D scalar field data..\n";
    }
    double S=0;
    double dT=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);

    for (int x1=starPt; x1<endPt; x1++){


     # pragma omp parallel for num_threads(num) private(dT)
      for(int y1=0; y1<Ny; y1++){
        Array<double,2> local_SF(Nr,q2-q1+1);
        Array<double,2> local_counter(Nr,q2-q1+1);
        Array<double,4> local_SF_Grid(2*Nx-1,2*Ny-1,2*Nz-1,q2-q1+1); //,local_SF_Grid_perp(2*Nx-1,2*Ny-1,Nz,q2-q1+1);
        Array<double,4> local_counter_Grid(2*Nx-1,2*Ny-1,2*Nz-1,q2-q1+1);
        local_SF=0;
        local_counter=0;
        local_counter_Grid=0;
        local_SF_Grid=0;
        for (int z1=0;z1<Nz;z1++){
          for (int x2=0; x2<Nx; x2++) {
              for(int y2=0; y2<Ny; y2++){
                for(int z2=0;z2<Nz;z2++){
                  int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((y2-y1),2)+pow((z2-z1),2)));
                  int r=(int)Lmag;
                  if (r==0){
                      dT=0;
                  }
                  else{
                      dT = T(x2,y2,z2)-T(x1,y1,z1);
                  }
                  for (int p=0;p<=q2-q1;p++){
                      local_SF(Lmag,p) += powInt(dT, q1+p);
                      local_counter(Lmag,p) += 1;
                      local_SF_Grid(Nx+(x2-x1-1), Ny+(y2-y1-1), Nz+(z2-z1),p) += pow(dT,q1+p);
                      local_counter_Grid(Nx+(x2-x1-1), Ny+(y2-y1-1), Nz+(z2-z1),p) += 1;

                  }
                }
              }
            }
          }
          # pragma omp critical
          {
             SF_Grid_Node += local_SF_Grid;
             SF_Node += local_SF;
             counter_Grid_Node += local_counter_Grid;
             counter_Node += local_counter;
          }
      }
  }
}

/**
 ********************************************************************************************************************************************
 * \brief   Function to calculate structure functions using 2D scalar field as function of (lx, lz) in addition to the structure functions
 *          as function of (l).
 *
 *          The following function computes both the nodal structure functions of 2D scalar field using four nested for-loops.
 *          The outer two for-loops correspond to the vector \f$\mathbf{r}\f$ and the inner two loops correspond to the vector \f$\mathbf{r+l}\f$. 
 *          The second for-loop is parallelized using OpenMP.
 * 
 * \param T is a 2D array representing the scalar field 
 * \param SF_Node is a 2D array containing the values of the nodal structure functions as function of (l) for a range of orders specified by the user.
 * \param counter_Node is a 2D array containing the numbers for dividing the values of SF_Node so as to get the average. 
 * \param SF_Grid_Node is a 3D array containing the values of the nodal structure functions as function of (lx,lz) for a range of orders specified by the user.
 * \param counter_Grid_Node is a 3D array containing the numbers for dividing the values of SF_Grid_pll_Node so as to get the average. 
 ********************************************************************************************************************************************
 */
void SF_scalar_2D(Array<double,2> T, 
                  Array<double,2>& SF_Node,
                  Array<double,2>& counter_Node,
                  Array<double,3>& SF_Grid_Node, 
                  Array<double,3>& counter_Grid_Node
              ) {
    
    if (rank_mpi==0){
        cout<<"\nComputing S(l) and S(lx, lz) using 2D scalar field data..\n";
    }
    double S=0;
    double dT=0;
    int starPt = rank_mpi*(Nx/P);
    int endPt = (rank_mpi+1)*(Nx/P);

    for (int x1=starPt; x1<endPt; x1++){


     # pragma omp parallel for num_threads(num) private(dT)
      for(int z1=0; z1<Nz; z1++){
        Array<double,2> local_SF(Nr,q2-q1+1);
        Array<double,2> local_counter(Nr,q2-q1+1);
        Array<double,3> local_SF_Grid(2*Nx-1,2*Nz-1,q2-q1+1); //,local_SF_Grid_perp(2*Nx-1,2*Ny-1,Nz,q2-q1+1);
        Array<double,3> local_counter_Grid(2*Nx-1,2*Nz-1,q2-q1+1);
        local_SF=0;
        local_counter=0;
        local_counter_Grid=0;
        local_SF_Grid=0;
        for (int x2=0;x2<Nx;x2++){
          for (int x2=0; x2<Nx; x2++) {
              for(int z2=0; z2<Nz; z2++){
                  int Lmag = ceil(sqrt(pow((x2-x1),2) + pow((z2-z1),2)));
                  int r=(int)Lmag;
                  if (r==0){
                      dT=0;
                  }
                  else{
                      dT = T(x2,z2)-T(x1,z1);
                  }
                  for (int p=0;p<=q2-q1;p++){
                      local_SF(Lmag,p) += powInt(dT, q1+p);
                      local_counter(Lmag,p) += 1;
                      local_SF_Grid(Nx+(x2-x1-1), Nz+(z2-z1-1), p) += pow(dT,q1+p);
                      local_counter_Grid(Nx+(x2-x1-1), Nz+(z2-z1-1), p) += 1;

                  } //for p
                } //for z2
              } //for x2
            
          # pragma omp critical
          {
             SF_Grid_Node += local_SF_Grid;
             SF_Node += local_SF;
             counter_Grid_Node += local_counter_Grid;
             counter_Node += local_counter;
          }
      }
      } //for z1
  } //for x1
}


