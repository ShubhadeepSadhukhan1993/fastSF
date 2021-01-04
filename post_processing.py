#############################################################################################################################################
 # fastSF
 # 
 # Copyright (C) 2020, Mahendra K. Verma
 #
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #     1. Redistributions of source code must retain the above copyright
 #        notice, this list of conditions and the following disclaimer.
 #     2. Redistributions in binary form must reproduce the above copyright
 #        notice, this list of conditions and the following disclaimer in the
 #        documentation and/or other materials provided with the distribution.
 #     3. Neither the name of the copyright holder nor the
 #        names of its contributors may be used to endorse or promote products
 #        derived from this software without specific prior written permission.
 # 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 # ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 # SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 ############################################################################################################################################
 ##
 ##! \file test.py
 #
 #   \brief Script to generate plots for the test cases.
 #
 #   \author Shashwat Bhattacharya, Shubhadeep Sadhukhan
 #   \date Feb 2020
 #   \copyright New BSD License
 #
 ############################################################################################################################################
##
import yaml
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker, colors

mpl.style.use('classic')

plt.rcParams['xtick.major.size'] = 4.2
plt.rcParams['xtick.major.width'] = 0.5
plt.rcParams['xtick.minor.size'] = 2.5
plt.rcParams['xtick.minor.width'] = 0.5
plt.rcParams['ytick.major.size'] = 4.2
plt.rcParams['ytick.major.width'] = 0.5
plt.rcParams['ytick.minor.size'] = 2.5
plt.rcParams['ytick.minor.width'] = 0.5
#plt.rcParams['axes.titlepad'] = 10



input_path=input("Enter the path of your in folder: ")
data_path=input("Enter the path of your out folder: ")
q=input("Enter the order: ")

with open(input_path+"/in/para.yaml", 'r') as stream:
	try:
		para=(yaml.safe_load(stream))
		Lxx=para['domain_dimension']['Lx']
		Lyy=para['domain_dimension']['Ly']
		Lzz=para['domain_dimension']['Lz']
		Nxx=para['grid']['Nx']
		Nyy=para['grid']['Ny']
		Nzz=para['grid']['Nz']
		scalar_switch=para['program']['scalar_switch']
		two_dim_switch=para['program']['2D_switch']
		longitudinal=para["program"]["Only_longitudinal"]
		
	except yaml.YAMLError as exc:
		print(exc)
        	

def hdf5_reader(filename,dataset):
	file_V1_read = h5py.File(filename, 'r')
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:,:,:]
	return V1

def hdf5_reader1D(filename,dataset):
	file_V1_read = h5py.File(filename, 'r')
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:]
	return V1


def hdf5_reader_plane(filename,dataset):
	file_V1_read = h5py.File(filename, 'r')
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:,:]
	return V1


def hdf5_reader_slice(filename,dataset, Ny):
	file_V1_read = h5py.File(filename)
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:,Ny,:]
	return V1
A=9.3
font = {'family' : 'serif', 'weight' : 'normal', 'size' : A}
plt.rc('font', **font)
B = 19





	
def plotSF_r_2D(data_path, q):
	SF = (hdf5_reader_plane(data_path+"/out/SF_Grid_pll.h5", "SF_Grid_pll"+str(q)))
	Nx, Nz = SF.shape
	
	Nr = int(np.ceil(np.sqrt((Nx-1)**2 + (Nz-1)**2)))+1

	r = np.zeros([Nr]) 
	for i in range(len(r)): #
		r[i]=np.sqrt(2)*i/(2*len(r))
	
	SF_r = np.zeros([Nr])
	
	counter = np.zeros([Nr])
	
	for x in range(Nx):
	    for z in range(Nz):
	        l = int(np.ceil(np.sqrt(x**2 + z**2)))
	        SF_r[l] = SF_r[l] + SF[x, z]
	        
	        counter[l] = counter[l] + 1	  
	  
	SF_r = SF_r/counter
	
	fig, axes = plt.subplots(figsize = (3.5, 2.57))
	
	axes.plot(r[1:len(r)], SF_r[1:len(r)], color='red', lw=1.5, label=r"$S_2^{u}(l)$") 
	
	 
	axes.set_xlabel('$l$')
	axes.set_ylabel('$S_q^{u}(l)$')
	axes.set_xscale('log')
	axes.set_yscale('log')
	fig.tight_layout()
	plt.savefig("SF_velocity_r2D.png", dpi=600)
	


def plotSF_r_3D(data_path, q):
	SF = (hdf5_reader(data_path+"/out/SF_Grid_pll.h5", "SF_Grid_pll"+str(q)))
	
	
	Nx, Ny, Nz = SF.shape
	
	Nr = int(np.ceil(np.sqrt((Nx-1)**2 + (Ny-1)**2 + (Nz-1)**2)))+1
	

	r = np.zeros([Nr]) 
	for i in range(len(r)): 
		r[i]=np.sqrt(3)*i/(2*len(r))
	
	SF_r = np.zeros([Nr])
	
	counter = np.zeros([Nr])
	
	for x in range(Nx):
	    for y in range(Ny):
	        for z in range(Nz):
	           l = int(np.ceil(np.sqrt(x**2 + y**2 + z**2)))
	           SF_r[l] = SF_r[l] + SF[x, y, z]
	           
	           counter[l] = counter[l] + 1	  
	
	SF_r = SF_r/counter
	
	fig, axes = plt.subplots(figsize = (3.5, 2.57))
	
	axes.plot(r[1:len(r)], SF_r[1:len(r)], color='red', lw=1.5, label=r"$S^{u}_{\parallel}(l)$")
	axes.set_xlabel('$l$')
	axes.set_ylabel('$S_q^{u}(l)$')
	axes.set_xscale('log')
	axes.set_yscale('log')
	fig.tight_layout()
	plt.savefig("SF_velocity_r3D.png", dpi=600)
		

	
def plot_SF2D_scalar(data_path, q):
    SF = (hdf5_reader_plane(data_path+"/out/SF_Grid_scalar.h5", "SF_Grid_scalar"+str(q)))
   
    
    Nlx, Nlz = SF.shape
    
    lx = np.linspace(0,Lxx/2.,Nlx)
    lz = np.linspace(0,Lzz/2.,Nlz)
    
    fig, axes = plt.subplots(1,1,figsize=(3.5,2.7),sharey=True)
    Lz,Lx=np.meshgrid(lz,lx)
    
    density = axes.contourf(lx, lz, np.transpose(SF), levels=np.linspace(SF.min(), SF.max(),50), cmap='jet')
    
    axes.set_xticks([0, Lxx/4., Lxx/2.])
    axes.set_yticks([0, Lzz/4., Lzz/2.])
    axes.set_xlabel('$l_x$')
    axes.set_ylabel('$l_z$')
    axes.tick_params(axis='x', which='major', pad=10)
    axes.tick_params(axis='y', which='major', pad=10)
    axes.set_title(r"$S^{\theta}$")

    
    axes.title.set_position([.5, 1.05])
    
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes,ticks=np.linspace(SF.min(), SF.max(),4))
    cb1.ax.tick_params(labelsize=A)
    
    fig.tight_layout()
    plt.savefig("SF_scalar2D.png", dpi=600)
    
   

	
def plot_SF2D_velocity(data_path, q):
    SFpll = (hdf5_reader_plane(data_path+"/out/SF_Grid_pll.h5", "SF_Grid_pll"+str(q)))
   
    
    Nlx, Nlz = SFpll.shape
    
    lx = np.linspace(0,Lxx/2.,Nlx)
    lz = np.linspace(0,Lzz/2.,Nlz)
    
    fig, axes = plt.subplots(1,1,figsize=(3.5,2.7),sharey=True)
  
    Lz,Lx=np.meshgrid(lz,lx)
   
    
    density = axes.contourf(lx, lz, np.transpose(SFpll), levels=np.linspace(SFpll.min(),SFpll.max(),50), cmap='jet')
    axes.set_xticks([0, Lxx/4., Lxx/2.])
    axes.set_yticks([0, Lzz/4., Lzz/2.])
    axes.set_xlabel('$l_x$')
    axes.set_ylabel('$l_z$')
    axes.tick_params(axis='x', which='major', pad=10)
    axes.tick_params(axis='y', which='major', pad=10)
    axes.set_title(r"$S^{u}_{\parallel}$")
    axes.title.set_position([.5, 1.05])
   
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes,ticks=np.linspace(SFpll.min(),SFpll.max(),4))#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb1.ax.tick_params(labelsize=A)
    
   
    fig.tight_layout()
    plt.savefig("SF_velocity2D_pll.png", dpi=600)
    
    if (longitudinal==False):
		SFperp = (hdf5_reader_plane(data_path+"/out/SF_Grid_perp.h5", "SF_Grid_perp"+str(q)))
		Nlx, Nlz = SFperp.shape
		lx = np.linspace(0,Lxx/2.,Nlx)
		lz = np.linspace(0,Lzz/2.,Nlz)
    
		fig, axes = plt.subplots(1,1,figsize=(3.5,2.7),sharey=True)
  
		Lz,Lx=np.meshgrid(lz,lx)
   
    
		density = axes.contourf(lx, lz, np.transpose(SFperp), levels=np.linspace(SFperp.min(),SFperp.max(),50), cmap='jet')
   
		axes.set_xticks([0, Lxx/4., Lxx/2.])
		axes.set_yticks([0, Lxx/4., Lzz/2.])
		axes.set_xlabel('$l_x$')
		axes.set_ylabel('$l_z$')
		axes.tick_params(axis='x', which='major', pad=10)
		axes.tick_params(axis='y', which='major', pad=10)
		axes.set_title(r"$S^{u}_{\perp}$")
		axes.title.set_position([.5, 1.05])
   
		cb1 = fig.colorbar(density, fraction=0.05, ax=axes,ticks=np.linspace(SFperp.min(),SFperp.max(),4))#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
		cb1.ax.tick_params(labelsize=A)
    
   
		fig.tight_layout()
		plt.savefig("SF_velocity2D_perp.png", dpi=600)
		
    
    
 

def plot_SF3D_scalar(data_path, q):
    SF = (hdf5_reader_slice(data_path+"/out/SF_Grid_scalar.h5", "SF_Grid_scalar"+str(q),-1))
   
    
    Nlx, Nlz= SF.shape

    lx = np.linspace(0, Lxx/2., Nlx)
    lz = np.linspace(0, Lzz/2., Nlz)
    
    
    fig, axes = plt.subplots(1,1,figsize=(3.5,2.7),sharey=True)
    
    Lz,Lx=np.meshgrid(lz,lx)
   

   
    
    
    density = axes.contourf(lx, lz, np.transpose(SF), levels=np.linspace(SF.min(),SF.max(),50), cmap='jet')
    
    
   
    axes.set_xticks([0, Lzz/4., Lxx/2.])
    axes.set_yticks([0, Lzz/4., Lzz/2.])
    axes.set_xlabel('$l_x$')
    axes.set_ylabel('$l_z$')
    axes.tick_params(axis='x', which='major', pad=10)
    axes.tick_params(axis='y', which='major', pad=10)
    axes.set_title(r"$S^{\theta}$")#, pad=10)

    
    axes.title.set_position([.5, 1.05])
    
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes,ticks=np.linspace(SF.min(),SF.max(),4))#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb1.ax.tick_params(labelsize=A)
    
    
    fig.tight_layout()
    plt.savefig("SF_scalar3D.png", dpi=600)
    
 
    
def plot_SF3D_velocity(data_path, q):
    SFpll = (hdf5_reader_slice(data_path+"/out/SF_Grid_pll.h5", "SF_Grid_pll"+str(q),-1))
    
    
    Nlx, Nlz= SFpll.shape

    lx = np.linspace(0, Lxx/2., Nlx)
    lz = np.linspace(0, Lzz/2., Nlz)
 
    Lz,Lx=np.meshgrid(lz,lx)
    fig, axes = plt.subplots(1,1,figsize=(3.5,2.7),sharey=True)
    
    density = axes.contourf(lx, lz, np.transpose(SFpll), levels=np.linspace(SFpll.min(),SFpll.max(),50), cmap='jet')
    
    axes.set_xticks([0, Lxx/4., Lxx/2.])
    axes.set_yticks([0, Lzz/4., Lzz/2.])
    axes.set_xlabel('$l_x$')
    axes.set_ylabel('$l_z$')
    axes.tick_params(axis='x', which='major', pad=10)
    axes.tick_params(axis='y', which='major', pad=10)
    axes.set_title(r"$S^{u}_{\parallel}$") #(l_x, l_y=0.5,l_z)$")#, pad=10)   
    axes.title.set_position([.5, 1.05])
    
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes,ticks=np.linspace(SFpll.min(),SFpll.max(),4))
    cb1.ax.tick_params(labelsize=A)
    
    
    fig.tight_layout()
    plt.savefig("SF_velocity3D_pll.png", dpi=600)
    if (longitudinal==False):
		SFperp = (hdf5_reader_slice(data_path+"/out/SF_Grid_perp.h5", "SF_Grid_perp"+str(q),-1))
		fig, axes = plt.subplots(1,1,figsize=(3.5,2.7),sharey=True)
		density = axes.contourf(lx, lz, np.transpose(SFperp), levels=np.linspace(SFperp.min(),SFperp.max(), 50), cmap='jet')
    
		axes.set_xticks([0, Lxx/4., Lxx/2.])
		axes.set_yticks([0, Lzz/4., Lzz/2.])
		axes.set_xlabel('$l_x$')
		axes.set_ylabel('$l_z$')
		axes.tick_params(axis='x', which='major', pad=10)
		axes.tick_params(axis='y', which='major', pad=10)
		axes.set_title(r"$S^{u}_{\perp}$") #(l_x, l_y=0.5,l_z)$")#, pad=10)   
		axes.title.set_position([.5, 1.05])
    
		cb1 = fig.colorbar(density, fraction=0.05, ax=axes,ticks=np.linspace(SFperp.min(),SFperp.max(),4))
		cb1.ax.tick_params(labelsize=A)
    
    
		fig.tight_layout()
		plt.savefig("SF_velocity3D_perp.png", dpi=600)
    
    
 





if scalar_switch:
	print ("Scalar")
	if two_dim_switch:
		print ("2D")
		plot_SF2D_scalar(data_path, q)
	else:
		print ("3D")
		plot_SF3D_scalar(data_path, q)
else:
	print ("Vector")
	if two_dim_switch:
		print ("2D")
		plot_SF2D_velocity(data_path, q)
		plotSF_r_2D(data_path,2)
	else:
		print ("3D")
		plot_SF3D_velocity(data_path, q)
		plotSF_r_3D(data_path,q)

plt.show()


