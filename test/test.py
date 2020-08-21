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
N=512


x=np.linspace(0,1,N)
y=np.linspace(0,1,N)
X,Y=np.meshgrid(x,y)




	
def plotSF_r_2D():
	SF2 = (hdf5_reader_plane("test_velocity_2D/out/SF_Grid_pll.h5", "SF_Grid_pll2"))
	SF3 = (hdf5_reader_plane("test_velocity_2D/out/SF_Grid_pll.h5", "SF_Grid_pll3"))
	
	Nx, Nz = SF2.shape
	
	#Nr = int(np.ceil(np.sqrt(Nx**2 + Nz**2)))
	Nr = int(np.ceil(np.sqrt((Nx-1)**2 + (Nz-1)**2)))+1

	r = np.zeros([Nr]) #
	for i in range(len(r)): #
		r[i]=np.sqrt(2)*i/(2*len(r))
	
	SF2_r = np.zeros([Nr])
	SF3_r = np.zeros([Nr])
	counter = np.zeros([Nr])
	
	for x in range(Nx):
	    for z in range(Nz):
	        l = int(np.ceil(np.sqrt(x**2 + z**2)))
	        SF2_r[l] = SF2_r[l] + SF2[x, z]
	        SF3_r[l] = SF3_r[l] + SF3[x, z]
	        counter[l] = counter[l] + 1	  
	  
	SF2_r = SF2_r/counter
	SF3_r = SF3_r/counter
	
	llim = 5
	ulim = 25
	fig, axes = plt.subplots(figsize = (3.5, 2.57))
	
	axes.plot(r[1:len(r)], SF2_r[1:len(r)], color='red', lw=1.5, label=r"$S_2^{u}(l)$") #
	axes.plot(r[1:len(r)], SF3_r[1:len(r)], color='green', lw=1.5, label=r"$S_3^{u}(l)$")
	
	
	s2 = r*r
	s3 = r*r*r
	
	axes.plot(r[llim:ulim], s2[llim:ulim], color='black', lw=1.5, linestyle= 'dashed')#, label=r"$S_q^{u}(l)$") #
	axes.plot(r[llim:ulim], s3[llim:ulim], color='black', lw=1.5, linestyle= 'dashed')#, label=r"$S_q^{u}(l)$")
	
	axes.text(0.4, 0.2, "$l^2$")
	axes.text(0.4, 0.035, "$l^3$")
	 
	axes.set_xlabel('$l$')
	axes.set_ylabel('$S_q^{u}(l)$')
	axes.set_xscale('log')
	axes.set_yscale('log')
	axes.set_xlim(0.1, 0.7)
	axes.set_ylim(2e-3, 1)
	axes.set_xticks([0.1, 0.4, 0.7])
	axes.set_xticklabels([0.1, 0.4, 0.7])
	axes.legend(scatterpoints=1, loc='lower right', prop={'size':0.95*A}, ncol = 1, frameon=False)
	fig.tight_layout()
	plt.savefig("SF_velocity_r2D.png", dpi=600)
	plt.show()


def plotSF_r_3D():
	SF2 = (hdf5_reader("test_velocity_3D/out/SF_Grid_pll.h5", "SF_Grid_pll2"))
	SF3 = (hdf5_reader("test_velocity_3D/out/SF_Grid_pll.h5", "SF_Grid_pll3"))
	
	Nx, Ny, Nz = SF2.shape
	
	#Nr = int(np.ceil(np.sqrt(Nx**2 + Ny**2 + Nz**2)))
	Nr = int(np.ceil(np.sqrt((Nx-1)**2 + (Ny-1)**2 + (Nz-1)**2)))+1
	

	r = np.zeros([Nr]) #
	for i in range(len(r)): #
		r[i]=np.sqrt(3)*i/(2*len(r))
	
	SF2_r = np.zeros([Nr])
	SF3_r = np.zeros([Nr])
	counter = np.zeros([Nr])
	
	for x in range(Nx):
	    for y in range(Ny):
	        for z in range(Nz):
	           l = int(np.ceil(np.sqrt(x**2 + y**2 + z**2)))
	           SF2_r[l] = SF2_r[l] + SF2[x, y, z]
	           SF3_r[l] = SF3_r[l] + SF3[x, y, z]
	           counter[l] = counter[l] + 1	  
	
	SF2_r = SF2_r/counter
	SF3_r = SF3_r/counter
	
	llim = 10
	ulim = 25
	fig, axes = plt.subplots(figsize = (3.5, 2.57))
	
	axes.plot(r[1:len(r)], SF2_r[1:len(r)], color='red', lw=1.5, label=r"$S_2^{u}(l)$") #
	axes.plot(r[1:len(r)], SF3_r[1:len(r)], color='green', lw=1.5, label=r"$S_3^{u}(l)$")
	
	
	s2 = r*r
	s3 = r*r*r
	
	axes.plot(r[llim:ulim], s2[llim:ulim], color='black', lw=1.5, linestyle= 'dashed')#, label=r"$S_q^{u}(l)$") #
	axes.plot(r[llim:ulim], s3[llim:ulim], color='black', lw=1.5, linestyle= 'dashed')#, label=r"$S_q^{u}(l)$")
	
	axes.text(0.4, 0.2, "$l^2$")
	axes.text(0.4, 0.035, "$l^3$")
	 
	axes.set_xlabel('$l$')
	axes.set_ylabel('$S_q^{u}(l)$')
	axes.set_xscale('log')
	axes.set_yscale('log')
	axes.set_xlim(0.15, 0.8)
	axes.set_ylim(2e-3, 1)
	axes.set_xticks([0.2, 0.4, 0.8])
	axes.set_xticklabels([0.2, 0.4, 0.8])
	axes.legend(scatterpoints=1, loc='lower right', prop={'size':0.95*A}, ncol = 1, frameon=False)
	fig.tight_layout()
	plt.savefig("SF_velocity_r3D.png", dpi=600)
	plt.show()	

	
def plot_SF2D_scalar():
    SF2 = (hdf5_reader_plane("test_scalar_2D/out/SF_Grid_scalar.h5", "SF_Grid_scalar2"))
    #SF3 = (hdf5_reader_plane("test_scalar_2D/out/SF_Grid_scalar3.h5", "SF_Grid_scalar3"))
    
    Nlx, Nlz = SF2.shape
    
    lx = np.linspace(0,0.5,Nlx)
    lz = np.linspace(0,0.5,Nlz)
    
    fig, axes = plt.subplots(1,2,figsize=(5,2.5),sharey=True)
    levels = []
    Lz,Lx=np.meshgrid(lz,lx)
    Z=(Lx+Lz)**2
    
    axes[1].contourf(Lx, Lz, Z, levels= np.linspace(0,1,50), cmap='jet')
    
    density = axes[0].contourf(lx, lz, np.transpose(SF2), levels=np.linspace(0,1,50), cmap='jet')
    #density = axes.pcolor(lx, lz, np.transpose(SF), cmap='jet', norm=colors.SymLogNorm(linthresh=1e-4, linscale=0.1, vmin=0, vmax=4.0))

    
    #axes[0].set_aspect(1)
    axes[0].set_xticks([0, 0.5])
    axes[0].set_yticks([0, 0.5])
    axes[0].set_xlabel('$l_x$')
    axes[0].set_ylabel('$l_z$')
    axes[0].tick_params(axis='x', which='major', pad=10)
    axes[0].tick_params(axis='y', which='major', pad=10)
    axes[0].set_title(r"$(\mathrm{a})$ $S_2^{\theta}(l_x,l_z)$")#, pad=10)

    axes[1].set_title(r"$(\mathrm{b})$ $(l_x + l_z)^2$")
    
    #axes[1].set_aspect(1)
    axes[1].set_xticks([0, 0.5])
    axes[1].set_yticks([0, 0.5])
    axes[1].set_xlabel('$l_x$')
    #axes[1].set_ylabel('$l_z$')
    axes[1].tick_params(axis='x', which='major', pad=10)
    axes[1].tick_params(axis='y', which='major', pad=10)
    axes[0].title.set_position([.5, 1.05])
    axes[1].title.set_position([.5, 1.05])
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes[0],ticks=[0, 0.5, 1])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb1.ax.tick_params(labelsize=A)
    
    cb2 = fig.colorbar(density, fraction=0.05, ax=axes[1],ticks=[0, 0.5, 1])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb2.ax.tick_params(labelsize=A)
    fig.tight_layout()
    plt.savefig("SF_scalar2D.png", dpi=600)
    
    plt.show()

	
def plot_SF2D_velocity():
    SF2 = (hdf5_reader_plane("test_velocity_2D/out/SF_Grid_pll.h5", "SF_Grid_pll2"))
    #SF3 = (hdf5_reader_plane("test_scalar_2D/out/SF_Grid_scalar3.h5", "SF_Grid_scalar3"))
    
    Nlx, Nlz = SF2.shape
    
    lx = np.linspace(0,0.5,Nlx)
    lz = np.linspace(0,0.5,Nlz)
    
    fig, axes = plt.subplots(1,2,figsize=(5,2.5),sharey=True)
    levels = []
    Lz,Lx=np.meshgrid(lz,lx)
    Z=(Lx**2+Lz**2)
    
    axes[1].contourf(Lx, Lz, Z, levels= np.linspace(0,0.5,50), cmap='jet')
    
    density = axes[0].contourf(lx, lz, np.transpose(SF2), levels=np.linspace(0,0.5,50), cmap='jet')
   
    axes[0].set_xticks([0, 0.25, 0.5])
    axes[0].set_yticks([0, 0.25, 0.5])
    axes[0].set_xlabel('$l_x$')
    axes[0].set_ylabel('$l_z$')
    axes[0].tick_params(axis='x', which='major', pad=10)
    axes[0].tick_params(axis='y', which='major', pad=10)
    axes[0].set_title(r"$(\mathrm{a})$ $S_2^{u}(l_x,l_z)$")#, pad=10)

    axes[1].set_title(r"$(\mathrm{b})$ $(l_x^2 + l_z^2)$")
    
    #axes[1].set_aspect(1)
    axes[1].set_xticks([0, 0.25, 0.5])
    axes[1].set_yticks([0, 0.25, 0.5])
    axes[1].set_xlabel('$l_x$')
    #axes[1].set_ylabel('$l_z$')
    axes[1].tick_params(axis='x', which='major', pad=10)
    axes[1].tick_params(axis='y', which='major', pad=10)
    axes[0].title.set_position([.5, 1.05])
    axes[1].title.set_position([.5, 1.05])
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes[0],ticks=[0, 0.25, 0.5])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb1.ax.tick_params(labelsize=A)
    
    cb2 = fig.colorbar(density, fraction=0.05, ax=axes[1],ticks=[0, 0.25, 0.5])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb2.ax.tick_params(labelsize=A)
    fig.tight_layout()
    plt.savefig("SF_velocity2D.png", dpi=600)
    
    plt.show()

def plot_SF3D_scalar():
    SF2 = (hdf5_reader_slice("test_scalar_3D/out/SF_Grid_scalar.h5", "SF_Grid_scalar2",-1))
    #SF3 = (hdf5_reader_slice("out/SF_Grid_scalar3.h5", "SF_Grid_scalar3",31))
    
    Nlx, Nlz= SF2.shape

    lx = np.linspace(0, 0.5, Nlx)
    lz = np.linspace(0, 0.5, Nlz)
    #lz = np.linspace(0,1,Nlz)
    
    fig, axes = plt.subplots(1,2,figsize=(5,2.5),sharey=True)
    
    Lz,Lx=np.meshgrid(lz,lx)
    Z=(Lx+Lz+0.5)**2

   
    axes[1].contourf(Lx, Lz, Z, levels= np.linspace(0,2.25,50), cmap='jet')
    
    density = axes[0].contourf(lx, lz, np.transpose(SF2), levels=np.linspace(0,2.25,50), cmap='jet')
    #density = axes.pcolor(lx, lz, np.transpose(SF), cmap='jet', norm=colors.SymLogNorm(linthresh=1e-4, linscale=0.1, vmin=0, vmax=4.0))

    
    #axes[0].set_aspect(1)
    axes[0].set_xticks([0, 0.25, 0.5])
    axes[0].set_yticks([0, 0.25, 0.5])
    axes[0].set_xlabel('$l_x$')
    axes[0].set_ylabel('$l_z$')
    axes[0].tick_params(axis='x', which='major', pad=10)
    axes[0].tick_params(axis='y', which='major', pad=10)
    axes[0].set_title(r"$(\mathrm{a})$ $S_2^{\theta}(l_x, l_y=0.5,l_z)$")#, pad=10)

    axes[1].set_title(r"$(\mathrm{b})$ $(l_x +0.5+ l_z)^2$")
    
    #axes[1].set_aspect(1)
    axes[1].set_xticks([0, 0.25, 0.5])
    axes[1].set_yticks([0, 0.25, 0.5])
    axes[1].set_xlabel('$l_x$')
    #axes[1].set_ylabel('$l_z$')
    axes[1].tick_params(axis='x', which='major', pad=10)
    axes[1].tick_params(axis='y', which='major', pad=10)
    axes[0].title.set_position([.5, 1.05])
    axes[1].title.set_position([.5, 1.05])
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes[0],ticks=[0, 0.75, 1.5, 2.25])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb1.ax.tick_params(labelsize=A)
    
    cb2 = fig.colorbar(density, fraction=0.05, ax=axes[1],ticks=[0, 0.75, 1.5, 2.25])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb2.ax.tick_params(labelsize=A)
    fig.tight_layout()
    plt.savefig("SF_scalar3D.png", dpi=600)
    
    plt.show()
    
def plot_SF3D_velocity():
    SF2 = (hdf5_reader_slice("test_velocity_3D/out/SF_Grid_pll.h5", "SF_Grid_pll2",-1))
    #SF3 = (hdf5_reader_slice("out/SF_Grid_scalar3.h5", "SF_Grid_scalar3",31))
    
    Nlx, Nlz= SF2.shape

    lx = np.linspace(0, 0.5, Nlx)
    lz = np.linspace(0, 0.5, Nlz)
    #lz = np.linspace(0,1,Nlz)
    
    fig, axes = plt.subplots(1,2,figsize=(5,2.5),sharey=True)
    
    Lz,Lx=np.meshgrid(lz,lx)
    Z=(Lx**2+Lz**2+0.5**2)

   
    axes[1].contourf(Lx, Lz, Z, levels= np.linspace(0,0.75,50), cmap='jet')
    
    density = axes[0].contourf(lx, lz, np.transpose(SF2), levels=np.linspace(0,0.75,50), cmap='jet')
    #density = axes.pcolor(lx, lz, np.transpose(SF), cmap='jet', norm=colors.SymLogNorm(linthresh=1e-4, linscale=0.1, vmin=0, vmax=4.0))

    
    #axes[0].set_aspect(1)
    axes[0].set_xticks([0, 0.25, 0.5])
    axes[0].set_yticks([0, 0.25, 0.5])
    axes[0].set_xlabel('$l_x$')
    axes[0].set_ylabel('$l_z$')
    axes[0].tick_params(axis='x', which='major', pad=10)
    axes[0].tick_params(axis='y', which='major', pad=10)
    axes[0].set_title(r"$(\mathrm{a})$ $S_2^{u}(l_x, l_y=0.5,l_z)$")#, pad=10)

    axes[1].set_title(r"$(\mathrm{b})$ $(l_x^2 +0.5^2+ l_z^2)$")
    
    #axes[1].set_aspect(1)
    axes[1].set_xticks([0, 0.25, 0.5])
    axes[1].set_yticks([0, 0.25, 0.5])
    axes[1].set_xlabel('$l_x$')
    #axes[1].set_ylabel('$l_z$')
    axes[1].tick_params(axis='x', which='major', pad=10)
    axes[1].tick_params(axis='y', which='major', pad=10)
    axes[0].title.set_position([.5, 1.05])
    axes[1].title.set_position([.5, 1.05])
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes[0],ticks=[0, 0.25,0.5,0.75])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb1.ax.tick_params(labelsize=A)
    
    cb2 = fig.colorbar(density, fraction=0.05, ax=axes[1],ticks=[0, 0.25, 0.5, 0.75])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb2.ax.tick_params(labelsize=A)
    fig.tight_layout()
    plt.savefig("SF_velocity3D.png", dpi=600)
    
    plt.show()


plotSF_r_2D()
plot_SF2D_scalar()
#plot_SF2D_velocity()

plotSF_r_3D()
plot_SF3D_scalar()
#plot_SF3D_velocity()