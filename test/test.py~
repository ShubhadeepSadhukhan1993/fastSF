import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker, colors

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
	file_V1_read = h5py.File(filename)
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:,:,:]
	return V1

def hdf5_reader1D(filename,dataset):
	file_V1_read = h5py.File(filename)
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:]
	return V1


def hdf5_reader_plane(filename,dataset):
	file_V1_read = h5py.File(filename)
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:,:]
	return V1



A=9.3
font = {'family' : 'serif', 'weight' : 'normal', 'size' : A}
plt.rc('font', **font)
B = 19
N=512


x=np.linspace(0,1,N)
y=np.linspace(0,1,N)
X,Y=np.meshgrid(x,y)


#print T[90,100,500]
#print T[:,255,255]


	
def plotSF_r2():
	S_r_Pll2 = hdf5_reader1D("test_velocity_2D/out/SF.h5","SF")[:,1] #
	S_r_Pll3 = hdf5_reader1D("test_velocity_2D/out/SF.h5","SF")[:,2] #
	
	#S_r_Pll_old2 = hdf5_reader1D("out/SF_old2.h5","SF")[:,1] #
	#S_r_Pll_old3 = hdf5_reader1D("out/SF_old3.h5","SF")[:,1] #
	r = np.zeros_like(S_r_Pll2) #
	for i in range(len(S_r_Pll2)): #
		r[i]=np.sqrt(2)*i/len(S_r_Pll2)     #
	
	llim = 10
	ulim = 25
	fig, axes = plt.subplots(figsize = (3.5, 2.57))
	
	axes.plot(r[1:len(r)], S_r_Pll2[1:len(r)], color='red', lw=1.5, label=r"$S_2^{u}(l)$") #
	axes.plot(r[1:len(r)], S_r_Pll3[1:len(r)], color='green', lw=1.5, label=r"$S_3^{u}(l)$")
	
	
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
	axes.set_xlim(0.15, 1)
	axes.set_ylim(2e-3, 1)
	axes.set_xticks([0.2, 0.4, 0.8])
	axes.set_xticklabels([0.2, 0.4, 0.8])
	axes.legend(scatterpoints=1, loc='lower right', prop={'size':0.95*A}, ncol = 1, frameon=False)
	fig.tight_layout()
	plt.savefig("SF_test.png", dpi=600)
	plt.show()

	

	
def plot_SF_density():
    SF = (hdf5_reader_plane("test_scalar_2D/out/SF_Grid_pll2.h5", "SF_Grid_pll2"))
    SF3 = (hdf5_reader_plane("test_scalar_2D/out/SF_Grid_pll3.h5", "SF_Grid_pll3"))
    
    Nlx, Nlz = SF.shape
    
    #for i in range(Nlx):
    #    for j in range(Nlz):
    #        if (SF[i,j]<1e-3):
    #            SF[i,j]=1e-3
    
    lx = np.linspace(-1,1,Nlx)
    lz = np.linspace(-1,1,Nlz)
    
    fig, axes = plt.subplots(1,2,figsize=(5,2.5),sharey=True)
    levels = []
    Lz,Lx=np.meshgrid(lz,lx)
    Z=(Lx+Lz)**2
    
    axes[1].contourf(Lx, Lz, Z, levels= np.linspace(0,4,50), cmap='jet')
    
    density = axes[0].contourf(lx, lz, np.transpose(SF), levels=np.linspace(0,4,50), cmap='jet')
    #density = axes.pcolor(lx, lz, np.transpose(SF), cmap='jet', norm=colors.SymLogNorm(linthresh=1e-4, linscale=0.1, vmin=0, vmax=4.0))

    
    #axes[0].set_aspect(1)
    axes[0].set_xticks([-1., 0, 1.0])
    axes[0].set_yticks([-1.0, 0, 1.0])
    axes[0].set_xlabel('$l_x$')
    axes[0].set_ylabel('$l_z$')
    axes[0].tick_params(axis='x', which='major', pad=10)
    axes[0].tick_params(axis='y', which='major', pad=10)
    axes[0].set_title(r"$(\mathrm{a})$ $S_2^{\theta}(l_x,l_z)$")#, pad=10)

    axes[1].set_title(r"$(\mathrm{b})$ $(l_x + l_z)^2$")
    
    #axes[1].set_aspect(1)
    axes[1].set_xticks([-1., 0, 1.0])
    axes[1].set_yticks([-1.0, 0, 1.0])
    axes[1].set_xlabel('$l_x$')
    #axes[1].set_ylabel('$l_z$')
    axes[1].tick_params(axis='x', which='major', pad=10)
    axes[1].tick_params(axis='y', which='major', pad=10)
    axes[0].title.set_position([.5, 1.05])
    axes[1].title.set_position([.5, 1.05])
    cb1 = fig.colorbar(density, fraction=0.05, ax=axes[0],ticks=[0,1,2,3,4])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb1.ax.tick_params(labelsize=A)
    
    cb2 = fig.colorbar(density, fraction=0.05, ax=axes[1],ticks=[0,1,2,3,4])#, ticks=[1e-4, 1e-2, 1e0]) ###### TICKS FOR THE COLORBARS ARE DEFINED HERE
    cb2.ax.tick_params(labelsize=A)
    fig.tight_layout()
    plt.savefig("SF_scalar.png", dpi=600)
    
    plt.show()


plotSF_r2()
plot_SF_density()
