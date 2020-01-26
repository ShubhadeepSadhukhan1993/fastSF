from __future__ import division
import h5py
from pylab import*
from matplotlib import ticker, colors
from mpl_toolkits.mplot3d import axes3d
from matplotlib import gridspec
from mayavi import mlab
from scipy.optimize import curve_fit

plt.rcParams['xtick.major.size'] = 9
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.major.size'] = 9
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['ytick.minor.width'] = 1


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


def hdf5_reader_plane(filename,dataset,N):
	file_V1_read = h5py.File(filename)
	dataset_V1_read = file_V1_read["/"+dataset]
	V1=dataset_V1_read[:,:,N]
	return V1
	
def f(x, c, alpha):
	return c*(x**alpha)


A=20
font = {'family' : 'serif', 'weight' : 'normal', 'size' : A}
plt.rc('font', **font)
B = 19
N=512


x=linspace(0,1,N)
y=linspace(0,1,N)
X,Y=meshgrid(x,y)


#print T[90,100,500]
#print T[:,255,255]


	
def plotSF_r2():
	S_r_Pll = hdf5_reader1D("out/SF.h5","SF")[:,2] #
	S_r_Pll_Shaheen = hdf5_reader1D("out/SF_Shaheen.h5","SF")[:,2] #
	#S_r_Pll_old2 = hdf5_reader1D("out/SF_old2.h5","SF")[:,1] #
	#S_r_Pll_old3 = hdf5_reader1D("out/SF_old3.h5","SF")[:,1] #
	r = zeros_like(S_r_Pll) #
	for i in range(len(S_r_Pll)): #
		r[i]=sqrt(3)*i/len(S_r_Pll)     #
	
	llim = 2
	ulim = 25
	fig, axes = plt.subplots(figsize = (7.5, 5.5))
	
	axes.plot(r[1:len(r)], -S_r_Pll[1:len(r)], color='blue', lw=3, label=r"$SF_{\pll}$") #
	axes.plot(r[1:len(r)], -S_r_Pll_Shaheen[1:len(r)], color='green', lw=2, label=r"$SF_{\pll}$") 
	#axes.plot(r2[1:len(r2)], S_r_Pll_old2[1:len(r2)], color='red', lw=1.5, label=r"$SF_{\pll}$")
	#axes.plot(r2[1:len(r2)], S_r_Pll_old3[1:len(r2)], color='black', lw=1, label=r"$SF_{\pll}$") 
	#popt, pcov = curve_fit(f, r2[llim:ulim], S_r_Pll2[llim:ulim])##
	#y = (popt[0]+1)*(r2**popt[1]) ##
	#axes.plot(r2[llim+1:ulim+11], y[llim+1:ulim+11], color='brown', lw=2) ##
	#print popt[1]
	axes.set_xlabel('$r$')
	axes.set_ylabel('$\mathrm{S_3^{(u)}}(r)$')
	axes.set_xscale('log')
	axes.set_yscale('log')
	axes.set_xlim(0.07, 2)
	#axes.set_ylim(1e-6, 1e-3)
	#axes.annotate(r"$\zeta_3 = %.2f$" %(popt[1]), fontsize=A*0.9,  xy=(0.01, 3.8e-5), xytext=(0.006, 0.00038), arrowprops=dict(facecolor='black', arrowstyle="->", linewidth = 1.5))
	#axes.set_ylim()
	fig.tight_layout()
	show()
	



#ploTDensity()
#plotVDens()
plotSF_r2()
#plotSFT()
#plotSF_r_T()
#plotSF_scalar()
