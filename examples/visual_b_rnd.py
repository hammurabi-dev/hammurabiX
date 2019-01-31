'''
snippet for visualising random/turbulent B fields

written by Jiaxin Wang based on original work from Tess Jaffe
'''
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
from matplotlib import cm
import sys

import hampyx as ham

# Read in a binary dump in either float or double format and convert to a numpy array of the given shape
def read_box( infile, inshape, dtype="double" ):
   	with open(infile, 'rb') as f:
  		#  C floats of 4 bytes aren't python floats which are really
 		#  doubles at 8 bytes.  So spcify f4.
		if dtype=="float":
   			data=np.fromstring(f.read(),dtype='<f4')
 		elif dtype=="double":
  			data=np.fromstring(f.read(),dtype='<f8')
 		else:
 			print("ERROR:  don't know how to read type %s." % dtype)
   			exit(1)
 	return np.reshape(data,inshape)

# 2D sample distribution plot
# 0th argument: type or turbulent B generator
# 1st argument: array of first sample points
# 2nd argument: array of second sample points
# 3rd argument: number of bins in one direction
def plot2d(_type,_xsample,_ysample,_bins,_tag):
	x_mean = np.mean(_xsample)
	y_mean = np.mean(_ysample)
	x_std = np.std(_xsample)
	y_std = np.std(_ysample)
	
	H,xedges,yedges = np.histogram2d(_xsample,_ysample,bins=_bins)
	x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,_bins))
    	y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((_bins,1))
    	X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
	pdf = (H*(x_bin_sizes*y_bin_sizes))
    	Z = pdf.T
	levels = np.linspace(0,max(Z.reshape(np.size(Z))),10)
	f2d,ax2d = plt.subplots()
	CS = ax2d.contourf(X, Y, Z,cmap=cm.Blues, levels=levels)
	ax2d.set_xlim((0,2*max(x_std,y_std)))
	ax2d.set_ylim((0,2*max(y_std,x_std)))
	ax2d.set_xlabel('$B_x$ ($\mu G$)',fontsize=15)
	ax2d.set_ylabel('$B_y$ ($\mu G$)',fontsize=15)
	xLimits, yLimits = ax2d.get_xlim(),ax2d.get_ylim() #plt.axes().get_xlim(), plt.axes().get_ylim()
	aspect =  (xLimits[1] - xLimits[0]) / (yLimits[1] - yLimits[0])
	ax2d.set_aspect(aspect) #plt.axes().set_aspect(aspect)
	proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) 
    		for pc in CS.collections]
	ax2d.legend(proxy, [_tag],fontsize=15)
	f2d.savefig(''.join((_type,'_2d_dist.pdf')))
	

# 2-point correlation plot
# 0th argument: type or turbulent B generator
# 1st argument: sample array in x direction
# 2nd argument: sample array in y direction
# 3rd argument: sample array in z direction
# 4th - 6th argument: grid size in x,y,z direction
def plotcorr(_type,_sx,_sy,_sz,_nx,_ny,_nz):
	x = np.arange(1,_nx/3)
	y = np.arange(1,_ny/3)
	z = np.arange(1,_nz/3)
	ex = list()
	ey = list()
	ez = list()
	#
	for offset in x:
		tmp = list()
		for i in range(0,_nx-offset):
			for j in range(0,_ny-offset):
				for k in range(0,_nz-offset):
					tmp.append(_sx[index(i+offset,j,k,_nx,_ny,_nz)]*_sx[index(i,j,k,_nx,_ny,_nz)]\
					+_sy[index(i+offset,j,k,_nx,_ny,_nz)]*_sy[index(i,j,k,_nx,_ny,_nz)]\
					+_sz[index(i+offset,j,k,_nx,_ny,_nz)]*_sz[index(i,j,k,_nx,_ny,_nz)])
		ex.append(np.mean(tmp))
	#
	for offset in y:
		tmp = list()
		for i in range(0,_nx-offset):
			for j in range(0,_ny-offset):
				for k in range(0,_nz-offset):
					tmp.append(_sx[index(i,j+offset,k,_nx,_ny,_nz)]*_sx[index(i,j,k,_nx,_ny,_nz)]\
					+_sy[index(i,j+offset,k,_nx,_ny,_nz)]*_sy[index(i,j,k,_nx,_ny,_nz)]\
					+_sz[index(i,j+offset,k,_nx,_ny,_nz)]*_sz[index(i,j,k,_nx,_ny,_nz)])
		ey.append(np.mean(tmp))
	#
	for offset in z:
		tmp = list()
		for i in range(0,_nx-offset):
			for j in range(0,_ny-offset):
				for k in range(0,_nz-offset):
					tmp.append(_sx[index(i,j,k+offset,_nx,_ny,_nz)]*_sx[index(i,j,k,_nx,_ny,_nz)]\
					+_sy[index(i,j,k+offset,_nx,_ny,_nz)]*_sy[index(i,j,k,_nx,_ny,_nz)]\
					+_sz[index(i,j,k+offset,_nx,_ny,_nz)]*_sz[index(i,j,k,_nx,_ny,_nz)])
		ez.append(np.mean(tmp))
	#
	fcr,axcr=plt.subplots()
	axcr.plot(x/(_nx+1.e-10),ex/ex[0],color='firebrick',linestyle='-',marker='D',linewidth=1.5,label="$\mathcal{E}(\hat{x})$")
	axcr.plot(y/(_ny+1.e-10),ey/ex[0],color='steelblue',linestyle='--',marker='o',linewidth=1.5,label="$\mathcal{E}(\hat{y})$")
	axcr.plot(z/(_nz+1.e-10),ez/ex[0],color='darkseagreen',linestyle=':',marker='v',linewidth=1.5,label="$\mathcal{E}(\hat{z})$")
	axcr.legend(loc=1,fontsize=15)
	axcr.set_ylabel("two point correlation",fontsize=15)
	axcr.set_xlabel("fractional displacement length",fontsize=15)
	axcr.set_xlim((0,0.3))
	axcr.set_ylim((-0.1,1.02))
	fcr.savefig(''.join((_type,'_2pc.pdf')))
	
	
# hammurabi X simulation setting
# 1st argument: random/turbulent field type
# 2st argument: level of anisotropy
# 3rd-5th argument: grid size in x,y,z direction
# 6th argument: grid points per unit length (kpc)
def simulator(_fieldtype,_anisotropy,_nx,_ny,_nz,_n2l):
	obj = ham.hampyx()
	obj.mod_par(['Grid','Box_GMF','nx'],{'value':str(_nx)})
	obj.mod_par(['Grid','Box_GMF','ny'],{'value':str(_ny)})
	obj.mod_par(['Grid','Box_GMF','nz'],{'value':str(_nz)})
	obj.mod_par(['Grid','Box_GMF','x_min'],{'value':str(-_nx/(_n2l*0.5))})
	obj.mod_par(['Grid','Box_GMF','x_max'],{'value':str(_nx/(_n2l*0.5))})
	obj.mod_par(['Grid','Box_GMF','y_min'],{'value':str(-_ny/(_n2l*0.5))})
	obj.mod_par(['Grid','Box_GMF','y_max'],{'value':str(_ny/(_n2l*0.5))})
	obj.mod_par(['Grid','Box_GMF','z_min'],{'value':str(-_nz/(_n2l*0.5))})
	obj.mod_par(['Grid','Box_GMF','z_max'],{'value':str(_nz/(_n2l*0.5))})
	
	obj.mod_par(['Fieldout','brnd_grid'],{'read':str(0),'write':str(1),'filename','brnd.bin'})
	obj.mod_par(['MagneticField','Regular'],{'cue':str(1),'type':'Test'}) #use "Test" regular B field setting
	obj.mod_par(['MagneticField','Regular','Test','b0'],{'value':str(2.0)})
	obj.mod_par(['MagneticField','Regular','Test','l0'],{'value':str(0.)}) # set field along x direction
	obj.mod_par(['MagneticField','Regular','Test','r'],{'value':str(0.)})
	obj.mod_par(['MagneticField','Random'],{'cue':str(1),'type':_fieldtype,'seed':str(0)}) #choose turublent/random field type

	obj.mod_par(['MagneticField','Random','Global'],{'type':'ES'})
	obj.mod_par(['MagneticField','Random','Global','rms'],{'value':str(1.)})
	obj.mod_par(['MagneticField','Random','Global','k0'],{'value':str(0.5)})
	obj.mod_par(['MagneticField','Random','Global','a0'],{'value':str(1.7)})
	obj.mod_par(['MagneticField','Random','Global','rho'],{'value':str(_anisotropy)})
	obj.mod_par(['MagneticField','Random','Global','r0'],{'value':str(1000)}) # no scaling
	obj.mod_par(['MagneticField','Random','Global','z0'],{'value':str(1000)}) # no scaling

	obj.mod_par(['MagneticField','Random','Local','pa0'],{'value':str(1.)})
	obj.mod_par(['MagneticField','Random','Local','pf0'],{'value':str(_anisotropy)})
	obj.mod_par(['MagneticField','Random','Local','ps0'],{'value':str(_anisotropy)})
	obj.mod_par(['MagneticField','Random','Local','aa0'],{'value':str(1.7)})
	obj.mod_par(['MagneticField','Random','Local','af0'],{'value':str(1.5)})
	obj.mod_par(['MagneticField','Random','Local','as0'],{'value':str(1.7)})
	obj.mod_par(['MagneticField','Random','Local','k0'],{'value':str(0.5)})
	obj.mod_par(['MagneticField','Random','Local','beta'],{'value':str(0.1)})
	obj.mod_par(['MagneticField','Random','Local','ma'],{'value':str(0.5)})
	#sys.exit(1)
	obj.call(True)

# auxiliray function
# calculate grid global index according to i,j,k position
def index(_i,_j,_k,_nx,_ny,_nz):
	return _i*_ny*_nz + _j*_nz + _k

# auxiliary function
# calculate divergence of vector field
def divergence(_i,_j,_k,_nx,_ny,_nz,_bx,_by,_bz):
	#idx = _i*_ny*_nz + _j*_nz + _k 
	div = (_bx[index(_i+1,_j,_k,_nx,_ny,_nz)] - _bx[index(_i-1,_j,_k,_nx,_ny,_nz)])/2 
	div += (_by[index(_i,_j+1,_k,_nx,_ny,_nz)] - _by[index(_i,_j-1,_k,_nx,_ny,_nz)])/2 
	div += (_bz[index(_i,_j,_k+1,_nx,_ny,_nz)] - _bz[index(_i,_j,_k-1,_nx,_ny,_nz)])/2
	return div
	
# get divergence info from a single run, invoked in multi_run_ana
def single_run_div(_type,_rho,_nx,_ny,_nz,_n2l):
	mG = 1.e-6 # CGS_unit, micro Gauss
	
	simulator(_type,_rho,_nx,_ny,_nz,_n2l) # setup and execute simulation
	brnd = read_box("brnd.bin",[_nx,_ny,_nz,3]) # get information from dumped binary file
	
	brndx = np.reshape(brnd[:,:,:,0],(_nx*_ny*_nz))/mG # B_x
	brndy = np.reshape(brnd[:,:,:,1],(_nx*_ny*_nz))/mG # B_y
	brndz = np.reshape(brnd[:,:,:,2],(_nx*_ny*_nz))/mG # B_z
	#brndp = np.sqrt(brndy**2 + brndz**2) # B_perp with respect to x direction
	
	div_b = list()
	for i in range(1,_nx-1):
		for j in range(1,_ny-1):
			for k in range(1,_nz-1):
				div_b.append(divergence(i,j,k,_nx,_ny,_nz,brndx,brndy,brndz))
	
	return np.std(div_b)/np.sqrt(np.var(brndx)+np.var(brndy)+np.var(brndz))
	
# print info from single run
# 1st argument: type of turublent B generator
# 2nd argument: level of anisotropy [0,1]
# 3rd-5th argument: grid size in each dimension
# 6th argument: grid points per unit length (1kpc)
# last argument: string tag for plot caption
def single_run_print(_type,_rho,_nx,_ny,_nz,_n2l):
	mG = 1.e-6 # CGS_unit, micro Gauss
	simulator(_type,_rho,_nx,_ny,_nz,_n2l) # setup and execute simulation
	brnd = read_box("brnd.bin",[_nx,_ny,_nz,3]) # get information from dumped binary file
	brndx = np.reshape(brnd[:,:,:,0],(_nx*_ny*_nz))/mG # B_x
	brndy = np.reshape(brnd[:,:,:,1],(_nx*_ny*_nz))/mG # B_y
	brndz = np.reshape(brnd[:,:,:,2],(_nx*_ny*_nz))/mG # B_z
	
	print 'x mean'
	print np.mean(brndx)
	print 'y mean'
	print np.mean(brndy)
	print 'z mean'
	print np.mean(brndz)
	print 'RMS'
	print np.sqrt(np.var(brndx)+np.var(brndy)+np.var(brndz))
	print 'x to y RMS ratio'
	print np.std(brndx)/np.std(brndy)
	print 'x to z RMS ratio'
	print np.std(brndx)/np.std(brndz)
	print 'y to z RMS ratio'
	print np.std(brndy)/np.std(brndz)
	
	div_b = list()
	for i in range(1,_nx-1):
		for j in range(1,_ny-1):
			for k in range(1,_nz-1):
				div_b.append(divergence(i,j,k,_nx,_ny,_nz,brndx,brndy,brndz))
	print 'divergence mean'
	print np.mean(div_b)
	print 'divergence std'
	print np.std(div_b)

# plot field anisotropy from single run
# 1st argument: type of turublent B generator
# 2nd argument: level of anisotropy [0,1]
# 3rd-5th argument: grid size in each dimension
# 6th argument: grid points per unit length (1kpc)
# last argument: string tag for plot caption
# _tag is a string tag for plot
def single_run_plot(_type,_rho,_nx,_ny,_nz,_n2l,_tag):
	mG = 1.e-6 # CGS_unit, micro Gauss
	simulator(_type,_rho,_nx,_ny,_nz,_n2l) # setup and execute simulation
	brnd = read_box("brnd.bin",[_nx,_ny,_nz,3]) # get information from dumped binary file
	brndx = np.reshape(brnd[:,:,:,0],(_nx*_ny*_nz))/mG # B_x
	brndy = np.reshape(brnd[:,:,:,1],(_nx*_ny*_nz))/mG # B_y
	brndz = np.reshape(brnd[:,:,:,2],(_nx*_ny*_nz))/mG # B_z
	# projected distribution plot
	plot2d(_type,brndx,brndy,100,_tag)
	# for anisotropy analysis
	plotcorr(_type,brndx,brndy,brndz,_nx,_ny,_nz)


# print results (to disk) from multiple runs
# 1st argument: type of turbulent B generator
def multi_run_print(_type):
	N = 5 # times of repeated execution with one setting
	x = np.arange(5,60,3) # change of grid points per unit length
	y = np.zeros((np.size(x),N))
	for j in range(0,N):
		for i in range(0,np.size(x)):
			y[i,j] = single_run_div(_type,0.5,200,200,200,x[i])
	dat = np.zeros((np.size(x),3))
	for i in range(0,np.size(x)):
		dat[i,0] = x[i]
		dat[i,1] = np.mean(y[i,:])
		dat[i,2] = np.std(y[i,:])
	np.savetxt(''.join((_type,'_div.txt')),dat)

# plot precision results from multiple runs
# we also provide a rough approximation on the shape of precision curves
def multi_run_plot():
	multi_run_print('Global')
	multi_run_print('Local')
	d1 = np.loadtxt('Local_div.txt')
	d2 = np.loadtxt('Global_div.txt')
	from scipy.optimize import curve_fit
	
	# fit curve
	def fit_func(x, a, b):
		return a*x**b
	
	params1 = curve_fit(fit_func, d1[:,0], d1[:,1])
	params2 = curve_fit(fit_func, d2[:,0], d2[:,1])
	
	[a1, b1] = params1[0]
	[a2, b2] = params2[0]
	
	print a1, b1
	print a2, b2
	plt.plot(np.arange(5,60,0.5),fit_func(np.arange(5,60,0.5),a1,b1),linestyle=':',color='firebrick',linewidth=2)
	plt.plot(np.arange(5,60,0.5),fit_func(np.arange(5,60,0.5),a2,b2),linestyle=':',color='steelblue',linewidth=2)
	plt.errorbar(d1[:,0],d1[:,1],yerr=d1[:,2],fmt='.',color='firebrick',label='local generator')
	plt.errorbar(d2[:,0]+1.5,d2[:,1],yerr=d2[:,2],fmt='d',color='steelblue',label='global generator')
	plt.xlabel('$N/L$ $(kpc^{-1})$',fontsize=15)
	plt.ylabel('$\sigma_{div}/RMS$',fontsize=15)
	plt.legend(loc=1,fontsize=15)
	plt.savefig('div.pdf')
	
	
# main routine
# check captions above for the meaning of each function
if __name__=="__main__":
	# do a single run and print info
	#single_run_print('Local',0.5,200,200,200,50)
	#single_run_print('Global',0.5,200,200,200,50)
	single_run_plot('Global',0.2,200,200,200,50,'$\\rho=0.2$')
	#single_run_plot('Local',0.9,200,200,200,50,'$\\rho=0.9$')
	single_run_plot('Local',0.2,200,200,200,50,'$P_{s,f}/P_A=0.2$')
	
	# do multiple runs and plot results
	# may take 1hr
	#multi_run_plot()



