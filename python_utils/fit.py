# -*- coding: utf-8 -*-
'''
@brief:
    	parameter constraining module
    	wrap up hammurabi executable by wrapper.py
    	we use PyMultiNest package as sampler

@author:
    	Jiaxin Wang, SISSA

@email:
    	jiwang@sissa.it

@notice:
    	1. use 'multinest_marginals.py out/' to observe results
    	2. check and modify specific settings/processes
    	in order to meet particular task
    	3. If you encounter "MKL FATAL ERROR: Cannot load neither libmkl_avx.so nor libmkl_def.so"
    	export LD_PRELOAD="/path/to/libmkl_core.so:/path/to/libmkl_sequential.so"
	this problem may happen when running python-C++ code
'''

import os
import numpy as np
import healpy as hp

import json
import pymultinest as pmt

import wrapper as wp
#import subprocess
#import xml.etree.ElementTree as et
#import random as rnd
#import tempfile as tf

# read WMAP-9yr K-band Q,U,QQ,UU
def read_wmap_file(nside):
	path = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))),'wmap_mock_iqu_n16.fits')
	Q = hp.read_map(path,verbose=False,field=1)
	U = hp.read_map(path,verbose=False,field=2)
	QQ = hp.read_map(path,verbose=False,field=4)
	UU = hp.read_map(path,verbose=False,field=5)
	return (Q,U,QQ,UU)


# pixel mask
def mask(nside,pix):
	l,b = hp.pix2ang(nside,pix,lonlat=True)
	cue = 0
   	# masking out low latitudes
	if b<60:
		cue = 1
	
    	# masking out NPS
	if b>0 and b<80:
		if l<50 or l>310:
			cue = 1
	
	return cue


# loglikelihood
def lnlike(theta,nside,Q,U,QQ,UU):
  	# seperate additional variance for likelihood
	V = 10.**theta[-1]
	
	hampy = wp.Wrapper()
	
	hampy.mod_par(keys=['Galaxy','MagneticField','Regular','WMAP','b0'],tag='value',attrib=str(theta[0]))
	hampy.mod_par(keys=['Galaxy','MagneticField','Regular','WMAP','psi0'],tag='value',attrib=str(theta[1]))
	hampy.mod_par(keys=['Galaxy','MagneticField','Regular','WMAP','psi1'],tag='value',attrib=str(theta[2]))
	hampy.mod_par(keys=['Galaxy','MagneticField','Regular','WMAP','chi0'],tag='value',attrib=str(theta[3]))
	hampy.mod_par(keys=['CRE','Analytic','alpha'],tag='value',attrib=str(theta[4]))

	hampy.call()
	
	loglike = 0.
	for i in range(0,12*nside**2):
		if mask(nside,i):
			continue
		else:
			# simulation in K, wmap in mK
			loglike += -0.5*(Q[i]-hampy.sim_map['Q'][i]*1000.)**2/(QQ[i]+V)
			loglike += -0.5*np.log( 2.*np.pi*(QQ[i]+V) )
			loglike += -0.5*(U[i]-hampy.sim_map['U'][i]*1000.)**2/(UU[i]+V)
			loglike += -0.5*np.log( 2.*np.pi*(UU[i]+V) )
	
	return loglike


# mcmc routine
def mult():
   	print 'INFO: analysis start with PyMultiNest'
    	nside = 16
    	# import observational maps
    	Q,U,QQ,UU = read_wmap_file(nside)
	
    	# initialize sampler
    	parameters = ['b_0','psi0','psi1','chi0','alpha','V']
    	n_params = len(parameters)
    	# define priors
    	def prior(cube,ndim,nparams):
        	#b0
        	cube[0] *= 10.
        	#psi0
        	cube[1] *= 50.
        	#psi1
        	cube[2] *= 5.
        	#chi0
        	cube[3] *= 50.
        	#alpha
        	cube[4] = cube[4]*3. + 1.
        	#V
        	cube[5] = cube[5]*8. - 12.
    
    	# define posterior
    	def loglike(cube,ndim,nparams):
        	theta = cube[0:n_params]
        	return lnlike(theta,nside,Q,U,QQ,UU)

    	# call PyMultiNest main function
    	path = os.path.join(os.path.abspath(os.path.dirname(os.path.realpath(__file__))),'out')
	if not os.path.isdir(path):
    		os.mkdir(path)
	
	json.dump(parameters,open('out/params.json','w'))
    	pmt.run(loglike,prior,n_params,outputfiles_basename='out/',resume=False,verbose=True)
    	
    	print 'INFO: analysis finished'


if __name__ == '__main__':
	mult()
