# -*- coding: utf-8 -*-
'''
this snippet serves as a simple precision test of hammurabiX
assuming simplest field modles:


contributed by Jiaxin Wang
'''
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import FormatStrFormatter
import hampyx as ham
import scipy.special as sp

# CGS units
pi = 3.14159265358979
kpc = 3.0856775806e+21
eV = 1.60217733e-12
GeV = 1.e+9*eV
GHz = 1.e9
mG = 1.e-6
q = 4.8032068e-10
mc2 = 0.51099907e-3*GeV
cc = 2.99792458e+10
mc = mc2/cc
kB = 1.380622e-16
h = 6.626075540e-27

# convert T_br into T_cmb
def t_conv(t,freq):
	p = (h*freq*GHz)/(kB*2.725)
	return t*(np.exp(p)-1.)**2/(p**2*np.exp(p))

# emissivity J_tot (in mK_cmb)
def J_tot(theta):
	b,r,l0,alpha,je,freq = theta
	# CRE normalization
	gamma10 = 10.*GeV/mc2+1.0
	beta10 = np.sqrt(1.-1./gamma10)
	rslt = je*4.*pi*mc*(gamma10**alpha)/((1.e+4)*GeV*beta10)
	# following Ribiki-Lightman eq(6.36)
	rslt *= np.sqrt(3)*(q**3)*b*mG/(mc2*(alpha+1))
	rslt *= sp.gamma(0.25*alpha+19./12.)*sp.gamma(0.25*alpha-1./12.)
	rslt *= (mc*2.*pi*freq*GHz/(3.*q*b*mG))**(0.5-0.5*alpha)
	# convert into mK_br
	rslt *= 1.e+3*cc*cc/(2.*kB*freq*GHz*freq*GHz)
	# convert into mK_cmb
	rslt = t_conv(rslt,freq)
	return rslt/(4.*pi)

# synchrotron Stokes I (in mK_cmb)
def i_th(theta,lon,lat):
	tmp = theta
	lon = (lon - tmp[2])*np.pi/180.
    	lat *= np.pi/180.
	# b_perp
    	tmp[0] *= np.sqrt(1+tmp[1]**2 - (np.cos(lat)*np.cos(lon) + tmp[1]*np.sin(lat))**2)
	return J_tot(tmp)

	
# emissivity J_pol (in mK_CMB)
def J_pol(theta):
	b,r,l0,alpha,je,freq = theta
	# CRE normalization
	gamma10 = 10.*GeV/mc2+1.0
	beta10 = np.sqrt(1.-1./gamma10)
	rslt = je*4.*pi*mc*(gamma10**alpha)/((1.e+4)*GeV*beta10)
	# following Ribiki-Lightman eq(6.36)
	rslt *= np.sqrt(3)*(q**3)*b*mG/(4.*mc2)
    	rslt *= sp.gamma(0.25*alpha+7./12.)*sp.gamma(0.25*alpha-1./12.)
     	rslt *= (mc*2.*pi*freq*GHz/(3.*q*b*mG))**(0.5-0.5*alpha)
  	# convert into mK_br
	rslt *= 1.e+3*cc*cc/(2.*kB*freq*GHz*freq*GHz)
 	# convert into mK_cmb
	rslt = t_conv(rslt,freq)
	return rslt/(4.*pi)	

# synchrotron Stokes Q (in mK_cmb)
def q_th(theta,lon,lat):
    	tmp = theta
	lon = (lon - tmp[2])*np.pi/180.
    	lat *= np.pi/180.
	# b_perp
	tmp[0] *= np.sqrt(1+tmp[1]**2 - (np.cos(lat)*np.cos(lon) + tmp[1]*np.sin(lat))**2)
	# cos(2X)
    	numerator = np.sin(lon)**2 - (tmp[1]*np.cos(lat)-np.sin(lat)*np.cos(lon))**2
    	denominator = np.sin(lon)**2 + (tmp[1]*np.cos(lat)-np.sin(lat)*np.cos(lon))**2
    	return J_pol(tmp)*numerator/denominator

# theoretical u
def u_th(theta,lon,lat):
	tmp = theta
	lon = (lon - tmp[2])*np.pi/180.
    	lat *= np.pi/180.
	# b_perp
    	tmp[0] *= np.sqrt(1+tmp[1]**2 - (np.cos(lat)*np.cos(lon) + tmp[1]*np.sin(lat))**2)
	# cos(2X)
    	numerator = 2.*np.sin(lon)*(tmp[1]*np.cos(lat)-np.sin(lat)*np.cos(lon))
    	denominator = np.sin(lon)**2 + (tmp[1]*np.cos(lat)-np.sin(lat)*np.cos(lon))**2
    	return J_pol(tmp)*numerator/denominator


# theoreitcal Fd
def fd_th(theta,lon,lat):
    	tmp = theta
    	lon = (lon - tmp[2])*np.pi/180.
    	lat *= np.pi/180.
	# b_parallel
	tmp[0] *= (np.cos(lat)*np.cos(lon) + tmp[1]*np.sin(lat))
    	return tmp[-1]*tmp[0]*mG*(-q**3/(2*pi*mc2**2))

# precision
def precision(_res):
	# observable controllers
	Nside = 8 # shouldn't affect precision
	Shell = 1 # shouldn't affect precision
	Res = _res # affecting precision
	radius = 4.0 # shouldn't affect precision
	
	# field controllers
	b0 = 6.0 #muG
	r = 0
	l0 = 0.0
	alpha = 3
	je = 0.25
	freq = 23 # affecting precision??
	ne = 0.01 #pccm
	
	# call hammurabiX wrapper
	obj = ham.hampyx()
	# assuming the xml file is not prepared
	obj.del_par(['Output','Sync'],'all')
	obj.add_par(['Output'],'Sync',{'cue':str(1),'freq':str(freq),'filename':'dumy'})
        obj.mod_par(['Output','DM'],{'cue':str(1)})
        obj.mod_par(['Output','Faraday'],{'cue':str(1)})
        # mute all field output/input
        obj.mod_par(['Fieldout','breg_grid'],{'read':str(0),'write',str(0)})
        obj.mod_par(['Fieldout','brnd_grid'],{'read':str(0),'write',str(0)})
        obj.mod_par(['Fieldout','fereg_grid'],{'read':str(0),'write',str(0)})
        obj.mod_par(['Fieldout','fernd_grid'],{'read':str(0),'write',str(0)})
        obj.mod_par(['Fieldout','cre_grid'],{'read':str(0),'write',str(0)})
	# calibrate simulation box
	obj.mod_par(['Grid','Shell','layer'],{'type':'auto'})
	obj.mod_par(['Grid','Shell','layer','auto','shell_num'],{'value':str(Shell)})
	obj.mod_par(['Grid','Shell','layer','auto','nside_min'],{'value':str(Nside)})
	obj.mod_par(['Grid','Shell','nside_sim'],{'value':str(Nside)})
	obj.mod_par(['Grid','Shell','ec_r_max'],{'value':str(radius)})
	obj.mod_par(['Grid','Shell','gc_r_max'],{'value':str(radius+9.)})
	obj.mod_par(['Grid','Shell','gc_z_max'],{'value':str(radius+1.)})
	obj.mod_par(['Grid','Shell','ec_r_res'],{'value':str(Res)})
        # fix GMF
        obj.mod_par(['MagneticField','Regular'],{'cue':str(1),'type':'Test'})
	obj.mod_par(['MagneticField','Regular','Test','b0'],{'value':str(b0)})
	obj.mod_par(['MagneticField','Regular','Test','l0'],{'value':str(l0)})
	obj.mod_par(['MagneticField','Regular','Test','r'],{'value':str(r)})
	obj.mod_par(['MagneticField','Random'],{'cue':str(0)})
	# fix FE
	obj.mod_par(['FreeElectron','Regular'],{'cue':str(1),'type':'Test'})
	obj.mod_par(['FreeElectron','Regular','Test','n0'],{'value':str(ne)})
	obj.mod_par(['FreeElectron','Regular','Test','r0'],{'value':str(radius)})
	obj.mod_par(['FreeElectron','Random'],{'cue':str(0)})
	# fix CRE
	obj.mod_par(['CRE'],{'type':'Test'})
	obj.mod_par(['CRE','Test','alpha'],{'value':str(alpha)})
	obj.mod_par(['CRE','Test','E0'],{'value':str(10.0)})
	obj.mod_par(['CRE','Test','j0'],{'value':str(je)})
	obj.mod_par(['CRE','Test','r0'],{'value':str(radius)})
	# call hammurabi executable
	obj.call()
	# (in mK_cmb)
	qsim = obj.sim_map['sync'][str(freq)]['Q']*1.e+3
    	usim = obj.sim_map['sync'][str(freq)]['U']*1.e+3
    	isim = obj.sim_map['sync'][str(freq)]['I']*1.e+3
    	fsim = obj.sim_map['fd']
	
	# get wmap 
    	#iwmap = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=0,hdu=1)
    	#iiwmap = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=0,hdu=2)
    	#for j in iiwmap:
    	#	j = 1.435**2/j
    	#iwmap = hp.ud_grade(iwmap,nside_out=Nside)
    	#iiwmap = hp.ud_grade(iiwmap,nside_out=Nside,power=2)
    	
	ith = np.zeros_like(isim)
	qth = np.zeros_like(qsim)
	uth = np.zeros_like(usim)
	fth = np.zeros_like(fsim)
    	for i in range(0,np.size(qth)):
        	l,b = hp.pix2ang(Nside,i,lonlat=True)
		ith[i] = i_th([b0,r,l0,alpha,je,freq],l,b)*radius*kpc
		qth[i] = q_th([b0,r,l0,alpha,je,freq],l,b)*radius*kpc
        	uth[i] = u_th([b0,r,l0,alpha,je,freq],l,b)*radius*kpc
        	fth[i] = fd_th([b0,r,l0,ne],l,b)*radius*kpc*1.e+4
	
	pith = np.sqrt(qth**2 + uth**2)
	pisim = np.sqrt(qsim**2+usim**2)	
	
	# get plot
	min_t = min(abs(isim-ith)/ith)
	min_p = min(abs(pisim-pith)/pith)
	
	max_t = max(abs(isim-ith)/ith)
	max_p = max(abs(pisim-pith)/pith)
	
	min_f = min(abs(fsim-fth)/fth)
	max_f = max(abs(fsim-fth)/fth)
	#hp.mollview((isim-ith)/ith)
	#plt.hist((isim-ith)/ith,bins=50,histtype='step',stacked=True,fill=False,label='$%E\pm%E$'%(mean_t,std_t))
	#plt.legend(loc=1)
	#plt.show()
	
	return (min_t,min_p,min_f,max_t,max_p,max_f)
	
# main
def main():
	from matplotlib.ticker import FuncFormatter
	num = 50
	x = np.linspace(0.1,100,num) #resolution pc
	z1_t = np.zeros_like(x)
	z2_t = np.zeros_like(x)
	z1_p = np.zeros_like(x)
	z2_p = np.zeros_like(x)
	z1_f = np.zeros_like(x)
	z2_f = np.zeros_like(x)
	for i in range(0,num):
		z1_t[i],z1_p[i],z1_f[i],z2_t[i],z2_p[i],z2_f[i] = precision(x[i]*0.001)
	
	fig, ax = plt.subplots()
	plt.plot(x,z2_t,'k-',label='total intensity')
	plt.plot(x,z2_p,'k-.',label='polarized intensity')
	plt.plot(x,z2_f,'k:',label='Faraday depth')
	def myfmt(x,pos):
		return '%0.2f' % (x*100)
	ax.yaxis.set_major_formatter(FuncFormatter(myfmt))
	plt.legend(loc=2,fontsize=13)
	#plt.ylim((-0.005,0.005))
	plt.ylabel('maximum relative error ($\%$)',fontsize=15)
	plt.xlabel('raidal integration resolution (pc)',fontsize=15)
	plt.savefig('precision.pdf')
	#plt.show()
		
if __name__ == '__main__':
	#precision(0.01)
	main()