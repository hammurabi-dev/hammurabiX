'''
@brief: this is a shabby code for prep mock data
	from fixed-parameter simulated observables

@notice:
	for captions please check my draft
'''

import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

# IQU pol. SMOOTHING
def smooth_iqu(maps,n_out,width):
  	alms = hp.map2alm(maps,lmax=3*n_out-1,pol=True)
    	return hp.alm2map(alms,n_out,pol=True,fwhm=width)

# masking
def mask_out(nside,ipix,lat_cut):
    	l,b = hp.pix2ang(nside,ipix,lonlat=True)
    	cue = 0
    	# nps
    	if(b>0. and b<80.):
        	if(l<50. or l>310.):
            		cue = 1
    	# latitude
    	if(b<lat_cut):
        	cue = 1
    	return cue


# pipeline
#--------------------------

t = hp.read_map('iqu_sync.fits',field=0)
q = hp.read_map('iqu_sync.fits',field=1)
u = hp.read_map('iqu_sync.fits',field=2)
i_map = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=0,hdu=1)
q_map = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=1,hdu=1)
u_map = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=2,hdu=1)

(i_smooth,q_smooth,u_smooth) = smooth_iqu((i_map,q_map,u_map),16,10*np.pi/180.)

# cut the maps and calculate RMS in masked maps
i_cut = list()
q_cut = list()
u_cut = list()
i_smooth_cut = list()
q_smooth_cut = list()
u_smooth_cut = list()
for i in range(0,np.size(q)):
	if(mask_out(16,i,60)):
		continue
	i_cut.append(t[i])
	q_cut.append(q[i])
	u_cut.append(u[i])
	i_smooth_cut.append(i_smooth[i])
	q_smooth_cut.append(q_smooth[i])
	u_smooth_cut.append(u_smooth[i])
t = np.multiply(t,np.std(i_smooth_cut)/np.std(i_cut))
q = np.multiply(q,np.std(q_smooth_cut)/np.std(q_cut))
u = np.multiply(u,np.std(u_smooth_cut)/np.std(u_cut))

#if need real turbulence in maps
q_map = hp.ud_grade(q_map,nside_out=16)
u_map = hp.ud_grade(u_map,nside_out=16)
q += q_map - q_smooth
u += u_map - u_smooth


#hp.mollview(t,norm='hist')
#plt.show()

ii = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=0,hdu=2)
qq = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=1,hdu=2)
uu = hp.read_map('wmap_band_iqumap_r9_9yr_K_v5.fits',field=3,hdu=2)
# downgrade variance
for j in range(0,np.size(i_map)):
	ii[j] = (1.435**2/ii[j])
	qq[j] = (1.435**2/qq[j])
	uu[j] = (1.435**2/uu[j])
ii = hp.ud_grade(ii,nside_out=16,power=2)
qq = hp.ud_grade(qq,nside_out=16,power=2)
uu = hp.ud_grade(uu,nside_out=16,power=2)


hp.write_map('wmap_mock_iqu_n16.fits',[t,q,u,ii,qq,uu])