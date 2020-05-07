"""
A simple mask map preparation script, healpy is required
"""

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

# be careful about mask map data type
DATATYPE = np.float64

def mask_map_prod(_nside,_clon,_clat,_sep):
    _c = np.pi/180
    def hav(_theta):
        return 0.5-0.5*np.cos(_theta)
    _tmp = np.ones(hp.nside2npix(_nside),dtype=DATATYPE)
    _count = 0
    for _ipix in range(len(_tmp)):
        _lon,_lat = hp.pix2ang(_nside,_ipix,lonlat=True)
        # iso-angle separation
        if((hav(np.fabs(_clat-_lat)*_c)+np.cos(_clat*_c)*np.cos(_lat*_c)*hav(np.fabs(_clon-_lon)*_c))>hav(_sep*_c)):
            _count = _count + 1
            _tmp[_ipix] = False
    print ('mask map sky-faction: %s' % str(1 - _count/(12*_nside*_nside)))
    return _tmp


if __name__=='__main__':
    example = mask_map_prod(32, 20, 30, 40)
    hp.mollview(example, min=0, max=1, title='example mask')
    plt.savefig('mask_illustration.pdf')
    example.tofile('mask.bin')
