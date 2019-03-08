"""
snippet for visualising random/turbulent B fields
"""

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
from matplotlib import cm
import sys

import hampyx as ham

matplotlib.use('Agg')


def read_box(infile, inshape, dtype='double'):
    """
    Read in a binary dump in either float or double format and convert to a numpy array of the given shape
    :param infile:
    :param inshape:
    :param dtype:
    :return:
    """
    with open(infile, 'rb') as f:
        #  C floats of 4 bytes aren't python floats which are really
        #  doubles at 8 bytes.  So spcify f4.
        if dtype == 'float':
            data = np.fromstring(f.read(), dtype='<f4')
        elif dtype == 'double':
            data = np.fromstring(f.read(), dtype='<f8')
        else:
            raise TypeError('unsupported type %s' % str(dtype))
    return np.reshape(data, inshape)


def plot2d(_type, _xsample, _ysample, _bins, _tag):
    """
    2D sample distribution plot
    :param _type: type or turbulent B generator
    :param _xsample: array of first sample points
    :param _ysample: array of second sample points
    :param _bins: number of bins in one direction
    :param _tag:
    :return:
    """
    x_std = np.std(_xsample)
    y_std = np.std(_ysample)
    
    h, xedges, yedges = np.histogram2d(_xsample, _ysample, bins=_bins)
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1, _bins))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((_bins, 1))
    x, y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    pdf = (h*(x_bin_sizes*y_bin_sizes))
    z = pdf.T
    levels = np.linspace(0, max(z.reshape(np.size(z))), 10)
    f2d, ax2d = plt.subplots()
    cs = ax2d.contourf(x, y, z, cmap=cm.Blues, levels=levels)
    ax2d.set_xlim((0, 2*max(x_std, y_std)))
    ax2d.set_ylim((0, 2*max(y_std, x_std)))
    ax2d.set_xlabel('$B_x$ ($\\mu G$)', fontsize=20)
    ax2d.set_ylabel('$B_y$ ($\\mu G$)', fontsize=20)
    xlimits, ylimits = ax2d.get_xlim(), ax2d.get_ylim()
    aspect = (xlimits[1] - xlimits[0]) / (ylimits[1] - ylimits[0])
    ax2d.set_aspect(aspect)
    proxy = [plt.Rectangle((0, 0), 1, 1, fc=pc.get_facecolor()[0])
             for pc in cs.collections]
    ax2d.legend(proxy, [_tag], fontsize=20)
    f2d.savefig(''.join((_type, '_2d_dist.pdf')))


def plotcorr(_type, _sx, _sy, _sz, _nx, _ny, _nz):
    """
    2-point correlation plot
    :param _type: type or turbulent B generator
    :param _sx: sample array in x direction
    :param _sy: sample array in y direction
    :param _sz: sample array in z direction
    :param _nx: grid size in x direction
    :param _ny: grid size in y direction
    :param _nz: grid size in z direction
    :return:
    """
    x = np.arange(1, _nx//3)
    y = np.arange(1, _ny//3)
    z = np.arange(1, _nz//3)
    ex = list()
    ey = list()
    ez = list()
    #
    for offset in x:
        tmp = list()
        for i in range(0, _nx-offset):
            for j in range(0, _ny-offset):
                for k in range(0, _nz-offset):
                    tmp.append(_sx[index(i+offset, j, k, _nx, _ny, _nz)]*_sx[index(i, j, k, _nx, _ny, _nz)]
                               + _sy[index(i+offset, j, k, _nx, _ny, _nz)]*_sy[index(i, j, k, _nx, _ny, _nz)]
                               + _sz[index(i+offset, j, k, _nx, _ny, _nz)]*_sz[index(i, j, k, _nx, _ny, _nz)])
        ex.append(np.mean(tmp))
    #
    for offset in y:
        tmp = list()
        for i in range(0, _nx-offset):
            for j in range(0, _ny-offset):
                for k in range(0, _nz-offset):
                    tmp.append(_sx[index(i, j+offset, k, _nx, _ny, _nz)]*_sx[index(i, j, k, _nx, _ny, _nz)]
                               + _sy[index(i, j+offset, k, _nx, _ny, _nz)]*_sy[index(i, j, k, _nx, _ny, _nz)]
                               + _sz[index(i, j+offset, k, _nx, _ny, _nz)]*_sz[index(i, j, k, _nx, _ny, _nz)])
        ey.append(np.mean(tmp))
    #
    for offset in z:
        tmp = list()
        for i in range(0, _nx-offset):
            for j in range(0, _ny-offset):
                for k in range(0, _nz-offset):
                    tmp.append(_sx[index(i, j, k+offset, _nx, _ny, _nz)]*_sx[index(i, j, k, _nx, _ny, _nz)]
                               + _sy[index(i, j, k+offset, _nx, _ny, _nz)]*_sy[index(i, j, k, _nx, _ny, _nz)]
                               + _sz[index(i, j, k+offset, _nx, _ny, _nz)]*_sz[index(i, j, k, _nx, _ny, _nz)])
        ez.append(np.mean(tmp))
    #
    fcr, axcr = plt.subplots()
    axcr.plot(x/(_nx+1.e-10), ex/ex[0], color='firebrick', linestyle='-', marker='D',
              linewidth=1.5, label='$\\mathcal{E}(\\hat{x})$')
    axcr.plot(y/(_ny+1.e-10), ey/ex[0], color='steelblue', linestyle='--', marker='o',
              linewidth=1.5, label='$\\mathcal{E}(\\hat{y})$')
    axcr.plot(z/(_nz+1.e-10), ez/ex[0], color='darkseagreen', linestyle=':', marker='v',
              linewidth=1.5, label='$\\mathcal{E}(\\hat{z})$')
    axcr.legend(loc=1, fontsize=20)
    axcr.set_ylabel('two point correlation', fontsize=20)
    axcr.set_xlabel('fractional displacement length', fontsize=20)
    axcr.set_xlim((0, 0.3))
    axcr.set_ylim((-0.1, 1.02))
    fcr.savefig(''.join((_type, '_2pc.pdf')))


def simulator(_fieldtype, _anisotropy, _nx, _ny, _nz, _n2l):
    """
    hammurabi X simulation setting
    :param _fieldtype: random/turbulent field type
    :param _anisotropy: level of anisotropy
    :param _nx: grid size in x,y,z direction
    :param _ny:
    :param _nz:
    :param _n2l: grid points per unit length (kpc)
    :return:
    """
    obj = ham.Hampyx()
    # assuming the xml file is not prepared
    # mute map outputs
    obj.del_par(['observable', 'sync'], 'all')
    obj.del_par(['observable', 'dm'], 'all')
    obj.del_par(['observable', 'faraday'], 'all')
    
    obj.mod_par(['grid', 'box_brnd', 'nx'], {'value': str(_nx)})
    obj.mod_par(['grid', 'box_brnd', 'ny'], {'value': str(_ny)})
    obj.mod_par(['grid', 'box_brnd', 'nz'], {'value': str(_nz)})
    obj.mod_par(['grid', 'box_brnd', 'x_min'], {'value': str(-_nx/(_n2l*0.5))})
    obj.mod_par(['grid', 'box_brnd', 'x_max'], {'value': str(_nx/(_n2l*0.5))})
    obj.mod_par(['grid', 'box_brnd', 'y_min'], {'value': str(-_ny/(_n2l*0.5))})
    obj.mod_par(['grid', 'box_brnd', 'y_max'], {'value': str(_ny/(_n2l*0.5))})
    obj.mod_par(['grid', 'box_brnd', 'z_min'], {'value': str(-_nz/(_n2l*0.5))})
    obj.mod_par(['grid', 'box_brnd', 'z_max'], {'value': str(_nz/(_n2l*0.5))})
    
    obj.mod_par(['fieldio', 'brnd'], {'read': str(0), 'write': str(1), 'filename': 'brnd.bin'})
    # use 'test' regular B field setting
    obj.mod_par(['magneticfield', 'regular'], {'cue': str(1), 'type': 'unif'})
    obj.mod_par(['magneticfield', 'regular', 'unif', 'b0'], {'value': str(2.0)})
    # set field along x direction
    obj.mod_par(['magneticfield', 'regular', 'unif', 'l0'], {'value': str(0.)})
    obj.mod_par(['magneticfield', 'regular', 'unif', 'r'], {'value': str(0.)})
    # choose turublent/random field type
    obj.mod_par(['magneticfield', 'random'], {'cue': str(1), 'type': _fieldtype, 'seed': str(0)})

    obj.mod_par(['magneticfield', 'random', 'global'], {'type': 'es'})
    obj.mod_par(['magneticfield', 'random', 'global', 'es', 'rms'], {'value': str(1.)})
    obj.mod_par(['magneticfield', 'random', 'global', 'es', 'k0'], {'value': str(0.5)})
    obj.mod_par(['magneticfield', 'random', 'global', 'es', 'a0'], {'value': str(1.7)})
    obj.mod_par(['magneticfield', 'random', 'global', 'es', 'rho'], {'value': str(_anisotropy)})
    obj.mod_par(['magneticfield', 'random', 'global', 'es', 'r0'], {'value': str(1000)})  # no scaling
    obj.mod_par(['magneticfield', 'random', 'global', 'es', 'z0'], {'value': str(1000)})  # no scaling

    obj.mod_par(['magneticfield', 'random', 'local'], {'type': 'mhd'})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'pa0'], {'value': str(1.)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'pf0'], {'value': str(_anisotropy)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'ps0'], {'value': str(_anisotropy)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'aa0'], {'value': str(1.7)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'af0'], {'value': str(1.5)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'as0'], {'value': str(1.7)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'k0'], {'value': str(0.5)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'beta'], {'value': str(0.1)})
    obj.mod_par(['magneticfield', 'random', 'local', 'mhd', 'ma'], {'value': str(0.5)})

    obj()


def index(_i, _j, _k, _nx, _ny, _nz):
    """
    auxiliray function
    calculate grid global index according to i,j,k position
    :param _i:
    :param _j:
    :param _k:
    :param _nx:
    :param _ny:
    :param _nz:
    :return:
    """
    return _i*_ny*_nz + _j*_nz + _k


def divergence(_i, _j, _k, _nx, _ny, _nz, _bx, _by, _bz):
    """
    auxiliary function
    calculate divergence of vector field
    :param _i:
    :param _j:
    :param _k:
    :param _nx:
    :param _ny:
    :param _nz:
    :param _bx:
    :param _by:
    :param _bz:
    :return:
    """
    div = (_bx[index(_i+1, _j, _k, _nx, _ny, _nz)] - _bx[index(_i-1, _j, _k, _nx, _ny, _nz)])/2
    div += (_by[index(_i, _j+1, _k, _nx, _ny, _nz)] - _by[index(_i, _j-1, _k, _nx, _ny, _nz)])/2
    div += (_bz[index(_i, _j, _k+1, _nx, _ny, _nz)] - _bz[index(_i, _j, _k-1, _nx, _ny, _nz)])/2
    return div


def single_run_div(_type, _rho, _nx, _ny, _nz, _n2l):
    """
    get divergence info from a single run, invoked in multi_run_ana
    :param _type:
    :param _rho:
    :param _nx:
    :param _ny:
    :param _nz:
    :param _n2l:
    :return:
    """
    mG = 1.e-6  # CGS_unit, micro Gauss
    
    simulator(_type, _rho, _nx, _ny, _nz, _n2l)  # setup and execute simulation
    brnd = read_box('brnd.bin', [_nx, _ny, _nz, 3])  # get information from dumped binary file
    
    brndx = np.reshape(brnd[:, :, :, 0], (_nx*_ny*_nz))/mG  # B_x
    brndy = np.reshape(brnd[:, :, :, 1], (_nx*_ny*_nz))/mG  # B_y
    brndz = np.reshape(brnd[:, :, :, 2], (_nx*_ny*_nz))/mG  # B_z
    """
    brndp = np.sqrt(brndy**2 + brndz**2)  # B_perp with respect to x direction
    """
    
    div_b = list()
    for i in range(1, _nx-1):
        for j in range(1, _ny-1):
            for k in range(1, _nz-1):
                div_b.append(divergence(i, j, k, _nx, _ny, _nz, brndx, brndy, brndz))
    
    return np.std(div_b)/np.sqrt(np.var(brndx)+np.var(brndy)+np.var(brndz))


def single_run_print(_type, _rho, _nx, _ny, _nz, _n2l):
    """
    print info from single run
    :param _type: type of turublent B generator
    :param _rho: level of anisotropy [0,1]
    :param _nx: grid size in each dimension
    :param _ny:
    :param _nz:
    :param _n2l: grid points per unit length (1kpc)
    :return:
    """
    mG = 1.e-6  # CGS_unit, micro Gauss
    simulator(_type, _rho, _nx, _ny, _nz, _n2l)  # setup and execute simulation
    brnd = read_box('brnd.bin', [_nx, _ny, _nz, 3])  # get information from dumped binary file
    brndx = np.reshape(brnd[:, :, :, 0], (_nx*_ny*_nz))/mG  # B_x
    brndy = np.reshape(brnd[:, :, :, 1], (_nx*_ny*_nz))/mG  # B_y
    brndz = np.reshape(brnd[:, :, :, 2], (_nx*_ny*_nz))/mG  # B_z
    
    print('x mean')
    print(np.mean(brndx))
    print('y mean')
    print(np.mean(brndy))
    print('z mean')
    print(np.mean(brndz))
    print('RMS')
    print(np.sqrt(np.var(brndx)+np.var(brndy)+np.var(brndz)))
    print('x to y RMS ratio')
    print(np.std(brndx)/np.std(brndy))
    print('x to z RMS ratio')
    print(np.std(brndx)/np.std(brndz))
    print('y to z RMS ratio')
    print(np.std(brndy)/np.std(brndz))
    
    div_b = list()
    for i in range(1, _nx-1):
        for j in range(1, _ny-1):
            for k in range(1, _nz-1):
                div_b.append(divergence(i, j, k, _nx, _ny, _nz, brndx, brndy, brndz))
    print('divergence mean')
    print(np.mean(div_b))
    print('divergence std')
    print(np.std(div_b))


def single_run_plot(_type, _rho, _nx, _ny, _nz, _n2l, _tag):
    """
    plot field anisotropy from single run
    :param _type: type of turublent B generator
    :param _rho: level of anisotropy [0,1]
    :param _nx: grid size in each dimension
    :param _ny:
    :param _nz:
    :param _n2l: grid points per unit length (1kpc)
    :param _tag: string tag for plot caption
    :return:
    """
    mG = 1.e-6  # CGS_unit, micro Gauss
    simulator(_type, _rho, _nx, _ny, _nz, _n2l)  # setup and execute simulation
    brnd = read_box('brnd.bin', [_nx, _ny, _nz, 3])  # get information from dumped binary file
    brndx = np.reshape(brnd[:, :, :, 0], (_nx*_ny*_nz))/mG  # B_x
    brndy = np.reshape(brnd[:, :, :, 1], (_nx*_ny*_nz))/mG  # B_y
    brndz = np.reshape(brnd[:, :, :, 2], (_nx*_ny*_nz))/mG  # B_z
    # projected distribution plot
    plot2d(_type, brndx, brndy, 100, _tag)
    # for anisotropy analysis
    plotcorr(_type, brndx, brndy, brndz, _nx, _ny, _nz)


def multi_run_print(_type):
    """
    print results (to disk) from multiple runs
    :param _type: type of turbulent B generator
    :return:
    """
    N = 5  # times of repeated execution with one setting
    x = np.arange(5, 60, 3)  # change of grid points per unit length
    y = np.zeros((np.size(x), N))
    for j in range(0, N):
        for i in range(0, np.size(x)):
            y[i, j] = single_run_div(_type, 0.5, 200, 200, 200, x[i])
    dat = np.zeros((np.size(x), 3))
    for i in range(0, np.size(x)):
        dat[i, 0] = x[i]
        dat[i, 1] = np.mean(y[i, :])
        dat[i, 2] = np.std(y[i, :])
    np.savetxt(''.join((_type, '_div.txt')), dat)


def multi_run_plot():
    """
    plot precision results from multiple runs
    we also provide a rough approximation on the shape of precision curves
    :return:
    """
    multi_run_print('global')
    multi_run_print('local')
    d1 = np.loadtxt('local_div.txt')
    d2 = np.loadtxt('global_div.txt')
    from scipy.optimize import curve_fit
    
    # fit curve
    def fit_func(x, a, b):
        return a*x**b
    
    params1 = curve_fit(fit_func, d1[:, 0], d1[:, 1])
    params2 = curve_fit(fit_func, d2[:, 0], d2[:, 1])
    
    [a1, b1] = params1[0]
    [a2, b2] = params2[0]
    
    print(a1, b1)
    print(a2, b2)
    fig, ax = plt.subplots(figsize=(9,9))
    ax.plot(np.arange(5, 60, 0.5), fit_func(np.arange(5, 60, 0.5), a1, b1),
             linestyle=':', color='firebrick', linewidth=2)
    ax.plot(np.arange(5, 60, 0.5), fit_func(np.arange(5, 60, 0.5), a2, b2),
             linestyle=':', color='steelblue', linewidth=2)
    ax.errorbar(d1[:, 0], d1[:, 1], yerr=d1[:, 2], fmt='.', color='firebrick', label='local generator')
    ax.errorbar(d2[:, 0]+1.5, d2[:, 1], yerr=d2[:, 2], fmt='d', color='steelblue', label='global generator')
    ax.tick_params(axis='both', which='major', labelsize='20')
    ax.xlabel('$N/L$ $(kpc^{-1})$', fontsize=20)
    ax.ylabel('$\\sigma_{div}/RMS$', fontsize=20)
    ax.legend(loc=1, fontsize=20)
    plt.savefig('div.pdf')


# main routine
# check captions above for the meaning of each function
if __name__ == '__main__':
    # do a single run and print info
    # single_run_print('local', 0.5, 200, 200, 200, 50)
    # single_run_print('global', 0.5, 200, 200, 200, 50)
    # single_run_plot('global', 0.2, 200, 200, 200, 50, '$\\rho=0.2$')
    # single_run_plot('local', 0.9, 200, 200, 200, 50, '$\\rho=0.9$')
    # single_run_plot('local', 0.2, 200, 200, 200, 50, '$P_{s,f}/P_A=0.2$')
    
    # do multiple runs and plot results may take 1hr
    multi_run_plot()
