#include <cmath>
#include <stdexcept>

#include <vec3.h>

#include <param.h>
#include <grid.h>
#include <brnd.h>
#include <breg.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <cassert>

vec3_t<double> Brnd::get_brnd (const vec3_t<double> &pos,
                               const Param *par,
                               const Grid_brnd *grid) const{
    if (par->grid_brnd.read_permission or par->grid_brnd.build_permission){
        return read_grid (pos,
                          par,
                          grid);
    }
    // if no specific random field model is called
    // base class will return zero vector
    else {
        return vec3_t<double> {0.,0.,0.};
    }
}

vec3_t<double> Brnd::read_grid (const vec3_t<double> &pos,
                                const Param *par,
                                const Grid_brnd *grid) const{
    double tmp {(par->grid_brnd.nx-1)*(pos.x-par->grid_brnd.x_min)/(par->grid_brnd.x_max-par->grid_brnd.x_min)};
    if (tmp<0 or tmp>par->grid_brnd.nx-1) { return vec3_t<double> {0.,0.,0.};}
    decltype(par->grid_brnd.nx) xl {(std::size_t)std::floor(tmp)};
    const double xd {tmp - xl};
    
    tmp = (par->grid_brnd.ny-1)*(pos.y-par->grid_brnd.y_min)/(par->grid_brnd.y_max-par->grid_brnd.y_min);
    if (tmp<0 or tmp>par->grid_brnd.ny-1) { return vec3_t<double> {0.,0.,0.};}
    decltype(par->grid_brnd.nx) yl {(std::size_t)std::floor(tmp)};
    const double yd {tmp - yl};
    
    tmp = (par->grid_brnd.nz-1)*(pos.z-par->grid_brnd.z_min)/(par->grid_brnd.z_max-par->grid_brnd.z_min);
    if (tmp<0 or tmp>par->grid_brnd.nz-1) { return vec3_t<double> {0.,0.,0.};}
    decltype(par->grid_brnd.nx) zl {(std::size_t)std::floor(tmp)};
    const double zd {tmp - zl};
    assert(xd>=0 and yd>=0 and zd>=0 and xd<=1 and yd<=1 and zd<=1);
    vec3_t<double> b_vec3;
    //trilinear interpolation
    if (xl+1<par->grid_brnd.nx and yl+1<par->grid_brnd.ny and zl+1<par->grid_brnd.nz){
        std::size_t idx1 {toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl,yl,zl)};
        std::size_t idx2 {toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl,yl,zl+1)};
        vec3_t<double> i1 {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        
        idx1 = toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl,yl+1,zl);
        idx2 = toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl,yl+1,zl+1);
        vec3_t<double> i2 {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        
        idx1 = toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl+1,yl,zl);
        idx2 = toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl+1,yl,zl+1);
        vec3_t<double> j1 {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        
        idx1 = toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl+1,yl+1,zl);
        idx2 = toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl+1,yl+1,zl+1);
        vec3_t<double> j2 {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        // interpolate along y direction, two interpolated vectors
        vec3_t<double> w1 {i1*(1.-yd) + i2*yd};
        vec3_t<double> w2 {j1*(1.-yd) + j2*yd};
        // interpolate along x direction
        b_vec3 = w1*(1.-xd) + w2*xd;
    }
    // on the boundary
    else {
        std::size_t idx {toolkit::index3d (par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,xl,yl,zl)};
        b_vec3 = vec3_t<double> {grid->bx[idx],grid->by[idx],grid->bz[idx]};
    }
    assert (b_vec3.Length()<1e+6*CGS_U_muGauss);
    return b_vec3;
}

void Brnd::write_grid (const Param *,
                       const Breg *,
                       const Grid_breg *,
                       Grid_brnd *) const{
    throw std::runtime_error("wrong inheritance");
}

// END
