#include <cmath>
#include <cassert>

#include <vec3.h>

#include <param.h>
#include <grid.h>
#include <breg.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

vec3_t<double> Breg::get_breg (const vec3_t<double> &pos,
                               const Param *par,
                               const Grid_breg *grid) const{
    if (par->grid_breg.read_permission){
        return read_grid (pos,
                          par,
                          grid);
    }
    else if (par->grid_breg.build_permission){
        return breg (pos,
                     par);
    }
    else {
        return vec3_t<double> {0.,0.,0.};
    }
}

// if no specified field model is built
// Breg object link directly here and return null field when invoked
vec3_t<double> Breg::breg (const vec3_t<double> &,
                           const Param *) const{
    return vec3_t<double> {0.,0.,0.};
}

vec3_t<double> Breg::read_grid (const vec3_t<double> &pos,
                                const Param *par,
                                const Grid_breg *grid) const{
    double tmp {(par->grid_breg.nx-1)*(pos.x-par->grid_breg.x_min)/(par->grid_breg.x_max-par->grid_breg.x_min)};
    if (tmp<0 or tmp>par->grid_breg.nx-1) {return vec3_t<double>(0.,0.,0.);}
    decltype(par->grid_breg.nx) xl {(std::size_t)floor(tmp)};
    const double xd {tmp - xl};
    tmp = (par->grid_breg.ny-1)*(pos.y-par->grid_breg.y_min)/(par->grid_breg.y_max-par->grid_breg.y_min);
    if (tmp<0 or tmp>par->grid_breg.ny-1) { return vec3_t<double>(0.,0.,0.);}
    decltype(par->grid_breg.nx) yl {(std::size_t)floor(tmp)};
    const double yd {tmp - yl};
    tmp = (par->grid_breg.nz-1)*(pos.z-par->grid_breg.z_min)/(par->grid_breg.z_max-par->grid_breg.z_min);
    if (tmp<0 or tmp>par->grid_breg.nz-1) { return vec3_t<double>(0.,0.,0.);}
    decltype(par->grid_breg.nx) zl {(std::size_t)floor(tmp)};
    const double zd {tmp - zl};
    assert(xd>=0 and yd>=0 and zd>=0 and xd<=1 and yd<=1 and zd<=1);
    vec3_t<double> b_vec3;
    // trilinear interpolation
    if (xl+1<par->grid_breg.nx and yl+1<par->grid_breg.ny and zl+1<par->grid_breg.nz) {
        // interpolate along z direction, there are four interpolated vectors
        std::size_t idx1 {toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl,yl,zl)};
        std::size_t idx2 {toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl,yl,zl+1)};
        vec3_t<double> i1 = {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        idx1 = toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl,yl+1,zl);
        idx2 = toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl,yl+1,zl+1);
        vec3_t<double> i2 {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        idx1 = toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl+1,yl,zl);
        idx2 = toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl+1,yl,zl+1);
        vec3_t<double> j1 {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        idx1 = toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl+1,yl+1,zl);
        idx2 = toolkit::index3d(par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl+1,yl+1,zl+1);
        vec3_t<double> j2 {grid->bx[idx1]*(1.-zd) + grid->bx[idx2]*zd,
            grid->by[idx1]*(1.-zd) + grid->by[idx2]*zd,
            grid->bz[idx1]*(1.-zd) + grid->bz[idx2]*zd};
        // interpolate along y direction, two interpolated vectors
        vec3_t<double> w1 {i1*(1.-yd) + i2*yd};
        vec3_t<double> w2 {j1*(1.-yd) + j2*yd};
        // interpolate along x direction
        b_vec3 = w1*(1.-xd) + w2*xd;
    }
    // no interpolation
    else {
        std::size_t idx {toolkit::index3d (par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,xl,yl,zl)};
        b_vec3 = vec3_t<double> {grid->bx[idx],grid->by[idx],grid->bz[idx]};
    }
    assert (b_vec3.Length()>1e+5*CGS_U_muGauss);
    return b_vec3;
}

void Breg::write_grid (const Param *par,
                       Grid_breg *grid) const{
    assert(par->grid_breg.write_permission);
    vec3_t<double> gc_pos, tmp_vec;
    double lx {par->grid_breg.x_max-par->grid_breg.x_min};
    double ly {par->grid_breg.y_max-par->grid_breg.y_min};
    double lz {par->grid_breg.z_max-par->grid_breg.z_min};
    for (decltype(par->grid_breg.nx) i=0;i!=par->grid_breg.nx;++i){
        gc_pos.x = i*lx/(par->grid_breg.nx-1) + par->grid_breg.x_min;
        for (decltype(par->grid_breg.ny) j=0;j!=par->grid_breg.ny;++j){
            gc_pos.y = j*ly/(par->grid_breg.ny-1) + par->grid_breg.y_min;
            for (decltype(par->grid_breg.nz) k=0;k!=par->grid_breg.nz;++k){
                gc_pos.z = k*lz/(par->grid_breg.nz-1) + par->grid_breg.z_min;
                tmp_vec = breg (gc_pos,par);
                std::size_t idx {toolkit::index3d (par->grid_breg.nx,par->grid_breg.ny,par->grid_breg.nz,i,j,k)};
                grid->bx[idx] = tmp_vec.x;
                grid->by[idx] = tmp_vec.y;
                grid->bz[idx] = tmp_vec.z;
            }
        }
    }
}

// END
