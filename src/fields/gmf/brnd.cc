#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include "param.h"
#include "grid.h"
#include "brnd.h"
#include "breg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

vec3_t<double> Brnd::get_brnd(const vec3_t<double> &pos, Grid_brnd *grid){
    if(grid->read_permission){
        return read_grid(pos,grid);
    }
    // if no specific random field model is called
    // base class will return zero vector
    else{
        return vec3_t<double> {0.,0.,0.};
    }
}

vec3_t<double> Brnd::read_grid(const vec3_t<double> &pos, Grid_brnd *grid){
    double tmp {(grid->nx-1)*(pos.x-grid->x_min)/(grid->x_max-grid->x_min)};
    if (tmp<0 or tmp>grid->nx-1) { return vec3_t<double> {0.,0.,0.};}
    decltype(grid->nx) xl {(std::size_t)floor(tmp)};
    const double xd {tmp - xl};
    
    tmp = (grid->ny-1)*(pos.y-grid->y_min)/(grid->y_max-grid->y_min);
    if (tmp<0 or tmp>grid->ny-1) { return vec3_t<double> {0.,0.,0.};}
    decltype(grid->nx) yl {(std::size_t)floor(tmp)};
    const double yd {tmp - yl};
    
    tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
    if (tmp<0 or tmp>grid->nz-1) { return vec3_t<double> {0.,0.,0.};}
    decltype(grid->nx) zl {(std::size_t)floor(tmp)};
    const double zd {tmp - zl};
#ifdef DEBUG
    if(xd<0 or yd<0 or zd<0 or xd>1 or yd>1 or zd>1){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG VALUE"<<endl;
        exit(1);
    }
#endif
    vec3_t<double> b_vec3;
    //trilinear interpolation
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        std::size_t idx1 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        std::size_t idx2 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1)};
        vec3_t<double> i1 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        vec3_t<double> i2 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        vec3_t<double> j1 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        vec3_t<double> j2 {grid->fftw_b_x[idx1]*(1.-zd) + grid->fftw_b_x[idx2]*zd,
            grid->fftw_b_y[idx1]*(1.-zd) + grid->fftw_b_y[idx2]*zd,
            grid->fftw_b_z[idx1]*(1.-zd) + grid->fftw_b_z[idx2]*zd};
        // interpolate along y direction, two interpolated vectors
        vec3_t<double> w1 {i1*(1.-yd) + i2*yd};
        vec3_t<double> w2 {j1*(1.-yd) + j2*yd};
        // interpolate along x direction
        b_vec3 = w1*(1.-xd) + w2*xd;
    }
    // on the boundary
    else {
        std::size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        b_vec3 = vec3_t<double> {grid->fftw_b_x[idx],grid->fftw_b_y[idx],grid->fftw_b_z[idx]};
    }
#ifdef DEBUG
    if (b_vec3.Length()>50.*CGS_U_muGauss) {
        cerr<<"WAR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<" too strong field at xl="<<xl<<", yl="<<yl<<", zl="<<zl<<endl
        <<" b_strength = "<< b_vec3.Length()/CGS_U_muGauss <<" microGauss"<<endl;
        exit(1);
    }
#endif
    return b_vec3;
}

void Brnd::write_grid(Param *, Breg *, Grid_breg *, Grid_brnd *){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// END
