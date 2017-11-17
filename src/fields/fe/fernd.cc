#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include "pond.h"
#include "grid.h"
#include "fernd.h"
#include "fereg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

double FErnd::get_fernd(const vec3_t<double> &pos,Grid_fernd *grid){
    if(grid->read_permission){
        return read_grid(pos,grid);
    }
    // if no specific turbulent field model is called
    // base class will return zero vector
    else{
        return 0.;
    }
}

double FErnd::read_grid(const vec3_t<double> &pos, Grid_fernd *grid){
    double tmp {(grid->nx-1)*(pos.x-grid->x_min)/(grid->x_max-grid->x_min)};
    if (tmp<0 or tmp>grid->nx-1) { return 0.;}
    decltype(grid->nx) xl {(unsigned int)floor(tmp)};
    const double xd {tmp - xl};
    
    tmp = (grid->ny-1)*(pos.y-grid->y_min)/(grid->y_max-grid->y_min);
    if (tmp<0 or tmp>grid->ny-1) { return 0.;}
    decltype(grid->nx) yl {(unsigned int)floor(tmp)};
    const double yd {tmp - yl};
    
    tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
    if (tmp<0 or tmp>grid->nz-1) { return 0.;}
    decltype(grid->nx) zl {(unsigned int)floor(tmp)};
    const double zd {tmp - zl};
#ifndef NDEBUG
    if(xd<0 or yd<0 or zd<0 or xd>1 or yd>1 or zd>1){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG VALUE"<<endl;
        exit(1);
    }
#endif
    double density;
    //trilinear interpolation
    if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
        //interpolate along z direction, there are four interpolated vectors
        std::size_t idx1 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        std::size_t idx2 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1)};
        double i1 {grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl+1,zl+1);
        double i2 {grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl,zl+1);
        double j1 {grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd};
        
        idx1 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl);
        idx2 = toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl+1,yl+1,zl+1);
        double j2 {grid->fftw_fe[idx1]*(1.-zd) + grid->fftw_fe[idx2]*zd};
        // interpolate along y direction, two interpolated vectors
        double w1 {i1*(1.-yd) + i2*yd};
        double w2 {j1*(1.-yd) + j2*yd};
        // interpolate along x direction
        density = w1*(1.-xd) + w2*xd;
    }
    // on the boundary
    else {
        std::size_t idx{toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        density = grid->fftw_fe[idx];
    }
    return density;
}

void FErnd::write_grid_iso(Pond *,Grid_fernd *){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// END
