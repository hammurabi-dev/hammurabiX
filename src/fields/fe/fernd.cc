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
#include "fe.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

double FErnd::get_fernd(const vec3 &pos,Grid_fernd *grid){
    if(grid->read_permission){
        return read_grid(pos,grid);
    }
    // if no specific turbulent field model is called
    // base class will return zero vector
    else{
        return 0.;
    }
}

// since we are using rms normalization
// p0 is hidden and not affecting anything
double FErnd::fe_spec(const double &k, Pond *par){
    //units fixing
    const double p0 {par->fernd_iso.rms}; //pccm
    const double k0 {par->fernd_iso.k0};
    const double a0 {par->fernd_iso.a0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // avoid nan
    if(k<=0.){
        return 0.;
    }
    // modified Giacalone-Jokipii model
    //const double P = 2*p0*pow(k/k0,a0)/(1.+pow(k/k0,a0-a1));
    // power law
    double P {0.};
    if(k>k0){
        P = p0*pow(k/k0,a0);
    }
    return P*unit;
}

double FErnd::rescal_fact(const vec3 &pos, Pond *par){
    const double r_cyl {(sqrt(pos.x*pos.x+pos.y*pos.y) - fabs(par->SunPosition.x))/CGS_U_kpc};
    const double z {(fabs(pos.z) - fabs(par->SunPosition.z))/CGS_U_kpc};
    const double r0 {par->fernd_scal.r0};
    const double z0 {par->fernd_scal.z0};
    if(r_cyl==0. or z==0) {return 1.;}
    else{
        return exp(-r_cyl/r0)*exp(-z/z0);
    }
}

double FErnd::read_grid(const vec3 &pos, Grid_fernd *grid){
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
        unsigned long int idx1 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
        unsigned long int idx2 {toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl+1)};
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
        unsigned long int idx{toolkit::Index3d(grid->nx,grid->ny,grid->nz,xl,yl,zl)};
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

void FErnd::write_grid_ani(Pond *,FE *,Grid_fe *,Grid_fernd *){
    cerr<<"WAR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// END
