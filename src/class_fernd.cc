#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include "class_pond.h"
#include "class_grid.h"
#include "class_fernd.h"
#include "class_fe.h"
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
    const double p0 {par->fernd_iso[0]}; //pccm
    const double k0 {par->fernd_iso[1]};
    const double a0 {par->fernd_iso[2]};
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
    const double r0 {par->fernd_scal[0]};
    const double z0 {par->fernd_scal[1]};
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

// isotropic turbulent field
double FErnd_iso::get_fernd(const vec3 &pos, Grid_fernd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid);
}

void FErnd_iso::write_grid_iso(Pond *par, Grid_fernd *grid){
    //PHASE I
    // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
    // initialize random seed
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    // from tools
    gsl_rng_set(r, toolkit::random_seed());
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in 1/kpc dimension
                double kx {i/(lx/CGS_U_kpc)};
                double ky {j/(ly/CGS_U_kpc)};
                double kz {l/(lz/CGS_U_kpc)};
                if(i>=grid->nx/2) kx -= grid->nx/(lx/CGS_U_kpc);
                if(j>=grid->ny/2) ky -= grid->ny/(ly/CGS_U_kpc);
                if(l>=grid->nz/2) kz -= grid->nz/(lz/CGS_U_kpc);
                const double k {sqrt(kx*kx + ky*ky + kz*kz)};
                // physical dk^3
                const double dk3 {pow(CGS_U_kpc,3.)/(lx*ly*lz)};
                // simpson's rule
                double element {2.*fe_spec(k,par)/3.};
                const double halfdk {0.5*sqrt(pow(CGS_U_kpc/lx,2)+pow(CGS_U_kpc/ly,2)+pow(CGS_U_kpc/lz,2))};
                element += fe_spec(k+halfdk,par)/6.;
                element += fe_spec(k-halfdk,par)/6.;
                // amplitude, dividing by two because equal allocation to Re and Im parts
                const double sigma {sqrt(element*dk3/2.0)};
                grid->fftw_fe_k[idx][0] = gsl_ran_gaussian(r,sigma);
                grid->fftw_fe_k[idx][1] = gsl_ran_gaussian(r,sigma);
                if(i==0 and j==0 and l==0){
                    grid->fftw_fe_k[idx][0] = 0.;
                    grid->fftw_fe_k[idx][1] = 0.;
                }
            }// l
        }// j
    }// i
    // free random memory
    gsl_rng_free(r);
    // execute DFT backward plan
    fftw_execute(grid->fftw_p);
    // PHASE II
    // RESCALING FIELD PROFILE IN REAL SPACE
    double fe_var {toolkit::Variance(grid->fftw_fe_k[0],grid->full_size)};
    for (decltype(grid->nx) i=0;i!=grid->nx;++i) {
        for (decltype(grid->ny) j=0;j!=grid->ny;++j) {
            for (decltype(grid->nz) l=0;l!=grid->nz;++l) {
                // get physical position
                vec3 pos {i*lx/(grid->nx-1) + grid->x_min,
                    j*ly/(grid->ny-1) + grid->y_min,
                    l*lz/(grid->nz-1) + grid->z_min};
                // get rescaling factor
                double ratio {sqrt(rescal_fact(pos,par))*par->fernd_iso[0]/sqrt(fe_var)};
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                grid->fftw_fe_k[idx][0] *= ratio;
            }
        }
    }
    // get real elements, use auxiliary function
    complex2real(grid->fftw_fe_k, grid->fftw_fe, grid->full_size);
    
#ifndef NDEBUG
    fe_var = toolkit::Variance(grid->fftw_fe, grid->full_size);
    cout<< "FERND: Numerical RMS: "<<sqrt(fe_var)<<" pccm"<<endl;
    //double fe_mean = toolkit::Mean(grid->fftw_fe,grid->full_size);
    //cout<<"FERND: Numrerical MEAN: "<<fe_mean<<"pccm"<<endl;
#endif
}


// get real components from fftw_complex arrays
void FErnd_iso::complex2real(const fftw_complex *input,double *output,const unsigned long int &size) {
    for(unsigned long int i=0;i!=size;++i){
        output[i] = input[i][0];
    }
}

// END
