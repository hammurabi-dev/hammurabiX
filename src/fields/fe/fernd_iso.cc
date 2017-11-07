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

// isotropic turbulent field
double FErnd_iso::get_fernd(const vec3_t<double> &pos, Grid_fernd *grid){
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
                vec3_t<double> pos {i*lx/(grid->nx-1) + grid->x_min,
                    j*ly/(grid->ny-1) + grid->y_min,
                    l*lz/(grid->nz-1) + grid->z_min};
                // get rescaling factor
                double ratio {sqrt(rescal_fact(pos,par))*par->fernd_iso.rms/sqrt(fe_var)};
                unsigned long int idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                grid->fftw_fe_k[idx][0] *= ratio;
            }
        }
    }
    // get real elements, use auxiliary function
    complex2real(grid->fftw_fe_k, grid->fftw_fe.get(), grid->full_size);
    
#ifndef NDEBUG
    fe_var = toolkit::Variance(grid->fftw_fe.get(), grid->full_size);
    cout<< "FERND: Numerical RMS: "<<sqrt(fe_var)<<" pccm"<<endl;
    //double fe_mean = toolkit::Mean(grid->fftw_fe.get(),grid->full_size);
    //cout<<"FERND: Numrerical MEAN: "<<fe_mean<<"pccm"<<endl;
#endif
}


// get real components from fftw_complex arrays
void FErnd_iso::complex2real(const fftw_complex *input,double *output,const unsigned long int &size) {
    for(unsigned long int i=0;i!=size;++i){
        output[i] = input[i][0];
    }
}
