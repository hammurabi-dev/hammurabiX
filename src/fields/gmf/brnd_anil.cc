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
#include "pond.h"
#include "grid.h"
#include "brnd.h"
#include "breg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"
// timer
//#include <sys/time.h>
//#include <sys/resource.h>
//#define RCPU_TIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + 1e-6 * ruse.ru_utime.tv_usec)

using namespace std;

//local anisotropic turbulent field
void Brnd_anil::anisotropy(double *tensor, Pond *par, Breg *breg, Grid_breg *gbreg){
    exit(1);
}

void Brnd_anil::write_grid_ani(Pond *par, Breg *breg, Grid_breg *gbreg, Grid_brnd *grid){
    exit(1);
    /*
     // initialize random seed
     gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
     gsl_rng_set(r, toolkit::random_seed());
     
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
     
     // cheating simpson evaluation
     const double halfdk {0.5*sqrt(pow(CGS_U_kpc/lx,2)+pow(CGS_U_kpc/ly,2)+pow(CGS_U_kpc/lz,2))};
     double element {2.*b_spec(k,par)/3.};
     element += b_spec(k+halfdk,par)/6.;
     element += b_spec(k-halfdk,par)/6.;
     
     // UNFINISHED
     }
     }// l
     }// j
     }// i
     // free random memory
     gsl_rng_free(r);
     
     // execute DFT plan
     fftw_execute(grid->fftw_px);
     fftw_execute(grid->fftw_py);
     fftw_execute(grid->fftw_pz);
     
     // get real elements, use auxiliary function
     complex2real(grid->fftw_b_kx, grid->fftw_b_x.get(), grid->full_size);
     complex2real(grid->fftw_b_ky, grid->fftw_b_y.get(), grid->full_size);
     complex2real(grid->fftw_b_kz, grid->fftw_b_z.get(), grid->full_size);
     
     #ifndef NDEBUG
     double b_var;
     b_var = toolkit::Variance(grid->fftw_b_x.get(), grid->full_size);
     b_var +=toolkit::Variance(grid->fftw_b_y.get(), grid->full_size);
     b_var +=toolkit::Variance(grid->fftw_b_z.get(), grid->full_size);
     cout<< "BRND: Numerical RMS: "<<sqrt(b_var)/CGS_U_muGauss<<" microG"<<endl;
     // check average divergence
     double div=0;
     for(decltype(grid->nx) i=2;i!=grid->nx-2;++i) {
     for(decltype(grid->ny) j=2;j!=grid->ny-2;++j) {
     for(decltype(grid->nz) k=2;k!=grid->nz-2;++k) {
     auto idxf = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i+1,j,k);
     auto idxb = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i-1,j,k);
     auto idyf = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j+1,k);
     auto idyb = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j-1,k);
     auto idzf = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k+1);
     auto idzb = toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,k-1);
     div += fabs( (grid->fftw_b_x[idxf]-grid->fftw_b_x[idxb]) + (grid->fftw_b_y[idyf]-grid->fftw_b_y[idyb]) + (grid->fftw_b_z[idzf]-grid->fftw_b_z[idzb]) );
     }
     }
     }
     cout<<"BRND: Averaged divergence: "<<div/(CGS_U_muGauss*grid->nx*grid->ny*grid->nz)<<" microG/grid"<<endl;
     #endif
     */
}
