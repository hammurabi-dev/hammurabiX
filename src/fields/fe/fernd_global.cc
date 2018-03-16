#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_integration.h>
#include <omp.h>
#include <param.h>
#include <grid.h>
#include <fernd.h>
#include <fereg.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

using namespace std;

// isotropic turbulent field
double FErnd_global::get_fernd(const vec3_t<double> &pos, Grid_fernd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid);
}

// since we are using rms normalization
// p0 is hidden and not affecting anything
double FErnd_global::spec(const double &k, Param *par){
    //units fixing
    const double p0 {par->fernd_global.rms}; //pccm
    const double kr {k/par->fernd_global.k0};
    const double a0 {par->fernd_global.a0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    double P {0.};
    if(kr>1){
        P = p0*pow(kr,a0);
    }
    return P*unit;
}

// galactic scaling of random field energy density
// set to 1 at observer's place
double FErnd_global::rescal(const vec3_t<double> &pos, Param *par){
    const double r_cyl {sqrt(pos.x*pos.x+pos.y*pos.y) - fabs(par->SunPosition.x)};
    const double z {fabs(pos.z) - fabs(par->SunPosition.z)};
    const double r0 {par->fernd_global.r0};
    const double z0 {par->fernd_global.z0};
    return exp(-r_cyl/r0)*exp(-z/z0);
}

void FErnd_global::write_grid_global(Param *par, Grid_fernd *grid){
    //PHASE I
    // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
    // initialize random seed
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed(par->fernd_seed));
    unique_ptr<double[]> gaussian_num = unique_ptr<double[]>(new double[2*grid->full_size]());
    for(decltype(grid->full_size)i=0;i<2*grid->full_size;++i)
        gaussian_num[i] = gsl_ran_gaussian(r,1);
    // free random memory
    gsl_rng_free(r);
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    // physical k in 1/kpc dimension
    // physical dk^3
    const double dk3 {CGS_U_kpc*CGS_U_kpc*CGS_U_kpc/(lx*ly*lz)};
    const double halfdk {0.5*sqrt( CGS_U_kpc*CGS_U_kpc/(lx*lx) + CGS_U_kpc*CGS_U_kpc/(ly*ly) + CGS_U_kpc*CGS_U_kpc/(lz*lz) )};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (decltype(grid->nx) i=0;i<grid->nx;++i) {
        double kx {CGS_U_kpc*i/lx};
        if(i>=grid->nx/2) kx -= CGS_U_kpc*grid->nx/lx;
        for (decltype(grid->ny) j=0;j<grid->ny;++j) {
            double ky {CGS_U_kpc*j/ly};
            if(j>=grid->ny/2) ky -= CGS_U_kpc*grid->ny/ly;
            for (decltype(grid->nz) l=0;l<grid->nz;++l) {
                std::size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                double kz {CGS_U_kpc*l/lz};
                if(l>=grid->nz/2) kz -= CGS_U_kpc*grid->nz/lz;
                const double k {sqrt(kx*kx + ky*ky + kz*kz)};
                
                // simpson's rule
                double element {2.*spec(k,par)/3.};
                element += spec(k+halfdk,par)/6.;
                element += spec(k-halfdk,par)/6.;
                // amplitude, dividing by two because equal allocation to Re and Im parts
                const double sigma {sqrt(0.5*element*dk3)};
                grid->fftw_fe_k[idx][0] = sigma*gaussian_num[2*idx];
                grid->fftw_fe_k[idx][1] = sigma*gaussian_num[2*idx+1];
            }// l
        }// j
    }// i
    grid->fftw_fe_k[0][0] = 0.;
    grid->fftw_fe_k[0][1] = 0.;
    // execute DFT backward plan
    fftw_execute(grid->fftw_p);
    // PHASE II
    // RESCALING FIELD PROFILE IN REAL SPACE
    double fe_var {toolkit::Variance(grid->fftw_fe_k[0],grid->full_size)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (decltype(grid->nx) i=0;i<grid->nx;++i) {
        vec3_t<double> pos {i*lx/(grid->nx-1) + grid->x_min,0,0};
        for (decltype(grid->ny) j=0;j<grid->ny;++j) {
            pos.y = j*ly/(grid->ny-1) + grid->y_min;
            for (decltype(grid->nz) l=0;l<grid->nz;++l) {
                // get physical position
                pos.z = l*lz/(grid->nz-1) + grid->z_min;
                // get rescaling factor
                double ratio {sqrt(rescal(pos,par))*par->fernd_global.rms/sqrt(fe_var)};
                std::size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                grid->fftw_fe_k[idx][0] *= ratio;
            }
        }
    }
    // get real elements, use auxiliary function
    toolkit::complex2real(grid->fftw_fe_k, grid->fftw_fe.get(), grid->full_size);
}
