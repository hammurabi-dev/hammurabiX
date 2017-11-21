#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <memory>
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

using namespace std;

// isotropic turbulent field
vec3_t<double> Brnd_iso::get_brnd(const vec3_t<double> &pos, Grid_brnd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid);
}

// since we are using rms normalization
// p0 is hidden and not affecting anything
double Brnd_iso::b_spec(const double &k, Pond *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_iso.rms*pow(CGS_U_muGauss,2)};
    const double k0 {par->brnd_iso.k0};
    const double a0 {par->brnd_iso.a0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // avoid nan
    if(k<=0.){
        return 0.;
    }
    // Giacalone-Jokipii
    //const double P = 2*p0/(1.+pow(k/k0,a0));
    // power law
    double P{0.};
    if(k>k0){
        P = p0/pow(k/k0,a0);
    }
    return P*unit;
}

// galactic scaling of random field energy density
// set to 1 at observer's place
double Brnd_iso::rescal_fact(const vec3_t<double> &pos, Pond *par){
    const double r_cyl {(sqrt(pos.x*pos.x+pos.y*pos.y) - fabs(par->SunPosition.x))/CGS_U_kpc};
    const double z {(fabs(pos.z) - fabs(par->SunPosition.z))/CGS_U_kpc};
    const double r0 {par->brnd_scal.r0};
    const double z0 {par->brnd_scal.z0};
    if(r_cyl==0. or z==0){return 1.;}
    else{
        return exp(-r_cyl/r0)*exp(-z/z0);
    }
}

// we tried to optimize the calc time cost for this function
void Brnd_iso::write_grid_iso(Pond *par, Grid_brnd *grid){
    // PHASE I
    // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
    // initialize random seed
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed(par->brnd_seed));
    // start Fourier space filling
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    // physical k in 1/kpc dimension
    // physical dk^3
    const double dk3 {CGS_U_kpc*CGS_U_kpc*CGS_U_kpc/(lx*ly*lz)};
    const double halfdk {0.5*sqrt( CGS_U_kpc*CGS_U_kpc/(lx*lx) + CGS_U_kpc*CGS_U_kpc/(ly*ly) + CGS_U_kpc*CGS_U_kpc/(lz*lz) )};
#pragma omp parallel for ordered schedule(static,1)
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
                double element {2.*b_spec(k,par)/3.};
                element += b_spec(k+halfdk,par)/6.;
                element += b_spec(k-halfdk,par)/6.;
                // amplitude, dividing by two because equal allocation to Re and Im parts
                const double sigma {sqrt(0.5*element*dk3)};
#pragma omp ordered
                {
                    grid->fftw_b_kx[idx][0] = gsl_ran_gaussian(r,sigma);
                    grid->fftw_b_ky[idx][0] = gsl_ran_gaussian(r,sigma);
                    grid->fftw_b_kz[idx][0] = gsl_ran_gaussian(r,sigma);
                    grid->fftw_b_kx[idx][1] = gsl_ran_gaussian(r,sigma);
                    grid->fftw_b_ky[idx][1] = gsl_ran_gaussian(r,sigma);
                    grid->fftw_b_kz[idx][1] = gsl_ran_gaussian(r,sigma);
                }
            }// l
        }// j
    }// i
    // fix the very 0th
    grid->fftw_b_kx[0][0] = 0.;
    grid->fftw_b_ky[0][0] = 0.;
    grid->fftw_b_kz[0][0] = 0.;
    grid->fftw_b_kx[0][1] = 0.;
    grid->fftw_b_ky[0][1] = 0.;
    grid->fftw_b_kz[0][1] = 0.;
    // free random memory
    gsl_rng_free(r);
    // no Hermiticity fixing, complex 2 complex
    // execute DFT backward plan
    fftw_execute(grid->fftw_px_bw);
    fftw_execute(grid->fftw_py_bw);
    fftw_execute(grid->fftw_pz_bw);
    // PHASE II
    // RESCALING FIELD PROFILE IN REAL SPACE
    double b_var {toolkit::Variance(grid->fftw_b_kx[0], grid->full_size)};
#pragma omp parallel for
    for (decltype(grid->nx) i=0;i<grid->nx;++i) {
        vec3_t<double> pos {i*lx/(grid->nx-1) + grid->x_min,0,0};
        for (decltype(grid->ny) j=0;j<grid->ny;++j) {
            pos.y = j*ly/(grid->ny-1) + grid->y_min;
            for (decltype(grid->nz) l=0;l<grid->nz;++l) {
                // get physical position
                pos.z = l*lz/(grid->nz-1) + grid->z_min;
                // get rescaling factor
                double ratio {sqrt(rescal_fact(pos,par))*par->brnd_iso.rms*CGS_U_muGauss/sqrt(3.*b_var)};
                // add anisotropic field to random one
                std::size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                grid->fftw_b_kx[idx][0] *= ratio;
                grid->fftw_b_ky[idx][0] *= ratio;
                grid->fftw_b_kz[idx][0] *= ratio;
                grid->fftw_b_kx[idx][1] *= ratio;
                grid->fftw_b_ky[idx][1] *= ratio;
                grid->fftw_b_kz[idx][1] *= ratio;
            } //l
        } //j
    } //i
    // execute DFT forward plan
    fftw_execute(grid->fftw_px_fw);
    fftw_execute(grid->fftw_py_fw);
    fftw_execute(grid->fftw_pz_fw);
    // PHASE III
    // RE-ORTHOGONALIZING IN FOURIER SPACE
    // Gram-Schmidt process
#pragma omp parallel for
    for (decltype(grid->nx) i=0;i<grid->nx;++i) {
        vec3_t<double> tmp_k {CGS_U_kpc*i/lx,0,0};
        if(i>=grid->nx/2) tmp_k.x -= CGS_U_kpc*grid->nx/lx;
        for (decltype(grid->ny) j=0;j<grid->ny;++j) {
            tmp_k.y = CGS_U_kpc*j/ly;
            if(j>=grid->ny/2) tmp_k.y -= CGS_U_kpc*grid->ny/ly;
            for (decltype(grid->nz) l=0;l<grid->nz;++l) {
                std::size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                // physical k in 1/kpc
                tmp_k.z = CGS_U_kpc*l/lz;
                if(l>=grid->nz/2) tmp_k.z -= CGS_U_kpc*grid->nz/lz;
                
                const vec3_t<double> tmp_b_re {grid->fftw_b_kx[idx][0],grid->fftw_b_ky[idx][0],grid->fftw_b_kz[idx][0]};
                const vec3_t<double> tmp_b_im {grid->fftw_b_kx[idx][1],grid->fftw_b_ky[idx][1],grid->fftw_b_kz[idx][1]};
                
                const vec3_t<double> free_b_re {gramschmidt(tmp_k,tmp_b_re)};
                const vec3_t<double> free_b_im {gramschmidt(tmp_k,tmp_b_im)};
                
                grid->fftw_b_kx[idx][0] = free_b_re.x;
                grid->fftw_b_ky[idx][0] = free_b_re.y;
                grid->fftw_b_kz[idx][0] = free_b_re.z;
                grid->fftw_b_kx[idx][1] = free_b_im.x;
                grid->fftw_b_ky[idx][1] = free_b_im.y;
                grid->fftw_b_kz[idx][1] = free_b_im.z;
            }// l
        }// j
    }// i
    // fixing the very 0th
    grid->fftw_b_kx[0][0] = 0.;
    grid->fftw_b_ky[0][0] = 0.;
    grid->fftw_b_kz[0][0] = 0.;
    grid->fftw_b_kx[0][1] = 0.;
    grid->fftw_b_ky[0][1] = 0.;
    grid->fftw_b_kz[0][1] = 0.;
    
    // execute DFT backward plan
    fftw_execute(grid->fftw_px_bw);
    fftw_execute(grid->fftw_py_bw);
    fftw_execute(grid->fftw_pz_bw);
    // get real elements, use auxiliary function
    complex2real(grid->fftw_b_kx, grid->fftw_b_x.get(), grid->full_size);
    complex2real(grid->fftw_b_ky, grid->fftw_b_y.get(), grid->full_size);
    complex2real(grid->fftw_b_kz, grid->fftw_b_z.get(), grid->full_size);
    // according to FFTW3 manual
    // transform forward followed by backword scale up array by nx*ny*nz
    double inv_grid_size = 1.0/grid->full_size;
#pragma omp parallel for
    for(std::size_t i=0;i<grid->full_size;++i){
        grid->fftw_b_x[i] *= inv_grid_size;
        grid->fftw_b_y[i] *= inv_grid_size;
        grid->fftw_b_z[i] *= inv_grid_size;
    }
    cout<<grid->fftw_b_x[30]/CGS_U_muGauss<<endl;
}

// get real components from fftw_complex arrays
void Brnd_iso::complex2real(const fftw_complex *input,double *output,const std::size_t &size) {
    for(std::size_t i=0;i!=size;++i){
        output[i] = input[i][0];
    }
}

// Gram-Schimdt, rewritten using Healpix vec3 library
// tiny error caused by machine is inevitable
vec3_t<double> Brnd_iso::gramschmidt(const vec3_t<double> &k,const vec3_t<double> &b){
    if(k.Length()==0 or b.Length()==0){
        return vec3_t<double> {0,0,0};
    }
    const double inv_k_mod = 1./k.SquaredLength();
    vec3_t<double> b_free {b.x - k.x*dotprod(k,b)*inv_k_mod,
        b.y - k.y*dotprod(k,b)*inv_k_mod,
        b.z - k.z*dotprod(k,b)*inv_k_mod};
    
    b_free = toolkit::versor(b_free)*b.Length();
    return b_free;
}
