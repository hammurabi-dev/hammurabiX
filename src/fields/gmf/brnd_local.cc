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
#include <param.h>
#include <grid.h>
#include <brnd.h>
#include <breg.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

using namespace std;

vec3_t<double> Brnd_local::get_brnd(const vec3_t<double> &pos, Grid_brnd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid);
}

void Brnd_local::write_grid(Param *par, Breg *breg, Grid_breg *gbreg, Grid_brnd *grid){
    // initialize random seed
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed(par->brnd_seed));
    unique_ptr<double[]> gaussian_num = unique_ptr<double[]>(new double[3*grid->full_size]());
    for(decltype(grid->full_size)i=0;i<3*grid->full_size;++i)
        gaussian_num[i] = gsl_ran_gaussian(r,1);
    unique_ptr<double[]> uniform_num = unique_ptr<double[]>(new double[2*grid->full_size]());
    for(decltype(grid->full_size)i=0;i<2*grid->full_size;++i)
        uniform_num[i] = gsl_rng_uniform(r);
    // free random memory
    gsl_rng_free(r);
    // start Fourier space filling
    double lx {grid->x_max-grid->x_min};
    double ly {grid->y_max-grid->y_min};
    double lz {grid->z_max-grid->z_min};
    const vec3_t<double> b {breg->get_breg(par->SunPosition,par,gbreg)};
    // physical dk^3
    const double dk3 {CGS_U_kpc*CGS_U_kpc*CGS_U_kpc/(lx*ly*lz)};
    const double halfdk {0.5*sqrt( CGS_U_kpc*CGS_U_kpc/(lx*lx) + CGS_U_kpc*CGS_U_kpc/(ly*ly) + CGS_U_kpc*CGS_U_kpc/(lz*lz) )};
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (decltype(grid->nx) i=0;i<grid->nx;++i) {
        // physical k in 1/kpc dimension
        vec3_t<double> k {CGS_U_kpc*i/lx,0,0};
        vec3_t<double> ep {0,0,0};
        vec3_t<double> em {0,0,0};
        if(i>=grid->nx/2) k.x -= CGS_U_kpc*grid->nx/lx;
        for (decltype(grid->ny) j=0;j<grid->ny;++j) {
            k.y = CGS_U_kpc*j/ly;
            if(j>=grid->ny/2) k.y -= CGS_U_kpc*grid->ny/ly;
            for (decltype(grid->nz) l=0;l<grid->nz;++l) {
                std::size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
                // FFT expects up to n/2 positive while n/2 to n negative
                k.z = CGS_U_kpc*l/lz;
                if(l>=grid->nz/2) k.z -= CGS_U_kpc*grid->nz/lz;
                const double ks {k.Length()};
                if(ks==0) continue;
                
                ep = eplus(b,k);
                em = eminus(b,k);
                // without Hermiticity fixing, all power multiplied by 2
                if(ep.Length()>1e-6){
                    double ang {cosa(b,k)};
                    double Pa {speca(ks,par)*0.66666667 + (speca(ks+halfdk,par) + speca(ks-halfdk,par))*0.16666667};
                    Pa *= fa(par->brnd_local.ma,ang)*dk3;
                    double Pf {specf(ks,par)*0.66666667 + (specf(ks+halfdk,par) + specf(ks-halfdk,par))*0.16666667};
                    Pf *= hf(par->brnd_local.beta,ang)*dk3;
                    double Ps {specs(ks,par)*0.66666667 + (specs(ks+halfdk,par) + specs(ks-halfdk,par))*0.16666667};
                    Ps *= fs(par->brnd_local.ma,ang)*hs(par->brnd_local.beta,ang)*dk3;
                    // complex angle
                    const double angp = 2*CGS_U_pi*uniform_num[2*idx];
                    const double angm = 2*CGS_U_pi*uniform_num[2*idx+1];
                    vec3_t<double> bkp {ep*gaussian_num[3*idx]*sqrt(2*Pa)};
                    vec3_t<double> bkm {em*(gaussian_num[3*idx+1]*sqrt(2.*Pf) + gaussian_num[3*idx+2]*sqrt(2.*Ps))};
                    grid->fftw_b_kx[idx][0] = bkp.x*cos(angp) + bkm.x*cos(angm);
                    grid->fftw_b_ky[idx][0] = bkp.y*cos(angp) + bkm.y*cos(angm);
                    grid->fftw_b_kz[idx][0] = bkp.z*cos(angp) + bkm.z*cos(angm);
                    grid->fftw_b_kx[idx][1] = bkp.x*sin(angp) + bkm.x*sin(angm);
                    grid->fftw_b_ky[idx][1] = bkp.y*sin(angp) + bkm.y*sin(angm);
                    grid->fftw_b_kz[idx][1] = bkp.z*sin(angp) + bkm.z*sin(angm);
                }
                else{
                    double Pf {specf(ks,par)*0.66666667 + (specf(ks+halfdk,par) + specf(ks-halfdk,par))*0.16666667};
                    Pf *= hf(par->brnd_local.beta,1)*dk3;
                    if(i==0 and j==0){
                        ep.x = k.x;
                        em.y = k.y;
                    }
                    else{
                        ep.x = k.x*k.z/(ks*sqrt(k.x*k.x+k.y*k.y));
                        ep.y = k.y*k.z/(ks*sqrt(k.x*k.x+k.y*k.y));
                        ep.z = -sqrt(k.x*k.x+k.y*k.y)/ks;
                        em.x = -k.y/sqrt(k.x*k.x+k.y*k.y);
                        em.y = k.x/sqrt(k.x*k.x+k.y*k.y);
                        em.z = 0.;
                    }
                    // complex angle
                    const double angp = 2*CGS_U_pi*uniform_num[2*idx];
                    const double angm = 2*CGS_U_pi*uniform_num[2*idx+1];
                    vec3_t<double> bkp {ep*gaussian_num[3*idx]*sqrt(Pf)};
                    vec3_t<double> bkm {ep*gaussian_num[3*idx+1]*sqrt(Pf)};
                    grid->fftw_b_kx[idx][0] = bkp.x*cos(angp) + bkm.x*cos(angm);
                    grid->fftw_b_ky[idx][0] = bkp.y*cos(angp) + bkm.y*cos(angm);
                    grid->fftw_b_kz[idx][0] = bkp.z*cos(angp) + bkm.z*cos(angm);
                    grid->fftw_b_kx[idx][1] = bkp.x*sin(angp) + bkm.x*sin(angm);
                    grid->fftw_b_ky[idx][1] = bkp.y*sin(angp) + bkm.y*sin(angm);
                    grid->fftw_b_kz[idx][1] = bkp.z*sin(angp) + bkm.z*sin(angm);
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
    // execute DFT backward plan
    fftw_execute(grid->fftw_px_bw);
    fftw_execute(grid->fftw_py_bw);
    fftw_execute(grid->fftw_pz_bw);
    // get real elements, use auxiliary function
    complex2real(grid->fftw_b_kx, grid->fftw_b_x.get(), grid->full_size);
    complex2real(grid->fftw_b_ky, grid->fftw_b_y.get(), grid->full_size);
    complex2real(grid->fftw_b_kz, grid->fftw_b_z.get(), grid->full_size);
#ifdef DEBUG
    const double b_var {toolkit::Variance(grid->fftw_b_x.get(),grid->full_size)+toolkit::Variance(grid->fftw_b_y.get(),grid->full_size)+toolkit::Variance(grid->fftw_b_z.get(),grid->full_size)};
    cout<< "BRND: Numerical RMS: "<<sqrt(b_var)/CGS_U_muGauss<<" microG"<<endl;
#endif
}

// PRIVATE FUNCTIONS FOR LOW-BETA SUB-ALFVENIC PLASMA
double Brnd_local::cosa(const vec3_t<double> &b,const vec3_t<double> &k){
    return dotprod(toolkit::versor(b),toolkit::versor(k));
}

vec3_t<double> Brnd_local::eplus(const vec3_t<double> &b,const vec3_t<double> &k){
    vec3_t<double> tmp {crossprod(toolkit::versor(k),toolkit::versor(b))};
    if(tmp.Length()<1e-6){
        return tmp;
    }
    else{
        return tmp/tmp.Length();
    }
}

vec3_t<double> Brnd_local::eminus(const vec3_t<double> &b,const vec3_t<double> &k){
    vec3_t<double> tmp {crossprod(crossprod(toolkit::versor(k),toolkit::versor(b)),toolkit::versor(k))};
    if(tmp.Length()<1e-6){
        return tmp;
    }
    else{
        return tmp/tmp.Length();
    }
}

double Brnd_local::dynamo(const double &beta,const double &cosa){
    return (1+0.5*beta)*(1+0.5*beta) - 2.*beta*cosa*cosa;
}

double Brnd_local::hs(const double &beta,const double &cosa){
    if(cosa<1e-6){
        return 0;
    }
    else{
        double dispersion {0.5*(1+0.5*beta)*( 1-sqrt(1-2*beta*cosa*cosa/((1+0.5*beta)*(1+0.5*beta))) )};
        double lambda {(1-sqrt(dynamo(beta,cosa))-0.5*beta)/(1+sqrt(dynamo(beta,cosa))+0.5*beta)};
        return lambda*lambda/( dispersion*(1./(cosa*cosa)+(lambda*lambda-1)) );
    }
}

double Brnd_local::hf(const double &beta,const double &cosa){
    if(cosa<1e-6){
        return 0;
    }
    else{
        double dispersion {0.5*(1+0.5*beta)*( 1+sqrt(1-2*beta*cosa*cosa/((1+0.5*beta)*(1+0.5*beta))) )};
        double lambda {(1-sqrt(dynamo(beta,cosa))+0.5*beta)/(1+sqrt(dynamo(beta,cosa))-0.5*beta)};
        return 1./(dispersion*(1+lambda*lambda*(1./(cosa*cosa) -1)));
    }
}

double Brnd_local::fa(const double &ma,const double &cosa){
    return exp( -pow(ma,-1.33333333)*cosa*cosa/pow(1-cosa*cosa,0.66666667) );
}

double Brnd_local::fs(const double &ma,const double &cosa){
    return exp( -pow(ma,-1.33333333)*cosa*cosa/pow(1-cosa*cosa,0.66666667) );
}

double Brnd_local::speca(const double &k,Param *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_local.pa0};
    const double k0 {par->brnd_local.k0};
    const double a0 {par->brnd_local.aa0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // avoid nan
    if(k<=0.){
        return 0.;
    }
    // power law
    double P{0.};
    if(k>k0){
        P = p0/pow(k/k0,a0);
    }
    return P*unit;
}

double Brnd_local::specf(const double &k,Param *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_local.pf0};
    const double k0 {par->brnd_local.k0};
    const double a0 {par->brnd_local.af0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // avoid nan
    if(k<=0.){
        return 0.;
    }
    // power law
    double P{0.};
    if(k>k0){
        P = p0/pow(k/k0,a0);
    }
    return P*unit;
}

double Brnd_local::specs(const double &k,Param *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_local.ps0};
    const double k0 {par->brnd_local.k0};
    const double a0 {par->brnd_local.as0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // avoid nan
    if(k<=0.){
        return 0.;
    }
    // power law
    double P{0.};
    if(k>k0){
        P = p0/pow(k/k0,a0);
    }
    return P*unit;
}

// get real components from fftw_complex arrays
void Brnd_local::complex2real(const fftw_complex *input,double *output,const std::size_t &size) {
    for(std::size_t i=0;i!=size;++i){
        output[i] = input[i][0];
    }
}
// END
