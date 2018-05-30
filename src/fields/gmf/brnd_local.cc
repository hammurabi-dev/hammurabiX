#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>
#include <omp.h>
#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <param.h>
#include <grid.h>
#include <brnd.h>
#include <breg.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>


vec3_t<double> Brnd_local::get_brnd(const vec3_t<double> &pos,
                                    Grid_brnd *grid){
    // interpolate written grid to given position
    // check if you have called ::write_grid
    return read_grid(pos,grid);
}

void Brnd_local::write_grid(Param *par,
                            Breg *breg,
                            Grid_breg *gbreg,
                            Grid_brnd *grid){
    // initialize random seed
#ifdef _OPENMP
    gsl_rng **threadvec = new gsl_rng *[omp_get_max_threads()];
    for (int b=0;b<omp_get_max_threads();++b){
        threadvec[b] = gsl_rng_alloc(gsl_rng_taus);
        gsl_rng_set(threadvec[b],b+toolkit::random_seed(par->brnd_seed));
    }
#else
    gsl_rng *r {gsl_rng_alloc(gsl_rng_taus)};
    gsl_rng_set(r, toolkit::random_seed(par->brnd_seed));
#endif
    // start Fourier space filling
    const double lx {grid->x_max-grid->x_min};
    const double ly {grid->y_max-grid->y_min};
    const double lz {grid->z_max-grid->z_min};
    const vec3_t<double> B {breg->get_breg(par->SunPosition,par,gbreg)};
    // physical dk^3
    const double dk3 {CGS_U_kpc*CGS_U_kpc*CGS_U_kpc/(lx*ly*lz)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static) // DO NOT CHANGE SCHEDULE TYPE
#endif
    for (decltype(grid->nx) i=0;i<grid->nx;++i) {
#ifdef _OPENMP
        auto seed_id = threadvec[omp_get_thread_num()];
#else
        auto seed_id = r;
#endif
        vec3_t<double> k {CGS_U_kpc*i/lx,0,0};
        if(i>=(grid->nx+1)/2) k.x -= CGS_U_kpc*grid->nx/lx;
        /**
         * it's better to calculate indeces manually
         * just for reference, how indeces are calculated
         * const size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
         */
        const std::size_t idx_lv1 {i*grid->ny*grid->nz};
        for (decltype(grid->ny) j=0;j<grid->ny;++j) {
            k.y = CGS_U_kpc*j/ly;
            if(j>=(grid->ny+1)/2) k.y -= CGS_U_kpc*grid->ny/ly;
            const std::size_t idx_lv2 {idx_lv1+j*grid->nz};
            for (decltype(grid->nz) l=0;l<grid->nz;++l) {
                /**
                 * the very 0th term is fixed to zero in allocation
                 */
                if(i==0 and j==0 and l==0) continue;
                k.z = CGS_U_kpc*l/lz;
                if(l>=(grid->nz+1)/2) k.z -= CGS_U_kpc*grid->nz/lz;
                const double ks {k.Length()};
                const std::size_t idx {idx_lv2+l};
                vec3_t<double> ep {eplus(B,k)};
                vec3_t<double> em {eminus(B,k)};
                /**
                 * since there is no specific rule about how to allocate spectrum power
                 * between Re and Im part in b+ and b-
                 * we multiply power by two and set Im parts to zero
                 * in the end of backward trans, we retrive Re part of bx,by,bz only
                 */
                if(ep.SquaredLength()>1e-6){
                    double ang {cosa(B,k)};
                    const double Pa {speca(ks,par)*fa(par->brnd_local.ma,ang)*dk3};
                    double Pf {specf(ks,par)*hf(par->brnd_local.beta,ang)*dk3};
                    double Ps {specs(ks,par)*fs(par->brnd_local.ma,ang)*hs(par->brnd_local.beta,ang)*dk3};
                    /**
                     * b+ is independent from b- in terms of power
                     * fast and slow modes are independent
                     */
                    const double Ap = gsl_ran_ugaussian(seed_id)*sqrt(2.*Pa);
                    const double Am = gsl_ran_ugaussian(seed_id)*sqrt(2.*Pf) + gsl_ran_ugaussian(seed_id)*sqrt(2.*Ps);
                    const vec3_t<double> bkp {ep*Ap};
                    const vec3_t<double> bkm {em*Am};
                    /**
                     * c0_R = bx_Re - by_Im
                     * c0_I = bx_Im + by_Re
                     * c1_R = by_Re - bz_Im
                     * c1_I = by_Im + bz_Re
                     * bx_R = bkp_Re.x+bkm_Re.x;
                     * etc.
                     * note that all Im parts are zero, c1_I = c0_R
                     */
                    grid->c0[idx][0] = bkp.x + bkm.x;
                    grid->c0[idx][1] = bkp.y + bkm.y;
                    grid->c1[idx][0] = bkp.y + bkm.y;
                    grid->c1[idx][1] = bkp.z + bkm.z;
                }
                else{
                    double Pf {specf(ks,par)*hf(par->brnd_local.beta,1)*dk3};
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
                    /**
                     * b+ and b- share power
                     */
                    const double Af = gsl_ran_ugaussian(seed_id)*sqrt(2.*Pf);
                    const double share = gsl_rng_uniform(seed_id);
                    const vec3_t<double> bkp {ep*Af*share};
                    const vec3_t<double> bkm {em*Af*(1.-share)};
                    grid->c0[idx][0] = bkp.x + bkm.x;
                    grid->c0[idx][1] = bkp.y + bkm.y;
                    grid->c1[idx][0] = grid->c0[idx][1];
                    grid->c1[idx][1] = bkp.z + bkm.z;
                }
            }// l
        }// j
    }// i
    // free random memory
#ifdef _OPENMP
    for (int b=0;b<omp_get_max_threads();++b)
        gsl_rng_free(threadvec[b]);
    delete [] threadvec;
#else
    gsl_rng_free(r);
#endif
    // execute DFT backward plan
    fftw_execute_dft(grid->plan_c0_bw,grid->c0,grid->c0);
    fftw_execute_dft(grid->plan_c1_bw,grid->c1,grid->c1);
    /**
     * now we need to get real parts manually
     * since we start in k-space real fields
     * the results in x-space are complex fields
     */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (decltype(grid->nx) i=0;i<grid->nx;++i) {
        decltype(grid->nx) i_sym {grid->nx-i};// apply Hermitian symmetry
        if(i==0) i_sym = i;
        /**
         * it's better to calculate indeces manually
         * just for reference, how indeces are calculated
         * const std::size_t idx {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i,j,l)};
         * const std::size_t idx_sym {toolkit::Index3d(grid->nx,grid->ny,grid->nz,i_sym,j_sym,l_sym)};
         */
        const std::size_t idx_lv1 {i*grid->ny*grid->nz};
        const std::size_t idx_sym_lv1 {i_sym*grid->ny*grid->nz};
        for (decltype(grid->ny) j=0;j<grid->ny;++j) {
            decltype(grid->ny) j_sym {grid->ny-j};// apply Hermitian symmetry
            if(j==0) j_sym = j;
            const std::size_t idx_lv2 {idx_lv1+j*grid->nz};
            const std::size_t idx_sym_lv2 {idx_sym_lv1+j_sym*grid->nz};
            for (decltype(grid->nz) l=0;l<grid->nz;++l) {
                decltype(grid->nz) l_sym {grid->nz-l};// apply Hermitian symmetry
                if(l==0) l_sym = l;
                const std::size_t idx {idx_lv2+l}; //k
                const std::size_t idx_sym {idx_sym_lv2+l_sym}; //-k
                /**
                 * reconstruct bx,by,bz from c0,c1,c*0,c*1
                 *
                 * c0(k) = bx(k) + i by(k)
                 * c*0(-k) = bx(k) - i by(k)
                 * c1(k) = by(k) + i bz(k)
                 * c1*1(-k) = by(k) - i bz(k)
                 */
                grid->bx[idx] = 0.5*(grid->c0[idx][0]+grid->c0[idx_sym][0]);
                grid->by[idx] = 0.5*(grid->c1[idx][0]+grid->c1[idx_sym][0]);
                grid->bz[idx] = 0.5*(grid->c1[idx_sym][1]+grid->c1[idx][1]);
            }
        }
    }
}

// PRIVATE FUNCTIONS FOR LOW-BETA SUB-ALFVENIC PLASMA
vec3_t<double> Brnd_local::eplus(const vec3_t<double> &b,
                                 const vec3_t<double> &k){
    vec3_t<double> tmp {crossprod(toolkit::versor(k),toolkit::versor(b))};
    if(tmp.SquaredLength()<1e-12){
        return tmp;
    }
    else{
        return tmp/tmp.Length();
    }
}

vec3_t<double> Brnd_local::eminus(const vec3_t<double> &b,
                                  const vec3_t<double> &k){
    vec3_t<double> tmp {crossprod(crossprod(toolkit::versor(k),toolkit::versor(b)),toolkit::versor(k))};
    if(tmp.SquaredLength()<1e-12){
        return tmp;
    }
    else{
        return tmp/tmp.Length();
    }
}

double Brnd_local::hs(const double &beta,
                      const double &cosa){
    if(cosa<1e-6){
        return 0;
    }
    else{
        double dispersion {0.5*(1+0.5*beta)*( 1-sqrt(1-2*beta*cosa*cosa/((1+0.5*beta)*(1+0.5*beta))) )};
        double lambda {(1-sqrt(dynamo(beta,cosa))-0.5*beta)/(1+sqrt(dynamo(beta,cosa))+0.5*beta)};
        return lambda*lambda/( dispersion*(1./(cosa*cosa)+(lambda*lambda-1)) );
    }
}

double Brnd_local::hf(const double &beta,
                      const double &cosa){
    if(cosa<1e-6){
        return 0;
    }
    else{
        double dispersion {0.5*(1+0.5*beta)*( 1+sqrt(1-2*beta*cosa*cosa/((1+0.5*beta)*(1+0.5*beta))) )};
        double lambda {(1-sqrt(dynamo(beta,cosa))+0.5*beta)/(1+sqrt(dynamo(beta,cosa))-0.5*beta)};
        return 1./(dispersion*(1+lambda*lambda*(1./(cosa*cosa) -1)));
    }
}

double Brnd_local::speca(const double &k,
                         Param *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_local.pa0};
    const double kr {k/par->brnd_local.k0};
    const double a0 {par->brnd_local.aa0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // power law
    double P {0.};
    if(kr>=1.){
        P = p0/pow(kr,a0);
    }
    return P*unit;
}

double Brnd_local::specf(const double &k,
                         Param *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_local.pf0};
    const double kr {k/par->brnd_local.k0};
    const double a0 {par->brnd_local.af0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // power law
    double P{0.};
    if(kr>=1.){
        P = p0/pow(kr,a0);
    }
    return P*unit;
}

double Brnd_local::specs(const double &k,
                         Param *par){
    //units fixing, wave vector in 1/kpc units
    const double p0 {par->brnd_local.ps0};
    const double kr {k/par->brnd_local.k0};
    const double a0 {par->brnd_local.as0};
    const double unit = 1./(4*CGS_U_pi*k*k);
    // power law
    double P {0.};
    if(kr>=1.){
        P = p0/pow(kr,a0);
    }
    return P*unit;
}
// END
