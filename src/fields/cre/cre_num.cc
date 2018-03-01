#include <iostream>
#include <memory>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <cre.h>
#include <param.h>
#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>
#include <ap_err.h>
using namespace std;

// numerical CRE flux
// use bilinear/trilinear interpolationi according to the dimension of CRE flux grid
double CRE_num::read_grid(const std::size_t &Eidx, const vec3_t<double> &pos,Grid_cre *grid){
    // if grid in spatial 2D
    if(grid->nr!=0){
        // sylindrical galactic centric position
        const double r {sqrt(pos.x*pos.x+pos.y*pos.y)};
        // bilinear interpolation
        // notice that lr is radius
        double tmp {(grid->nr-1)*(r/grid->r_max)};
        if(tmp<0 or tmp>grid->nr-1) {return 0.;}
        decltype(grid->nr) rl {(std::size_t)floor(tmp)};
        const double rd {tmp - rl};
        tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
        if(tmp<0 or tmp>grid->nz-1) {return 0.;}
        decltype(grid->nr) zl {(std::size_t)floor(tmp)};
        const double zd {tmp - zl};
#ifdef DEBUG
        if(rd<0 or zd<0 or rd>1 or zd>1){
            ap_err("wrong value");
            exit(1);
        }
#endif
        double cre;
        if(rl+1<grid->nr and zl+1<grid->nz){
            std::size_t idx1 {toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl)};
            std::size_t idx2 {toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl+1,zl)};
            double i1 {grid->cre_flux[idx1]*(1-rd) + grid->cre_flux[idx2]*rd};
            idx1 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl+1);
            idx2 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl+1,zl+1);
            double i2 {grid->cre_flux[idx1]*(1-rd) + grid->cre_flux[idx2]*rd};
            cre = (i1*(1-zd)+i2*zd);
        }
        else{
            std::size_t idx1 {toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl)};
            cre = grid->cre_flux[idx1];
        }
#ifdef DEBUG
        if(cre<0){
            ap_err("negative CRE flux");
            exit(1);
        }
#endif
        return cre;
    }
    // if grid in spatial 3D
    else if(grid->nx!=0){
        //trilinear interpolation
        double tmp {(grid->nx-1)*(pos.x-grid->x_min)/(grid->x_max-grid->x_min)};
        if (tmp<0 or tmp>grid->nx-1) { return 0.;}
        decltype(grid->nx) xl {(std::size_t)floor(tmp)};
        const double xd {tmp - xl};
        tmp = (grid->ny-1)*(pos.y-grid->y_min)/(grid->y_max-grid->y_min);
        if (tmp<0 or tmp>grid->ny-1) { return 0.;}
        decltype(grid->nx) yl {(std::size_t)floor(tmp)};
        const double yd {tmp - yl};
        tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
        if (tmp<0 or tmp>grid->nz-1) { return 0.;}
        decltype(grid->nx) zl {(std::size_t)floor(tmp)};
        const double zd {tmp - zl};
#ifdef DEBUG
        if(xd<0 or yd<0 or zd<0 or xd>1 or yd>1 or zd>1){
            ap_err("wrong value");
            exit(1);
        }
#endif
        double cre;
        if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
            std::size_t idx1 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl)};
            std::size_t idx2 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl+1)};
            const double i1 {grid->cre_flux[idx1]*(1.-zd) + grid->cre_flux[idx2]*zd};
            idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl+1,zl);
            idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl+1,zl+1);
            const double i2 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
            idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl,zl);
            idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl,zl+1);
            const double j1 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
            idx1 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl+1,zl);
            idx2 = toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl+1,yl+1,zl+1);
            const double j2 {grid->cre_flux[idx1]*(1-zd) + grid->cre_flux[idx2]*zd};
            const double w1 {i1*(1-yd)+i2*yd};
            const double w2 {j1*(1-yd)+j2*yd};
            cre = (w1*(1-xd)+w2*xd);
            
        }
        else {
            std::size_t idx1 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl)};
            cre = grid->cre_flux[idx1];
        }
#ifdef DEBUG
        if(cre<0){
            ap_err("negative CRE flux");
            exit(1);
        }
#endif
        return cre;
    }
    else{
        ap_err("wrong CRE grid dimension");
        exit(1);
    }
}

// J_tot(\nu)
double CRE_num::get_emissivity_t(const vec3_t<double> &pos,Param *par,Grid_cre *grid,const double &Bper){
    double J {0.};
#ifdef DEBUG
    if(!grid->read_permission){
        ap_err("no input");
        exit(1);
    }
#endif
    // allocate energy grid
    unique_ptr<double[]> KE = unique_ptr<double[]> (new double[grid->nE] {0.});
    // we need F(x[E]) and G(x[E]) in spectral integration
    unique_ptr<double[]> x = unique_ptr<double[]> (new double[grid->nE] {0.});
    unique_ptr<double[]> beta = unique_ptr<double[]> (new double[grid->nE] {0.});
    // consts used in loop, using cgs units
    const double x_fact {(2.*CGS_U_MEC*CGS_U_MEC2*CGS_U_MEC2*2.*CGS_U_pi*par->sim_freq)/(3.*CGS_U_qe*Bper)};
    // KE in cgs units
    for(decltype(grid->nE) i=0;i!=grid->nE;++i){
        KE[i] = grid->E_min*exp(i*grid->E_fact);
        x[i] = x_fact/(KE[i]*KE[i]);
        beta[i] = sqrt(1-CGS_U_MEC2/KE[i]);
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and density n
    // ref: "Cosmic rays n' particle physics", A3
    const double fore_factor {4.*CGS_U_pi*sqrt(3.)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*abs(Bper)/(CGS_U_MEC2*CGS_U_C_light*CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec)};
    for(decltype(grid->nE) i=0;i!=grid->nE-1;++i){
        const double xv {(x[i+1]+x[i])/2.};
        // avoid underflow in gsl functions
        if(xv>100) {continue;}
        const double dE {fabs(KE[i+1]-KE[i])};
        // we put beta here
        const double de {(read_grid(i+1,pos,grid)/beta[i+1]+read_grid(i,pos,grid)/beta[i])/2.};
#ifdef DEBUG
        if(de<0){
            ap_err("negative CRE density");
            exit(1);
        }
#endif
        J += gsl_sf_synchrotron_1(xv)*de*dE;
    }
    return fore_factor*J/(4.*CGS_U_pi);
}

// J_pol(\nu)
double CRE_num::get_emissivity_p(const vec3_t<double> &pos,Param *par,Grid_cre *grid,const double &Bper){
    double J {0.};
#ifdef DEBUG
    if(!grid->read_permission){
        ap_err("no input");
        exit(1);
    }
#endif
    // allocate energy grid
    unique_ptr<double[]> KE = unique_ptr<double[]> (new double[grid->nE] {0.});
    // we need F(x[E]) and G(x[E]) in spectral integration
    unique_ptr<double[]> x = unique_ptr<double[]> (new double[grid->nE] {0.});
    unique_ptr<double[]> beta = unique_ptr<double[]> (new double[grid->nE] {0.});
    // consts used in loop, using cgs units
    const double x_fact {(2.*CGS_U_MEC*CGS_U_MEC2*CGS_U_MEC2*2.*CGS_U_pi*par->sim_freq)/(3.*CGS_U_qe*Bper)};
    // KE in cgs units
    for(decltype(grid->nE) i=0;i!=grid->nE;++i){
        KE[i] = grid->E_min*exp(i*grid->E_fact);
        x[i] = x_fact/(KE[i]*KE[i]);
        beta[i] = sqrt(1-CGS_U_MEC2/KE[i]);
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and density n
    // ref: "Cosmic rays n' particle physics", A3
    const double fore_factor {4.*CGS_U_pi*sqrt(3.)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*abs(Bper)/(CGS_U_MEC2*CGS_U_C_light*CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec)};
    for(decltype(grid->nE) i=0;i!=grid->nE-1;++i){
        const double xv {(x[i+1]+x[i])/2.};
        // avoid underflow in gsl functions
        if(xv>100) {continue;}
        const double dE {fabs(KE[i+1]-KE[i])};
        // we put beta here
        const double de {(read_grid(i+1,pos,grid)/beta[i+1]+read_grid(i,pos,grid)/beta[i])/2.};
#ifdef DEBUG
        if(de<0){
            ap_err("negative CRE density");
            exit(1);
        }
#endif
        J += gsl_sf_synchrotron_2(xv)*de*dE;
    }
    return fore_factor*J/(4.*CGS_U_pi);
}
