#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include "class_cre.h"
#include "class_pond.h"
#include "class_grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

double CRE::get_emissivity(const vec3 &,Pond *,Grid_cre *,const double &,const bool &){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

double CRE::read_grid(const unsigned int &,const vec3 &,Grid_cre *){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

void CRE::write_grid(Pond *,Grid_cre *){
    cerr<<"ERR:"<<__FILE__
    <<" : in function "<<__func__<<endl
    <<" at line "<<__LINE__<<endl
    <<"DYNAMIC BINDING FAILURE"<<endl;
    exit(1);
}

// analytical CRE flux
// give values to spectral index and norm factor, in cgs units
// analytical CRE integration use N(\gamma)
void CRE_ana::flux_param(const vec3 &pos,Pond *par,double &index,double &norm){
    // units
    const double alpha {par->creana[0]};
    const double beta {par->creana[1]};
    const double theta {par->creana[2]};
    const double hr {par->creana[3]*CGS_U_kpc};
    const double hz {par->creana[4]*CGS_U_kpc};
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->creana[5]};
    const double R0 {sqrt(par->SunPosition.x*par->SunPosition.x+par->SunPosition.y*par->SunPosition.y)};
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)};
    const double z {fabs(pos.z)};
    const double cre_gamma_10 {10.*CGS_U_GeV/CGS_U_MEC2};
    const double cre_beta_10 {sqrt(1.-1./cre_gamma_10)};
    // from PHI(E) to N(\gamma) convertion
    const double unit_factor {(4.*CGS_U_pi*CGS_U_MEC)/(CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec*cre_beta_10)};
    // MODEL DEPENDENT PARAMETERS
    const double norm_factor {je*pow(cre_gamma_10,alpha-beta*R0)*exp(R0/hr)};
    const double scal_factor {exp(-r/hr)*(1./pow(cosh(z/hz),2.))};
    // this is changeable by users
    index = -alpha+beta*r+theta*z;
    
    norm = norm_factor*scal_factor*unit_factor;
}

// analytical modeling use N(\gamma) while flux is PHI(E)
// En in CGS units, return in [GeV m^2 s Sr]^-1
double CRE_ana::flux(const vec3 &pos,Pond *par,const double &En){
    // units
    const double alpha {par->creana[0]};
    const double beta {par->creana[1]};
    const double theta {par->creana[2]};
    const double hr {par->creana[3]*CGS_U_kpc};
    const double hz {par->creana[4]*CGS_U_kpc};
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->creana[5]};
    const double R0 {sqrt(par->SunPosition.x*par->SunPosition.x+par->SunPosition.y*par->SunPosition.y)};
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)};
    const double z {fabs(pos.z)};
    const double gamma {En/CGS_U_MEC2};
    const double cre_gamma_10 {10.*CGS_U_GeV/CGS_U_MEC2};
    // converting from N to PHI
    const double unit_factor {sqrt((1.-1./gamma)/(1.-1./cre_gamma_10))};
    // MODEL DEPENDENT PARAMETERS
    // CRE flux normalizaton factor at earth, model dependent
    const double norm_factor {je*pow(cre_gamma_10,alpha-beta*R0)*exp(R0/hr)};
    const double scal_factor {exp(-r/hr)*(1./pow(cosh(z/hz),2.))};
    const double index {-alpha+beta*r+theta*z};
    
    return norm_factor*scal_factor*unit_factor*pow(gamma,index);
}

double CRE_ana::get_emissivity(const vec3 &pos,Pond *par,Grid_cre *grid,const double &Bper,const bool &cue){
    if(grid->read_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG MODULE"<<endl;
        exit(1);
    }
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    double index, norm;
    flux_param(pos,par,index,norm);
    // coefficients which do not attend integration
    norm *= sqrt(3)*pow(CGS_U_qe,3)*fabs(Bper)/(4*CGS_U_pi*CGS_U_MEC2);
    // synchrotron integration
    const double A {sqrt(2.*CGS_U_MEC*2.*CGS_U_pi*par->sim_freq)/sqrt(3.*CGS_U_qe*fabs(Bper))};
    const double mu {-(3.+index)/2.0};
    
    if(cue){
        double J {pow(A,index+1)*pow(2,mu+1)*gsl_sf_gamma(0.5*mu+7./3.)*gsl_sf_gamma(0.5*mu+2./3.)/(mu+2.)};
        return norm*J/(4.*CGS_U_pi);
    }
    else{
        double J {pow(A,index+1)*pow(2,mu)*gsl_sf_gamma(0.5*mu+4./3.)*gsl_sf_gamma(0.5*mu+2./3.)};
        return norm*J/(4.*CGS_U_pi);
    }
    /* the last 4pi comes from solid-angle integration/deviation,
     check eq(6.16) in Ribiki-Lightman's where Power is defined,
     we need isotropic power which means we need a 1/4pi factor!
     */
}

// writing out CRE DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
void CRE_ana::write_grid(Pond *par, Grid_cre *grid){
    if(!grid->write_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NO PERMISSION"<<endl;
        exit(1);
    }
    vec3 gc_pos;
    // 2D grid
    if(grid->nr!=0){
        double lr {grid->r_max};
        double lz {grid->z_max-grid->z_min};
        for(decltype(grid->nE) i=0;i!=grid->nE;++i){
            for(decltype(grid->nr) j=0;j!=grid->nr;++j){
                for(decltype(grid->nz) k=0;k!=grid->nz;++k){
                    // on y=0 2D slide
                    gc_pos = vec3 {lr*j/(grid->nr-1),
                        0,
                        lz*k/(grid->nz-1) + grid->z_min};
                    double E {grid->E_min*exp(i*grid->E_fact)};
                    unsigned long int idx {toolkit::Index3d(grid->nE,grid->nr,grid->nz,i,j,k)};
                    grid->cre_flux[idx] = flux(gc_pos,par,E);
                }
            }
        }
    }
    // 3D grid
    else if(grid->nx!=0){
        double lx {grid->x_max-grid->x_min};
        double ly {grid->y_max-grid->y_min};
        double lz {grid->z_max-grid->z_min};
        for(decltype(grid->nE) i=0;i!=grid->nE;++i){
            for(decltype(grid->nx) j=0;j!=grid->nx;++j){
                for(decltype(grid->ny) k=0;k!=grid->ny;++k){
                    for(decltype(grid->nz) m=0;m!=grid->nz;++m){
                        gc_pos = vec3 {lx*j/(grid->nx-1) + grid->x_min,
                            ly*k/(grid->ny-1) + grid->y_min,
                            lz*m/(grid->nz-1) + grid->z_min};
                        double E {grid->E_min*exp(i*grid->E_fact)};
                        unsigned long int idx {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,i,j,k,m)};
                        grid->cre_flux[idx] = flux(gc_pos,par,E);
                    }
                }
            }
        }
    }
    
}

// numerical CRE flux
// use bilinear/trilinear interpolationi according to the dimension of CRE flux grid
double CRE_num::read_grid(const unsigned int &Eidx, const vec3 &pos,Grid_cre *grid){
    // if grid in spatial 2D
    if(grid->nr!=0){
        // sylindrical galactic centric position
        const double r {sqrt(pos.x*pos.x+pos.y*pos.y)};
        // bilinear interpolation
        // notice that lr is radius
        double tmp {(grid->nr-1)*(r/grid->r_max)};
        if(tmp<0 or tmp>grid->nr-1) {return 0.;}
        decltype(grid->nr) rl {(unsigned int)floor(tmp)};
        const double rd {tmp - rl};
        tmp = (grid->nz-1)*(pos.z-grid->z_min)/(grid->z_max-grid->z_min);
        if(tmp<0 or tmp>grid->nz-1) {return 0.;}
        decltype(grid->nr) zl {(unsigned int)floor(tmp)};
        const double zd {tmp - zl};
#ifndef NDEBUG
        if(rd<0 or zd<0 or rd>1 or zd>1){
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"WRONG VALUE: "<<endl;
            exit(1);
        }
#endif
        double cre;
        if(rl+1<grid->nr and zl+1<grid->nz){
            unsigned long int idx1 {toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl)};
            unsigned long int idx2 {toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl+1,zl)};
            double i1 {grid->cre_flux[idx1]*(1-rd) + grid->cre_flux[idx2]*rd};
            idx1 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl+1);
            idx2 = toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl+1,zl+1);
            double i2 {grid->cre_flux[idx1]*(1-rd) + grid->cre_flux[idx2]*rd};
            cre = (i1*(1-zd)+i2*zd);
        }
        else{
            unsigned long int idx1 {toolkit::Index3d(grid->nE,grid->nr,grid->nz,Eidx,rl,zl)};
            cre = grid->cre_flux[idx1];
        }
#ifndef NDEBUG
        if(cre<0){
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"NEGATIVE CRE FLUX"<<endl;
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
            <<"WRONG VALUE: "<<endl;
            exit(1);
        }
#endif
        double cre;
        if (xl+1<grid->nx and yl+1<grid->ny and zl+1<grid->nz) {
            unsigned long int idx1 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl)};
            unsigned long int idx2 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl+1)};
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
            unsigned long int idx1 {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,Eidx,xl,yl,zl)};
            cre = grid->cre_flux[idx1];
        }
#ifndef NDEBUG
        if(cre<0){
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"NEGATIVE CRE FLUX"<<endl;
            exit(1);
        }
#endif
        return cre;
    }
    else{
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG CRE GRID DIMENSION"<<endl;
        exit(1);
    }
}

double CRE_num::get_emissivity(const vec3 &pos,Pond *par,Grid_cre *grid,const double &Bper,const bool &cue){
    double J {0.};
    if(!grid->read_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NO INPUT"<<endl;
        exit(1);
    }
    // allocate energy grid
    double KE[grid->nE] {0.};
    // we need F(x[E]) and G(x[E]) in spectral integration
    double x[grid->nE] {0.};
    double beta[grid->nE] {0.};
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
    const double fore_factor {2.*sqrt(3.)*pow(CGS_U_qe,3.)*abs(Bper)/(CGS_U_MEC2*CGS_U_C_light*CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec)};
    for(decltype(grid->nE) i=0;i!=grid->nE-1;++i){
        const double xv {(x[i+1]+x[i])/2.};
        // avoid underflow in gsl functions
        if(xv>100) {continue;}
        const double dE {fabs(KE[i+1]-KE[i])};
        // we put beta here
        const double de {(read_grid(i+1,pos,grid)/beta[i+1]+read_grid(i,pos,grid)/beta[i])/2.};
#ifndef NDEBUG
        if(de<0){
            cerr<<"ERR:"<<__FILE__
            <<" : in function "<<__func__<<endl
            <<" at line "<<__LINE__<<endl
            <<"NEGATIVE CRE DENSITY"<<endl;
            exit(1);
        }
#endif
        if(cue){
            J += gsl_sf_synchrotron_1(xv)*de*dE;
        }
        else{
            J += gsl_sf_synchrotron_2(xv)*de*dE;
        }
    }
    return fore_factor*J/(4.*CGS_U_pi);
}

// END
