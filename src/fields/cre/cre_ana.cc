#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include "cre.h"
#include "pond.h"
#include "grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

// analytical CRE flux
// give values to spectral index and norm factor, in cgs units
// analytical CRE integration use N(\gamma)
void CRE_ana::flux_param(const vec3 &pos,Pond *par,double &index,double &norm){
    // units
    const double alpha {par->cre_ana.alpha};
    const double beta {par->cre_ana.beta};
    const double theta {par->cre_ana.theta};
    const double hr {par->cre_ana.hr*CGS_U_kpc};
    const double hz {par->cre_ana.hz*CGS_U_kpc};
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->cre_ana.je};
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
    const double alpha {par->cre_ana.alpha};
    const double beta {par->cre_ana.beta};
    const double theta {par->cre_ana.theta};
    const double hr {par->cre_ana.hr*CGS_U_kpc};
    const double hz {par->cre_ana.hz*CGS_U_kpc};
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->cre_ana.je};
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
    vec3 gc_pos {0.,0.,0.};
    // 2D grid
    if(grid->nr!=0){
        double lr {grid->r_max};
        double lz {grid->z_max-grid->z_min};
        for(decltype(grid->nE) i=0;i!=grid->nE;++i){
            double E {grid->E_min*exp(i*grid->E_fact)};
            for(decltype(grid->nr) j=0;j!=grid->nr;++j){
                gc_pos.x = lr*j/(grid->nr-1);
                for(decltype(grid->nz) k=0;k!=grid->nz;++k){
                    // on y=0 2D slide
                    gc_pos.z = lz*k/(grid->nz-1) + grid->z_min;
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
            double E {grid->E_min*exp(i*grid->E_fact)};
            for(decltype(grid->nx) j=0;j!=grid->nx;++j){
                gc_pos.x = lx*j/(grid->nx-1) + grid->x_min;
                for(decltype(grid->ny) k=0;k!=grid->ny;++k){
                    gc_pos.y = ly*k/(grid->ny-1) + grid->y_min;
                    for(decltype(grid->nz) m=0;m!=grid->nz;++m){
                        gc_pos.z = lz*m/(grid->nz-1) + grid->z_min;
                        unsigned long int idx {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,i,j,k,m)};
                        grid->cre_flux[idx] = flux(gc_pos,par,E);
                    }
                }
            }
        }
    }
    
}
