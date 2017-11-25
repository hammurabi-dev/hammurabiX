#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include "cre.h"
#include "param.h"
#include "grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

// analytical CRE flux
// give values to spectral index and norm factor, in cgs units
// analytical CRE integration use N(\gamma)
void CRE_verify::flux_param(const vec3_t<double> &/*pos*/,Param *par,double &index,double &norm){
    // units
    const double alpha {par->cre_verify.alpha};
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->cre_verify.je};
    const double cre_gamma_10 {10.*CGS_U_GeV/CGS_U_MEC2};
    const double cre_beta_10 {sqrt(1.-1./cre_gamma_10)};
    // from PHI(E) to N(\gamma) convertion
    const double unit_factor {(4.*CGS_U_pi*CGS_U_MEC)/(CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec*cre_beta_10)};
    const double norm_factor {je*pow(cre_gamma_10,alpha)};
    index = -alpha;
    norm = norm_factor*unit_factor;
}

// analytical modeling use N(\gamma) while flux is PHI(E)
// En in CGS units, return in [GeV m^2 s Sr]^-1
double CRE_verify::flux(const vec3_t<double> &pos,Param *par,const double &En){
    if((pos-par->SunPosition).Length() > par->cre_verify.r0){
        return 0.;
    }
    // units
    const double alpha {par->cre_verify.alpha};
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->cre_verify.je};
    const double gamma {En/CGS_U_MEC2};
    const double cre_gamma_10 {10.*CGS_U_GeV/CGS_U_MEC2};
    // converting from N to PHI
    const double unit_factor {sqrt((gamma-1.)*cre_gamma_10/((cre_gamma_10-1.)*gamma))};
    const double norm_factor {je*pow(cre_gamma_10,alpha)};
    return norm_factor*unit_factor*pow(gamma,-alpha);
}

// J_tot(\nu)
double CRE_verify::get_emissivity_t(const vec3_t<double> &pos,Param *par,Grid_cre *grid,const double &Bper){
#ifndef NDEBUG
    if(grid->read_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG MODULE"<<endl;
        exit(1);
    }
#endif
    if((pos-par->SunPosition).Length() > par->cre_verify.r0){
        return 0.;
    }
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    double index, norm;
    flux_param(pos,par,index,norm);
    // coefficients which do not attend integration
    norm *= sqrt(3)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*fabs(Bper)/(2.*CGS_U_MEC2);
    // synchrotron integration
    const double A {4.*CGS_U_MEC*CGS_U_pi*par->sim_freq/(3.*CGS_U_qe*fabs(Bper))};
    const double mu {-0.5*(3.+index)};
    return norm*( pow(A,0.5*(index+1))*pow(2,mu+1)*gsl_sf_gamma(0.5*mu+7./3.)*gsl_sf_gamma(0.5*mu+2./3.)/(mu+2.) )/(4.*CGS_U_pi);
    /* the last 4pi comes from solid-angle integration/deviation,
     check eq(6.16) in Ribiki-Lightman's where Power is defined,
     we need isotropic power which means we need a 1/4pi factor!
     */
}

// J_pol(\nu)
double CRE_verify::get_emissivity_p(const vec3_t<double> &pos,Param *par,Grid_cre *grid,const double &Bper){
#ifndef NDEBUG
    if(grid->read_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"WRONG MODULE"<<endl;
        exit(1);
    }
#endif
    if((pos-par->SunPosition).Length() > par->cre_verify.r0){
        return 0.;
    }
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    double index, norm;
    flux_param(pos,par,index,norm);
    // coefficients which do not attend integration
    norm *= sqrt(3)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*fabs(Bper)/(2.*CGS_U_MEC2);
    // synchrotron integration
    const double A {4.*CGS_U_MEC*CGS_U_pi*par->sim_freq/(3.*CGS_U_qe*fabs(Bper))};
    const double mu {-0.5*(3.+index)};
    return norm*( pow(A,0.5*(index+1))*pow(2,mu)*gsl_sf_gamma(0.5*mu+4./3.)*gsl_sf_gamma(0.5*mu+2./3.) )/(4.*CGS_U_pi);
    /* the last 4pi comes from solid-angle integration/deviation,
     check eq(6.16) in Ribiki-Lightman's where Power is defined,
     we need isotropic power which means we need a 1/4pi factor!
     */
}

// writing out CRE DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
void CRE_verify::write_grid(Param *par, Grid_cre *grid){
#ifndef NDEBUG
    if(!grid->write_permission){
        cerr<<"ERR:"<<__FILE__
        <<" : in function "<<__func__<<endl
        <<" at line "<<__LINE__<<endl
        <<"NO PERMISSION"<<endl;
        exit(1);
    }
#endif
    vec3_t<double> gc_pos {0.,0.,0.};
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
                    std::size_t idx {toolkit::Index3d(grid->nE,grid->nr,grid->nz,i,j,k)};
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
                        std::size_t idx {toolkit::Index4d(grid->nE,grid->nx,grid->ny,grid->nz,i,j,k,m)};
                        grid->cre_flux[idx] = flux(gc_pos,par,E);
                    }
                }
            }
        }
    }
    
}
