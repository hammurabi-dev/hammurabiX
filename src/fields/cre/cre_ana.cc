#include <iostream>
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
#include <cassert>


// CRE flux spatial rescaling
double CRE_ana::rescal(const vec3_t<double> &pos,
                       Param *par){
    const double r0 {par->cre_ana.r0};
    const double z0 {par->cre_ana.z0};
    const double R0 {sqrt(par->SunPosition.x*par->SunPosition.x+par->SunPosition.y*par->SunPosition.y)};
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)};
    return exp((R0-r)/r0)*(1./(cosh(pos.z/z0)*cosh(pos.z/z0)));
}

// CRE spectral index
double CRE_ana::flux_idx(const vec3_t<double> &pos,
                         Param *par){
    const double alpha {par->cre_ana.alpha};
    const double beta {par->cre_ana.beta};
    const double theta {par->cre_ana.theta};
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)};
    const double z {fabs(pos.z)};
    return -alpha+beta*r+theta*z;
}

// analytical CRE flux normalization factor at E0
// analytical CRE spectral integrations use N(\gamma)
double CRE_ana::flux_norm(const vec3_t<double> &pos,
                          Param *par){
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->cre_ana.j0};
    const double gamma0 {par->cre_ana.E0/CGS_U_MEC2+1};
    const double beta0 {sqrt(1.-1./gamma0)};
    // from PHI(E) to N(\gamma) convertion
    const double unit {(4.*CGS_U_pi*CGS_U_MEC)/(CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec*beta0)};
    const double norm {je*pow(gamma0,-flux_idx(par->SunPosition,par))};
    
    return norm*unit*rescal(pos,par);
}

// analytical modelings use N(\gamma) while flux is PHI(E)
// En in CGS units, return in [GeV m^2 s Sr]^-1
double CRE_ana::flux(const vec3_t<double> &pos,
                     Param *par,
                     const double &En){
    // units
    // je is in [GeV m^2 s sr]^-1 units
    const double je {par->cre_ana.j0};
    const double gamma {En/CGS_U_MEC2};
    const double gamma0 {par->cre_ana.E0/CGS_U_MEC2+1};
    // converting from N to PHI
    const double unit {sqrt((1.-1./gamma)/(1.-1./gamma0))};
    const double norm {je*pow(gamma0,-flux_idx(par->SunPosition,par))};
    
    return norm*unit*pow(gamma,flux_idx(pos,par))*rescal(pos,par);
}

// J_tot(\nu)
double CRE_ana::get_emissivity_t(const vec3_t<double> &pos,
                                 Param *par,
                                 Grid_cre */*grid*/,
                                 const double &Bper){
    //assert(!grid->read_permission);
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    const double index {flux_idx(pos,par)};
    // coefficients which do not attend integration
    const double norm {flux_norm(pos,par)*sqrt(3)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*fabs(Bper)/(2.*CGS_U_MEC2)};
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
double CRE_ana::get_emissivity_p(const vec3_t<double> &pos,
                                 Param *par,
                                 Grid_cre */*grid*/,
                                 const double &Bper){
    //assert(!grid->read_permission);
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    const double index {flux_idx(pos,par)};
    // coefficients which do not attend integration
    const double norm {flux_norm(pos,par)*sqrt(3)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*fabs(Bper)/(2.*CGS_U_MEC2)};
    // synchrotron integration
    const double A {4.*CGS_U_MEC*CGS_U_pi*par->sim_freq/(3.*CGS_U_qe*fabs(Bper))};
    const double mu {-0.5*(3.+index)};
    return norm*( pow(A,0.5*(index+1))*pow(2,mu)*gsl_sf_gamma(0.5*mu+4./3.)*gsl_sf_gamma(0.5*mu+2./3.) )/(4.*CGS_U_pi);
    /* the last 4pi comes from solid-angle integration/deviation,
     check eq(6.16) in Ribiki-Lightman's where Power is defined,
     we need isotropic power which means we need a 1/4pi factor!
     */
}
