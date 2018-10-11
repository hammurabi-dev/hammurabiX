#include <memory>
#include <vector>
#include <array>
#include <cmath>
#include <cassert>

#include <vec3.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_sf_gamma.h>

#include <cre.h>
#include <param.h>
#include <grid.h>
#include <cgs_units_file.h>

// numerical CRE flux
// J_tot(\nu)
double CRE_num::get_emissivity_t (const vec3_t<double> &pos,
                                  const Param *par,
                                  const Grid_cre *grid,
                                  const double &Bper) const{
    double J {0.};
    assert(grid->read_permission);
    // allocate energy grid
    std::unique_ptr<double[]> KE = std::make_unique<double[]>(grid->nE);
    // we need F(x[E]) and G(x[E]) in spectral integration
    std::unique_ptr<double[]> x = std::make_unique<double[]>(grid->nE);
    std::unique_ptr<double[]> beta = std::make_unique<double[]>(grid->nE);
    // consts used in loop, using cgs units
    const double x_fact {(2.*CGS_U_MEC*CGS_U_MEC2*CGS_U_MEC2*2.*CGS_U_pi*par->sim_freq)/(3.*CGS_U_qe*Bper)};
    // KE in cgs units
    for (decltype(grid->nE) i=0;i!=grid->nE;++i){
        KE[i] = grid->E_min*std::exp(i*grid->E_fact);
        x[i] = x_fact/(KE[i]*KE[i]);
        beta[i] = std::sqrt(1-CGS_U_MEC2/KE[i]);
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and density n
    // ref: "Cosmic rays n' particle physics", A3
    const double fore_factor {4.*CGS_U_pi*std::sqrt(3.)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*std::fabs(Bper)/(CGS_U_MEC2*CGS_U_C_light*CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec)};
    for (decltype(grid->nE) i=0;i!=grid->nE-1;++i){
        const double xv {(x[i+1]+x[i])/2.};
        // avoid underflow in gsl functions
        if(xv>100) {continue;}
        const double dE {std::fabs(KE[i+1]-KE[i])};
        // we put beta here
        const double de {(read_grid(i+1,pos,grid)/beta[i+1]+read_grid(i,pos,grid)/beta[i])/2.};
        assert(de>=0);
        J += gsl_sf_synchrotron_1(xv)*de*dE;
    }
    return fore_factor*J/(4.*CGS_U_pi);
}

// J_pol(\nu)
double CRE_num::get_emissivity_p (const vec3_t<double> &pos,
                                  const Param *par,
                                  const Grid_cre *grid,
                                  const double &Bper) const{
    double J {0.};
    assert (grid->read_permission);
    // allocate energy grid
    std::unique_ptr<double[]> KE = std::make_unique<double[]>(grid->nE);
    // we need F(x[E]) and G(x[E]) in spectral integration
    std::unique_ptr<double[]> x = std::make_unique<double[]>(grid->nE);
    std::unique_ptr<double[]> beta = std::make_unique<double[]>(grid->nE);
    // consts used in loop, using cgs units
    const double x_fact {(2.*CGS_U_MEC*CGS_U_MEC2*CGS_U_MEC2*2.*CGS_U_pi*par->sim_freq)/(3.*CGS_U_qe*Bper)};
    // KE in cgs units
    for (decltype(grid->nE) i=0;i!=grid->nE;++i){
        KE[i] = grid->E_min*std::exp(i*grid->E_fact);
        x[i] = x_fact/(KE[i]*KE[i]);
        beta[i] = std::sqrt(1-CGS_U_MEC2/KE[i]);
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and density n
    // ref: "Cosmic rays n' particle physics", A3
    const double fore_factor {4.*CGS_U_pi*std::sqrt(3.)*(CGS_U_qe*CGS_U_qe*CGS_U_qe)*abs(Bper)/(CGS_U_MEC2*CGS_U_C_light*CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec)};
    for (decltype(grid->nE) i=0;i!=grid->nE-1;++i){
        const double xv {(x[i+1]+x[i])/2.};
        // avoid underflow in gsl functions
        if(xv>100) {continue;}
        const double dE {std::fabs(KE[i+1]-KE[i])};
        // we put beta here
        const double de {(read_grid(i+1,pos,grid)/beta[i+1]+read_grid(i,pos,grid)/beta[i])/2.};
        assert(de>=0);
        J += gsl_sf_synchrotron_2(xv)*de*dE;
    }
    return fore_factor*J/(4.*CGS_U_pi);
}

// END
