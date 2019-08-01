#include <cassert>
#include <cmath>

#include <gsl/gsl_sf_gamma.h>

#include <cgs_units_file.h>
#include <crefield.h>
#include <grid.h>
#include <hamvec.h>
#include <param.h>

// CRE flux spatial reprofiling
double CRE_unif::spatial_profile(const hamvec<3, double> &pos,
                                 const Param *par) const {
  if ((pos - par->observer).length() > par->cre_unif.r0) {
    return 0;
  } else
    return 1;
}

// CRE spectral index
double CRE_unif::flux_idx(const hamvec<3, double> & /*pos*/,
                          const Param *par) const {
  return -(par->cre_unif.alpha);
}

// analytical CRE flux
// give values to spectral index and norm factor, in cgs units
// analytical CRE integration use N(\gamma)
double CRE_unif::flux_norm(const hamvec<3, double> &pos,
                           const Param *par) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  const double gamma0{par->cre_unif.E0 / CGS_U_MEC2 + 1};
  const double beta0{std::sqrt(1. - 1. / gamma0)};
  // from PHI(E) to N(\gamma) convertion
  const double unit{
      (4. * CGS_U_pi * CGS_U_MEC) /
      (CGS_U_GeV * 100. * CGS_U_cm * 100. * CGS_U_cm * CGS_U_sec * beta0)};
  const double norm{par->cre_unif.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * unit * spatial_profile(pos, par);
}

// analytical modeling use N(\gamma) while flux is PHI(E)
// En in CGS units, return in [GeV m^2 s Sr]^-1
double CRE_unif::flux(const hamvec<3, double> &pos, const Param *par,
                      const double &En) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  const double gamma{En / CGS_U_MEC2};
  const double gamma0{par->cre_unif.E0 / CGS_U_MEC2 + 1};
  // converting from N to PHI
  const double unit{std::sqrt((1. - 1. / gamma) / (1. - 1. / gamma0))};
  const double norm{par->cre_unif.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * unit * std::pow(gamma, flux_idx(pos, par)) *
         spatial_profile(pos, par);
}

// J_tot(\nu)
double CRE_unif::read_emissivity_t(const hamvec<3, double> &pos,
                                   const Param *par, const Grid_cre * /*grid*/,
                                   const double &Bper) const {
  // assert(!par->grid_cre.read_permission);
  // allocating values to index, norm according to user defined model
  // user may consider building derived class from CRE_ana
  const double index{flux_idx(pos, par)};
  // coefficients which do not attend integration
  const double norm{flux_norm(pos, par) * std::sqrt(3) *
                    (CGS_U_qe * CGS_U_qe * CGS_U_qe) * std::fabs(Bper) /
                    (2. * CGS_U_MEC2)};
  // synchrotron integration
  const double A{4. * CGS_U_MEC * CGS_U_pi *
                 par->grid_obs.sim_sync_freq.back() /
                 (3. * CGS_U_qe * std::fabs(Bper))};
  const double mu{-0.5 * (3. + index)};
  return norm *
         (std::pow(A, 0.5 * (index + 1)) * std::pow(2, mu + 1) *
          gsl_sf_gamma(0.5 * mu + 7. / 3.) * gsl_sf_gamma(0.5 * mu + 2. / 3.) /
          (mu + 2.)) /
         (4. * CGS_U_pi);
  // the last 4pi comes from solid-angle integration/deviation,
  // check eq(6.16) in Ribiki-Lightman's where Power is defined,
  // we need isotropic power which means we need a 1/4pi factor!
}

// J_pol(\nu)
double CRE_unif::read_emissivity_p(const hamvec<3, double> &pos,
                                   const Param *par, const Grid_cre * /*grid*/,
                                   const double &Bper) const {
  // assert(!par->grid_cre.read_permission);
  // allocating values to index, norm according to user defined model
  // user may consider building derived class from CRE_ana
  const double index{flux_idx(pos, par)};
  // coefficients which do not attend integration
  const double norm{flux_norm(pos, par) * std::sqrt(3) *
                    (CGS_U_qe * CGS_U_qe * CGS_U_qe) * std::fabs(Bper) /
                    (2. * CGS_U_MEC2)};
  // synchrotron integration
  const double A{4. * CGS_U_MEC * CGS_U_pi *
                 par->grid_obs.sim_sync_freq.back() /
                 (3. * CGS_U_qe * std::fabs(Bper))};
  const double mu{-0.5 * (3. + index)};
  return norm *
         (std::pow(A, 0.5 * (index + 1)) * std::pow(2, mu) *
          gsl_sf_gamma(0.5 * mu + 4. / 3.) * gsl_sf_gamma(0.5 * mu + 2. / 3.)) /
         (4. * CGS_U_pi);
  // the last 4pi comes from solid-angle integration/deviation,
  // check eq(6.16) in Ribiki-Lightman's where Power is defined,
  // we need isotropic power which means we need a 1/4pi factor!
}
