#include <cassert>
#include <cmath>

#include <cgs_units.h>
#include <crefield.h>
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
  const double gamma0{par->cre_unif.E0 / cgs_mec2 + 1};
  const double beta0{std::sqrt(1. - 1. / gamma0)};
  // from PHI(E) to N(\gamma) convertion
  const double unit{(4. * cgs_pi * cgs_mec) / (cgs_GeV * 100. * cgs_cm * 100. *
                                               cgs_cm * cgs_sec * beta0)};
  const double norm{par->cre_unif.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * unit * spatial_profile(pos, par);
}

// analytical modeling use N(\gamma) while flux is PHI(E)
// En in CGS units, return in [GeV m^2 s Sr]^-1
double CRE_unif::write_field(const hamvec<3, double> &pos, const double &En,
                             const Param *par) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  const double gamma{En / cgs_mec2};
  const double gamma0{par->cre_unif.E0 / cgs_mec2 + 1};
  // converting from N to PHI
  const double unit{std::sqrt((1. - 1. / gamma) / (1. - 1. / gamma0))};
  const double norm{par->cre_unif.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * unit * std::pow(gamma, flux_idx(pos, par)) *
         spatial_profile(pos, par);
}
