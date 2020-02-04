#include <cassert>
#include <cmath>

#include <cgsunits.h>
#include <crefield.h>
#include <hamvec.h>
#include <param.h>

// CRE flux spatial profiling (dimensionless)
double CRE_unif::spatial_profile(const hamvec<3, double> &pos,
                                 const Param *par) const {
  if ((pos - par->observer).length() > par->cre_unif.r0) {
    return 0;
  } else
    return 1;
}

// CRE flux spectral index
// for N(\gamma) = flux_norm * \gamma^{flux_idx}
double CRE_unif::flux_idx(const hamvec<3, double> & /*pos*/,
                          const Param *par) const {
  return -(par->cre_unif.alpha);
}

// CRE flux norm factor in cgs units
// the so called flux norm is actually normalization factor for N(\gamma)
// which satisfies N(\gamma) = N_0 * \gamma^(flux_idx)
double CRE_unif::flux_norm(const hamvec<3, double> &pos,
                           const Param *par) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  // note that the j0 has already been defined by per steradian
  // which is the conventional units for PHI(E)
  const double gamma0{par->cre_unif.E0 / cgs::mec2};
  const double beta0{std::sqrt(1. - 1. / gamma0)};
  // N(\gamma) = (4\pi mc/\beta) * PHI(E) by definition
  // the extra 'unit' coefficients come from the definition of 'j0'
  const double unit{4. * cgs::pi * cgs::mec /
                    (cgs::GeV * cgs::m * cgs::m * cgs::sec * beta0)};
  const double norm{par->cre_unif.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * unit * spatial_profile(pos, par);
}

// return the actual PHI(E) according to the definition of N(\gamma)
// PHI(E) = (\beta/4\pi mc) N(\gamma)
// En in CGS units, return in units [GeV m^2 s Sr]^-1
double CRE_unif::write_field(const hamvec<3, double> &pos, const double &En,
                             const Param *par) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  // note that the j0 has already been defined by per steradian
  const double gamma{En / cgs::mec2};
  const double gamma0{par->cre_unif.E0 / cgs::mec2};
  const double beta_ratio{std::sqrt((1. - 1. / gamma) / (1. - 1. / gamma0))};
  const double norm{par->cre_unif.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * beta_ratio * std::pow(gamma, flux_idx(pos, par)) *
         spatial_profile(pos, par);
}
