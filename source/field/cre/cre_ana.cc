#include <cassert>
#include <cmath>

#include <cgs_units.h>
#include <crefield.h>
#include <hamvec.h>
#include <param.h>

// CRE flux spatial profiling (dimensionless)
double CRE_ana::spatial_profile(const hamvec<3, double> &pos,
                                const Param *par) const {
  const double R0{std::sqrt(par->observer[0] * par->observer[0] +
                            par->observer[1] * par->observer[1])};
  const double r{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  return std::exp((R0 - r) / par->cre_ana.r0) *
         (1. /
          (cosh(pos[2] / par->cre_ana.z0) * cosh(pos[2] / par->cre_ana.z0)));
}

// CRE spectral index
// for N(\gamma) = flux_norm * \gamma^{flux_idx}
double CRE_ana::flux_idx(const hamvec<3, double> &pos, const Param *par) const {
  const double r{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]) / cgs_kpc};
  const double z{std::fabs(pos[2]) / cgs_kpc};
  return -par->cre_ana.alpha + par->cre_ana.beta * r + par->cre_ana.theta * z;
}

// CRE flux norm factor in cgs units
// the so called flux norm is actually normalization N_0 for N(\gamma)
// which satisfies N(\gamma) = N_0 * \gamma^(flux_idx)
double CRE_ana::flux_norm(const hamvec<3, double> &pos,
                          const Param *par) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  // note that the j0 has already been defined by per steradian
  // which is the conventional units for PHI(E)
  const double gamma0{par->cre_ana.E0 / cgs_mec2};
  const double beta0{std::sqrt(1. - 1. / gamma0)};
  // N(\gamma) = (4\pi mc/\beta) * PHI(E) by definition
  // the extra 'unit' coefficients come from the definition of 'j0'
  const double unit{4. * cgs_pi * cgs_mec /
                    (cgs_GeV * cgs_m * cgs_m * cgs_sec * beta0)};
  const double norm{par->cre_ana.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * unit * spatial_profile(pos, par);
}

// return the actual PHI(E) according to the definition of N(\gamma)
// PHI(E) = (\beta/4\pi mc) N(\gamma)
// En in CGS units, return in units [GeV m^2 s Sr]^-1
double CRE_ana::write_field(const hamvec<3, double> &pos, const double &En,
                            const Param *par) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  // note that the j0 has already been defined by per steradian
  const double gamma{En / cgs_mec2};
  const double gamma0{par->cre_ana.E0 / cgs_mec2};
  const double beta_ratio{std::sqrt((1. - 1. / gamma) / (1. - 1. / gamma0))};
  const double norm{par->cre_ana.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};
  return norm * beta_ratio * std::pow(gamma, flux_idx(pos, par)) *
         spatial_profile(pos, par);
}
