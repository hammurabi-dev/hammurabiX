#include <cassert>
#include <cmath>

#include <cgs_units.h>
#include <crefield.h>
#include <hamvec.h>
#include <param.h>

// CRE flux spatial rescaling
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
double CRE_ana::flux_idx(const hamvec<3, double> &pos, const Param *par) const {
  const double r{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]) / cgs_kpc};
  const double z{std::fabs(pos[2]) / cgs_kpc};
  return -par->cre_ana.alpha + par->cre_ana.beta * r + par->cre_ana.theta * z;
}

// analytical CRE flux normalization factor at E0
// analytical CRE spectral integrations use N(\gamma)
double CRE_ana::flux_norm(const hamvec<3, double> &pos,
                          const Param *par) const {
  // j0 is in [GeV m^2 s sr]^-1 units
  const double gamma0{par->cre_ana.E0 / cgs_mec2 + 1};
  const double beta0{std::sqrt(1. - 1. / gamma0)};
  // from PHI(E) to N(\gamma) convertion
  const double unit{(4. * cgs_pi * cgs_mec) / (cgs_GeV * 100. * cgs_cm * 100. *
                                               cgs_cm * cgs_sec * beta0)};
  const double norm{par->cre_ana.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};

  return norm * unit * spatial_profile(pos, par);
}

// analytical modelings use N(\gamma) while flux is PHI(E)
// En in CGS units, return in [GeV m^2 s Sr]^-1
double CRE_ana::write_field(const hamvec<3, double> &pos, const double &En,
                            const Param *par) const {
  // units
  // j0 is in [GeV m^2 s sr]^-1 units
  const double gamma{En / cgs_mec2};
  const double gamma0{par->cre_ana.E0 / cgs_mec2 + 1};
  // converting from N to PHI
  const double unit{std::sqrt((1. - 1. / gamma) / (1. - 1. / gamma0))};
  const double norm{par->cre_ana.j0 *
                    std::pow(gamma0, -flux_idx(par->observer, par))};

  return norm * unit * std::pow(gamma, flux_idx(pos, par)) *
         spatial_profile(pos, par);
}
