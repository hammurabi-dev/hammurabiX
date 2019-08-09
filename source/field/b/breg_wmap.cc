#include <cmath>

#include <bfield.h>
#include <cgs_units.h>
#include <hamvec.h>
#include <param.h>

hamvec<3, double> Breg_wmap::write_field(const hamvec<3, double> &pos,
                                         const Param *par) const {
  const double r{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  if (r > (20. * cgs_kpc) or r < (3. * cgs_kpc)) {
    return hamvec<3, double>{0., 0., 0.};
  }
  const double cos_p{pos[0] / r};
  const double sin_p{pos[1] / r};
  const double psi{par->breg_wmap.psi0 +
                   par->breg_wmap.psi1 * std::log(r / (8. * cgs_kpc))};
  const double chi{par->breg_wmap.chi0 * std::tanh(pos[2] / (1. * cgs_kpc))};
  return hamvec<3, double>{
      par->breg_wmap.b0 * (std::sin(psi) * std::cos(chi) * cos_p -
                           std::cos(psi) * std::cos(chi) * sin_p),
      par->breg_wmap.b0 * (std::sin(psi) * std::cos(chi) * sin_p +
                           std::cos(psi) * std::cos(chi) * cos_p),
      par->breg_wmap.b0 * std::sin(chi)};
}
