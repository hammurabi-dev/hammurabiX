#include <cmath>

#include <bfield.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>

Hamvec<3, ham_float> Breg_lsa::write_field(const Hamvec<3, ham_float> &pos,
                                           const Param *par) const {
  const ham_float r{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  if (r > (20. * cgs::kpc) or r < (3. * cgs::kpc)) {
    return Hamvec<3, ham_float>{0., 0., 0.};
  }
  const ham_float cos_p{pos[0] / r};
  const ham_float sin_p{pos[1] / r};
  const ham_float psi{par->breg_lsa.psi0 +
                      par->breg_lsa.psi1 * std::log(r / (8. * cgs::kpc))};
  const ham_float chi{par->breg_lsa.chi0 * std::tanh(pos[2] / (1. * cgs::kpc))};
  return Hamvec<3, ham_float>{
      par->breg_lsa.b0 * (std::sin(psi) * std::cos(chi) * cos_p -
                          std::cos(psi) * std::cos(chi) * sin_p),
      par->breg_lsa.b0 * (std::sin(psi) * std::cos(chi) * sin_p +
                          std::cos(psi) * std::cos(chi) * cos_p),
      par->breg_lsa.b0 * std::sin(chi)};
}
