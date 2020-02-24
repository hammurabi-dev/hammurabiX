#include <bfield.h>
#include <cmath>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>

Hamvec<3, ham_float> Breg_unif::write_field(const Hamvec<3, ham_float> &,
                                            const Param *par) const {
  return Hamvec<3, ham_float>{par->breg_unif.bp * std::cos(par->breg_unif.l0),
                              par->breg_unif.bp * std::sin(par->breg_unif.l0),
                              par->breg_unif.bv};
}
