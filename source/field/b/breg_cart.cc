#include <bfield.h>
#include <cmath>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>

Hamvec<3, ham_float> Breg_cart::write_field(const Hamvec<3, ham_float> &,
                                            const Param *par) const {
  return Hamvec<3, ham_float>{par->breg_cart.bx,
                              par->breg_cart.by,
                              par->breg_cart.bz};
}
