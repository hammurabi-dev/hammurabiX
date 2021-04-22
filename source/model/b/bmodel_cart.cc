#include <bfield.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>

Hamvec<3, ham_float> Bmodel_cart::write_field(const Hamvec<3, ham_float> &,
                                              const Param *par) const {
  return Hamvec<3, ham_float>{par->bmodel_cart.bx, par->bmodel_cart.by,
                              par->bmodel_cart.bz};
}
