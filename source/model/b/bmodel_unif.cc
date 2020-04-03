#include <bfield.h>
#include <cmath>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>

Hamvec<3, ham_float> Bmodel_unif::read_model(const Hamvec<3, ham_float> &,
                                             const Param *par) const {
  return Hamvec<3, ham_float>{
      par->bmodel_unif.bp * std::cos(par->bmodel_unif.l0),
      par->bmodel_unif.bp * std::sin(par->bmodel_unif.l0), par->bmodel_unif.bv};
}
