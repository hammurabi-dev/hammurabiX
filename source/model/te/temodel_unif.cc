#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>

ham_float TEmodel_unif::read_model(const Hamvec<3, ham_float> &pos,
                                   const Param *par) const {
  if ((pos - par->observer).length() > par->temodel_unif.r0) {
    return 0.;
  }
  return par->temodel_unif.n0;
}
