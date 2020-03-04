#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>

ham_float TEreg_unif::write_field(const Hamvec<3, ham_float> &pos,
                                  const Param *par) const {
  if ((pos - par->observer).length() > par->tereg_unif.r0) {
    return 0.;
  }
  return par->tereg_unif.n0;
}
