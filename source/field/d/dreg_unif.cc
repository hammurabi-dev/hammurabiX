#include <dfield.h>
#include <hamvec.h>
#include <param.h>

double Dreg_unif::write_field(const hamvec<3, double> &pos,
                              const Param *par) const {
  if ((pos - par->observer).length() > par->dreg_unif.r0) {
    return 0.;
  }
  return par->dreg_unif.n0;
}
