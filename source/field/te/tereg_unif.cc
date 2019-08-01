#include <hamvec.h>
#include <param.h>
#include <tefield.h>

double TEreg_unif::write_field(const hamvec<3, double> &pos,
                               const Param *par) const {
  if ((pos - par->observer).length() > par->tereg_unif.r0) {
    return 0.;
  }
  return par->tereg_unif.n0;
}
