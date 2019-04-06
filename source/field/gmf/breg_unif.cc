#include <cmath>

#include <hvec.h>

#include <breg.h>
#include <param.h>

hvec<3, double> Breg_unif::gmf(const hvec<3, double> &,
                               const Param *par) const {
  return hvec<3, double>{par->breg_unif.b0 * std::cos(par->breg_unif.l0),
                         par->breg_unif.b0 * std::sin(par->breg_unif.l0),
                         par->breg_unif.b0 * par->breg_unif.r};
}

// END
