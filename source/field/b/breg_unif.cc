#include <bfield.h>
#include <cmath>
#include <hamvec.h>
#include <param.h>

hamvec<3, double> Breg_unif::write_field(const hamvec<3, double> &,
                                         const Param *par) const {
  return hamvec<3, double>{par->breg_unif.bp * std::cos(par->breg_unif.l0),
                           par->breg_unif.bp * std::sin(par->breg_unif.l0),
                           par->breg_unif.bv};
}
