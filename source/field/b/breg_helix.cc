#include <cmath>
#include <vector>

#include <bfield.h>
#include <grid.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>

// Create helical magnetic field from rmin to rmax with form (bx cos(phi), by
// sin(phi), bz phi)
Hamvec<3, ham_float> Breg_helix::write_field(const Hamvec<3, ham_float> &pos,
                                             const Param *par) const {

  const ham_float r{sqrt(pos[0] * pos[0] +
                         pos[1] * pos[1])}; // radius in cylindrical coordinates
  const ham_float phi{atan2(pos[1], pos[0]) +
                      M_PI / 2.0}; // azimuthal angle in cylindrical coordinates
  const ham_float rmin{par->breg_helix.r_min};
  const ham_float rmax{par->breg_helix.r_max};

  if ((r > rmin) && (r < rmax)) {
    return Hamvec<3, ham_float>{par->breg_helix.bx * std::cos(phi),
                                par->breg_helix.by * std::sin(phi),
                                par->breg_helix.bz};
  } else {
    return Hamvec<3, ham_float>{0.0, 0.0, 0.0};
  }
}
