#include <cassert>
#include <cmath>
#include <stdexcept>

#include <bfield.h>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

Hamvec<3, ham_float> Brnd::read_field(const Hamvec<3, ham_float> &pos,
                                      const Param *par,
                                      const Grid_brnd *grid) const {
  if (par->grid_brnd.read_permission or par->grid_brnd.build_permission) {
    return read_grid(pos, par, grid);
  }
  // if no specific random field model is called
  // base class will return zero vector
  else {
    return Hamvec<3, ham_float>{0., 0., 0.};
  }
}

Hamvec<3, ham_float> Brnd::read_grid(const Hamvec<3, ham_float> &pos,
                                     const Param *par,
                                     const Grid_brnd *grid) const {
  ham_float tmp{(par->grid_brnd.nx - 1) * (pos[0] - par->grid_brnd.x_min) /
                (par->grid_brnd.x_max - par->grid_brnd.x_min)};
  if (tmp <= 0 or tmp >= par->grid_brnd.nx - 1) {
    return Hamvec<3, ham_float>{0., 0., 0.};
  }
  decltype(par->grid_brnd.nx) xl{(ham_uint)std::floor(tmp)};
  const ham_float xd{tmp - xl};
  tmp = (par->grid_brnd.ny - 1) * (pos[1] - par->grid_brnd.y_min) /
        (par->grid_brnd.y_max - par->grid_brnd.y_min);
  if (tmp <= 0 or tmp >= par->grid_brnd.ny - 1) {
    return Hamvec<3, ham_float>{0., 0., 0.};
  }
  decltype(par->grid_brnd.nx) yl{(ham_uint)std::floor(tmp)};
  const ham_float yd{tmp - yl};
  tmp = (par->grid_brnd.nz - 1) * (pos[2] - par->grid_brnd.z_min) /
        (par->grid_brnd.z_max - par->grid_brnd.z_min);
  if (tmp <= 0 or tmp >= par->grid_brnd.nz - 1) {
    return Hamvec<3, ham_float>{0., 0., 0.};
  }
  decltype(par->grid_brnd.nx) zl{(ham_uint)std::floor(tmp)};
  const ham_float zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  // linear interpolation
  ham_uint idx1{toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                                 par->grid_brnd.nz, xl, yl, zl)};
  ham_uint idx2{toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                                 par->grid_brnd.nz, xl, yl, zl + 1)};
  const Hamvec<3, ham_float> i1{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl, yl + 1, zl + 1);
  const Hamvec<3, ham_float> i2{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl, zl + 1);
  const Hamvec<3, ham_float> j1{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl + 1, zl + 1);
  const Hamvec<3, ham_float> j2{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  // interpolate along y direction, two interpolated vectors
  const Hamvec<3, ham_float> w1{i1 * (1. - yd) + i2 * yd};
  const Hamvec<3, ham_float> w2{j1 * (1. - yd) + j2 * yd};
  // interpolate along x direction
  return w1 * (1. - xd) + w2 * xd;
}

void Brnd::write_grid(const Param *, const Breg *, const Grid_breg *,
                      Grid_brnd *) const {
  throw std::runtime_error("wrong inheritance");
}
