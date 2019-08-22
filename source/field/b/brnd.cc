#include <cassert>
#include <cmath>
#include <stdexcept>

#include <bfield.h>
#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

hamvec<3, double> Brnd::read_field(const hamvec<3, double> &pos,
                                   const Param *par,
                                   const Grid_brnd *grid) const {
  if (par->grid_brnd.read_permission or par->grid_brnd.build_permission) {
    return read_grid(pos, par, grid);
  }
  // if no specific random field model is called
  // base class will return zero vector
  else {
    return hamvec<3, double>{0., 0., 0.};
  }
}

hamvec<3, double> Brnd::read_grid(const hamvec<3, double> &pos,
                                  const Param *par,
                                  const Grid_brnd *grid) const {
  double tmp{(par->grid_brnd.nx - 1) * (pos[0] - par->grid_brnd.x_min) /
             (par->grid_brnd.x_max - par->grid_brnd.x_min)};
  if (tmp <= 0 or tmp >= par->grid_brnd.nx - 1) {
    return hamvec<3, double>{0., 0., 0.};
  }
  decltype(par->grid_brnd.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd{tmp - xl};
  tmp = (par->grid_brnd.ny - 1) * (pos[1] - par->grid_brnd.y_min) /
        (par->grid_brnd.y_max - par->grid_brnd.y_min);
  if (tmp <= 0 or tmp >= par->grid_brnd.ny - 1) {
    return hamvec<3, double>{0., 0., 0.};
  }
  decltype(par->grid_brnd.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_brnd.nz - 1) * (pos[2] - par->grid_brnd.z_min) /
        (par->grid_brnd.z_max - par->grid_brnd.z_min);
  if (tmp <= 0 or tmp >= par->grid_brnd.nz - 1) {
    return hamvec<3, double>{0., 0., 0.};
  }
  decltype(par->grid_brnd.nx) zl{(std::size_t)std::floor(tmp)};
  const double zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  // linear interpolation
  std::size_t idx1{toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                                    par->grid_brnd.nz, xl, yl, zl)};
  std::size_t idx2{toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                                    par->grid_brnd.nz, xl, yl, zl + 1)};
  const hamvec<3, double> i1{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl, yl + 1, zl + 1);
  const hamvec<3, double> i2{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl, zl + 1);
  const hamvec<3, double> j1{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                          par->grid_brnd.nz, xl + 1, yl + 1, zl + 1);
  const hamvec<3, double> j2{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  // interpolate along y direction, two interpolated vectors
  const hamvec<3, double> w1{i1 * (1. - yd) + i2 * yd};
  const hamvec<3, double> w2{j1 * (1. - yd) + j2 * yd};
  // interpolate along x direction
  return w1 * (1. - xd) + w2 * xd;
}

void Brnd::write_grid(const Param *, const Breg *, const Grid_breg *,
                      Grid_brnd *) const {
  throw std::runtime_error("wrong inheritance");
}
