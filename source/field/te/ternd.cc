#include <cassert>
#include <cmath>
#include <stdexcept>

#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>
#include <toolkit.h>

double TErnd::read_field(const hamvec<3, double> &pos, const Param *par,
                         const Grid_ternd *grid) const {
  if (par->grid_ternd.read_permission or par->grid_ternd.build_permission) {
    return read_grid(pos, par, grid);
  } else {
    return 0.;
  }
}

double TErnd::read_grid(const hamvec<3, double> &pos, const Param *par,
                        const Grid_ternd *grid) const {
  double tmp{(par->grid_ternd.nx - 1) * (pos[0] - par->grid_ternd.x_min) /
             (par->grid_ternd.x_max - par->grid_ternd.x_min)};
  if (tmp <= 0 or tmp >= par->grid_ternd.nx - 1) {
    return 0.;
  }
  decltype(par->grid_ternd.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd{tmp - xl};
  tmp = (par->grid_ternd.ny - 1) * (pos[1] - par->grid_ternd.y_min) /
        (par->grid_ternd.y_max - par->grid_ternd.y_min);
  if (tmp <= 0 or tmp >= par->grid_ternd.ny - 1) {
    return 0.;
  }
  decltype(par->grid_ternd.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_ternd.nz - 1) * (pos[2] - par->grid_ternd.z_min) /
        (par->grid_ternd.z_max - par->grid_ternd.z_min);
  if (tmp <= 0 or tmp >= par->grid_ternd.nz - 1) {
    return 0.;
  }
  decltype(par->grid_ternd.nx) zl{(std::size_t)std::floor(tmp)};
  const double zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  // linear interpolation
  std::size_t idx1{toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                                    par->grid_ternd.nz, xl, yl, zl)};
  std::size_t idx2{toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                                    par->grid_ternd.nz, xl, yl, zl + 1)};
  const double i1{grid->te[idx1] * (1. - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                          par->grid_ternd.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                          par->grid_ternd.nz, xl, yl + 1, zl + 1);
  const double i2{grid->te[idx1] * (1. - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                          par->grid_ternd.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                          par->grid_ternd.nz, xl + 1, yl, zl + 1);
  const double j1{grid->te[idx1] * (1. - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                          par->grid_ternd.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                          par->grid_ternd.nz, xl + 1, yl + 1, zl + 1);
  const double j2{grid->te[idx1] * (1. - zd) + grid->te[idx2] * zd};
  // interpolate along y direction, two interpolated vectors
  const double w1{i1 * (1. - yd) + i2 * yd};
  const double w2{j1 * (1. - yd) + j2 * yd};
  // interpolate along x direction
  return w1 * (1. - xd) + w2 * xd;
}

void TErnd::write_grid(const Param *, const TEreg *, const Grid_tereg *,
                       Grid_ternd *) const {
  throw std::runtime_error("wrong inheritance");
}
