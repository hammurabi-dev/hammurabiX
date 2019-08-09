#include <cassert>
#include <cmath>
#include <stdexcept>

#include <dfield.h>
#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

double Drnd::read_field(const hamvec<3, double> &pos, const Param *par,
                        const Grid_drnd *grid) const {
  if (par->grid_drnd.read_permission or par->grid_drnd.build_permission) {
    return read_grid(pos, par, grid);
  } else {
    return 0.;
  }
}

double Drnd::read_grid(const hamvec<3, double> &pos, const Param *par,
                       const Grid_drnd *grid) const {
  double tmp{(par->grid_drnd.nx - 1) * (pos[0] - par->grid_drnd.x_min) /
             (par->grid_drnd.x_max - par->grid_drnd.x_min)};
  if (tmp <= 0 or tmp >= par->grid_drnd.nx - 1) {
    return 0.;
  }
  decltype(par->grid_drnd.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd{tmp - xl};
  tmp = (par->grid_drnd.ny - 1) * (pos[1] - par->grid_drnd.y_min) /
        (par->grid_drnd.y_max - par->grid_drnd.y_min);
  if (tmp <= 0 or tmp >= par->grid_drnd.ny - 1) {
    return 0.;
  }
  decltype(par->grid_drnd.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_drnd.nz - 1) * (pos[2] - par->grid_drnd.z_min) /
        (par->grid_drnd.z_max - par->grid_drnd.z_min);
  if (tmp <= 0 or tmp >= par->grid_drnd.nz - 1) {
    return 0.;
  }
  decltype(par->grid_drnd.nx) zl{(std::size_t)std::floor(tmp)};
  const double zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  // trilinear interpolation
  std::size_t idx1{toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                                    par->grid_drnd.nz, xl, yl, zl)};
  std::size_t idx2{toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                                    par->grid_drnd.nz, xl, yl, zl + 1)};
  const double i1{grid->d[idx1] * (1. - zd) + grid->d[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                          par->grid_drnd.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                          par->grid_drnd.nz, xl, yl + 1, zl + 1);
  const double i2{grid->d[idx1] * (1. - zd) + grid->d[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                          par->grid_drnd.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                          par->grid_drnd.nz, xl + 1, yl, zl + 1);
  const double j1{grid->d[idx1] * (1. - zd) + grid->d[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                          par->grid_drnd.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_drnd.nx, par->grid_drnd.ny,
                          par->grid_drnd.nz, xl + 1, yl + 1, zl + 1);
  const double j2{grid->d[idx1] * (1. - zd) + grid->d[idx2] * zd};
  // interpolate along y direction, two interpolated vectors
  const double w1{i1 * (1. - yd) + i2 * yd};
  const double w2{j1 * (1. - yd) + j2 * yd};
  // interpolate along x direction
  return w1 * (1. - xd) + w2 * xd;
}

void Drnd::write_grid(const Param *, Grid_drnd *) const {
  throw std::runtime_error("wrong inheritance");
}
