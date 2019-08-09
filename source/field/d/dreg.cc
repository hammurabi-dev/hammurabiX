#include <cassert>
#include <cmath>

#include <dfield.h>
#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

double Dreg::read_field(const hamvec<3, double> &pos, const Param *par,
                        const Grid_dreg *grid) const {
  if (par->grid_dreg.read_permission) {
    return read_grid(pos, par, grid);
  } else if (par->grid_dreg.build_permission) {
    return write_field(pos, par);
  } else {
    return 0.;
  }
}

double Dreg::write_field(const hamvec<3, double> &, const Param *) const {
  return 0.;
}

double Dreg::read_grid(const hamvec<3, double> &pos, const Param *par,
                       const Grid_dreg *grid) const {
  double tmp{(par->grid_dreg.nx - 1) * (pos[0] - par->grid_dreg.x_min) /
             (par->grid_dreg.x_max - par->grid_dreg.x_min)};
  if (tmp <= 0 or tmp >= par->grid_dreg.nx - 1) {
    return 0.;
  }
  decltype(par->grid_dreg.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd{tmp - xl};
  tmp = (par->grid_dreg.ny - 1) * (pos[1] - par->grid_dreg.y_min) /
        (par->grid_dreg.y_max - par->grid_dreg.y_min);
  if (tmp <= 0 or tmp >= par->grid_dreg.ny - 1) {
    return 0.;
  }
  decltype(par->grid_dreg.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_dreg.nz - 1) * (pos[2] - par->grid_dreg.z_min) /
        (par->grid_dreg.z_max - par->grid_dreg.z_min);
  if (tmp <= 0 or tmp >= par->grid_dreg.nz - 1) {
    return 0.;
  }
  decltype(par->grid_dreg.nx) zl{(std::size_t)std::floor(tmp)};
  const double zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  std::size_t idx1{toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                                    par->grid_dreg.nz, xl, yl, zl)};
  std::size_t idx2{toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                                    par->grid_dreg.nz, xl, yl, zl + 1)};
  const double i1{grid->d[idx1] * (1. - zd) + grid->d[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                          par->grid_dreg.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                          par->grid_dreg.nz, xl, yl + 1, zl + 1);
  const double i2{grid->d[idx1] * (1 - zd) + grid->d[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                          par->grid_dreg.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                          par->grid_dreg.nz, xl + 1, yl, zl + 1);
  const double j1{grid->d[idx1] * (1 - zd) + grid->d[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                          par->grid_dreg.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_dreg.nx, par->grid_dreg.ny,
                          par->grid_dreg.nz, xl + 1, yl + 1, zl + 1);
  const double j2{grid->d[idx1] * (1 - zd) + grid->d[idx2] * zd};
  const double w1{i1 * (1 - yd) + i2 * yd};
  const double w2{j1 * (1 - yd) + j2 * yd};
  return w1 * (1 - xd) + w2 * xd;
}

void Dreg::write_grid(const Param *par, Grid_dreg *grid) const {
  assert(par->grid_dreg.write_permission);
  hamvec<3, double> gc_pos;
  double lx{par->grid_dreg.x_max - par->grid_dreg.x_min};
  double ly{par->grid_dreg.y_max - par->grid_dreg.y_min};
  double lz{par->grid_dreg.z_max - par->grid_dreg.z_min};
  for (decltype(par->grid_dreg.nx) i = 0; i != par->grid_dreg.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_dreg.nx - 1) + par->grid_dreg.x_min;
    const std::size_t idx_lv1{i * par->grid_dreg.ny * par->grid_dreg.nz};
    for (decltype(par->grid_dreg.ny) j = 0; j != par->grid_dreg.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_dreg.ny - 1) + par->grid_dreg.y_min;
      const std::size_t idx_lv2{idx_lv1 + j * par->grid_dreg.nz};
      for (decltype(par->grid_dreg.nz) k = 0; k != par->grid_dreg.nz; ++k) {
        const std::size_t idx{idx_lv2 + k};
        gc_pos[2] = lz * k / (par->grid_dreg.nz - 1) + par->grid_dreg.z_min;
        grid->d[idx] = write_field(gc_pos, par);
      }
    }
  }
}
