#include <cassert>
#include <cmath>

#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>
#include <toolkit.h>

double TEreg::read_field(const hamvec<3, double> &pos, const Param *par,
                         const Grid_tereg *grid) const {
  if (par->grid_tereg.read_permission) {
    return read_grid(pos, par, grid);
  } else if (par->grid_tereg.build_permission) {
    return write_field(pos, par);
  } else {
    return 0.;
  }
}

double TEreg::write_field(const hamvec<3, double> &, const Param *) const {
  return 0.;
}

double TEreg::read_grid(const hamvec<3, double> &pos, const Param *par,
                        const Grid_tereg *grid) const {
  double tmp{(par->grid_tereg.nx - 1) * (pos[0] - par->grid_tereg.x_min) /
             (par->grid_tereg.x_max - par->grid_tereg.x_min)};
  if (tmp <= 0 or tmp >= par->grid_tereg.nx - 1) {
    return 0.;
  }
  decltype(par->grid_tereg.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd{tmp - xl};
  tmp = (par->grid_tereg.ny - 1) * (pos[1] - par->grid_tereg.y_min) /
        (par->grid_tereg.y_max - par->grid_tereg.y_min);
  if (tmp <= 0 or tmp >= par->grid_tereg.ny - 1) {
    return 0.;
  }
  decltype(par->grid_tereg.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_tereg.nz - 1) * (pos[2] - par->grid_tereg.z_min) /
        (par->grid_tereg.z_max - par->grid_tereg.z_min);
  if (tmp <= 0 or tmp >= par->grid_tereg.nz - 1) {
    return 0.;
  }
  decltype(par->grid_tereg.nx) zl{(std::size_t)std::floor(tmp)};
  const double zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  std::size_t idx1{toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                                    par->grid_tereg.nz, xl, yl, zl)};
  std::size_t idx2{toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                                    par->grid_tereg.nz, xl, yl, zl + 1)};
  const double i1{grid->te[idx1] * (1. - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl, yl + 1, zl + 1);
  const double i2{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl, zl + 1);
  const double j1{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl + 1, zl + 1);
  const double j2{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  const double w1{i1 * (1 - yd) + i2 * yd};
  const double w2{j1 * (1 - yd) + j2 * yd};
  return w1 * (1 - xd) + w2 * xd;
}

void TEreg::write_grid(const Param *par, Grid_tereg *grid) const {
  assert(par->grid_tereg.write_permission);
  hamvec<3, double> gc_pos;
  double lx{par->grid_tereg.x_max - par->grid_tereg.x_min};
  double ly{par->grid_tereg.y_max - par->grid_tereg.y_min};
  double lz{par->grid_tereg.z_max - par->grid_tereg.z_min};
  for (decltype(par->grid_tereg.nx) i = 0; i != par->grid_tereg.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_tereg.nx - 1) + par->grid_tereg.x_min;
    const std::size_t idx_lv1{i * par->grid_tereg.ny * par->grid_tereg.nz};
    for (decltype(par->grid_tereg.ny) j = 0; j != par->grid_tereg.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_tereg.ny - 1) + par->grid_tereg.y_min;
      const std::size_t idx_lv2{idx_lv1 + j * par->grid_tereg.nz};
      for (decltype(par->grid_tereg.nz) k = 0; k != par->grid_tereg.nz; ++k) {
        const std::size_t idx{idx_lv2 + k};
        gc_pos[2] = lz * k / (par->grid_tereg.nz - 1) + par->grid_tereg.z_min;
        grid->te[idx] = write_field(gc_pos, par);
      }
    }
  }
}
