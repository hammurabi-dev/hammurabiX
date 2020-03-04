#include <cassert>
#include <cmath>

#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>
#include <toolkit.h>

ham_float TEreg::read_field(const Hamvec<3, ham_float> &pos, const Param *par,
                            const Grid_tereg *grid) const {
  if (par->grid_tereg.read_permission) {
    return read_grid(pos, par, grid);
  } else if (par->grid_tereg.build_permission) {
    return write_field(pos, par);
  } else {
    return 0.;
  }
}

ham_float TEreg::write_field(const Hamvec<3, ham_float> &,
                             const Param *) const {
  return 0.;
}

ham_float TEreg::read_grid(const Hamvec<3, ham_float> &pos, const Param *par,
                           const Grid_tereg *grid) const {
  ham_float tmp{(par->grid_tereg.nx - 1) * (pos[0] - par->grid_tereg.x_min) /
                (par->grid_tereg.x_max - par->grid_tereg.x_min)};
  if (tmp <= 0 or tmp >= par->grid_tereg.nx - 1) {
    return 0.;
  }
  decltype(par->grid_tereg.nx) xl{(ham_uint)std::floor(tmp)};
  const ham_float xd{tmp - xl};
  tmp = (par->grid_tereg.ny - 1) * (pos[1] - par->grid_tereg.y_min) /
        (par->grid_tereg.y_max - par->grid_tereg.y_min);
  if (tmp <= 0 or tmp >= par->grid_tereg.ny - 1) {
    return 0.;
  }
  decltype(par->grid_tereg.nx) yl{(ham_uint)std::floor(tmp)};
  const ham_float yd{tmp - yl};
  tmp = (par->grid_tereg.nz - 1) * (pos[2] - par->grid_tereg.z_min) /
        (par->grid_tereg.z_max - par->grid_tereg.z_min);
  if (tmp <= 0 or tmp >= par->grid_tereg.nz - 1) {
    return 0.;
  }
  decltype(par->grid_tereg.nx) zl{(ham_uint)std::floor(tmp)};
  const ham_float zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  ham_uint idx1{toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                                 par->grid_tereg.nz, xl, yl, zl)};
  ham_uint idx2{toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                                 par->grid_tereg.nz, xl, yl, zl + 1)};
  const ham_float i1{grid->te[idx1] * (1. - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl, yl + 1, zl + 1);
  const ham_float i2{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl, zl + 1);
  const ham_float j1{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                          par->grid_tereg.nz, xl + 1, yl + 1, zl + 1);
  const ham_float j2{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  const ham_float w1{i1 * (1 - yd) + i2 * yd};
  const ham_float w2{j1 * (1 - yd) + j2 * yd};
  return w1 * (1 - xd) + w2 * xd;
}

void TEreg::write_grid(const Param *par, Grid_tereg *grid) const {
  assert(par->grid_tereg.write_permission);
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_tereg.x_max - par->grid_tereg.x_min};
  ham_float ly{par->grid_tereg.y_max - par->grid_tereg.y_min};
  ham_float lz{par->grid_tereg.z_max - par->grid_tereg.z_min};
  for (decltype(par->grid_tereg.nx) i = 0; i != par->grid_tereg.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_tereg.nx - 1) + par->grid_tereg.x_min;
    const ham_uint idx_lv1{i * par->grid_tereg.ny * par->grid_tereg.nz};
    for (decltype(par->grid_tereg.ny) j = 0; j != par->grid_tereg.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_tereg.ny - 1) + par->grid_tereg.y_min;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_tereg.nz};
      for (decltype(par->grid_tereg.nz) k = 0; k != par->grid_tereg.nz; ++k) {
        const ham_uint idx{idx_lv2 + k};
        gc_pos[2] = lz * k / (par->grid_tereg.nz - 1) + par->grid_tereg.z_min;
        grid->te[idx] = write_field(gc_pos, par);
      }
    }
  }
}
