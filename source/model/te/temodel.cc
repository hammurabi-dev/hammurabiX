#include <cassert>
#include <cmath>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>
#include <toolkit.h>

ham_float TEmodel::read_model(const Hamvec<3, ham_float> &,
                              const Param *) const {
  return 0.;
}

ham_float TEmodel::read_grid(const ham_uint &idx, const Grid_te *grid) const {
  return grid->te[idx];
}

ham_float TEmodel::read_grid(const Hamvec<3, ham_float> &pos, const Param *par,
                             const Grid_te *grid) const {
  ham_float tmp{(par->grid_te.nx - 1) * (pos[0] - par->grid_te.x_min) /
                (par->grid_te.x_max - par->grid_te.x_min)};
  if (tmp <= 0 or tmp >= par->grid_te.nx - 1) {
    return 0.;
  }
  decltype(par->grid_te.nx) xl{(ham_uint)std::floor(tmp)};
  const ham_float xd{tmp - xl};
  tmp = (par->grid_te.ny - 1) * (pos[1] - par->grid_te.y_min) /
        (par->grid_te.y_max - par->grid_te.y_min);
  if (tmp <= 0 or tmp >= par->grid_te.ny - 1) {
    return 0.;
  }
  decltype(par->grid_te.nx) yl{(ham_uint)std::floor(tmp)};
  const ham_float yd{tmp - yl};
  tmp = (par->grid_te.nz - 1) * (pos[2] - par->grid_te.z_min) /
        (par->grid_te.z_max - par->grid_te.z_min);
  if (tmp <= 0 or tmp >= par->grid_te.nz - 1) {
    return 0.;
  }
  decltype(par->grid_te.nx) zl{(ham_uint)std::floor(tmp)};
  const ham_float zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  ham_uint idx1{toolkit::index3d(par->grid_te.nx, par->grid_te.ny,
                                 par->grid_te.nz, xl, yl, zl)};
  ham_uint idx2{toolkit::index3d(par->grid_te.nx, par->grid_te.ny,
                                 par->grid_te.nz, xl, yl, zl + 1)};
  const ham_float i1{grid->te[idx1] * (1. - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_te.nx, par->grid_te.ny, par->grid_te.nz, xl,
                          yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_te.nx, par->grid_te.ny, par->grid_te.nz, xl,
                          yl + 1, zl + 1);
  const ham_float i2{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_te.nx, par->grid_te.ny, par->grid_te.nz,
                          xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_te.nx, par->grid_te.ny, par->grid_te.nz,
                          xl + 1, yl, zl + 1);
  const ham_float j1{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_te.nx, par->grid_te.ny, par->grid_te.nz,
                          xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_te.nx, par->grid_te.ny, par->grid_te.nz,
                          xl + 1, yl + 1, zl + 1);
  const ham_float j2{grid->te[idx1] * (1 - zd) + grid->te[idx2] * zd};
  const ham_float w1{i1 * (1 - yd) + i2 * yd};
  const ham_float w2{j1 * (1 - yd) + j2 * yd};
  return w1 * (1 - xd) + w2 * xd;
}

void TEmodel::write_grid(const Param *par, Grid_te *grid) const {
  assert(par->grid_te.build_permission or par->grid_te.write_permission or
         par->grid_te.read_permission);
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_te.x_max - par->grid_te.x_min};
  ham_float ly{par->grid_te.y_max - par->grid_te.y_min};
  ham_float lz{par->grid_te.z_max - par->grid_te.z_min};
  for (decltype(par->grid_te.nx) i = 0; i != par->grid_te.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_te.nx - 1) + par->grid_te.x_min;
    const ham_uint idx_lv1{i * par->grid_te.ny * par->grid_te.nz};
    for (decltype(par->grid_te.ny) j = 0; j != par->grid_te.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_te.ny - 1) + par->grid_te.y_min;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_te.nz};
      for (decltype(par->grid_te.nz) k = 0; k != par->grid_te.nz; ++k) {
        const ham_uint idx{idx_lv2 + k};
        gc_pos[2] = lz * k / (par->grid_te.nz - 1) + par->grid_te.z_min;
        grid->te[idx] = read_model(gc_pos, par);
      }
    }
  }
}
