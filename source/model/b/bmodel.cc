#include <bfield.h>
#include <cassert>
#include <cmath>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

Hamvec<3, ham_float> Bmodel::read_model(const Hamvec<3, ham_float> &,
                                        const Param *) const {
  return Hamvec<3, ham_float>(0., 0., 0.);
}

Hamvec<3, ham_float> Bmodel::read_grid(const ham_uint &idx,
                                       const Grid_b *grid) const {
  return Hamvec<3, ham_float>(grid->bx[idx], grid->by[idx], grid->bz[idx]);
}

Hamvec<3, ham_float> Bmodel::read_grid(const Hamvec<3, ham_float> &pos,
                                       const Param *par,
                                       const Grid_b *grid) const {
  // x position in bin-space
  ham_float tmp{(par->grid_b.nx - 1) * (pos[0] - par->grid_b.x_min) /
                (par->grid_b.x_max - par->grid_b.x_min)};
  // outside box, surface excluded
  if (tmp <= 0 or tmp >= par->grid_b.nx - 1) {
    return Hamvec<3, ham_float>(0., 0., 0.);
  }
  // current bin id
  ham_uint xl{(ham_uint)floor(tmp)};
  // bin-space distance to lower surface at current bin
  // must be within (0,1)
  const ham_float xd{tmp - xl};
  tmp = (par->grid_b.ny - 1) * (pos[1] - par->grid_b.y_min) /
        (par->grid_b.y_max - par->grid_b.y_min);
  if (tmp <= 0 or tmp >= par->grid_b.ny - 1) {
    return Hamvec<3, ham_float>(0., 0., 0.);
  }
  decltype(par->grid_b.nx) yl{(ham_uint)floor(tmp)};
  const ham_float yd{tmp - yl};
  tmp = (par->grid_b.nz - 1) * (pos[2] - par->grid_b.z_min) /
        (par->grid_b.z_max - par->grid_b.z_min);
  if (tmp <= 0 or tmp >= par->grid_b.nz - 1) {
    return Hamvec<3, ham_float>(0., 0., 0.);
  }
  ham_uint zl{(ham_uint)floor(tmp)};
  const ham_float zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  // linear interpolation
  // interpolate along z direction, there are four interpolated vectors
  ham_uint idx1{toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz,
                                 xl, yl, zl)};
  ham_uint idx2{toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz,
                                 xl, yl, zl + 1)};
  const Hamvec<3, ham_float> i1{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz, xl,
                          yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz, xl,
                          yl + 1, zl + 1);
  const Hamvec<3, ham_float> i2{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz,
                          xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz,
                          xl + 1, yl, zl + 1);
  const Hamvec<3, ham_float> j1{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz,
                          xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz,
                          xl + 1, yl + 1, zl + 1);
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

void Bmodel::write_grid(const Param *par, Grid_b *grid) const {
  assert(par->grid_b.build_permission or par->grid_b.write_permission or
         par->grid_b.read_permission);
  Hamvec<3, ham_float> gc_pos, tmp_vec;
  ham_float lx{par->grid_b.x_max - par->grid_b.x_min};
  ham_float ly{par->grid_b.y_max - par->grid_b.y_min};
  ham_float lz{par->grid_b.z_max - par->grid_b.z_min};
  for (ham_uint i = 0; i != par->grid_b.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_b.nx - 1) + par->grid_b.x_min;
    const ham_uint idx_lv1{i * par->grid_b.ny * par->grid_b.nz};
    for (ham_uint j = 0; j != par->grid_b.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_b.ny - 1) + par->grid_b.y_min;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_b.nz};
      for (ham_uint k = 0; k != par->grid_b.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_b.nz - 1) + par->grid_b.z_min;
        const ham_uint idx{idx_lv2 + k};
        tmp_vec = read_model(gc_pos, par);
        grid->bx[idx] = tmp_vec[0];
        grid->by[idx] = tmp_vec[1];
        grid->bz[idx] = tmp_vec[2];
      }
    }
  }
}
