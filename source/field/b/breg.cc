#include <cassert>
#include <cmath>

#include <bfield.h>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

Hamvec<3, ham_float> Breg::read_field(const Hamvec<3, ham_float> &pos,
                                      const Param *par,
                                      const Grid_breg *grid) const {
  if (par->grid_breg.read_permission) {
    return read_grid(pos, par, grid);
  } else if (par->grid_breg.build_permission) {
    return write_field(pos, par);
  } else {
    return Hamvec<3, ham_float>{0., 0., 0.};
  }
}

Hamvec<3, ham_float> Breg::write_field(const Hamvec<3, ham_float> &,
                                       const Param *) const {
  return Hamvec<3, ham_float>{0., 0., 0.};
}

Hamvec<3, ham_float> Breg::read_grid(const Hamvec<3, ham_float> &pos,
                                     const Param *par,
                                     const Grid_breg *grid) const {
  // x position in bin-space
  ham_float tmp{(par->grid_breg.nx - 1) * (pos[0] - par->grid_breg.x_min) /
                (par->grid_breg.x_max - par->grid_breg.x_min)};
  // outside box, surface excluded
  if (tmp <= 0 or tmp >= par->grid_breg.nx - 1) {
    return Hamvec<3, ham_float>(0., 0., 0.);
  }
  // current bin id
  decltype(par->grid_breg.nx) xl{(ham_uint)floor(tmp)};
  // bin-space distance to lower surface at current bin
  // must be within (0,1)
  const ham_float xd{tmp - xl};
  tmp = (par->grid_breg.ny - 1) * (pos[1] - par->grid_breg.y_min) /
        (par->grid_breg.y_max - par->grid_breg.y_min);
  if (tmp <= 0 or tmp >= par->grid_breg.ny - 1) {
    return Hamvec<3, ham_float>(0., 0., 0.);
  }
  decltype(par->grid_breg.nx) yl{(ham_uint)floor(tmp)};
  const ham_float yd{tmp - yl};
  tmp = (par->grid_breg.nz - 1) * (pos[2] - par->grid_breg.z_min) /
        (par->grid_breg.z_max - par->grid_breg.z_min);
  if (tmp <= 0 or tmp >= par->grid_breg.nz - 1) {
    return Hamvec<3, ham_float>(0., 0., 0.);
  }
  decltype(par->grid_breg.nx) zl{(ham_uint)floor(tmp)};
  const ham_float zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  // linear interpolation
  // interpolate along z direction, there are four interpolated vectors
  ham_uint idx1{toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                                 par->grid_breg.nz, xl, yl, zl)};
  ham_uint idx2{toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                                 par->grid_breg.nz, xl, yl, zl + 1)};
  const Hamvec<3, ham_float> i1{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl, yl + 1, zl + 1);
  const Hamvec<3, ham_float> i2{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl, zl + 1);
  const Hamvec<3, ham_float> j1{
      grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
      grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
      grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl + 1, zl + 1);
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

void Breg::write_grid(const Param *par, Grid_breg *grid) const {
  assert(par->grid_breg.write_permission);
  Hamvec<3, ham_float> gc_pos, tmp_vec;
  ham_float lx{par->grid_breg.x_max - par->grid_breg.x_min};
  ham_float ly{par->grid_breg.y_max - par->grid_breg.y_min};
  ham_float lz{par->grid_breg.z_max - par->grid_breg.z_min};
  for (decltype(par->grid_breg.nx) i = 0; i != par->grid_breg.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_breg.nx - 1) + par->grid_breg.x_min;
    const ham_uint idx_lv1{i * par->grid_breg.ny * par->grid_breg.nz};
    for (decltype(par->grid_breg.ny) j = 0; j != par->grid_breg.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_breg.ny - 1) + par->grid_breg.y_min;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_breg.nz};
      for (decltype(par->grid_breg.nz) k = 0; k != par->grid_breg.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_breg.nz - 1) + par->grid_breg.z_min;
        const ham_uint idx{idx_lv2 + k};
        tmp_vec = write_field(gc_pos, par);
        grid->bx[idx] = tmp_vec[0];
        grid->by[idx] = tmp_vec[1];
        grid->bz[idx] = tmp_vec[2];
      }
    }
  }
}
