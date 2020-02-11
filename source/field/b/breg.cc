#include <cassert>
#include <cmath>
#include <cstddef> // for std::size_t

#include <bfield.h>
#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

hamvec<3, double> Breg::read_field(const hamvec<3, double> &pos,
                                   const Param *par,
                                   const Grid_breg *grid) const {
  if (par->grid_breg.read_permission) {
    return read_grid(pos, par, grid);
  } else if (par->grid_breg.build_permission) {
    return write_field(pos, par);
  } else {
    return hamvec<3, double>{0., 0., 0.};
  }
}

hamvec<3, double> Breg::write_field(const hamvec<3, double> &,
                                    const Param *) const {
  return hamvec<3, double>{0., 0., 0.};
}

hamvec<3, double> Breg::read_grid(const hamvec<3, double> &pos,
                                  const Param *par,
                                  const Grid_breg *grid) const {
  // x position in bin-space
  double tmp{(par->grid_breg.nx - 1) * (pos[0] - par->grid_breg.x_min) /
             (par->grid_breg.x_max - par->grid_breg.x_min)};
  // outside box, surface excluded
  if (tmp <= 0 or tmp >= par->grid_breg.nx - 1) {
    return hamvec<3, double>(0., 0., 0.);
  }
  // current bin id
  decltype(par->grid_breg.nx) xl{(std::size_t)floor(tmp)};
  // bin-space distance to lower surface at current bin
  // must be within (0,1)
  const double xd{tmp - xl};
  tmp = (par->grid_breg.ny - 1) * (pos[1] - par->grid_breg.y_min) /
        (par->grid_breg.y_max - par->grid_breg.y_min);
  if (tmp <= 0 or tmp >= par->grid_breg.ny - 1) {
    return hamvec<3, double>(0., 0., 0.);
  }
  decltype(par->grid_breg.nx) yl{(std::size_t)floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_breg.nz - 1) * (pos[2] - par->grid_breg.z_min) /
        (par->grid_breg.z_max - par->grid_breg.z_min);
  if (tmp <= 0 or tmp >= par->grid_breg.nz - 1) {
    return hamvec<3, double>(0., 0., 0.);
  }
  decltype(par->grid_breg.nx) zl{(std::size_t)floor(tmp)};
  const double zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  // linear interpolation
  // interpolate along z direction, there are four interpolated vectors
  std::size_t idx1{toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                                    par->grid_breg.nz, xl, yl, zl)};
  std::size_t idx2{toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                                    par->grid_breg.nz, xl, yl, zl + 1)};
  const hamvec<3, double> i1{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl, yl + 1, zl + 1);
  const hamvec<3, double> i2{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl, zl);
  idx2 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl, zl + 1);
  const hamvec<3, double> j1{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  idx1 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl + 1, zl);
  idx2 = toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                          par->grid_breg.nz, xl + 1, yl + 1, zl + 1);
  const hamvec<3, double> j2{grid->bx[idx1] * (1. - zd) + grid->bx[idx2] * zd,
                             grid->by[idx1] * (1. - zd) + grid->by[idx2] * zd,
                             grid->bz[idx1] * (1. - zd) + grid->bz[idx2] * zd};
  // interpolate along y direction, two interpolated vectors
  const hamvec<3, double> w1{i1 * (1. - yd) + i2 * yd};
  const hamvec<3, double> w2{j1 * (1. - yd) + j2 * yd};
  // interpolate along x direction
  return w1 * (1. - xd) + w2 * xd;
}

void Breg::write_grid(const Param *par, Grid_breg *grid) const {
  assert(par->grid_breg.write_permission);
  hamvec<3, double> gc_pos, tmp_vec;
  double lx{par->grid_breg.x_max - par->grid_breg.x_min};
  double ly{par->grid_breg.y_max - par->grid_breg.y_min};
  double lz{par->grid_breg.z_max - par->grid_breg.z_min};
  for (decltype(par->grid_breg.nx) i = 0; i != par->grid_breg.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_breg.nx - 1) + par->grid_breg.x_min;
    const std::size_t idx_lv1{i * par->grid_breg.ny * par->grid_breg.nz};
    for (decltype(par->grid_breg.ny) j = 0; j != par->grid_breg.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_breg.ny - 1) + par->grid_breg.y_min;
      const std::size_t idx_lv2{idx_lv1 + j * par->grid_breg.nz};
      for (decltype(par->grid_breg.nz) k = 0; k != par->grid_breg.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_breg.nz - 1) + par->grid_breg.z_min;
        const std::size_t idx{idx_lv2 + k};
        tmp_vec = write_field(gc_pos, par);
        grid->bx[idx] = tmp_vec[0];
        grid->by[idx] = tmp_vec[1];
        grid->bz[idx] = tmp_vec[2];
      }
    }
  }
}
