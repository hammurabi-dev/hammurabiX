#include <cassert>
#include <cmath>
#include <stdexcept>

#include <crefield.h>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

ham_float CREfield::read_field(const Hamvec<3, ham_float> &pos,
                               const ham_float &E, const Param *par,
                               const Grid_cre *grid) const {
  if (par->grid_cre.read_permission) {
    return read_grid(pos, E, par, grid);
  } else if (par->grid_cre.build_permission) {
    return write_field(pos, E, par);
  } else {
    return 0.;
  }
}

ham_float CREfield::write_field(const Hamvec<3, ham_float> &, const ham_float &,
                                const Param *) const {
  return 0.;
}

ham_float CREfield::flux_norm(const Hamvec<3, ham_float> &,
                              const Param *) const {
  throw std::runtime_error("wrong inheritance");
}

ham_float CREfield::flux_idx(const Hamvec<3, ham_float> &,
                             const Param *) const {
  throw std::runtime_error("wrong inheritance");
}

ham_float CREfield::spatial_profile(const Hamvec<3, ham_float> &,
                                    const Param *) const {
  throw std::runtime_error("wrong inheritance");
}

// linear interpolate in log(E) frame
// with tri-linear interpolation in spatial frame
// En in CGS units, return in [GeV m^2 s Sr]^-1
ham_float CREfield::read_grid(const Hamvec<3, ham_float> &pos,
                              const ham_float &En, const Param *par,
                              const Grid_cre *grid) const {
  // linear interpolation in log(E)
  ham_float tmp{std::log(En / par->grid_cre.E_min) /
                std::log(par->grid_cre.E_max / par->grid_cre.E_min)};
  if (tmp <= 0 or tmp >= par->grid_cre.nE - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nE) El{(ham_uint)std::floor(tmp)};
  const ham_float Ed{tmp - El};
  // linear interpolation in the spatial domain
  tmp = (par->grid_cre.nx - 1) * (pos[0] - par->grid_cre.x_min) /
        (par->grid_cre.x_max - par->grid_cre.x_min);
  if (tmp <= 0 or tmp >= par->grid_cre.nx - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) xl{(ham_uint)std::floor(tmp)};
  const ham_float xd{tmp - xl};
  tmp = (par->grid_cre.ny - 1) * (pos[1] - par->grid_cre.y_min) /
        (par->grid_cre.y_max - par->grid_cre.y_min);
  if (tmp <= 0 or tmp >= par->grid_cre.ny - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) yl{(ham_uint)std::floor(tmp)};
  const ham_float yd{tmp - yl};
  tmp = (par->grid_cre.nz - 1) * (pos[2] - par->grid_cre.z_min) /
        (par->grid_cre.z_max - par->grid_cre.z_min);
  if (tmp <= 0 or tmp >= par->grid_cre.nz - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) zl{(ham_uint)std::floor(tmp)};
  const ham_float zd{tmp - zl};
  assert(Ed >= 0 and Ed < 1 and xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and
         yd < 1 and zd < 1);
  // @ El
  ham_uint idx1{toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx,
                                 par->grid_cre.ny, par->grid_cre.nz, xl, yl, zl,
                                 El)};
  ham_uint idx2{toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx,
                                 par->grid_cre.ny, par->grid_cre.nz, xl, yl,
                                 zl + 1, El)};
  const ham_float i1{grid->cre_flux[idx1] * (1. - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl + 1, zl, El);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl + 1, zl + 1, El);
  const ham_float i2{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl, zl, El);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl, zl + 1, El);
  const ham_float j1{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl + 1, zl, El);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl + 1, zl + 1, El);
  const ham_float j2{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  const ham_float w1{i1 * (1 - yd) + i2 * yd};
  const ham_float w2{j1 * (1 - yd) + j2 * yd};
  const ham_float q1{w1 * (1 - xd) + w2 * xd};
  // @ El+1
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl, zl, El + 1);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl, zl + 1, El + 1);
  const ham_float i3{grid->cre_flux[idx1] * (1. - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl + 1, zl, El + 1);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl + 1, zl + 1, El + 1);
  const ham_float i4{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl, zl, El + 1);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl, zl + 1, El + 1);
  const ham_float j3{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl + 1, zl, El + 1);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl + 1, zl + 1, El + 1);
  const ham_float j4{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  const ham_float w3{i3 * (1 - yd) + i4 * yd};
  const ham_float w4{j3 * (1 - yd) + j4 * yd};
  const ham_float q2{w3 * (1 - xd) + w4 * xd};
  // linear interpolate in log(E)
  return q1 * (1 - Ed) + q2 * Ed;
}

ham_float CREfield::read_grid_num(const Hamvec<3, ham_float> &pos,
                                  const ham_uint &Eidx, const Param *par,
                                  const Grid_cre *grid) const {
  // linear interpolation
  ham_float tmp{(par->grid_cre.nx - 1) * (pos[0] - par->grid_cre.x_min) /
                (par->grid_cre.x_max - par->grid_cre.x_min)};
  if (tmp <= 0 or tmp >= par->grid_cre.nx - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) xl{(ham_uint)std::floor(tmp)};
  const ham_float xd{tmp - xl};
  tmp = (par->grid_cre.ny - 1) * (pos[1] - par->grid_cre.y_min) /
        (par->grid_cre.y_max - par->grid_cre.y_min);
  if (tmp <= 0 or tmp >= par->grid_cre.ny - 1) {
    return 0.;
  }
  decltype(par->grid_cre.ny) yl{(ham_uint)std::floor(tmp)};
  const ham_float yd{tmp - yl};
  tmp = (par->grid_cre.nz - 1) * (pos[2] - par->grid_cre.z_min) /
        (par->grid_cre.z_max - par->grid_cre.z_min);
  if (tmp <= 0 or tmp >= par->grid_cre.nz - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nz) zl{(ham_uint)std::floor(tmp)};
  const ham_float zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  ham_uint idx1{toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny,
                                 par->grid_cre.nz, par->grid_cre.nE, xl, yl, zl,
                                 Eidx)};
  ham_uint idx2{toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny,
                                 par->grid_cre.nz, par->grid_cre.nE, xl, yl,
                                 zl + 1, Eidx)};
  const ham_float i1{grid->cre_flux[idx1] * (1. - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl + 1, zl, Eidx);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl, yl + 1, zl + 1, Eidx);
  const ham_float i2{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl, zl, Eidx);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl, zl + 1, Eidx);
  const ham_float j1{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl + 1, zl, Eidx);
  idx2 = toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny, par->grid_cre.nz,
                          par->grid_cre.nE, xl + 1, yl + 1, zl + 1, Eidx);
  const ham_float j2{grid->cre_flux[idx1] * (1 - zd) +
                     grid->cre_flux[idx2] * zd};
  const ham_float w1{i1 * (1 - yd) + i2 * yd};
  const ham_float w2{j1 * (1 - yd) + j2 * yd};
  return w1 * (1 - xd) + w2 * xd;
}

// writing CRE DIFFERENTIAL density flux, in [GeV m^2 s sr]^-1
void CREfield::write_grid(const Param *par, Grid_cre *grid) const {
  assert(par->grid_cre.write_permission);
  Hamvec<3, ham_float> gc_pos;
  ham_float E;
  ham_float lx{par->grid_cre.x_max - par->grid_cre.x_min};
  ham_float ly{par->grid_cre.y_max - par->grid_cre.y_min};
  ham_float lz{par->grid_cre.z_max - par->grid_cre.z_min};
  for (decltype(par->grid_cre.nx) i = 0; i != par->grid_cre.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_cre.nx - 1) + par->grid_cre.x_min;
    const ham_uint idx_lv1{i * par->grid_cre.ny * par->grid_cre.nz *
                           par->grid_cre.nE};
    for (decltype(par->grid_cre.ny) j = 0; j != par->grid_cre.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_cre.ny - 1) + par->grid_cre.y_min;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_cre.nz * par->grid_cre.nE};
      for (decltype(par->grid_cre.nz) k = 0; k != par->grid_cre.nz; ++k) {
        gc_pos[2] = lz * k / (par->grid_cre.nz - 1) + par->grid_cre.z_min;
        const ham_uint idx_lv3{idx_lv2 + k * par->grid_cre.nE};
        for (decltype(par->grid_cre.nE) m = 0; m != par->grid_cre.nE; ++m) {
          E = par->grid_cre.E_min * std::exp(m * par->grid_cre.E_fact);
          const ham_uint idx{idx_lv3 + m};
          grid->cre_flux[idx] = write_field(gc_pos, E, par);
        }
      }
    }
  }
}
