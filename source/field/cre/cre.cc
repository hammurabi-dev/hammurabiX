#include <cassert>
#include <cmath>
#include <stdexcept>

#include <crefield.h>
#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

double CRE::read_field(const hamvec<3, double> &pos, const double &E,
                       const Param *par, const Grid_cre *grid) const {
  if (par->grid_cre.read_permission) {
    return read_grid(pos, E, par, grid);
  } else if (par->grid_cre.build_permission) {
    return write_field(pos, E, par);
  } else {
    return 0.;
  }
}

double CRE::write_field(const hamvec<3, double> &, const double &,
                        const Param *) const {
  return 0.;
}

double CRE::flux_norm(const hamvec<3, double> &, const Param *) const {
  throw std::runtime_error("wrong inheritance");
}

double CRE::flux_idx(const hamvec<3, double> &, const Param *) const {
  throw std::runtime_error("wrong inheritance");
}

double CRE::spatial_profile(const hamvec<3, double> &, const Param *) const {
  throw std::runtime_error("wrong inheritance");
}

// linear interpolate in log(E) frame
// with tri-linear interpolation in spatial frame
// En in CGS units, return in [GeV m^2 s Sr]^-1
double CRE::read_grid(const hamvec<3, double> &pos, const double &En,
                      const Param *par, const Grid_cre *grid) const {
  // linear interpolation in log(E)
  double tmp{std::log(En / par->grid_cre.E_min) /
             std::log(par->grid_cre.E_max / par->grid_cre.E_min)};
  if (tmp <= 0 or tmp >= par->grid_cre.nE - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nE) El{(std::size_t)std::floor(tmp)};
  const double Ed{tmp - El};
  // linear interpolation in the spatial domain
  tmp = (par->grid_cre.nx - 1) * (pos[0] - par->grid_cre.x_min) /
        (par->grid_cre.x_max - par->grid_cre.x_min);
  if (tmp <= 0 or tmp >= par->grid_cre.nx - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd{tmp - xl};
  tmp = (par->grid_cre.ny - 1) * (pos[1] - par->grid_cre.y_min) /
        (par->grid_cre.y_max - par->grid_cre.y_min);
  if (tmp <= 0 or tmp >= par->grid_cre.ny - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_cre.nz - 1) * (pos[2] - par->grid_cre.z_min) /
        (par->grid_cre.z_max - par->grid_cre.z_min);
  if (tmp <= 0 or tmp >= par->grid_cre.nz - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) zl{(std::size_t)std::floor(tmp)};
  const double zd{tmp - zl};
  assert(Ed >= 0 and Ed < 1 and xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and
         yd < 1 and zd < 1);
  // @ El
  std::size_t idx1{toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx,
                                    par->grid_cre.ny, par->grid_cre.nz, El, xl,
                                    yl, zl)};
  std::size_t idx2{toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx,
                                    par->grid_cre.ny, par->grid_cre.nz, El, xl,
                                    yl, zl + 1)};
  const double i1{grid->cre_flux[idx1] * (1. - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El, xl, yl + 1, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El, xl, yl + 1, zl + 1);
  const double i2{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El, xl + 1, yl, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El, xl + 1, yl, zl + 1);
  const double j1{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El, xl + 1, yl + 1, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El, xl + 1, yl + 1, zl + 1);
  const double j2{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  const double w1{i1 * (1 - yd) + i2 * yd};
  const double w2{j1 * (1 - yd) + j2 * yd};
  const double q1{w1 * (1 - xd) + w2 * xd};
  // @ El+1
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl, yl, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl, yl, zl + 1);
  const double i3{grid->cre_flux[idx1] * (1. - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl, yl + 1, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl, yl + 1, zl + 1);
  const double i4{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl + 1, yl, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl + 1, yl, zl + 1);
  const double j3{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl + 1, yl + 1, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, El + 1, xl + 1, yl + 1, zl + 1);
  const double j4{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  const double w3{i3 * (1 - yd) + i4 * yd};
  const double w4{j3 * (1 - yd) + j4 * yd};
  const double q2{w3 * (1 - xd) + w4 * xd};
  // linear interpolate in log(E)
  return q1 * (1 - Ed) + q2 * Ed;
}

double CRE::read_grid_num(const hamvec<3, double> &pos, const std::size_t &Eidx,
                          const Param *par, const Grid_cre *grid) const {
  // linear interpolation
  double tmp{(par->grid_cre.nx - 1) * (pos[0] - par->grid_cre.x_min) /
             (par->grid_cre.x_max - par->grid_cre.x_min)};
  if (tmp <= 0 or tmp >= par->grid_cre.nx - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd{tmp - xl};
  tmp = (par->grid_cre.ny - 1) * (pos[1] - par->grid_cre.y_min) /
        (par->grid_cre.y_max - par->grid_cre.y_min);
  if (tmp <= 0 or tmp >= par->grid_cre.ny - 1) {
    return 0.;
  }
  decltype(par->grid_cre.ny) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_cre.nz - 1) * (pos[2] - par->grid_cre.z_min) /
        (par->grid_cre.z_max - par->grid_cre.z_min);
  if (tmp <= 0 or tmp >= par->grid_cre.nz - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nz) zl{(std::size_t)std::floor(tmp)};
  const double zd{tmp - zl};
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd < 1 and yd < 1 and zd < 1);
  std::size_t idx1{toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx,
                                    par->grid_cre.ny, par->grid_cre.nz, Eidx,
                                    xl, yl, zl)};
  std::size_t idx2{toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx,
                                    par->grid_cre.ny, par->grid_cre.nz, Eidx,
                                    xl, yl, zl + 1)};
  const double i1{grid->cre_flux[idx1] * (1. - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, Eidx, xl, yl + 1, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, Eidx, xl, yl + 1, zl + 1);
  const double i2{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, Eidx, xl + 1, yl, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, Eidx, xl + 1, yl, zl + 1);
  const double j1{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  idx1 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, Eidx, xl + 1, yl + 1, zl);
  idx2 = toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx, par->grid_cre.ny,
                          par->grid_cre.nz, Eidx, xl + 1, yl + 1, zl + 1);
  const double j2{grid->cre_flux[idx1] * (1 - zd) + grid->cre_flux[idx2] * zd};
  const double w1{i1 * (1 - yd) + i2 * yd};
  const double w2{j1 * (1 - yd) + j2 * yd};
  return w1 * (1 - xd) + w2 * xd;
}

// writing CRE DIFFERENTIAL density flux, in [GeV m^2 s sr]^-1
void CRE::write_grid(const Param *par, Grid_cre *grid) const {
  assert(par->grid_cre.write_permission);
  hamvec<3, double> gc_pos;
  double E;
  double lx{par->grid_cre.x_max - par->grid_cre.x_min};
  double ly{par->grid_cre.y_max - par->grid_cre.y_min};
  double lz{par->grid_cre.z_max - par->grid_cre.z_min};
  for (decltype(par->grid_cre.nE) i = 0; i != par->grid_cre.nE; ++i) {
    E = par->grid_cre.E_min * std::exp(i * par->grid_cre.E_fact);
    const std::size_t idx_lv1{i * par->grid_cre.nx * par->grid_cre.ny *
                              par->grid_cre.nz};
    for (decltype(par->grid_cre.nx) j = 0; j != par->grid_cre.nx; ++j) {
      gc_pos[0] = lx * j / (par->grid_cre.nx - 1) + par->grid_cre.x_min;
      const std::size_t idx_lv2{idx_lv1 +
                                j * par->grid_cre.ny * par->grid_cre.nz};
      for (decltype(par->grid_cre.ny) k = 0; k != par->grid_cre.ny; ++k) {
        gc_pos[1] = ly * k / (par->grid_cre.ny - 1) + par->grid_cre.y_min;
        const std::size_t idx_lv3{idx_lv2 + k * par->grid_cre.nz};
        for (decltype(par->grid_cre.nz) m = 0; m != par->grid_cre.nz; ++m) {
          gc_pos[2] = lz * m / (par->grid_cre.nz - 1) + par->grid_cre.z_min;
          const std::size_t idx{idx_lv3 + m};
          grid->cre_flux[idx] = write_field(gc_pos, E, par);
        }
      }
    }
  }
}
