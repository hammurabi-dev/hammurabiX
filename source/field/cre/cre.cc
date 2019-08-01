#include <cassert>
#include <cmath>
#include <stdexcept>

#include <crefield.h>
#include <grid.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

double CRE::flux(const hamvec<3, double> &, const Param *,
                 const double &) const {
  throw std::runtime_error("wrong inheritance");
}

double CRE::read_emissivity_t(const hamvec<3, double> &, const Param *,
                              const Grid_cre *, const double &) const {
  throw std::runtime_error("wrong inheritance");
}

double CRE::read_emissivity_p(const hamvec<3, double> &, const Param *,
                              const Grid_cre *, const double &) const {
  throw std::runtime_error("wrong inheritance");
}

double CRE::read_grid(const std::size_t &Eidx, const hamvec<3, double> &pos,
                      const Param *par, const Grid_cre *grid) const {
  // trilinear interpolation
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
  decltype(par->grid_cre.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd{tmp - yl};
  tmp = (par->grid_cre.nz - 1) * (pos[2] - par->grid_cre.z_min) /
        (par->grid_cre.z_max - par->grid_cre.z_min);
  if (tmp <= 0 or tmp >= par->grid_cre.nz - 1) {
    return 0.;
  }
  decltype(par->grid_cre.nx) zl{(std::size_t)std::floor(tmp)};
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
  const double cre{w1 * (1 - xd) + w2 * xd};
  return cre;
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
          grid->cre_flux[idx] = flux(gc_pos, par, E);
        }
      }
    }
  }
}
