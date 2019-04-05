#include <cassert>
#include <cmath>
#include <hvec.h>
#include <omp.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <cgs_units_file.h>
#include <fereg.h>
#include <grid.h>
#include <namespace_toolkit.h>
#include <param.h>

double FEreg::get_density(const hvec<3, double> &pos, const Param *par,
                          const Grid_fereg *grid) const {
  if (par->grid_fereg.read_permission) {
    return read_grid(pos, par, grid);
  } else if (par->grid_fereg.build_permission) {
    return density(pos, par);
  } else {
    return 0.;
  }
}

// not recommended to use without enough computing source
// recommend to use this once (replace density in write_grid) if no
// free parameters in FE
double FEreg::density_blur(const hvec<3, double> &pos, const Param *par) const {
  double ne_blur{0.};
  // sampling point number
  std::size_t step{1000};
  // gaussian blur scale
  double blur_scale_x{(par->grid_fereg.x_max - par->grid_fereg.x_min) /
                      (par->grid_fereg.nx * CGS_U_kpc)};
  double blur_scale_y{(par->grid_fereg.y_max - par->grid_fereg.y_min) /
                      (par->grid_fereg.ny * CGS_U_kpc)};
  double blur_scale_z{(par->grid_fereg.z_max - par->grid_fereg.z_min) /
                      (par->grid_fereg.nz * CGS_U_kpc)};
  // sample position
  hvec<3, double> pos_s;
  gsl_rng *r{gsl_rng_alloc(gsl_rng_taus)};
  gsl_rng_set(r, toolkit::random_seed(par->fernd_seed));
#pragma omp parallel for ordered schedule(static, 1) reduction(+ : ne_blur)
  for (decltype(step) i = 0; i < step; ++i) {
#pragma omp ordered
    {
      pos_s =
          pos + hvec<3, double>{
                    gsl_ran_gaussian(r, (blur_scale_x / 2.355)) * CGS_U_kpc,
                    gsl_ran_gaussian(r, (blur_scale_y / 2.355)) * CGS_U_kpc,
                    gsl_ran_gaussian(r, (blur_scale_z / 2.355)) * CGS_U_kpc};
    }
    ne_blur += density(pos_s, par);
  }
  gsl_rng_free(r);
  return ne_blur / step;
}

// if no specified field model is built
// FEreg object link directly here and return null field when invoked
double FEreg::density(const hvec<3, double> &, const Param *) const {
  return 0.;
}

double FEreg::read_grid(const hvec<3, double> &pos, const Param *par,
                        const Grid_fereg *grid) const {
  double tmp{(par->grid_fereg.nx - 1) * (pos[0] - par->grid_fereg.x_min) /
             (par->grid_fereg.x_max - par->grid_fereg.x_min)};
  if (tmp < 1 or tmp > par->grid_fereg.nx - 1) {
    return 0.;
  }
  decltype(par->grid_fereg.nx) xl{(std::size_t)std::floor(tmp)};
  const double xd = tmp - xl;

  tmp = (par->grid_fereg.ny - 1) * (pos[1] - par->grid_fereg.y_min) /
        (par->grid_fereg.y_max - par->grid_fereg.y_min);
  if (tmp < 1 or tmp > par->grid_fereg.ny - 1) {
    return 0.;
  }
  decltype(par->grid_fereg.nx) yl{(std::size_t)std::floor(tmp)};
  const double yd = tmp - yl;

  tmp = (par->grid_fereg.nz - 1) * (pos[2] - par->grid_fereg.z_min) /
        (par->grid_fereg.z_max - par->grid_fereg.z_min);
  if (tmp < 1 or tmp > par->grid_fereg.nz - 1) {
    return 0.;
  }
  decltype(par->grid_fereg.nx) zl{(std::size_t)std::floor(tmp)};
  const double zd = tmp - zl;
  assert(xd >= 0 and yd >= 0 and zd >= 0 and xd <= 1 and yd <= 1 and zd <= 1);
  double fe;
  if (xl + 1 < par->grid_fereg.nx and yl + 1 < par->grid_fereg.ny and
      zl + 1 < par->grid_fereg.nz) {
    std::size_t idx1{toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                                      par->grid_fereg.nz, xl, yl, zl)};
    std::size_t idx2{toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                                      par->grid_fereg.nz, xl, yl, zl + 1)};
    double i1{grid->fe[idx1] * (1. - zd) + grid->fe[idx2] * zd};
    idx1 = toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                            par->grid_fereg.nz, xl, yl + 1, zl);
    idx2 = toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                            par->grid_fereg.nz, xl, yl + 1, zl + 1);
    double i2{grid->fe[idx1] * (1 - zd) + grid->fe[idx2] * zd};
    idx1 = toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                            par->grid_fereg.nz, xl + 1, yl, zl);
    idx2 = toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                            par->grid_fereg.nz, xl + 1, yl, zl + 1);
    double j1{grid->fe[idx1] * (1 - zd) + grid->fe[idx2] * zd};
    idx1 = toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                            par->grid_fereg.nz, xl + 1, yl + 1, zl);
    idx2 = toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                            par->grid_fereg.nz, xl + 1, yl + 1, zl + 1);
    double j2{grid->fe[idx1] * (1 - zd) + grid->fe[idx2] * zd};
    double w1{i1 * (1 - yd) + i2 * yd};
    double w2{j1 * (1 - yd) + j2 * yd};
    fe = (w1 * (1 - xd) + w2 * xd);
  } else {
    std::size_t idx1{toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                                      par->grid_fereg.nz, xl, yl, zl)};
    fe = grid->fe[idx1];
  }
  assert(fe >= 0);
  return fe;
}

void FEreg::write_grid(const Param *par, Grid_fereg *grid) const {
  assert(par->grid_fereg.write_permission);
  hvec<3, double> gc_pos;
  double lx{par->grid_fereg.x_max - par->grid_fereg.x_min};
  double ly{par->grid_fereg.y_max - par->grid_fereg.y_min};
  double lz{par->grid_fereg.z_max - par->grid_fereg.z_min};
  for (decltype(par->grid_fereg.nx) i = 0; i != par->grid_fereg.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_fereg.nx - 1) + par->grid_fereg.x_min;
    for (decltype(par->grid_fereg.ny) j = 0; j != par->grid_fereg.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_fereg.ny - 1) + par->grid_fereg.y_min;
      for (decltype(par->grid_fereg.nz) k = 0; k != par->grid_fereg.nz; ++k) {
        std::size_t idx{toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                                         par->grid_fereg.nz, i, j, k)};
        gc_pos[2] = lz * k / (par->grid_fereg.nz - 1) + par->grid_fereg.z_min;
        // two solutions
        // grid->fe[idx] = density_blur(gc_pos, par, grid);
        grid->fe[idx] = density(gc_pos, par);
      }
    }
  }
}

// END
