#include <cmath>
#include <omp.h>

#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <hvec.h>

#include <cgs_units_file.h>
#include <fereg.h>
#include <fernd.h>
#include <grid.h>
#include <namespace_toolkit.h>
#include <param.h>

// since we are using rms normalization
// p0 is hidden and not affecting anything
double FErnd_dft::spec(const double &k, const Param *par) const {
  if (k < par->fernd_dft.k0)
    return 0.;
  else
    return par->fernd_dft.rms *
           std::pow(k / par->fernd_dft.k0, par->fernd_dft.a0) /
           (4 * CGS_U_pi * k * k);
}

// galactic scaling of random field energy density
// set to 1 at observer's place
double FErnd_dft::rescal(const hvec<3, double> &pos, const Param *par) const {
  const double r_cyl{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]) -
                     std::fabs(par->observer[0])};
  const double z{std::fabs(pos[2]) - std::fabs(par->observer[2])};
  return std::exp(-r_cyl / par->fernd_dft.r0) *
         std::exp(-z / par->fernd_dft.z0);
}

void FErnd_dft::write_grid(const Param *par, Grid_fernd *grid) const {
  // PHASE I
  // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
  // initialize random seed
#ifdef _OPENMP
  gsl_rng **threadvec = new gsl_rng *[omp_get_max_threads()];
  for (int b = 0; b < omp_get_max_threads(); ++b) {
    threadvec[b] = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(threadvec[b], b + toolkit::random_seed(par->brnd_seed));
  }
#else
  gsl_rng *r{gsl_rng_alloc(gsl_rng_taus)};
  gsl_rng_set(r, toolkit::random_seed(par->brnd_seed));
#endif
  const double lx{par->grid_fernd.x_max - par->grid_fernd.x_min};
  const double ly{par->grid_fernd.y_max - par->grid_fernd.y_min};
  const double lz{par->grid_fernd.z_max - par->grid_fernd.z_min};
  // physical k in 1/kpc dimension
  // physical dk^3
  const double dk3{CGS_U_kpc * CGS_U_kpc * CGS_U_kpc / (lx * ly * lz)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_fernd.nx) i = 0; i < par->grid_fernd.nx; ++i) {
#ifdef _OPENMP
    auto seed_id = threadvec[omp_get_thread_num()];
#else
    auto seed_id = r;
#endif
    double kx{CGS_U_kpc * i / lx};
    if (i >= (par->grid_fernd.nx + 1) / 2)
      kx -= CGS_U_kpc * par->grid_fernd.nx / lx;
    /**
     * it's better to calculate indeces manually
     * just for reference, how indeces are calculated
     * const size_t idx
     * {toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,i,j,l)};
     */
    const size_t idx_lv1{i * par->grid_fernd.ny * par->grid_fernd.nz};
    for (decltype(par->grid_fernd.ny) j = 0; j < par->grid_fernd.ny; ++j) {
      double ky{CGS_U_kpc * j / ly};
      if (j >= (par->grid_fernd.ny + 1) / 2)
        ky -= CGS_U_kpc * par->grid_fernd.ny / ly;
      const size_t idx_lv2{idx_lv1 + j * par->grid_fernd.nz};
      for (decltype(par->grid_fernd.nz) l = 0; l < par->grid_fernd.nz; ++l) {
        /**
         * the very 0th term is fixed to zero in allocation
         */
        if (i == 0 and j == 0 and l == 0)
          continue;
        double kz{CGS_U_kpc * l / lz};
        if (l >= (par->grid_fernd.nz + 1) / 2)
          kz -= CGS_U_kpc * par->grid_fernd.nz / lz;
        const double ks{std::sqrt(kx * kx + ky * ky + kz * kz)};
        const size_t idx{idx_lv2 + l};
        /**
         * since we drop Im part after DFT
         * P ~ fe_Re^2 ~ fe_Im^2
         */
        const double sigma{std::sqrt(spec(ks, par) * dk3)};
        grid->fe_k[idx][0] = sigma * gsl_ran_ugaussian(seed_id);
        grid->fe_k[idx][1] = sigma * gsl_ran_ugaussian(seed_id);
      } // l
    }   // j
  }     // i
        // ks=0 should be automatically addressed in P(k)
        // free random memory
#ifdef _OPENMP
  for (int b = 0; b < omp_get_max_threads(); ++b)
    gsl_rng_free(threadvec[b]);
  delete[] threadvec;
#else
  gsl_rng_free(r);
#endif
  // execute DFT backward plan
  fftw_execute_dft(grid->plan_fe_bw, grid->fe_k, grid->fe_k);
  // PHASE II
  // RESCALING FIELD PROFILE IN REAL SPACE
  // 1/sqrt(fe_var)
  const double fe_var_invsq{
      1. /
      std::sqrt(toolkit::variance(grid->fe_k[0], par->grid_fernd.full_size))};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_fernd.nx) i = 0; i < par->grid_fernd.nx; ++i) {
    hvec<3, double> pos{
        i * lx / (par->grid_fernd.nx - 1) + par->grid_fernd.x_min, 0, 0};
    /**
     * it's better to calculate indeces manually
     * just for reference, how indeces are calculated
     * const size_t idx
     * {toolkit::index3d(par->grid_fernd.nx,par->grid_fernd.ny,par->grid_fernd.nz,i,j,l)};
     */
    const size_t idx_lv1{i * par->grid_fernd.ny * par->grid_fernd.nz};
    for (decltype(par->grid_fernd.ny) j = 0; j < par->grid_fernd.ny; ++j) {
      pos[1] = j * ly / (par->grid_fernd.ny - 1) + par->grid_fernd.y_min;
      const size_t idx_lv2{idx_lv1 + j * par->grid_fernd.nz};
      for (decltype(par->grid_fernd.nz) l = 0; l < par->grid_fernd.nz; ++l) {
        // get physical position
        pos[2] = l * lz / (par->grid_fernd.nz - 1) + par->grid_fernd.z_min;
        // get rescaling factor
        double ratio{std::sqrt(rescal(pos, par)) * par->fernd_dft.rms *
                     fe_var_invsq};
        const size_t idx{idx_lv2 + l};
        // manually pass back rescaled Re part
        grid->fe[idx] = grid->fe_k[idx][0] * ratio;
      }
    }
  }
}
