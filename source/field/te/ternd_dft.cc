#include <cassert>
#include <cmath>
#include <omp.h>

#include <fftw3.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <grid.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>
#include <toolkit.h>

ham_float TErnd_dft::spectrum(const ham_float &k, const Param *par) const {
  if (k < par->ternd_dft.k0)
    return 0.;
  else
    return par->ternd_dft.rms *
           std::pow(k / par->ternd_dft.k0, par->ternd_dft.a0) /
           (4 * cgs::pi * k * k);
}

ham_float TErnd_dft::spatial_profile(const Hamvec<3, ham_float> &pos,
                                     const Param *par) const {
  const ham_float r_cyl{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]) -
                        std::fabs(par->observer[0])};
  const ham_float z{std::fabs(pos[2]) - std::fabs(par->observer[2])};
  return std::exp(-r_cyl / par->ternd_dft.r0) *
         std::exp(-z / par->ternd_dft.z0);
}

void TErnd_dft::write_grid(const Param *par, const TEreg *, const Grid_tereg *,
                           Grid_ternd *grid) const {
  // STEP I
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
  const ham_float lx{par->grid_ternd.x_max - par->grid_ternd.x_min};
  const ham_float ly{par->grid_ternd.y_max - par->grid_ternd.y_min};
  const ham_float lz{par->grid_ternd.z_max - par->grid_ternd.z_min};
  // physical k in 1/kpc dimension
  // physical dk^3
  const ham_float dk3{cgs::kpc * cgs::kpc * cgs::kpc / (lx * ly * lz)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_ternd.nx) i = 0; i < par->grid_ternd.nx; ++i) {
#ifdef _OPENMP
    auto seed_id = threadvec[omp_get_thread_num()];
#else
    auto seed_id = r;
#endif
    ham_float kx{cgs::kpc * i / lx};
    if (i >= (par->grid_ternd.nx + 1) / 2)
      kx -= cgs::kpc * par->grid_ternd.nx / lx;
    const size_t idx_lv1{i * par->grid_ternd.ny * par->grid_ternd.nz};
    for (decltype(par->grid_ternd.ny) j = 0; j < par->grid_ternd.ny; ++j) {
      ham_float ky{cgs::kpc * j / ly};
      if (j >= (par->grid_ternd.ny + 1) / 2)
        ky -= cgs::kpc * par->grid_ternd.ny / ly;
      const size_t idx_lv2{idx_lv1 + j * par->grid_ternd.nz};
      for (decltype(par->grid_ternd.nz) l = 0; l < par->grid_ternd.nz; ++l) {
        if (i == 0 and j == 0 and l == 0)
          continue;
        ham_float kz{cgs::kpc * l / lz};
        if (l >= (par->grid_ternd.nz + 1) / 2)
          kz -= cgs::kpc * par->grid_ternd.nz / lz;
        const ham_float ks{std::sqrt(kx * kx + ky * ky + kz * kz)};
        const size_t idx{idx_lv2 + l};
        // since we drop Im part after DFT
        // P ~ te_Re^2 ~ te_Im^2
        const ham_float sigma{std::sqrt(spectrum(ks, par) * dk3)};
        grid->te_k[idx][0] = sigma * gsl_ran_ugaussian(seed_id);
        grid->te_k[idx][1] = sigma * gsl_ran_ugaussian(seed_id);
      } // l
    }   // j
  }     // i
#ifdef _OPENMP
  for (int b = 0; b < omp_get_max_threads(); ++b)
    gsl_rng_free(threadvec[b]);
  delete[] threadvec;
#else
  gsl_rng_free(r);
#endif
  // execute DFT backward plan
  fftw_execute_dft(grid->plan_te_bw, grid->te_k, grid->te_k);
  // STEP II
  // RESCALING FIELD PROFILE IN REAL SPACE
  // 1/sqrt(te_var)
  const ham_float te_var_invsq{
      1. /
      std::sqrt(toolkit::variance(grid->te_k[0], par->grid_ternd.full_size))};
  assert(std::isfinite(te_var_invsq));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_ternd.nx) i = 0; i < par->grid_ternd.nx; ++i) {
    Hamvec<3, ham_float> pos{
        i * lx / (par->grid_ternd.nx - 1) + par->grid_ternd.x_min, 0, 0};
    const size_t idx_lv1{i * par->grid_ternd.ny * par->grid_ternd.nz};
    for (decltype(par->grid_ternd.ny) j = 0; j < par->grid_ternd.ny; ++j) {
      pos[1] = j * ly / (par->grid_ternd.ny - 1) + par->grid_ternd.y_min;
      const size_t idx_lv2{idx_lv1 + j * par->grid_ternd.nz};
      for (decltype(par->grid_ternd.nz) l = 0; l < par->grid_ternd.nz; ++l) {
        // get physical position
        pos[2] = l * lz / (par->grid_ternd.nz - 1) + par->grid_ternd.z_min;
        // get reprofiling factor
        ham_float ratio{std::sqrt(spatial_profile(pos, par)) *
                        par->ternd_dft.rms * te_var_invsq};
        const size_t idx{idx_lv2 + l};
        // manually pass back reprofiled Re part
        grid->te[idx] = grid->te_k[idx][0] * ratio;
      }
    }
  }
}
