#include <bfield.h>
#include <cassert>
#include <cmath>
#include <fftw3.h>
#include <grid.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <omp.h>
#include <param.h>
#include <stdexcept>
#include <toolkit.h>

// since we are using rms normalization
// p0 is hidden and not affecting anything
ham_float Bmodel_es::spectrum(const ham_float &k, const Param *par) const {
  // units fixing, wave vector in 1/kpc units
  const ham_float p0{par->bmodel_es.rms * par->bmodel_es.rms};
  const ham_float unit{1. / (4 * cgs::pi * k * k)};
  // power laws
  const ham_float band1{ham_float(k < par->bmodel_es.k1)};
  const ham_float band2{ham_float(k > par->bmodel_es.k1) *
                        ham_float(k < par->bmodel_es.k0)};
  const ham_float band3{ham_float(k > par->bmodel_es.k0)};
  const ham_float P{
      band1 *
          std::pow(par->bmodel_es.k0 / par->bmodel_es.k1, par->bmodel_es.a1) *
          std::pow(k / par->bmodel_es.k1, 6.0) +
      band2 / std::pow(k / par->bmodel_es.k0, par->bmodel_es.a1) +
      band3 / std::pow(k / par->bmodel_es.k0, par->bmodel_es.a0)};
  return P * p0 * unit;
}

// galactic scaling of random field energy density
// set to 1 at observer's place
ham_float Bmodel_es::spatial_profile(const Hamvec<3, ham_float> &pos,
                                     const Param *par) const {
  const ham_float r_cyl{std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]) -
                        std::fabs(par->observer[0])};
  const ham_float z{std::fabs(pos[2]) - std::fabs(par->observer[2])};
  return std::exp(-r_cyl / par->bmodel_es.r0) *
         std::exp(-z / par->bmodel_es.z0);
}

ham_float Bmodel_es::spatial_profile(const ham_uint &, const Param *,
                                     const Grid_b *) const {
  throw std::runtime_error("un-implemented");
  return 0.;
}

Hamvec<3, ham_float>
Bmodel_es::anisotropy_direction(const Hamvec<3, ham_float> &,
                                const Param *) const {
  throw std::runtime_error("un-implemented");
  return Hamvec<3, ham_float>(0., 0., 0.);
}

Hamvec<3, ham_float> Bmodel_es::anisotropy_direction(const ham_uint &idx,
                                                     const Param *,
                                                     const Grid_b *grid) const {
  return (read_grid(idx, grid)).versor();
}

ham_float Bmodel_es::anisotropy_ratio(const Hamvec<3, ham_float> &,
                                      const Param *par) const {
  // the simplest case, const.
  return par->bmodel_es.rho;
}

ham_float Bmodel_es::anisotropy_ratio(const ham_uint &, const Param *,
                                      const Grid_b *) const {
  throw std::runtime_error("un-implemented");
  return 0.;
}

// Gram-Schimdt, rewritten using Healpix vec3 library
// tiny error caused by machine is inevitable
Hamvec<3, ham_float>
Bmodel_es::gramschmidt(const Hamvec<3, ham_float> &k,
                       const Hamvec<3, ham_float> &b) const {
  if (k.lengthsq() == 0 or b.lengthsq() == 0)
    return b;
  const ham_float inv_k_mod{1. / k.lengthsq()};
  // multiply \sqrt(3) for preserving spectral power statistically
  return Hamvec<3, ham_float>{
      cgs::sqrt3 * (b[0] - k[0] * k.dotprod(b) * inv_k_mod),
      cgs::sqrt3 * (b[1] - k[1] * k.dotprod(b) * inv_k_mod),
      cgs::sqrt3 * (b[2] - k[2] * k.dotprod(b) * inv_k_mod)};
}

void Bmodel_es::write_grid(const Param *par, Grid_b *grid) const {
  // STEP I
  // GENERATE GAUSSIAN RANDOM FROM SPECTRUM
  // initialize random seed
#ifdef _OPENMP
  gsl_rng **threadvec = new gsl_rng *[omp_get_max_threads()];
  for (int b = 0; b < omp_get_max_threads(); ++b) {
    threadvec[b] = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(threadvec[b], b + toolkit::random_seed(par->bmodel_es.seed));
  }
#else
  gsl_rng *r{gsl_rng_alloc(gsl_rng_taus)};
  gsl_rng_set(r, toolkit::random_seed(par->bmodel_es.seed));
#endif
  // start Fourier space filling, physical k in 1/kpc dimension
  const ham_float lx{par->grid_b.x_max - par->grid_b.x_min};
  const ham_float ly{par->grid_b.y_max - par->grid_b.y_min};
  const ham_float lz{par->grid_b.z_max - par->grid_b.z_min};
  // physical dk^3
  const ham_float dk3{cgs::kpc * cgs::kpc * cgs::kpc / (lx * ly * lz)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (ham_uint i = 0; i < par->grid_b.nx; ++i) {
#ifdef _OPENMP
    auto seed_id = threadvec[omp_get_thread_num()];
#else
    auto seed_id = r;
#endif
    ham_float kx{cgs::kpc * i / lx};
    if (i >= (par->grid_b.nx + 1) / 2)
      kx -= cgs::kpc * par->grid_b.nx / lx;
    // it's faster to calculate indeces manually
    const ham_uint idx_lv1{i * par->grid_b.ny * par->grid_b.nz};
    for (ham_uint j = 0; j < par->grid_b.ny; ++j) {
      ham_float ky{cgs::kpc * j / ly};
      if (j >= (par->grid_b.ny + 1) / 2)
        ky -= cgs::kpc * par->grid_b.ny / ly;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_b.nz};
      for (ham_uint l = 0; l < par->grid_b.nz; ++l) {
        // 0th term is fixed to zero in allocation
        if (i == 0 and j == 0 and l == 0)
          continue;
        ham_float kz{cgs::kpc * l / lz};
        if (l >= (par->grid_b.nz + 1) / 2)
          kz -= cgs::kpc * par->grid_b.nz / lz;
        const ham_float ks{std::sqrt(kx * kx + ky * ky + kz * kz)};
        const ham_uint idx{idx_lv2 + l};
        // turbulent power is shared in following pattern
        // P ~ (bx^2 + by^2 + bz^2)
        // c0^2 ~ c1^2 ~ (bx^2 + by^2) ~ P*2/3
        // as renormalization comes in PHASE II,
        // 1/3, P0 in spectrum, dk3 are numerically redundant
        // while useful for precision check
        const ham_float sigma{std::sqrt(spectrum(ks, par) * dk3 / 3.)};
        grid->c0[idx][0] = sigma * gsl_ran_ugaussian(seed_id);
        grid->c0[idx][1] = sigma * gsl_ran_ugaussian(seed_id);
        grid->c1[idx][0] = sigma * gsl_ran_ugaussian(seed_id);
        grid->c1[idx][1] = sigma * gsl_ran_ugaussian(seed_id);
      } // l
    }   // j
  }     // i
  // ks=0 should be automatically addressed in P(k)
  // execute DFT backward plan
  fftw_execute_dft(grid->plan_c0_bw, grid->c0, grid->c0);
  fftw_execute_dft(grid->plan_c1_bw, grid->c1, grid->c1);
  // free random memory
#ifdef _OPENMP
  for (int b = 0; b < omp_get_max_threads(); ++b)
    gsl_rng_free(threadvec[b]);
  delete[] threadvec;
#else
  // free random memory
  gsl_rng_free(r);
#endif
  // STEP II
  // RESCALING FIELD PROFILE IN REAL SPACE
  // 1./std::sqrt(3*bi_var)
  // after 1st Fourier transformation, c0_R = bx, c0_I = by, c1_I = bz
  const ham_float b_var_invsq{
      1. /
      std::sqrt(3. * toolkit::variance(grid->c0[0], par->grid_b.full_size))};
  assert(std::isfinite(b_var_invsq));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (ham_uint i = 0; i < par->grid_b.nx; ++i) {
    Hamvec<3, ham_float> pos{i * lx / (par->grid_b.nx - 1) + par->grid_b.x_min,
                             0, 0};
    const ham_uint idx_lv1{i * par->grid_b.ny * par->grid_b.nz};
    for (ham_uint j = 0; j < par->grid_b.ny; ++j) {
      pos[1] = j * ly / (par->grid_b.ny - 1) + par->grid_b.y_min;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_b.nz};
      for (ham_uint l = 0; l < par->grid_b.nz; ++l) {
        // get physical position
        pos[2] = l * lz / (par->grid_b.nz - 1) + par->grid_b.z_min;
        // get reprofiling factor
        ham_float ratio{std::sqrt(spatial_profile(pos, par)) *
                        par->bmodel_es.rms * b_var_invsq};
        const ham_uint idx{idx_lv2 + l};
        // assemble b_Re
        // after 1st Fourier transformation, c0_R = bx, c0_I = by, c1_I = bz
        Hamvec<3, ham_float> b_re{grid->c0[idx][0] * ratio,
                                  grid->c0[idx][1] * ratio,
                                  grid->c1[idx][1] * ratio};
        // impose anisotropy
        Hamvec<3, ham_float> H_versor = anisotropy_direction(idx, par, grid);
        const ham_float rho{anisotropy_ratio(pos, par)};
        assert(rho >= 0.);
        const ham_float rho2 = rho * rho;
        const ham_float rhonorm =
            1. / std::sqrt(rho2 * cgs::onethird + cgs::twothirds / rho2);
        // zero regular field, no prefered anisotropy
        if (H_versor.lengthsq() < 1e-10)
          continue;
        Hamvec<3, ham_float> b_re_par{H_versor * H_versor.dotprod(b_re)};
        Hamvec<3, ham_float> b_re_perp{b_re - b_re_par};
        b_re = (b_re_par * rho + b_re_perp / rho) * rhonorm;
        // push b_re back to c0 and c1
        grid->c0[idx][0] = b_re[0];
        grid->c0[idx][1] = b_re[1];
        grid->c1[idx][0] = b_re[1];
        grid->c1[idx][1] = b_re[2];
      } // l
    }   // j
  }     // i
  // execute DFT forward plan
  fftw_execute_dft(grid->plan_c0_fw, grid->c0, grid->c0);
  fftw_execute_dft(grid->plan_c1_fw, grid->c1, grid->c1);
  // STEP III
  // RE-ORTHOGONALIZING IN FOURIER SPACE
  // Gram-Schmidt process
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (ham_uint i = 0; i < par->grid_b.nx; ++i) {
    // apply Hermitian symmetry
    ham_uint i_sym{par->grid_b.nx - i};
    if (i == 0)
      i_sym = i;
    Hamvec<3, ham_float> tmp_k{cgs::kpc * i / lx, 0, 0};
    // it's better to calculate indeces manually
    // just for reference, how indeces are calculated
    // const ham_uint idx
    // {toolkit::index3d(par->grid_b.nx,par->grid_b.ny,par->grid_b.nz,i,j,l)};
    // const ham_uint idx_sym
    // {toolkit::index3d(par->grid_b.nx,par->grid_b.ny,par->grid_b.nz,i_sym,j_sym,l_sym)};
    const ham_uint idx_lv1{i * par->grid_b.ny * par->grid_b.nz};
    const ham_uint idx_sym_lv1{i_sym * par->grid_b.ny * par->grid_b.nz};
    if (i >= (par->grid_b.nx + 1) / 2)
      tmp_k[0] -= cgs::kpc * par->grid_b.nx / lx;
    for (ham_uint j = 0; j < par->grid_b.ny; ++j) {
      // apply Hermitian symmetry
      ham_uint j_sym{par->grid_b.ny - j};
      if (j == 0)
        j_sym = j;
      tmp_k[1] = cgs::kpc * j / ly;
      if (j >= (par->grid_b.ny + 1) / 2)
        tmp_k[1] -= cgs::kpc * par->grid_b.ny / ly;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_b.nz};
      const ham_uint idx_sym_lv2{idx_sym_lv1 + j_sym * par->grid_b.nz};
      for (ham_uint l = 0; l < par->grid_b.nz; ++l) {
        // apply Hermitian symmetry
        ham_uint l_sym{par->grid_b.nz - l};
        if (l == 0)
          l_sym = l;
        tmp_k[2] = cgs::kpc * l / lz;
        if (l >= (par->grid_b.nz + 1) / 2)
          tmp_k[2] -= cgs::kpc * par->grid_b.nz / lz;
        const ham_uint idx{idx_lv2 + l};             // k
        const ham_uint idx_sym{idx_sym_lv2 + l_sym}; //-k
        // reconstruct bx,by,bz from c0,c1,c*0,c*1
        // c0(k) = bx(k) + i by(k)
        // c*0(-k) = bx(k) - i by(k)
        // c1(k) = by(k) + i bz(k)
        // c1*1(-k) = by(k) - i bz(k)
        const Hamvec<3, ham_float> tmp_b_re{
            0.5 * (grid->c0[idx][0] + grid->c0[idx_sym][0]),
            0.5 * (grid->c1[idx][0] + grid->c1[idx_sym][0]),
            0.5 * (grid->c1[idx_sym][1] + grid->c1[idx][1])};
        // Gram-Schmidt process
        const Hamvec<3, ham_float> free_b_re{gramschmidt(tmp_k, tmp_b_re)};
        // reassemble c0,c1 from bx,by,bz
        // c0(k) = bx(k) + i by(k)
        // c1(k) = by(k) + i bz(k)
        // we take only the real part of b, multiply it by sqrt(2)
        // cause after G-S process, conjugate symmetry might have been destroied
        // sqrt(2) preserve the total spectral power
        grid->c0[idx][0] = cgs::sqrt2 * free_b_re[0];
        grid->c0[idx][1] = cgs::sqrt2 * free_b_re[1];
        grid->c1[idx][0] = cgs::sqrt2 * free_b_re[1];
        grid->c1[idx][1] = cgs::sqrt2 * free_b_re[2];
      } // l
    }   // j
  }     // i
  // execute DFT backward plan
  fftw_execute_dft(grid->plan_c0_bw, grid->c0, grid->c0);
  fftw_execute_dft(grid->plan_c1_bw, grid->c1, grid->c1);
  // according to FFTW convention
  // transform forward followed by backword scale up array by nx*ny*nz
  const ham_float inv_grid_size = 1.0 / par->grid_b.full_size;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (ham_uint i = 0; i < par->grid_b.nx; ++i) {
    // apply Hermitian symmetry
    ham_uint i_sym{par->grid_b.nx - i};
    if (i == 0)
      i_sym = i;
    // it's faster to calculate indeces manually
    const ham_uint idx_lv1{i * par->grid_b.ny * par->grid_b.nz};
    const ham_uint idx_sym_lv1{i_sym * par->grid_b.ny * par->grid_b.nz};
    for (ham_uint j = 0; j < par->grid_b.ny; ++j) {
      // apply Hermitian symmetry
      ham_uint j_sym{par->grid_b.ny - j};
      if (j == 0)
        j_sym = j;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_b.nz};
      const ham_uint idx_sym_lv2{idx_sym_lv1 + j_sym * par->grid_b.nz};
      for (ham_uint l = 0; l < par->grid_b.nz; ++l) {
        // apply Hermitian symmetry
        ham_uint l_sym{par->grid_b.nz - l};
        if (l == 0)
          l_sym = l;
        const ham_uint idx{idx_lv2 + l};             // q
        const ham_uint idx_sym{idx_sym_lv2 + l_sym}; //-q
        // reconstruct bx,by,bz from c0,c1,c*0,c*1
        // c0(q) = bx(q) + i by(q)
        // c*0(-q) = bx(q) - i by(q)
        // c1(q) = by(q) + i bz(q)
        // c1*1(-q) = by(q) - i bz(q)
        // grid has already been filled, we add new component instead of reset
        grid->bx[idx] +=
            0.5 * (grid->c0[idx][0] + grid->c0[idx_sym][0]) * inv_grid_size;
        grid->by[idx] +=
            0.5 * (grid->c1[idx][0] + grid->c1[idx_sym][0]) * inv_grid_size;
        grid->bz[idx] +=
            0.5 * (grid->c1[idx_sym][1] + grid->c1[idx][1]) * inv_grid_size;
      }
    }
  }
}
