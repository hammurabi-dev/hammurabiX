#include <cassert>
#include <cmath>
#include <omp.h>

#include <bfield.h>
#include <fftw3.h>
#include <grid.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

// global anisotropic turbulent field
Hamvec<3, ham_float>
Brnd_jf12::anisotropy_direction(const Hamvec<3, ham_float> &pos,
                                const Param *par, const Breg *breg,
                                const Grid_breg *gbreg) const {
  return (breg->read_field(pos, par, gbreg)).versor();
}

// global anisotropic turbulent field
ham_float Brnd_jf12::anisotropy_ratio(const Hamvec<3, ham_float> &,
                                      const Param *par, const Breg *,
                                      const Grid_breg *) const {
  // the simplest case, const.
  return par->brnd_jf12.rho;
}

// since we are using rms normalization
// p0 is hidden and not affecting anything
ham_float Brnd_jf12::spectrum(const ham_float &k, const Param *par) const {
  // units fixing, wave vector in 1/kpc units
  const ham_float p0{par->brnd_jf12.rms * par->brnd_jf12.rms};
  const ham_float unit = 1. / (4 * cgs::pi * k * k);
  // power laws
  const ham_float band1{ham_float(k < par->brnd_jf12.k1)};
  const ham_float band2{ham_float(k > par->brnd_jf12.k1) *
                        ham_float(k < par->brnd_jf12.k0)};
  const ham_float band3{ham_float(k > par->brnd_jf12.k0)};
  const ham_float P =
      band1 *
          std::pow(par->brnd_jf12.k0 / par->brnd_jf12.k1, par->brnd_jf12.a1) *
          std::pow(k / par->brnd_jf12.k1, 6.0) +
      band2 / std::pow(k / par->brnd_jf12.k0, par->brnd_jf12.a1) +
      band3 / std::pow(k / par->brnd_jf12.k0, par->brnd_jf12.a0);
  return P * p0 * unit;
}

// galactic scaling of random field energy density
// set to 1 at observer's place
ham_float Brnd_jf12::spatial_profile(const Hamvec<3, ham_float> &pos,
                                     const Param *par) const {

  const ham_float b0_1 = par->brnd_jf12.b0_1;
  const ham_float b0_2 = par->brnd_jf12.b0_2;
  const ham_float b0_3 = par->brnd_jf12.b0_3;
  const ham_float b0_4 = par->brnd_jf12.b0_4;
  const ham_float b0_5 = par->brnd_jf12.b0_5;
  const ham_float b0_6 = par->brnd_jf12.b0_6;
  const ham_float b0_7 = par->brnd_jf12.b0_7;
  const ham_float b0_8 = par->brnd_jf12.b0_8;

  const ham_float b0_int = par->brnd_jf12.b0_int;
  const ham_float z0_spiral = par->brnd_jf12.z0_spiral;
  const ham_float b0_halo = par->brnd_jf12.b0_halo;
  const ham_float r0_halo = par->brnd_jf12.r0_halo;
  const ham_float z0_halo = par->brnd_jf12.z0_halo;

  const ham_float Rmax = 20. * cgs::kpc;
  const ham_float rho_GC = 1. * cgs::kpc;
  const ham_float r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  const ham_float rho{
      sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2])};
  const ham_float z{pos[2]};
  const ham_float phi{atan2(pos[1], pos[0])};

  const ham_float rc_B[8] = {
      5.1, 6.3,  7.1,  8.3,
      9.8, 11.4, 12.7, 15.5}; // neg x crossings of spiral arms
  const ham_float inc = 11.5; // inclination, in degrees

  const ham_float b_arms[8] = {b0_1, b0_2, b0_3, b0_4, b0_5, b0_6, b0_7, b0_8};

  ham_float scaling_disk = 0.0 * cgs::muGauss;
  ham_float scaling_halo = 0.0 * cgs::muGauss;

  // boundaries outside which B is zero, not sure if this works?
  if (r > Rmax || rho < rho_GC) {
    return 0.0;
  }
  if (r < 5. * cgs::kpc) {
    scaling_disk = b0_int;
  } else {
    ham_float r_negx =
        r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi - M_PI));
    if (r_negx > rc_B[7] * cgs::kpc) {
      r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + M_PI));
    }
    if (r_negx > rc_B[7] * cgs::kpc) {
      r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + 3 * M_PI));
    }
    for (int i = 7; i >= 0; i--) {
      if (r_negx < rc_B[i] * cgs::kpc) {
        scaling_disk = b_arms[i] * (5. * cgs::kpc) / r;
      }
    } // "region 8,7,6,..,2"
  }

  scaling_disk = scaling_disk * exp(-0.5 * z * z / (z0_spiral * z0_spiral));
  scaling_halo =
      b0_halo * exp(-std::fabs(r / r0_halo)) * exp(-std::fabs(z / z0_halo));

  return (scaling_disk * scaling_disk + scaling_halo * scaling_halo) /
         (cgs::muGauss * cgs::muGauss);
}

// Gram-Schimdt, rewritten using Healpix vec3 library
// tiny error caused by machine is inevitable
Hamvec<3, ham_float>
Brnd_jf12::gramschmidt(const Hamvec<3, ham_float> &k,
                       const Hamvec<3, ham_float> &b) const {
  if (k.lengthsq() == 0 or b.lengthsq() == 0) {
    return b;
  }
  const ham_float inv_k_mod = 1. / k.lengthsq();
  // multiply \sqrt(3) for preserving spectral power statistically
  return Hamvec<3, ham_float>{
      1.73205081 * (b[0] - k[0] * k.dotprod(b) * inv_k_mod),
      1.73205081 * (b[1] - k[1] * k.dotprod(b) * inv_k_mod),
      1.73205081 * (b[2] - k[2] * k.dotprod(b) * inv_k_mod)};
}

void Brnd_jf12::write_grid(const Param *par, const Breg *breg,
                           const Grid_breg *gbreg, Grid_brnd *grid) const {
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
  // start Fourier space filling, physical k in 1/kpc dimension
  const ham_float lx{par->grid_brnd.x_max - par->grid_brnd.x_min};
  const ham_float ly{par->grid_brnd.y_max - par->grid_brnd.y_min};
  const ham_float lz{par->grid_brnd.z_max - par->grid_brnd.z_min};
  // physical dk^3
  const ham_float dk3{cgs::kpc * cgs::kpc * cgs::kpc / (lx * ly * lz)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_brnd.nx) i = 0; i < par->grid_brnd.nx; ++i) {
#ifdef _OPENMP
    auto seed_id = threadvec[omp_get_thread_num()];
#else
    auto seed_id = r;
#endif
    ham_float kx{cgs::kpc * i / lx};
    if (i >= (par->grid_brnd.nx + 1) / 2)
      kx -= cgs::kpc * par->grid_brnd.nx / lx;
    // it's faster to calculate indeces manually
    const ham_uint idx_lv1{i * par->grid_brnd.ny * par->grid_brnd.nz};
    for (decltype(par->grid_brnd.ny) j = 0; j < par->grid_brnd.ny; ++j) {
      ham_float ky{cgs::kpc * j / ly};
      if (j >= (par->grid_brnd.ny + 1) / 2)
        ky -= cgs::kpc * par->grid_brnd.ny / ly;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_brnd.nz};
      for (decltype(par->grid_brnd.nz) l = 0; l < par->grid_brnd.nz; ++l) {
        // 0th term is fixed to zero in allocation
        if (i == 0 and j == 0 and l == 0)
          continue;
        ham_float kz{cgs::kpc * l / lz};
        if (l >= (par->grid_brnd.nz + 1) / 2)
          kz -= cgs::kpc * par->grid_brnd.nz / lz;
        const ham_float ks{std::sqrt(kx * kx + ky * ky + kz * kz)};
        const ham_uint idx{idx_lv2 + l};
        // turbulent power is shared in following pattern
        // P ~ (bx^2 + by^2 + bz^2)
        // c0^2 ~ c1^2 ~ (bx^2 + by^2) ~ P*2/3
        // as renormalization comes in PHASE II,
        // 1/3, P0 in spectrum, dk3 are numerically redundant
        // while useful for precision check
        const ham_float sigma{std::sqrt(0.33333333 * spectrum(ks, par) * dk3)};
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
      std::sqrt(3. * toolkit::variance(grid->c0[0], par->grid_brnd.full_size))};
  assert(std::isfinite(b_var_invsq));
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_brnd.nx) i = 0; i < par->grid_brnd.nx; ++i) {
    Hamvec<3, ham_float> pos{
        i * lx / (par->grid_brnd.nx - 1) + par->grid_brnd.x_min, 0, 0};
    const ham_uint idx_lv1{i * par->grid_brnd.ny * par->grid_brnd.nz};
    for (decltype(par->grid_brnd.ny) j = 0; j < par->grid_brnd.ny; ++j) {
      pos[1] = j * ly / (par->grid_brnd.ny - 1) + par->grid_brnd.y_min;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_brnd.nz};
      for (decltype(par->grid_brnd.nz) l = 0; l < par->grid_brnd.nz; ++l) {
        // get physical position
        pos[2] = l * lz / (par->grid_brnd.nz - 1) + par->grid_brnd.z_min;
        // get reprofiling factor
        ham_float ratio{std::sqrt(spatial_profile(pos, par)) *
                        par->brnd_jf12.rms * b_var_invsq};
        const ham_uint idx{idx_lv2 + l};
        // assemble b_Re
        // after 1st Fourier transformation, c0_R = bx, c0_I = by, c1_I = bz
        Hamvec<3, ham_float> b_re{grid->c0[idx][0] * ratio,
                                  grid->c0[idx][1] * ratio,
                                  grid->c1[idx][1] * ratio};
        // impose anisotropy
        Hamvec<3, ham_float> H_versor =
            anisotropy_direction(pos, par, breg, gbreg);
        const ham_float rho{anisotropy_ratio(pos, par, breg, gbreg)};
        assert(rho >= 0.);
        const ham_float rho2 = rho * rho;
        const ham_float rhonorm =
            1. / std::sqrt(0.33333333 * rho2 + 0.66666667 / rho2);
        if (H_versor.lengthsq() <
            1e-10) // zero regular field, no prefered anisotropy
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
  for (decltype(par->grid_brnd.nx) i = 0; i < par->grid_brnd.nx; ++i) {
    decltype(par->grid_brnd.nx) i_sym{par->grid_brnd.nx -
                                      i}; // apply Hermitian symmetry
    if (i == 0)
      i_sym = i;
    Hamvec<3, ham_float> tmp_k{cgs::kpc * i / lx, 0, 0};
    // it's better to calculate indeces manually
    // just for reference, how indeces are calculated
    // const ham_uint idx
    // {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i,j,l)};
    // const ham_uint idx_sym
    // {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i_sym,j_sym,l_sym)};
    const ham_uint idx_lv1{i * par->grid_brnd.ny * par->grid_brnd.nz};
    const ham_uint idx_sym_lv1{i_sym * par->grid_brnd.ny * par->grid_brnd.nz};
    if (i >= (par->grid_brnd.nx + 1) / 2)
      tmp_k[0] -= cgs::kpc * par->grid_brnd.nx / lx;
    for (decltype(par->grid_brnd.ny) j = 0; j < par->grid_brnd.ny; ++j) {
      decltype(par->grid_brnd.ny) j_sym{par->grid_brnd.ny -
                                        j}; // apply Hermitian symmetry
      if (j == 0)
        j_sym = j;
      tmp_k[1] = cgs::kpc * j / ly;
      if (j >= (par->grid_brnd.ny + 1) / 2)
        tmp_k[1] -= cgs::kpc * par->grid_brnd.ny / ly;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_brnd.nz};
      const ham_uint idx_sym_lv2{idx_sym_lv1 + j_sym * par->grid_brnd.nz};
      for (decltype(par->grid_brnd.nz) l = 0; l < par->grid_brnd.nz; ++l) {
        decltype(par->grid_brnd.nz) l_sym{par->grid_brnd.nz -
                                          l}; // apply Hermitian symmetry
        if (l == 0)
          l_sym = l;
        tmp_k[2] = cgs::kpc * l / lz;
        if (l >= (par->grid_brnd.nz + 1) / 2)
          tmp_k[2] -= cgs::kpc * par->grid_brnd.nz / lz;
        const ham_uint idx{idx_lv2 + l};             // k
        const ham_uint idx_sym{idx_sym_lv2 + l_sym}; //-k
        // reconstruct bx,by,bz from c0,c1,c*0,c*1, keep real parts
        const Hamvec<3, ham_float> tmp_b_re{
            0.5 * (grid->c0[idx][0] + grid->c0[idx_sym][0]),
            0.5 * (grid->c1[idx][0] + grid->c1[idx_sym][0]),
            0.5 * (grid->c1[idx][1] + grid->c1[idx_sym][1])};
        // Gram-Schmidt process
        const Hamvec<3, ham_float> free_b_re{gramschmidt(tmp_k, tmp_b_re)};
        // reassemble c0,c1 from bx,by,bz
        // c0(k) = bx(k) + i by(k)
        // c1(k) = by(k) + i bz(k)
        // we take only the real part of b, multiply it by sqrt(2)
        // cause after G-S process, conjugate symmetry might have been destroied
        // sqrt(2) preserve the total spectral power
        grid->bx[idx] = cgs::sqrtwo * free_b_re[0];
        grid->by[idx] = cgs::sqrtwo * free_b_re[1];
        grid->bz[idx] = cgs::sqrtwo * free_b_re[2];
      } // l
    }   // j
  }     // i
  // re-assign the field back, avoid thread crash
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_brnd.full_size) idx = 0;
       idx < par->grid_brnd.full_size; ++idx) {
    grid->c0[idx][0] = grid->bx[idx];
    grid->c0[idx][1] = grid->by[idx];
    grid->c1[idx][0] = grid->c0[idx][1];
    grid->c1[idx][1] = grid->bz[idx];
  }
  // execute DFT backward plan
  fftw_execute_dft(grid->plan_c0_bw, grid->c0, grid->c0);
  fftw_execute_dft(grid->plan_c1_bw, grid->c1, grid->c1);
  // according to FFTW convention
  // transform forward followed by backword scale up array by nx*ny*nz
  const ham_float inv_grid_size = 1.0 / par->grid_brnd.full_size;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_brnd.nx) i = 0; i < par->grid_brnd.nx; ++i) {
    decltype(par->grid_brnd.nx) i_sym{par->grid_brnd.nx -
                                      i}; // apply Hermitian symmetry
    if (i == 0)
      i_sym = i;
    // it's faster to calculate indeces manually
    const ham_uint idx_lv1{i * par->grid_brnd.ny * par->grid_brnd.nz};
    const ham_uint idx_sym_lv1{i_sym * par->grid_brnd.ny * par->grid_brnd.nz};
    for (decltype(par->grid_brnd.ny) j = 0; j < par->grid_brnd.ny; ++j) {
      decltype(par->grid_brnd.ny) j_sym{par->grid_brnd.ny -
                                        j}; // apply Hermitian symmetry
      if (j == 0)
        j_sym = j;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_brnd.nz};
      const ham_uint idx_sym_lv2{idx_sym_lv1 + j_sym * par->grid_brnd.nz};
      for (decltype(par->grid_brnd.nz) l = 0; l < par->grid_brnd.nz; ++l) {
        decltype(par->grid_brnd.nz) l_sym{par->grid_brnd.nz -
                                          l}; // apply Hermitian symmetry
        if (l == 0)
          l_sym = l;
        const ham_uint idx{idx_lv2 + l};             // q
        const ham_uint idx_sym{idx_sym_lv2 + l_sym}; //-q
        // reconstruct bx,by,bz from c0,c1,c*0,c*1
        // c0(q) = bx(q) + i by(q)
        // c*0(-q) = bx(q) - i by(q)
        // c1(q) = by(q) + i bz(q)
        // c1*1(-q) = by(q) - i bz(q)
        grid->bx[idx] =
            0.5 * (grid->c0[idx][0] + grid->c0[idx_sym][0]) * inv_grid_size;
        grid->by[idx] =
            0.5 * (grid->c1[idx][0] + grid->c1[idx_sym][0]) * inv_grid_size;
        grid->bz[idx] =
            0.5 * (grid->c1[idx_sym][1] + grid->c1[idx][1]) * inv_grid_size;
      }
    }
  }
}
