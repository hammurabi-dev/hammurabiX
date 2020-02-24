#include <cassert>
#include <cmath>
#include <omp.h>

#include <bfield.h>
#include <fftw3.h>
#include <grid.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>
#include <toolkit.h>

void Brnd_mhd::write_grid(const Param *par, const Breg *breg,
                          const Grid_breg *gbreg, Grid_brnd *grid) const {
  // initialize random seed
#ifdef _OPENMP
  gsl_rng **threadvec = new gsl_rng *[omp_get_max_threads()];
  for (ham_int b = 0; b < omp_get_max_threads(); ++b) {
    threadvec[b] = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(threadvec[b], b + toolkit::random_seed(par->brnd_seed));
  }
#else
  gsl_rng *r{gsl_rng_alloc(gsl_rng_taus)};
  gsl_rng_set(r, toolkit::random_seed(par->brnd_seed));
#endif
  // start Fourier space filling
  const ham_float lx{par->grid_brnd.x_max - par->grid_brnd.x_min};
  const ham_float ly{par->grid_brnd.y_max - par->grid_brnd.y_min};
  const ham_float lz{par->grid_brnd.z_max - par->grid_brnd.z_min};
  const Hamvec<3, ham_float> B{breg->read_field(par->observer, par, gbreg)};
  // physical dk^3
  const ham_float dk3{cgs::kpc * cgs::kpc * cgs::kpc / (lx * ly * lz)};
#ifdef _OPENMP
#pragma omp parallel for schedule(static) // DO NOT CHANGE SCHEDULE TYPE
#endif
  for (decltype(par->grid_brnd.nx) i = 0; i < par->grid_brnd.nx; ++i) {
#ifdef _OPENMP
    auto seed_id = threadvec[omp_get_thread_num()];
#else
    auto seed_id = r;
#endif
    Hamvec<3, ham_float> k{cgs::kpc * i / lx, 0, 0};
    if (i >= (par->grid_brnd.nx + 1) / 2)
      k[0] -= cgs::kpc * par->grid_brnd.nx / lx;
    // it's better to calculate indeces manually
    // just for reference, how indeces are calculated
    // const size_t idx
    // {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i,j,l)};
    const ham_uint idx_lv1{i * par->grid_brnd.ny * par->grid_brnd.nz};
    for (decltype(par->grid_brnd.ny) j = 0; j < par->grid_brnd.ny; ++j) {
      k[1] = cgs::kpc * j / ly;
      if (j >= (par->grid_brnd.ny + 1) / 2)
        k[1] -= cgs::kpc * par->grid_brnd.ny / ly;
      const ham_uint idx_lv2{idx_lv1 + j * par->grid_brnd.nz};
      for (decltype(par->grid_brnd.nz) l = 0; l < par->grid_brnd.nz; ++l) {
        // the very 0th term is fixed to zero in allocation
        if (i == 0 and j == 0 and l == 0)
          continue;
        k[2] = cgs::kpc * l / lz;
        if (l >= (par->grid_brnd.nz + 1) / 2)
          k[2] -= cgs::kpc * par->grid_brnd.nz / lz;
        const ham_float ks{k.length()};
        const ham_uint idx{idx_lv2 + l};
        Hamvec<3, ham_float> ep{e_plus(B, k)};
        Hamvec<3, ham_float> em{e_minus(B, k)};
        // since there is no specific rule about how to allocate spectrum power
        // between Re and Im part in b+ and b-
        // we multiply power by two and set Im parts to zero
        // in the end of backward trans, we retrive Re part of bx,by,bz only
        if (ep.lengthsq() > 1e-6) {
          ham_float ang{cosine(B, k)};
          const ham_float Pa{spectrum_a(ks, par) * F_a(par->brnd_mhd.ma, ang) *
                             dk3};
          ham_float Pf{spectrum_f(ks, par) * h_f(par->brnd_mhd.beta, ang) *
                       dk3};
          ham_float Ps{spectrum_s(ks, par) * F_s(par->brnd_mhd.ma, ang) *
                       h_s(par->brnd_mhd.beta, ang) * dk3};
          // b+ is independent from b- in terms of power
          // fast and slow modes are independent
          const ham_float Ap = gsl_ran_ugaussian(seed_id) * std::sqrt(2. * Pa);
          const ham_float Am = gsl_ran_ugaussian(seed_id) * std::sqrt(2. * Pf) +
                               gsl_ran_ugaussian(seed_id) * std::sqrt(2. * Ps);
          const Hamvec<3, ham_float> bkp{ep * Ap};
          const Hamvec<3, ham_float> bkm{em * Am};
          // c0_R = bx_Re - by_Im
          // c0_I = bx_Im + by_Re
          // c1_R = by_Re - bz_Im
          // c1_I = by_Im + bz_Re
          // bx_R = bkp_Re[0]+bkm_Re[0];
          // etc.
          // note that all Im parts are zero, c1_R = c0_I
          grid->c0[idx][0] = bkp[0] + bkm[0];  // bx_Re
          grid->c0[idx][1] = bkp[1] + bkm[1];  // by_Re
          grid->c1[idx][0] = grid->c0[idx][1]; // by_Re
          grid->c1[idx][1] = bkp[2] + bkm[2];  // bz_Re
        } else {
          ham_float Pf{spectrum_f(ks, par) * h_f(par->brnd_mhd.beta, 1) * dk3};
          if (i == 0 and j == 0) {
            ep[0] = k[0];
            em[1] = k[1];
          } else {
            ep[0] = k[0] * k[2] / (ks * std::sqrt(k[0] * k[0] + k[1] * k[1]));
            ep[1] = k[1] * k[2] / (ks * std::sqrt(k[0] * k[0] + k[1] * k[1]));
            ep[2] = -std::sqrt(k[0] * k[0] + k[1] * k[1]) / ks;
            em[0] = -k[1] / std::sqrt(k[0] * k[0] + k[1] * k[1]);
            em[1] = k[0] / std::sqrt(k[0] * k[0] + k[1] * k[1]);
            em[2] = 0.;
          }
          // b+ and b- share power
          const ham_float Af = gsl_ran_ugaussian(seed_id) * std::sqrt(2. * Pf);
          const ham_float share = gsl_rng_uniform(seed_id);
          const Hamvec<3, ham_float> bkp{ep * Af * share};
          const Hamvec<3, ham_float> bkm{em * Af * (1. - share)};
          grid->c0[idx][0] = bkp[0] + bkm[0];
          grid->c0[idx][1] = bkp[1] + bkm[1];
          grid->c1[idx][0] = grid->c0[idx][1];
          grid->c1[idx][1] = bkp[2] + bkm[2];
        }
      } // l
    }   // j
  }     // i
  // free random memory
#ifdef _OPENMP
  for (ham_int b = 0; b < omp_get_max_threads(); ++b)
    gsl_rng_free(threadvec[b]);
  delete[] threadvec;
#else
  gsl_rng_free(r);
#endif
  // execute DFT backward plan
  fftw_execute_dft(grid->plan_c0_bw, grid->c0, grid->c0);
  fftw_execute_dft(grid->plan_c1_bw, grid->c1, grid->c1);
  // now we need to get real parts manually
  // since we start in k-space real fields
  // the results in x-space are complex fields
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (decltype(par->grid_brnd.nx) i = 0; i < par->grid_brnd.nx; ++i) {
    decltype(par->grid_brnd.nx) i_sym{par->grid_brnd.nx -
                                      i}; // apply Hermitian symmetry
    if (i == 0)
      i_sym = i;
    // it's better to calculate indeces manually
    // just for reference, how indeces are calculated
    // const ham_uint idx
    // {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i,j,l)};
    // const ham_uint idx_sym
    // {toolkit::index3d(par->grid_brnd.nx,par->grid_brnd.ny,par->grid_brnd.nz,i_sym,j_sym,l_sym)};
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
        grid->bx[idx] = 0.5 * (grid->c0[idx][0] + grid->c0[idx_sym][0]);
        grid->by[idx] = 0.5 * (grid->c1[idx][0] + grid->c1[idx_sym][0]);
        grid->bz[idx] = 0.5 * (grid->c1[idx_sym][1] + grid->c1[idx][1]);
      }
    }
  }
}

Hamvec<3, ham_float> Brnd_mhd::e_plus(const Hamvec<3, ham_float> &b,
                                      const Hamvec<3, ham_float> &k) const {
  Hamvec<3, ham_float> tmp{(k.versor()).crossprod(b.versor())};
  if (tmp.lengthsq() < 1e-12) {
    return tmp;
  } else {
    return tmp / tmp.length();
  }
}

Hamvec<3, ham_float> Brnd_mhd::e_minus(const Hamvec<3, ham_float> &b,
                                       const Hamvec<3, ham_float> &k) const {
  Hamvec<3, ham_float> tmp{
      ((k.versor()).crossprod(b.versor())).crossprod(k.versor())};
  if (tmp.lengthsq() < 1e-12) {
    return tmp;
  } else {
    return tmp / tmp.length();
  }
}

ham_float Brnd_mhd::h_s(const ham_float &beta, const ham_float &cosa) const {
  const ham_float sqrtD{std::sqrt(dynamo(beta, cosa))};
  const ham_float dpp{1 + sqrtD + 0.5 * beta};
  const ham_float dmm{1 - sqrtD - 0.5 * beta};
  const ham_float dmp{1 - sqrtD + 0.5 * beta};
  const ham_float demon{
      1 / (dmp * (cosa * cosa + dpp * dpp * (1 - cosa * cosa) / (dmm * dmm)))};
  if (std::isfinite(demon))
    return 2. * cosa * cosa * demon * (demon > 0);
  return 0;
}

ham_float Brnd_mhd::h_f(const ham_float &beta, const ham_float &cosa) const {
  if (cosa == 0)
    return 0;
  const ham_float sqrtD{std::sqrt(dynamo(beta, cosa))};
  const ham_float dpp{1 + sqrtD + 0.5 * beta};
  const ham_float dpm{1 + sqrtD - 0.5 * beta};
  const ham_float dmp{(1 - sqrtD + 0.5 * beta)};
  const ham_float demon{
      1 / (dpp * (cosa * cosa + dmp * dmp * (1 - cosa * cosa) / (dpm * dpm)))};
  if (std::isfinite(demon))
    return 2. * cosa * cosa * demon;
  return 0;
}

ham_float Brnd_mhd::spectrum_a(const ham_float &k, const Param *par) const {
  // units fixing, wave vector in 1/kpc units
  const ham_float unit = 0.25 / (cgs::pi * k * k);
  // power laws
  const ham_float band1{ham_float(k < par->brnd_mhd.k1)};
  const ham_float band2{ham_float(k > par->brnd_mhd.k1) *
                        ham_float(k < par->brnd_mhd.k0)};
  const ham_float band3{ham_float(k > par->brnd_mhd.k0)};
  const ham_float P =
      band1 * std::pow(par->brnd_mhd.k0 / par->brnd_mhd.k1, par->brnd_mhd.a1) *
          std::pow(k / par->brnd_mhd.k1, 6.0) +
      band2 / std::pow(k / par->brnd_mhd.k0, par->brnd_mhd.a1) +
      band3 / std::pow(k / par->brnd_mhd.k0, par->brnd_mhd.aa0);
  return P * par->brnd_mhd.pa0 * unit;
}

ham_float Brnd_mhd::spectrum_f(const ham_float &k, const Param *par) const {
  // units fixing, wave vector in 1/kpc units
  const ham_float unit = 0.25 / (cgs::pi * k * k);
  // power laws
  const ham_float band1{ham_float(k < par->brnd_mhd.k1)};
  const ham_float band2{ham_float(k > par->brnd_mhd.k1) *
                        ham_float(k < par->brnd_mhd.k0)};
  const ham_float band3{ham_float(k > par->brnd_mhd.k0)};
  const ham_float P =
      band1 * std::pow(par->brnd_mhd.k0 / par->brnd_mhd.k1, par->brnd_mhd.a1) *
          std::pow(k / par->brnd_mhd.k1, 6.0) +
      band2 / std::pow(k / par->brnd_mhd.k0, par->brnd_mhd.a1) +
      band3 / std::pow(k / par->brnd_mhd.k0, par->brnd_mhd.af0);
  return P * par->brnd_mhd.pf0 * unit;
}

ham_float Brnd_mhd::spectrum_s(const ham_float &k, const Param *par) const {
  // units fixing, wave vector in 1/kpc units
  const ham_float unit = 0.25 / (cgs::pi * k * k);
  // power laws
  const ham_float band1{ham_float(k < par->brnd_mhd.k1)};
  const ham_float band2{ham_float(k > par->brnd_mhd.k1) *
                        ham_float(k < par->brnd_mhd.k0)};
  const ham_float band3{ham_float(k > par->brnd_mhd.k0)};
  const ham_float P =
      band1 * std::pow(par->brnd_mhd.k0 / par->brnd_mhd.k1, par->brnd_mhd.a1) *
          std::pow(k / par->brnd_mhd.k1, 6.0) +
      band2 / std::pow(k / par->brnd_mhd.k0, par->brnd_mhd.a1) +
      band3 / std::pow(k / par->brnd_mhd.k0, par->brnd_mhd.as0);
  return P * par->brnd_mhd.ps0 * unit;
}
