#include <cassert>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

#include <fitshandle.h>
#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <pointing.h>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_synchrotron.h>

#include <bfield.h>
#include <cgs_units.h>
#include <crefield.h>
#include <grid.h>
#include <hamvec.h>
#include <integrator.h>
#include <param.h>
#include <tefield.h>
#include <hampixma.h>

void Integrator::write_grid(const Breg *breg, const Brnd *brnd,
                            const TEreg *tereg, const TErnd *ternd,
                            const CREfield *cre, const Grid_breg *gbreg,
                            const Grid_brnd *gbrnd, const Grid_tereg *gtereg,
                            const Grid_ternd *gternd, const Grid_cre *gcre,
                            Grid_obs *gobs, const Param *par) const {
  auto mask = std::make_unique<hampixma>();
  if (par->grid_obs.do_mask) {
    mask->import(par);
  }
  if (par->grid_obs.do_dm) {
    gobs->dm_map->fill(0.);
  }
  if (par->grid_obs.do_sync.back()) {
    gobs->is_map->fill(0.);
    gobs->qs_map->fill(0.);
    gobs->us_map->fill(0.);
  }
  if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
    gobs->fd_map->fill(0.);
  }
  auto shell_ref = std::make_unique<struct_shell>();
  // loop through shells
  for (decltype(par->grid_obs.total_shell) current_shell = 0;
       current_shell != par->grid_obs.total_shell; ++current_shell) {
    // get current shell nside & npix
    const std::size_t current_nside{
        par->grid_obs.nside_shell[current_shell]};
    const std::size_t current_npix{12 * current_nside * current_nside};
    // get current mask
    if (par->grid_obs.do_mask) {
      mask->duplicate(current_nside);
    }
    // prepare temporary maps for current shell
    if (par->grid_obs.do_dm) {
      gobs->tmp_dm_map->SetNside(current_nside, RING);
      gobs->tmp_dm_map->fill(0.);
    }
    if (par->grid_obs.do_sync.back()) {
      gobs->tmp_is_map->SetNside(current_nside, RING);
      gobs->tmp_qs_map->SetNside(current_nside, RING);
      gobs->tmp_us_map->SetNside(current_nside, RING);
      gobs->tmp_is_map->fill(0.);
      gobs->tmp_qs_map->fill(0.);
      gobs->tmp_us_map->fill(0.);
    }
    if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
      gobs->tmp_fd_map->SetNside(current_nside, RING);
      gobs->tmp_fd_map->fill(0.);
    }
    // setting for radial_integration
    // call auxiliary function assemble_shell_ref
    assemble_shell_ref(shell_ref.get(), par, current_shell);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t ipix = 0; ipix < current_npix; ++ipix) {
      auto observables = std::make_unique<struct_observables>();
      observables->is = 0.;
      observables->qs = 0.;
      observables->us = 0.;
      observables->dm = 0.;
      // check pixel masking
      if ((not par->grid_obs.do_mask) or mask->info(current_nside,ipix) == 1.0) {
        // remember to complete logic for ptg assignment!
        // and for caching Faraday depth and/or optical depth
        // make serious tests after changing this part!
        pointing ptg;
        if (par->grid_obs.do_dm) {
          ptg = gobs->tmp_dm_map->pix2ang(ipix);
        }
        if (par->grid_obs.do_fd) {
          ptg = gobs->tmp_fd_map->pix2ang(ipix);
          // cache Faraday rotation from inner shells
          observables->fd = gobs->fd_map->interpolated_value(ptg);
        } else if (par->grid_obs.do_sync.back()) {
          ptg = gobs->tmp_is_map->pix2ang(ipix);
          // cache Faraday rotation from inner shells
          observables->fd = gobs->fd_map->interpolated_value(ptg);
        }
        // core function!
        radial_integration(shell_ref.get(), ptg, observables.get(), breg, brnd,
                           tereg, ternd, cre, gbreg, gbrnd, gtereg, gternd, gcre,
                           par);
      }
      // collect from pixels
      if (par->grid_obs.do_dm) {
        (*gobs->tmp_dm_map)[ipix] = observables->dm;
      }
      if (par->grid_obs.do_sync.back()) {
        (*gobs->tmp_is_map)[ipix] =
            temp_convert(observables->is, par->grid_obs.sim_sync_freq.back());
        (*gobs->tmp_qs_map)[ipix] =
            temp_convert(observables->qs, par->grid_obs.sim_sync_freq.back());
        (*gobs->tmp_us_map)[ipix] =
            temp_convert(observables->us, par->grid_obs.sim_sync_freq.back());
      }
      if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
        (*gobs->tmp_fd_map)[ipix] = observables->fd;
      }
    }
    // accumulating new shell map to sim map
    if (par->grid_obs.do_dm) {
      std::size_t npix_dm = gobs->dm_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t ipix = 0; ipix < npix_dm; ++ipix) {
        pointing ptg{gobs->dm_map->pix2ang(ipix)};
        (*gobs->dm_map)[ipix] += gobs->tmp_dm_map->interpolated_value(ptg);
      }
    }
    if (par->grid_obs.do_sync.back()) {
      std::size_t npix_sync = gobs->is_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t ipix = 0; ipix < npix_sync; ++ipix) {
        pointing ptg = {gobs->is_map->pix2ang(ipix)};
        (*gobs->is_map)[ipix] += gobs->tmp_is_map->interpolated_value(ptg);
        (*gobs->qs_map)[ipix] += gobs->tmp_qs_map->interpolated_value(ptg);
        (*gobs->us_map)[ipix] += gobs->tmp_us_map->interpolated_value(ptg);
      }
    }
    if (par->grid_obs.do_fd) {
      std::size_t npix_fd = gobs->fd_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t ipix = 0; ipix < npix_fd; ++ipix) {
        pointing ptg = {gobs->fd_map->pix2ang(ipix)};
        (*gobs->fd_map)[ipix] += gobs->tmp_fd_map->interpolated_value(ptg);
      }
    } // end shell accumulation
  }   // end shell iteration
}

void Integrator::radial_integration(
    const struct_shell *shell_ref, const pointing &ptg_in,
    struct_observables *pixobs, const Breg *breg, const Brnd *brnd,
    const TEreg *tereg, const TErnd *ternd, const CREfield *cre,
    const Grid_breg *gbreg, const Grid_brnd *gbrnd, const Grid_tereg *gtereg,
    const Grid_ternd *gternd, const Grid_cre *gcre, const Param *par) const {
  // pass in fd, zero others
  double inner_shells_fd{0.};
  if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
    inner_shells_fd = pixobs->fd;
  }
  pixobs->dm = 0.;
  pixobs->fd = 0.;
  pixobs->is = 0.;
  pixobs->us = 0.;
  pixobs->qs = 0.;
  // angular position
  const double THE{ptg_in.theta};
  const double PHI{ptg_in.phi};
  double lambda_square = 0., i2bt_sync = 0., fd_forefactor = 0.;
  if (par->grid_obs.do_sync.back()) {
    // for calculating synchrotron emission
    lambda_square = (cgs_c_light / par->grid_obs.sim_sync_freq.back()) *
                    (cgs_c_light / par->grid_obs.sim_sync_freq.back());
    // convert sync intensity(freq) to brightness temperature, Rayleigh-Jeans
    // law
    i2bt_sync = cgs_c_light * cgs_c_light /
                (2. * cgs_kB * par->grid_obs.sim_sync_freq.back() *
                 par->grid_obs.sim_sync_freq.back());
  }
  if (par->grid_obs.do_fd) {
    fd_forefactor =
        -(cgs_qe * cgs_qe * cgs_qe) / (2. * cgs_pi * cgs_mec2 * cgs_mec2);
  }
  // radial accumulation
  for (decltype(shell_ref->step) looper = 0; looper < shell_ref->step;
       ++looper) {
    // ec and gc position
    hamvec<3, double> oc_pos{los_versor(THE, PHI) * shell_ref->dist[looper]};
    hamvec<3, double> pos{oc_pos + par->observer};
    // check LoS depth limit
    if (check_simulation_lower_limit(pos.length(), par->grid_obs.gc_r_min)) {
      continue;
    }
    if (check_simulation_upper_limit(pos.length(), par->grid_obs.gc_r_max)) {
      continue;
    }
    if (check_simulation_lower_limit(pos[2], par->grid_obs.gc_z_min)) {
      continue;
    }
    if (check_simulation_upper_limit(pos[2], par->grid_obs.gc_z_max)) {
      continue;
    }
    // regular magnetic field
    hamvec<3, double> B_vec{breg->read_field(pos, par, gbreg)};
    // add random magnetic field
    B_vec += brnd->read_field(pos, par, gbrnd);
    const double B_par{los_parproj(B_vec, THE, PHI)};
    assert(std::isfinite(B_par));
    // be aware of un-resolved random B_per in calculating emissivity
    const double B_per{los_perproj(B_vec, THE, PHI)};
    assert(std::isfinite(B_per));
    // thermal electron field
    double te{tereg->read_field(pos, par, gtereg)};
    // add random thermal electron field
    te += ternd->read_field(pos, par, gternd);
    // to avoid negative value
    te *= double(te > 0.);
    assert(std::isfinite(te));
    // dispersion measure
    if (par->grid_obs.do_dm) {
      pixobs->dm += te * shell_ref->delta_d;
    }
    // Faraday depth
    if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
      pixobs->fd += te * B_par * fd_forefactor * shell_ref->delta_d;
    }
    // Synchrotron emission
    if (par->grid_obs.do_sync.back()) {
      const double Jtot{sync_emissivity_t(pos, par, cre, gcre, B_per) *
                        shell_ref->delta_d * i2bt_sync};
      // J_pol receives no contribution from unresolved random field
      const double Jpol{sync_emissivity_p(pos, par, cre, gcre, B_per) *
                        shell_ref->delta_d * i2bt_sync};
      assert(Jtot < 1e30 and Jpol < 1e30 and Jtot >= 0 and Jpol >= 0);
      pixobs->is += Jtot;
      // intrinsic polarization angle, following IAU definition
      const double qui{(inner_shells_fd + pixobs->fd) * lambda_square +
                       sync_ipa(B_vec, THE, PHI)};
      assert(std::isfinite(qui));
      pixobs->qs += std::cos(2. * qui) * Jpol;
      pixobs->us += std::sin(2. * qui) * Jpol;
    }
  }
}

// assembling shell_ref structure
void Integrator::assemble_shell_ref(struct_shell *target, const Param *par,
                                    const std::size_t &shell_num) const {
  target->shell_num = shell_num;
  target->d_start = par->grid_obs.radii_shell[shell_num];
  target->d_stop = par->grid_obs.radii_shell[shell_num+1];
  target->delta_d = par->grid_obs.oc_r_res;
  target->step = floor(
      (target->d_stop / target->delta_d - target->d_start / target->delta_d));
  // get rid of error in the previous step
  target->delta_d = (target->d_stop - target->d_start) / (target->step);
  target->dist.clear(); // clean cache
  for (std::size_t i = 0; i < target->step; ++i) {
    target->dist.push_back(target->d_start + (i + 0.5) * target->delta_d);
  }
#ifdef VERBOSE
  std::cout << "shell reference: " << std::endl
            << "shell No. " << target->shell_num << std::endl
            << "resolution " << target->delta_d / cgs_kpc << " kpc" << std::endl
            << "d_start " << target->d_start / cgs_kpc << " kpc" << std::endl
            << "d_stop " << target->d_stop / cgs_kpc << " kpc" << std::endl
            << "steps " << target->step << std::endl;
#endif
}

// calculate synchrotron emission intrinsic polarization angle
double Integrator::sync_ipa(const hamvec<3, double> &input,
                            const double &the_ec, const double &phi_ec) const {
  const hamvec<3, double> sph_unit_v_the(
      cos(the_ec) * cos(phi_ec), cos(the_ec) * sin(phi_ec), -sin(the_ec));
  const hamvec<3, double> sph_unit_v_phi(-sin(phi_ec), cos(phi_ec), 0.);
  // IAU convention
  return atan2(-sph_unit_v_the.dotprod(input), -sph_unit_v_phi.dotprod(input));
}

// Carteisan unit vector of given LoS direction
hamvec<3, double> Integrator::los_versor(const double &the_los,
                                         const double &phi_los) const {
  return hamvec<3, double>{cos(phi_los) * sin(the_los),
                           sin(phi_los) * sin(the_los), cos(the_los)};
}

// perpendicular projection of a vector wrt LoS direction
double Integrator::los_perproj(const hamvec<3, double> &input,
                               const double &the_los,
                               const double &phi_los) const {
  return (input.crossprod(los_versor(the_los, phi_los))).length();
}

// (signed) parallel projection of a vector wrt LoS direction
double Integrator::los_parproj(const hamvec<3, double> &input,
                               const double &the_los,
                               const double &phi_los) const {
  return input.dotprod(los_versor(the_los, phi_los));
}

// converting brightness temp into thermal temp with T_0 = 2.725K
double Integrator::temp_convert(const double &temp_br,
                                const double &freq) const {
  const double p{cgs_h_planck * freq / (cgs_kB * 2.725)};
  return temp_br * (exp(p) - 1.) * (exp(p) - 1.) / (p * p * exp(p));
}

// cre synchrotron J_tot(\nu)
double Integrator::sync_emissivity_t(const hamvec<3, double> &pos,
                                     const Param *par, const CREfield *cre,
                                     const Grid_cre *grid,
                                     const double &Bper) const {
  double J{0};
  // calculate from grid
  if (par->grid_cre.read_permission) {
    // allocate energy grid
    std::unique_ptr<double[]> KE = std::make_unique<double[]>(par->grid_cre.nE);
    // we need F(x[E]) in spectral integration
    std::unique_ptr<double[]> x = std::make_unique<double[]>(par->grid_cre.nE);
    std::unique_ptr<double[]> beta =
        std::make_unique<double[]>(par->grid_cre.nE);
    // consts used for converting E to x, using cgs units
    const double x_fact{4. * cgs_pi * par->grid_obs.sim_sync_freq.back() *
                        cgs_mec * cgs_mec2 * cgs_mec2 /
                        (3. * cgs_qe * std::fabs(Bper))};
    // KE, x, beta arrays in cgs units
    for (decltype(par->grid_cre.nE) i = 0; i != par->grid_cre.nE; ++i) {
      KE[i] = par->grid_cre.E_min * std::exp(i * par->grid_cre.E_fact);
      x[i] = x_fact / (KE[i] * KE[i]);
      beta[i] = std::sqrt(1 - cgs_mec2 / KE[i]);
    }
    // spectral integral
    for (decltype(par->grid_cre.nE) i = 0; i != par->grid_cre.nE - 1; ++i) {
      const double xv{0.5 * (x[i + 1] + x[i])};
      // avoid underflow in gsl functions
      if (xv > 100) {
        continue;
      }
      const double dE{std::fabs(KE[i + 1] - KE[i])};
      // we put beta here, midpoint rule
      const double flux{
          0.5 * (cre->read_grid_num(pos, i + 1, par, grid) / beta[i + 1] +
                 cre->read_grid_num(pos, i, par, grid) / beta[i])};
      assert(flux >= 0);
      J += gsl_sf_synchrotron_1(xv) * flux * dE;
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and
    // density n ref: "Cosmic rays n' particle physics", A3
    const double fore_factor{
        1.73205081 * cgs_qe * cgs_qe * cgs_qe * std::fabs(Bper) /
        (cgs_mec2 * cgs_c_light * cgs_GeV * cgs_m * cgs_m * cgs_sec)};
    J *= fore_factor;
  }
  // calculate from N(\gamma) with local constant spectral index
  else {
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    const double index{cre->flux_idx(pos, par)};
    // coefficients which do not attend integration
    const double norm{cre->flux_norm(pos, par) * 1.73205081 *
                      (cgs_qe * cgs_qe * cgs_qe) * std::fabs(Bper) /
                      (4. * cgs_pi * cgs_mec2 * (1. - index))};
    // synchrotron integration
    const double A{2. * cgs_pi * par->grid_obs.sim_sync_freq.back() * cgs_mec /
                   (3. * cgs_qe * std::fabs(Bper))};
    J = norm * (std::pow(A, 0.5 * (index + 1)) *
                gsl_sf_gamma(-0.25 * index + 19. / 12.) *
                gsl_sf_gamma(-0.25 * index - 1. / 12.));
  }
  return J;
}

// cre synchrotron J_pol(\nu)
double Integrator::sync_emissivity_p(const hamvec<3, double> &pos,
                                     const Param *par, const CREfield *cre,
                                     const Grid_cre *grid,
                                     const double &Bper) const {
  double J{0};
  // calculate from grid
  if (par->grid_cre.read_permission) {
    // allocate energy grid
    std::unique_ptr<double[]> KE = std::make_unique<double[]>(par->grid_cre.nE);
    // we need F(x[E]) in spectral integration
    std::unique_ptr<double[]> x = std::make_unique<double[]>(par->grid_cre.nE);
    std::unique_ptr<double[]> beta =
        std::make_unique<double[]>(par->grid_cre.nE);
    // consts used for converting E to x, using cgs units
    const double x_fact{(2. * cgs_mec * cgs_mec2 * cgs_mec2 * 2. * cgs_pi *
                         par->grid_obs.sim_sync_freq.back()) /
                        (3. * cgs_qe * std::fabs(Bper))};
    // KE, x, beta arrays in cgs units
    for (decltype(par->grid_cre.nE) i = 0; i != par->grid_cre.nE; ++i) {
      KE[i] = par->grid_cre.E_min * std::exp(i * par->grid_cre.E_fact);
      x[i] = x_fact / (KE[i] * KE[i]);
      beta[i] = std::sqrt(1 - cgs_mec2 / KE[i]);
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and
    // density n ref: "Cosmic rays n' particle physics", A3
    const double fore_factor{
        1.73205081 * cgs_qe * cgs_qe * cgs_qe * std::fabs(Bper) /
        (cgs_mec2 * cgs_c_light * cgs_GeV * cgs_m * cgs_m * cgs_sec)};
    // spectral integral
    for (decltype(par->grid_cre.nE) i = 0; i != par->grid_cre.nE - 1; ++i) {
      const double xv{0.5 * (x[i + 1] + x[i])};
      // avoid underflow in gsl functions
      if (xv > 100) {
        continue;
      }
      const double dE{std::fabs(KE[i + 1] - KE[i])};
      // we put beta here, midpoint rule
      const double flux{
          0.5 * (cre->read_grid_num(pos, i + 1, par, grid) / beta[i + 1] +
                 cre->read_grid_num(pos, i, par, grid) / beta[i])};
      assert(flux >= 0);
      J += gsl_sf_synchrotron_2(xv) * flux * dE;
    }
    J *= fore_factor;
  }
  // calculate from N(\gamma) with local constant spectral index
  else {
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    const double index{cre->flux_idx(pos, par)};
    // coefficients which do not attend integration
    const double norm{cre->flux_norm(pos, par) * 1.73205081 *
                      (cgs_qe * cgs_qe * cgs_qe) * std::fabs(Bper) /
                      (16. * cgs_pi * cgs_mec2)};
    // synchrotron integration
    const double A{2. * cgs_pi * par->grid_obs.sim_sync_freq.back() * cgs_mec /
                   (3. * cgs_qe * std::fabs(Bper))};
    J = norm * (std::pow(A, 0.5 * (index + 1)) *
                gsl_sf_gamma(-0.25 * index + 7. / 12.) *
                gsl_sf_gamma(-0.25 * index - 1. / 12.));
  }
  // the last 4pi comes from solid-angle integration/deviation,
  // check eq(6.16) in Ribiki-Lightman's where Power is defined,
  // we need isotropic power which means we need a 1/4pi factor!
  return J;
}
