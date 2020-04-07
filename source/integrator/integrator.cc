#include <bfield.h>
#include <cassert>
#include <cmath>
#include <crefield.h>
#include <grid.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <hamdis.h>
#include <hamp.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <integrator.h>
#include <iostream>
#include <omp.h>
#include <param.h>
#include <tefield.h>
#include <timer.h>
#include <vector>

void Integrator::write_grid(const Param *par) const {
  auto shell_ref = std::make_unique<struct_shell>();
  // loop through shells
  for (decltype(par->grid_obs.total_shell) current_shell = 0;
       current_shell != par->grid_obs.total_shell; ++current_shell) {
    // get current shell nside & npix
    const ham_uint current_nside{par->grid_obs.nside_shell[current_shell]};
    const ham_uint current_npix{12 * current_nside * current_nside};
    // get current mask
    if (par->grid_obs.do_mask) {
      grids.obs->mask_map->duplicate(current_nside);
    }
    // prepare temporary maps for current shell
    if (par->grid_obs.do_dm) {
      grids.obs->tmp_dm_map->reset(current_nside);
    }
    if (par->grid_obs.do_sync.back()) {
      grids.obs->tmp_is_map->reset(current_nside);
      grids.obs->tmp_qs_map->reset(current_nside);
      grids.obs->tmp_us_map->reset(current_nside);
    }
    if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
      grids.obs->tmp_fd_map->reset(current_nside);
    }
    // setting for radial_integration
    // call auxiliary function assemble_shell_ref
    assemble_shell_ref(shell_ref.get(), par, current_shell);
#ifndef NTIMING
    auto tmr = std::make_unique<Timer>();
    tmr->start("pixLoS");
#endif
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (ham_uint ipix = 0; ipix < current_npix; ++ipix) {
      auto observables = std::make_unique<struct_observables>();
      observables->is = 0.;
      observables->qs = 0.;
      observables->us = 0.;
      observables->dm = 0.;
      // check pixel masking
      if ((not par->grid_obs.do_mask) or
          grids.obs->mask_map->data(current_nside, ipix) == 1.0) {
        // remember to complete logic for ptg assignment!
        // and for caching Faraday depth and/or optical depth
        // make serious tests after changing this part!
        Hamp ptg;
        if (par->grid_obs.do_dm) {
          ptg = grids.obs->tmp_dm_map->pointing(ipix);
        }
        if (par->grid_obs.do_fd) {
          ptg = grids.obs->tmp_fd_map->pointing(ipix);
          // cache Faraday rotation from inner shells
          observables->fd = grids.obs->fd_map->interpolate(ptg);
        } else if (par->grid_obs.do_sync.back()) {
          ptg = grids.obs->tmp_is_map->pointing(ipix);
          // cache Faraday rotation from inner shells
          observables->fd = grids.obs->fd_map->interpolate(ptg);
        }
        // core function!
        radial_integration(shell_ref.get(), ptg, observables.get(), par);
      }
      // collect from pixels
      if (par->grid_obs.do_dm) {
        grids.obs->tmp_dm_map->data(ipix, observables->dm);
      }
      if (par->grid_obs.do_sync.back()) {
        grids.obs->tmp_is_map->data(
            ipix,
            temp_convert(observables->is, par->grid_obs.sim_sync_freq.back()));
        grids.obs->tmp_qs_map->data(
            ipix,
            temp_convert(observables->qs, par->grid_obs.sim_sync_freq.back()));
        grids.obs->tmp_us_map->data(
            ipix,
            temp_convert(observables->us, par->grid_obs.sim_sync_freq.back()));
      }
      if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
        grids.obs->tmp_fd_map->data(ipix, observables->fd);
      }
    }
#ifndef NTIMING
    tmr->stop("pixLoS");
    tmr->print();
#endif
    // accumulating new shell map to sim map
    if (par->grid_obs.do_dm) {
      grids.obs->dm_map->accumulate(*(grids.obs->tmp_dm_map));
    }
    if (par->grid_obs.do_sync.back()) {
      grids.obs->is_map->accumulate(*(grids.obs->tmp_is_map));
      grids.obs->qs_map->accumulate(*(grids.obs->tmp_qs_map));
      grids.obs->us_map->accumulate(*(grids.obs->tmp_us_map));
    }
    if (par->grid_obs.do_fd) {
      grids.obs->fd_map->accumulate(*(grids.obs->tmp_fd_map));
    } // end shell accumulation
  }   // end shell iteration
}

void Integrator::radial_integration(const struct_shell *shell_ref,
                                    const Hamp &ptg_in,
                                    struct_observables *pixobs,
                                    const Param *par) const {
  // pass in fd, zero others
  ham_float inner_shells_fd{0.};
  if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
    inner_shells_fd = pixobs->fd;
  }
  pixobs->dm = 0.;
  pixobs->fd = 0.;
  pixobs->is = 0.;
  pixobs->us = 0.;
  pixobs->qs = 0.;
  // angular position
  const ham_float theta{ptg_in.theta()};
  const ham_float phi{ptg_in.phi()};
  // pre-calculated LoS versor
  const Hamvec<3, ham_float> los_direction{los_versor(theta, phi)};
  // emission related constants
  ham_float lambda_square = 0., i2bt_sync = 0., fd_forefactor = 0.;
  if (par->grid_obs.do_sync.back()) {
    // for calculating synchrotron emission
    lambda_square = (cgs::c_light / par->grid_obs.sim_sync_freq.back()) *
                    (cgs::c_light / par->grid_obs.sim_sync_freq.back());
    // convert sync intensity(freq) to brightness temperature, Rayleigh-Jeans
    // law
    i2bt_sync = cgs::c_light * cgs::c_light /
                (2. * cgs::kB * par->grid_obs.sim_sync_freq.back() *
                 par->grid_obs.sim_sync_freq.back());
  }
  if (par->grid_obs.do_fd) {
    fd_forefactor =
        -(cgs::qe * cgs::qe * cgs::qe) / (2. * cgs::pi * cgs::mec2 * cgs::mec2);
  }
  // radial accumulation
  for (decltype(shell_ref->step) looper = 0; looper < shell_ref->step;
       ++looper) {
    // ec and gc position
    Hamvec<3, ham_float> oc_pos{los_direction * shell_ref->dist[looper]};
    Hamvec<3, ham_float> pos{oc_pos + par->observer};
    // check LoS depth limit
    const bool gc_r_min_flag{
        check_simulation_lower_limit(pos.length(), par->grid_obs.gc_r_min)};
    const bool gc_r_max_flag{
        check_simulation_upper_limit(pos.length(), par->grid_obs.gc_r_max)};
    const bool gc_z_min_flag{
        check_simulation_lower_limit(pos[2], par->grid_obs.gc_z_min)};
    const bool gc_z_max_flag{
        check_simulation_upper_limit(pos[2], par->grid_obs.gc_z_max)};
    if (gc_r_min_flag or gc_r_max_flag or gc_z_min_flag or gc_z_max_flag)
      continue;
    // regular magnetic field
    Hamvec<3, ham_float> b_vec{fields.b->read_field(pos, par, grids.b)};
    const ham_float b_par{los_parproj(b_vec, los_direction)};
    assert(std::isfinite(b_par));
    // be aware of un-resolved random b_per in calculating emissivity
    const ham_float b_per{los_perproj(b_vec, los_direction)};
    assert(std::isfinite(b_per));
    // thermal electron field
    ham_float te{fields.te->read_field(pos, par, grids.te)};
    // to avoid negative value
    te *= ham_float(te > 0.);
    assert(std::isfinite(te));
    // dispersion measure
    if (par->grid_obs.do_dm) {
      pixobs->dm += te * shell_ref->delta_d;
    }
    // Faraday depth
    if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
      pixobs->fd += te * b_par * fd_forefactor * shell_ref->delta_d;
    }
    // Synchrotron emission
    if (par->grid_obs.do_sync.back()) {
      const ham_float Jtot{sync_emissivity_t(pos, par, b_per) *
                           shell_ref->delta_d * i2bt_sync};
      // J_pol receives no contribution from unresolved random field
      const ham_float Jpol{sync_emissivity_p(pos, par, b_per) *
                           shell_ref->delta_d * i2bt_sync};
      assert(Jtot < 1e30 and Jpol < 1e30 and Jtot >= 0 and Jpol >= 0);
      pixobs->is += Jtot;
      // intrinsic polarization angle, following IAU definition
      const ham_float qui{(inner_shells_fd + pixobs->fd) * lambda_square +
                          sync_ipa(b_vec, theta, phi)};
      assert(std::isfinite(qui));
      pixobs->qs += std::cos(2. * qui) * Jpol;
      pixobs->us += std::sin(2. * qui) * Jpol;
    }
  }
}

// assembling shell_ref structure
void Integrator::assemble_shell_ref(struct_shell *target, const Param *par,
                                    const ham_uint &shell_num) const {
  target->shell_num = shell_num;
  target->d_start = par->grid_obs.radii_shell[shell_num];
  target->d_stop = par->grid_obs.radii_shell[shell_num + 1];
  target->delta_d = par->grid_obs.oc_r_res;
  target->step = floor(
      (target->d_stop / target->delta_d - target->d_start / target->delta_d));
  // get rid of error in the previous step
  target->delta_d = (target->d_stop - target->d_start) / (target->step);
  target->dist.clear(); // clean cache
  for (ham_uint i = 0; i < target->step; ++i) {
    target->dist.push_back(target->d_start + (i + 0.5) * target->delta_d);
  }
#ifdef VERBOSE
  std::cout << "shell reference: " << std::endl
            << "shell No. " << target->shell_num << std::endl
            << "resolution " << target->delta_d / cgs::kpc << " kpc"
            << std::endl
            << "d_start " << target->d_start / cgs::kpc << " kpc" << std::endl
            << "d_stop " << target->d_stop / cgs::kpc << " kpc" << std::endl
            << "steps " << target->step << std::endl;
#endif
}

// calculate synchrotron emission intrinsic polarization angle
ham_float Integrator::sync_ipa(const Hamvec<3, ham_float> &input,
                               const ham_float &the_ec,
                               const ham_float &phi_ec) const {
  const Hamvec<3, ham_float> sph_unit_v_the(
      cos(the_ec) * cos(phi_ec), cos(the_ec) * sin(phi_ec), -sin(the_ec));
  const Hamvec<3, ham_float> sph_unit_v_phi(-sin(phi_ec), cos(phi_ec), 0.);
  // IAU convention
  return atan2(-sph_unit_v_the.dotprod(input), -sph_unit_v_phi.dotprod(input));
}

// Carteisan unit vector of given LoS direction
Hamvec<3, ham_float> Integrator::los_versor(const ham_float &the_los,
                                            const ham_float &phi_los) const {
  return Hamvec<3, ham_float>{cos(phi_los) * sin(the_los),
                              sin(phi_los) * sin(the_los), cos(the_los)};
}

// converting brightness temp into thermal temp with T_0 = 2.725K
ham_float Integrator::temp_convert(const ham_float &temp_br,
                                   const ham_float &freq) const {
  const ham_float p{cgs::h_planck * freq / (cgs::kB * 2.725)};
  return temp_br * (exp(p) - 1.) * (exp(p) - 1.) / (p * p * exp(p));
}

// cre synchrotron J_tot(\nu)
ham_float Integrator::sync_emissivity_t(const Hamvec<3, ham_float> &pos,
                                        const Param *par,
                                        const ham_float &b_per) const {
  ham_float J{0};
  // calculate from grid
  if (par->grid_cre.read_permission or par->grid_cre.build_permission) {
    // allocate energy grid
    std::unique_ptr<ham_float[]> ke =
        std::make_unique<ham_float[]>(par->grid_cre.ne);
    // we need F(x[E]) in spectral integration
    std::unique_ptr<ham_float[]> x =
        std::make_unique<ham_float[]>(par->grid_cre.ne);
    std::unique_ptr<ham_float[]> beta =
        std::make_unique<ham_float[]>(par->grid_cre.ne);
    // consts used for converting E to x, using cgs units
    const ham_float x_fact{4. * cgs::pi * par->grid_obs.sim_sync_freq.back() *
                           cgs::mec * cgs::mec2 * cgs::mec2 /
                           (3. * cgs::qe * b_per)};
    // ke, x, beta arrays in cgs units
    for (decltype(par->grid_cre.ne) i = 0; i != par->grid_cre.ne; ++i) {
      ke[i] = par->grid_cre.e_min * std::exp(i * par->grid_cre.e_fact);
      x[i] = x_fact / (ke[i] * ke[i]);
      beta[i] = std::sqrt(1. - cgs::mec2 / ke[i]);
    }
    // spectral integral
    for (decltype(par->grid_cre.ne) i = 0; i != par->grid_cre.ne - 1; ++i) {
      const ham_float xv{0.5 * (x[i + 1] + x[i])};
      // avoid underflow in gsl functions
      if (xv > 100.) {
        continue;
      }
      const ham_float de{std::fabs(ke[i + 1] - ke[i])};
      // we put beta here, midpoint rule
      const ham_float flux{
          0.5 *
          (fields.cre->read_field(pos, i + 1, par, grids.cre) / beta[i + 1] +
           fields.cre->read_field(pos, i, par, grids.cre) / beta[i])};
      assert(flux >= 0);
      J += gsl_sf_synchrotron_1(xv) * flux * de;
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and
    // density n ref: "Cosmic rays n' particle physics", A3
    const ham_float fore_factor{
        cgs::sqrt3 * cgs::qe * cgs::qe * cgs::qe * b_per /
        (cgs::mec2 * cgs::c_light * cgs::GeV * cgs::m * cgs::m * cgs::sec)};
    J *= fore_factor;
  }
  // calculate from N(\gamma) with local constant spectral index
  else {
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    const ham_float index{fields.cre->flux_idx(pos, par)};
    // coefficients which do not attend integration
    const ham_float norm{fields.cre->flux_norm(pos, par) * cgs::sqrt3 *
                         (cgs::qe * cgs::qe * cgs::qe) * b_per /
                         (4. * cgs::pi * cgs::mec2 * (1. - index))};
    // synchrotron integration
    const ham_float A{2. * cgs::pi * par->grid_obs.sim_sync_freq.back() *
                      cgs::mec / (3. * cgs::qe * b_per)};
    J = norm * (std::pow(A, 0.5 * (index + 1.)) *
                gsl_sf_gamma(-0.25 * index + 19. / 12.) *
                gsl_sf_gamma(-0.25 * index - 1. / 12.));
  }
  return J;
}

// cre synchrotron J_pol(\nu)
ham_float Integrator::sync_emissivity_p(const Hamvec<3, ham_float> &pos,
                                        const Param *par,
                                        const ham_float &b_per) const {
  ham_float J{0};
  // calculate from grid
  if (par->grid_cre.read_permission or par->grid_cre.build_permission) {
    // allocate energy grid
    std::unique_ptr<ham_float[]> ke =
        std::make_unique<ham_float[]>(par->grid_cre.ne);
    // we need F(x[E]) in spectral integration
    std::unique_ptr<ham_float[]> x =
        std::make_unique<ham_float[]>(par->grid_cre.ne);
    std::unique_ptr<ham_float[]> beta =
        std::make_unique<ham_float[]>(par->grid_cre.ne);
    // consts used for converting E to x, using cgs units
    const ham_float x_fact{(2. * cgs::mec * cgs::mec2 * cgs::mec2 * 2. *
                            cgs::pi * par->grid_obs.sim_sync_freq.back()) /
                           (3. * cgs::qe * b_per)};
    // ke, x, beta arrays in cgs units
    for (decltype(par->grid_cre.ne) i = 0; i != par->grid_cre.ne; ++i) {
      ke[i] = par->grid_cre.e_min * std::exp(i * par->grid_cre.e_fact);
      x[i] = x_fact / (ke[i] * ke[i]);
      beta[i] = std::sqrt(1 - cgs::mec2 / ke[i]);
    }
    // do energy spectrum integration at given position
    // unit_factor for DIFFERENTIAL density flux, [GeV m^2 s sr]^-1
    // n(E,pos) = \phi(E,pos)*(4\pi/\beta*c), the relatin between flux \phi and
    // density n ref: "Cosmic rays n' particle physics", A3
    const ham_float fore_factor{
        cgs::sqrt3 * cgs::qe * cgs::qe * cgs::qe * b_per /
        (cgs::mec2 * cgs::c_light * cgs::GeV * cgs::m * cgs::m * cgs::sec)};
    // spectral integral
    for (decltype(par->grid_cre.ne) i = 0; i != par->grid_cre.ne - 1; ++i) {
      const ham_float xv{0.5 * (x[i + 1] + x[i])};
      // avoid underflow in gsl functions
      if (xv > 100.) {
        continue;
      }
      const ham_float de{std::fabs(ke[i + 1] - ke[i])};
      // we put beta here, midpoint rule
      const ham_float flux{
          0.5 *
          (fields.cre->read_field(pos, i + 1, par, grids.cre) / beta[i + 1] +
           fields.cre->read_field(pos, i, par, grids.cre) / beta[i])};
      assert(flux >= 0);
      J += gsl_sf_synchrotron_2(xv) * flux * de;
    }
    J *= fore_factor;
  }
  // calculate from N(\gamma) with local constant spectral index
  else {
    // allocating values to index, norm according to user defined model
    // user may consider building derived class from CRE_ana
    const ham_float index{fields.cre->flux_idx(pos, par)};
    // coefficients which do not attend integration
    const ham_float norm{fields.cre->flux_norm(pos, par) * cgs::sqrt3 *
                         (cgs::qe * cgs::qe * cgs::qe) * b_per /
                         (16. * cgs::pi * cgs::mec2)};
    // synchrotron integration
    const ham_float A{2. * cgs::pi * par->grid_obs.sim_sync_freq.back() *
                      cgs::mec / (3. * cgs::qe * b_per)};
    J = norm * (std::pow(A, 0.5 * (index + 1.)) *
                gsl_sf_gamma(-0.25 * index + 7. / 12.) *
                gsl_sf_gamma(-0.25 * index - 1. / 12.));
  }
  // the last 4pi comes from solid-angle integration/deviation,
  // check eq(6.16) in Ribiki-Lightman's where Power is defined,
  // we need isotropic power which means we need a 1/4pi factor!
  return J;
}
