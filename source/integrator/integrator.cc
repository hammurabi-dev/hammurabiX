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

#include <bfield.h>
#include <cgs_units_file.h>
#include <crefield.h>
#include <grid.h>
#include <hamvec.h>
#include <integrator.h>
#include <param.h>
#include <tefield.h>

void Integrator::write_grid(Breg *breg, Brnd *brnd, TEreg *tereg, TErnd *ternd,
                            CRE *cre, Grid_breg *gbreg, Grid_brnd *gbrnd,
                            Grid_tereg *gtereg, Grid_ternd *gternd,
                            Grid_cre *gcre, Grid_obs *gobs,
                            const Param *par) const {
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
  for (decltype(par->grid_obs.total_shell) current_shell = 1;
       current_shell != (par->grid_obs.total_shell + 1); ++current_shell) {
    // get current shell nside & npix
    const std::size_t current_nside{
        par->grid_obs.nside_shell[current_shell - 1]};
    const std::size_t current_npix{12 * current_nside * current_nside};
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
    this->assemble_shell_ref(shell_ref.get(), par, current_shell);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t ipix = 0; ipix < current_npix; ++ipix) {
      struct_observables observables;
      observables.is = 0.;
      observables.qs = 0.;
      observables.us = 0.;
      observables.dm = 0.;
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
        observables.fd = gobs->fd_map->interpolated_value(ptg);
      } else if (par->grid_obs.do_sync.back()) {
        ptg = gobs->tmp_is_map->pix2ang(ipix);
        // cache Faraday rotation from inner shells
        observables.fd = gobs->fd_map->interpolated_value(ptg);
      }
      // check angular direction boundaries
      if (this->check_simulation_lower_limit(0.5 * CGS_U_pi - ptg.theta,
                                             par->grid_obs.oc_lat_min)) {
        continue;
      }
      if (this->check_simulation_upper_limit(0.5 * CGS_U_pi - ptg.theta,
                                             par->grid_obs.oc_lat_max)) {
        continue;
      }
      if (this->check_simulation_lower_limit(ptg.phi,
                                             par->grid_obs.oc_lon_min)) {
        continue;
      }
      if (this->check_simulation_upper_limit(ptg.phi,
                                             par->grid_obs.oc_lon_max)) {
        continue;
      }
      // core function!
      this->radial_integration(shell_ref.get(), ptg, observables, breg, brnd,
                               tereg, ternd, cre, gbreg, gbrnd, gtereg, gternd,
                               gcre, par);
      // collect from pixels
      if (par->grid_obs.do_dm) {
        (*gobs->tmp_dm_map)[ipix] = observables.dm;
      }
      if (par->grid_obs.do_sync.back()) {
        (*gobs->tmp_is_map)[ipix] = this->temp_convert(
            observables.is, par->grid_obs.sim_sync_freq.back());
        (*gobs->tmp_qs_map)[ipix] = this->temp_convert(
            observables.qs, par->grid_obs.sim_sync_freq.back());
        (*gobs->tmp_us_map)[ipix] = this->temp_convert(
            observables.us, par->grid_obs.sim_sync_freq.back());
      }
      if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
        (*gobs->tmp_fd_map)[ipix] = observables.fd;
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

void Integrator::radial_integration(struct_shell *shell_ref, pointing &ptg_in,
                                    struct_observables &pixobs, Breg *breg,
                                    Brnd *brnd, TEreg *tereg, TErnd *ternd,
                                    CRE *cre, Grid_breg *gbreg,
                                    Grid_brnd *gbrnd, Grid_tereg *gtereg,
                                    Grid_ternd *gternd, Grid_cre *gcre,
                                    const Param *par) const {
  // pass in fd, zero others
  double inner_shells_fd{0.};
  if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
    inner_shells_fd = pixobs.fd;
  }
  pixobs.dm = 0.;
  pixobs.fd = 0.;
  pixobs.is = 0.;
  pixobs.us = 0.;
  pixobs.qs = 0.;
  // angular position
  const double THE{ptg_in.theta};
  const double PHI{ptg_in.phi};
  double lambda_square = 0., i2bt_sync = 0., fd_forefactor = 0.;
  if (par->grid_obs.do_sync.back()) {
    // for calculating synchrotron emission
    lambda_square = (CGS_U_C_light / par->grid_obs.sim_sync_freq.back()) *
                    (CGS_U_C_light / par->grid_obs.sim_sync_freq.back());
    // convert sync intensity(freq) to brightness temperature, Rayleigh-Jeans
    // law
    i2bt_sync = CGS_U_C_light * CGS_U_C_light /
                (2. * CGS_U_kB * par->grid_obs.sim_sync_freq.back() *
                 par->grid_obs.sim_sync_freq.back());
  }
  if (par->grid_obs.do_fd) {
    fd_forefactor = -(CGS_U_qe * CGS_U_qe * CGS_U_qe) /
                    (2. * CGS_U_pi * CGS_U_MEC2 * CGS_U_MEC2);
  }
  // radial accumulation
  for (decltype(shell_ref->step) looper = 0; looper < shell_ref->step;
       ++looper) {
    // ec and gc position
    hamvec<3, double> oc_pos{this->los_versor(THE, PHI) *
                             shell_ref->dist[looper]};
    hamvec<3, double> pos{oc_pos + par->observer};
    // check LoS depth limit
    if (this->check_simulation_lower_limit(pos.length(),
                                           par->grid_obs.gc_r_min)) {
      continue;
    }
    if (this->check_simulation_upper_limit(pos.length(),
                                           par->grid_obs.gc_r_max)) {
      continue;
    }
    if (this->check_simulation_lower_limit(pos[2], par->grid_obs.gc_z_min)) {
      continue;
    }
    if (this->check_simulation_upper_limit(pos[2], par->grid_obs.gc_z_max)) {
      continue;
    }
    // regular magnetic field
    hamvec<3, double> B_vec{breg->read_field(pos, par, gbreg)};
    // add random magnetic field
    B_vec += brnd->read_field(pos, par, gbrnd);
    const double B_par{this->los_parproj(B_vec, THE, PHI)};
    assert(std::isfinite(B_par));
    // be aware of un-resolved random B_per in calculating emissivity
    const double B_per{this->los_perproj(B_vec, THE, PHI)};
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
      pixobs.dm += te * shell_ref->delta_d;
    }
    // Faraday depth
    if (par->grid_obs.do_fd or par->grid_obs.do_sync.back()) {
      pixobs.fd += te * B_par * fd_forefactor * shell_ref->delta_d;
    }
    // Synchrotron emission
    if (par->grid_obs.do_sync.back()) {
      const double Jtot{cre->read_emissivity_t(pos, par, gcre, B_per) *
                        shell_ref->delta_d * i2bt_sync};
      // J_pol receives no contribution from unresolved random field
      const double Jpol{cre->read_emissivity_p(pos, par, gcre, B_per) *
                        shell_ref->delta_d * i2bt_sync};
      assert(Jtot < 1e30 and Jpol < 1e30 and Jtot >= 0 and Jpol >= 0);
      pixobs.is += Jtot;
      // intrinsic polarization angle, following IAU definition
      const double qui{(inner_shells_fd + pixobs.fd) * lambda_square +
                       this->sync_ipa(B_vec, THE, PHI)};
      assert(std::isfinite(qui));
      pixobs.qs += std::cos(2. * qui) * Jpol;
      pixobs.us += std::sin(2. * qui) * Jpol;
    }
  }
}

// assembling shell_ref structure
void Integrator::assemble_shell_ref(struct_shell *target, const Param *par,
                                    const std::size_t &shell_num) const {
  target->shell_num = shell_num;
  target->d_start = par->grid_obs.radii_shell[shell_num - 1];
  target->d_stop = par->grid_obs.radii_shell[shell_num];
  target->delta_d = par->grid_obs.oc_r_res;
  target->step = floor(
      (target->d_stop / target->delta_d - target->d_start / target->delta_d));
  // get rid of error in the previous step
  target->delta_d = (target->d_stop - target->d_start) / (target->step);
  for (std::size_t i = 0; i < target->step; ++i) {
    target->dist.push_back(target->d_start + i * 0.5 * target->delta_d);
  }
#ifdef VERBOSE
  std::cout << "shell reference: " << std::endl
            << "shell No. " << target->shell_num << std::endl
            << "resolution " << target->delta_d / CGS_U_kpc << " kpc"
            << std::endl
            << "d_start " << target->d_start / CGS_U_kpc << " kpc" << std::endl
            << "d_stop " << target->d_stop / CGS_U_kpc << " kpc" << std::endl
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
  return (input.crossprod(this->los_versor(the_los, phi_los))).length();
}

// (signed) parallel projection of a vector wrt LoS direction
double Integrator::los_parproj(const hamvec<3, double> &input,
                               const double &the_los,
                               const double &phi_los) const {
  return input.dotprod(this->los_versor(the_los, phi_los));
}

// converting brightness temp into thermal temp with T_0 = 2.725K
double Integrator::temp_convert(const double &temp_br,
                                const double &freq) const {
  const double p{CGS_U_h_planck * freq / (CGS_U_kB * 2.725)};
  return temp_br * (exp(p) - 1.) * (exp(p) - 1.) / (p * p * exp(p));
}
