#include <cassert>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <vector>

#include <breg.h>
#include <brnd.h>
#include <cgs_units_file.h>
#include <cre.h>
#include <fereg.h>
#include <fernd.h>
#include <fitshandle.h>
#include <grid.h>
#include <healpix_base.h>
#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <hvec.h>
#include <integrator.h>
#include <namespace_toolkit.h>
#include <param.h>
#include <pointing.h>

void Integrator::write_grid(Breg *breg, Brnd *brnd, FEreg *fereg, FErnd *fernd,
                            CRE *cre, Grid_breg *gbreg, Grid_brnd *gbrnd,
                            Grid_fereg *gfereg, Grid_fernd *gfernd,
                            Grid_cre *gcre, Grid_int *gint,
                            const Param *par) const {
  if (par->grid_int.do_dm) {
    gint->dm_map->fill(0.);
  }
  if (par->grid_int.do_sync.back()) {
    gint->Is_map->fill(0.);
    gint->Qs_map->fill(0.);
    gint->Us_map->fill(0.);
  }
  if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
    gint->fd_map->fill(0.);
  }
  auto shell_ref = std::make_unique<struct_shell>();
  // loop through shells
  for (decltype(par->grid_int.total_shell) current_shell = 1;
       current_shell != (par->grid_int.total_shell + 1); ++current_shell) {
    // fetch current shell nside & npix
    const std::size_t current_nside{
        par->grid_int.nside_shell[current_shell - 1]};
    const std::size_t current_npix{12 * current_nside * current_nside};
    // prepare temporary maps for current shell
    if (par->grid_int.do_dm) {
      gint->tmp_dm_map->SetNside(current_nside, RING);
      gint->tmp_dm_map->fill(0.);
    }
    if (par->grid_int.do_sync.back()) {
      gint->tmp_Is_map->SetNside(current_nside, RING);
      gint->tmp_Qs_map->SetNside(current_nside, RING);
      gint->tmp_Us_map->SetNside(current_nside, RING);
      gint->tmp_Is_map->fill(0.);
      gint->tmp_Qs_map->fill(0.);
      gint->tmp_Us_map->fill(0.);
    }
    if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
      gint->tmp_fd_map->SetNside(current_nside, RING);
      gint->tmp_fd_map->fill(0.);
    }
    // setting for radial_integration
    // call auxiliary function assemble_shell_ref
    assemble_shell_ref(shell_ref.get(), par, current_shell);
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (std::size_t ipix = 0; ipix < current_npix; ++ipix) {
      struct_observables observables;
      observables.Is = 0.;
      observables.Qs = 0.;
      observables.Us = 0.;
      observables.dm = 0.;
      // remember to complete logic for ptg assignment!
      pointing ptg;
      if (par->grid_int.do_dm) {
        ptg = gint->tmp_dm_map->pix2ang(ipix);
      } else if (par->grid_int.do_fd) {
        ptg = gint->tmp_fd_map->pix2ang(ipix);
        // accumulate Faraday rotation
        observables.fd = gint->fd_map->interpolated_value(ptg);
      } else if (par->grid_int.do_sync.back()) {
        ptg = gint->tmp_Is_map->pix2ang(ipix);
        // accumulate Faraday rotation
        observables.fd = gint->fd_map->interpolated_value(ptg);
      }
      // check angular direction boundaries
      if (check_simulation_lower_limit(0.5 * CGS_U_pi - ptg.theta,
                                       par->grid_int.lat_min)) {
        continue;
      }
      if (check_simulation_upper_limit(0.5 * CGS_U_pi - ptg.theta,
                                       par->grid_int.lat_max)) {
        continue;
      }
      if (check_simulation_lower_limit(ptg.phi, par->grid_int.lon_min)) {
        continue;
      }
      if (check_simulation_upper_limit(ptg.phi, par->grid_int.lon_max)) {
        continue;
      }
      // core function!
      radial_integration(shell_ref.get(), ptg, observables, breg, brnd, fereg,
                         fernd, cre, gbreg, gbrnd, gfereg, gfernd, gcre, par);
      // assembling new shell
      if (par->grid_int.do_dm) {
        (*gint->tmp_dm_map)[ipix] += observables.dm;
      }
      if (par->grid_int.do_sync.back()) {
        (*gint->tmp_Is_map)[ipix] += toolkit::temp_convert(
            observables.Is, par->grid_int.sim_sync_freq.back());
        (*gint->tmp_Qs_map)[ipix] += toolkit::temp_convert(
            observables.Qs, par->grid_int.sim_sync_freq.back());
        (*gint->tmp_Us_map)[ipix] += toolkit::temp_convert(
            observables.Us, par->grid_int.sim_sync_freq.back());
      }
      // accumulate Faraday rotation
      if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
        (*gint->tmp_fd_map)[ipix] += observables.fd;
      }
    }
    // accumulating new shell map to sim map
    if (par->grid_int.do_dm) {
      std::size_t npix_dm = gint->dm_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t ipix = 0; ipix < npix_dm; ++ipix) {
        pointing ptg{gint->dm_map->pix2ang(ipix)};
        (*gint->dm_map)[ipix] += gint->tmp_dm_map->interpolated_value(ptg);
      }
    }
    if (par->grid_int.do_sync.back()) {
      std::size_t npix_sync = gint->Is_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t ipix = 0; ipix < npix_sync; ++ipix) {
        pointing ptg = {gint->Is_map->pix2ang(ipix)};
        (*gint->Is_map)[ipix] += gint->tmp_Is_map->interpolated_value(ptg);
        (*gint->Qs_map)[ipix] += gint->tmp_Qs_map->interpolated_value(ptg);
        (*gint->Us_map)[ipix] += gint->tmp_Us_map->interpolated_value(ptg);
      }
    }
    if (par->grid_int.do_fd) {
      std::size_t npix_fd = gint->fd_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (std::size_t ipix = 0; ipix < npix_fd; ++ipix) {
        pointing ptg = {gint->fd_map->pix2ang(ipix)};
        (*gint->fd_map)[ipix] += gint->tmp_fd_map->interpolated_value(ptg);
      }
    } // end shell accumulation
  }   // end shell iteration
}

void Integrator::radial_integration(struct_shell *shell_ref, pointing &ptg_in,
                                    struct_observables &pixobs, Breg *breg,
                                    Brnd *brnd, FEreg *fereg, FErnd *fernd,
                                    CRE *cre, Grid_breg *gbreg,
                                    Grid_brnd *gbrnd, Grid_fereg *gfereg,
                                    Grid_fernd *gfernd, Grid_cre *gcre,
                                    const Param *par) const {
  // pass in fd, zero others
  double inner_shells_fd{0.};
  if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
    inner_shells_fd = pixobs.fd;
  }
  pixobs.dm = 0.;
  pixobs.fd = 0.;
  pixobs.Is = 0.;
  pixobs.Us = 0.;
  pixobs.Qs = 0.;
  // angular position
  const double THE{ptg_in.theta};
  const double PHI{ptg_in.phi};
  double lambda_square = 0., i2bt_sync = 0., fd_forefactor = 0.;
  if (par->grid_int.do_sync.back()) {
    // for calculating synchrotron emission
    lambda_square = (CGS_U_C_light / par->grid_int.sim_sync_freq.back()) *
                    (CGS_U_C_light / par->grid_int.sim_sync_freq.back());
    // convert sync intensity(freq) to brightness temperature, Rayleigh-Jeans
    // law
    i2bt_sync = CGS_U_C_light * CGS_U_C_light /
                (2. * CGS_U_kB * par->grid_int.sim_sync_freq.back() *
                 par->grid_int.sim_sync_freq.back());
  }
  if (par->grid_int.do_fd) {
    fd_forefactor = -(CGS_U_qe * CGS_U_qe * CGS_U_qe) /
                    (2. * CGS_U_pi * CGS_U_MEC2 * CGS_U_MEC2);
  }
  // radial accumulation
  for (decltype(shell_ref->step) looper = 0; looper < shell_ref->step;
       ++looper) {
    // ec and gc position
    hvec<3, double> ec_pos{toolkit::los_versor(THE, PHI) *
                           shell_ref->dist[looper]};
    hvec<3, double> pos{ec_pos + par->observer};
    // check LoS depth limit
    if (check_simulation_lower_limit(pos.length(), par->grid_int.gc_r_min)) {
      continue;
    }
    if (check_simulation_upper_limit(pos.length(), par->grid_int.gc_r_max)) {
      continue;
    }
    if (check_simulation_lower_limit(pos[2], par->grid_int.gc_z_min)) {
      continue;
    }
    if (check_simulation_upper_limit(pos[2], par->grid_int.gc_z_max)) {
      continue;
    }
    // B field
    hvec<3, double> B_vec{breg->get_vector(pos, par, gbreg)};
    // add random field
    B_vec += brnd->get_vector(pos, par, gbrnd);
    const double B_par{toolkit::par2los(B_vec, THE, PHI)};
    assert(std::isfinite(B_par));
    // be aware of un-resolved random B_per in calculating emissivity
    const double B_per{toolkit::perp2los(B_vec, THE, PHI)};
    assert(std::isfinite(B_per));
    // FE field
    double te{fereg->get_density(pos, par, gfereg)};
    // add random field
    te += fernd->get_density(pos, par, gfernd);
    // to avoid negative value
    te *= double(te > 0.);
    assert(std::isfinite(te));
    // DM
    if (par->grid_int.do_dm) {
      pixobs.dm += te * shell_ref->delta_d;
    }
    // Faraday depth
    if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
      pixobs.fd += te * B_par * fd_forefactor * shell_ref->delta_d;
    }
    // Synchrotron Emission
    if (par->grid_int.do_sync.back()) {
      const double Jtot{cre->get_emissivity_t(pos, par, gcre, B_per) *
                        shell_ref->delta_d * i2bt_sync};
      // J_pol receives no contribution from unresolved random field
      const double Jpol{cre->get_emissivity_p(pos, par, gcre, B_per) *
                        shell_ref->delta_d * i2bt_sync};
      assert(Jtot < 1e30 and Jpol < 1e30 and Jtot >= 0 and Jpol >= 0);
      pixobs.Is += Jtot;
      // intrinsic polarization angle, following IAU definition
      const double qui{(inner_shells_fd + pixobs.fd) * lambda_square +
                       toolkit::intr_pol_ang(B_vec, THE, PHI)};
      assert(std::isfinite(qui));
      pixobs.Qs += cos(2. * qui) * Jpol;
      pixobs.Us += sin(2. * qui) * Jpol;
    }
  } // precalc
} // end of radial_integrate

// TOOLS
//---------------------------------------------------------

// assembling shell_ref structure
void Integrator::assemble_shell_ref(struct_shell *target, const Param *par,
                                    const std::size_t &shell_num) const {
  target->shell_num = shell_num;
  target->d_start = par->grid_int.radii_shell[shell_num - 1];
  target->d_stop = par->grid_int.radii_shell[shell_num];
  target->delta_d = par->grid_int.radial_res;
  target->step = floor((target->d_stop / target->delta_d -
                        target->d_start / target->delta_d));
  // get rid of error in the previous step
  target->delta_d = (target->d_stop - target->d_start) / (target->step);
  for (std::size_t i = 0; i < target->step; ++i) {
    target->dist.push_back(target->d_start + i * 0.5 * target->delta_d);
  }
}

// END ALL
