#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <cassert>

#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include <healpix_base.h>
#include <healpix_map.h>
#include <pointing.h>
#include <hvec.h>
#include <breg.h>
#include <brnd.h>
#include <cre.h>
#include <integrator.h>
#include <fereg.h>
#include <fernd.h>
#include <param.h>
#include <grid.h>
#include <cgs_units_file.h>
#include <namespace_toolkit.h>

void Integrator::write_grid (Breg *breg,
                             Brnd *brnd,
                             FEreg *fereg,
                             FErnd *fernd,
                             CRE *cre,
                             Grid_breg *gbreg,
                             Grid_brnd *gbrnd,
                             Grid_fereg *gfereg,
                             Grid_fernd *gfernd,
                             Grid_cre *gcre,
                             Grid_int *gint,
                             const Param *par) const{
    if (par->grid_int.do_dm) {
        gint->tmp_dm_map->fill (0.);
    }
    if (par->grid_int.do_sync.back()) {
        gint->tmp_Is_map->fill (0.);
        gint->tmp_Qs_map->fill (0.);
        gint->tmp_Us_map->fill (0.);
    }
    if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
        gint->tmp_fd_map->fill (0.);
    }
    auto shell_ref = std::make_unique<struct_shell>();
    const std::size_t npix_sim {12*par->grid_int.nside_sim*par->grid_int.nside_sim};
    for (decltype(par->grid_int.total_shell) current_shell=1;current_shell!=(par->grid_int.total_shell+1);++current_shell) {
        // setting for radial_integration
        // call auxiliary function assemble_shell_ref
        // use unique_ptr to avoid stack overflow
        assemble_shell_ref (shell_ref.get(),
                            par,
                            current_shell);
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (std::size_t ipix=0;ipix<npix_sim;++ipix) {
            struct_observables observables;
            observables.Is = 0.;
            observables.Qs = 0.;
            observables.Us = 0.;
            observables.dm = 0.;
            // accumulate Faraday rotation
            observables.fd = (*gint->tmp_fd_map)[ipix];
            // remember to complete logic for ptg assignment!
            pointing ptg;
            if (par->grid_int.do_dm) {
                ptg = gint->tmp_dm_map->pix2ang (ipix);
            }
            else if (par->grid_int.do_fd){
                ptg = gint->tmp_fd_map->pix2ang (ipix);
            }
            else if (par->grid_int.do_sync.back()) {
                ptg = gint->tmp_Is_map->pix2ang (ipix);
            }
            // check angular direction boundaries
            if (check_simulation_lower_limit (0.5*CGS_U_pi-ptg.theta,par->grid_int.lat_min)) {continue;}
            if (check_simulation_upper_limit (0.5*CGS_U_pi-ptg.theta,par->grid_int.lat_max)) {continue;}
            if (check_simulation_lower_limit (ptg.phi,par->grid_int.lon_min)) {continue;}
            if (check_simulation_upper_limit (ptg.phi,par->grid_int.lon_max)) {continue;}
            // core function!
            radial_integration (shell_ref.get(),
                                ptg,
                                observables,
                                breg,
                                brnd,
                                fereg,
                                fernd,
                                cre,
                                gbreg,
                                gbrnd,
                                gfereg,
                                gfernd,
                                gcre,
                                par);
            // assembling new shell
            if (par->grid_int.do_dm) {
                (*gint->tmp_dm_map)[ipix] += observables.dm;
            }
            if (par->grid_int.do_sync.back()) {
                (*gint->tmp_Is_map)[ipix] += toolkit::temp_convert (observables.Is,par->grid_int.sim_sync_freq.back());
                (*gint->tmp_Qs_map)[ipix] += toolkit::temp_convert (observables.Qs,par->grid_int.sim_sync_freq.back());
                (*gint->tmp_Us_map)[ipix] += toolkit::temp_convert (observables.Us,par->grid_int.sim_sync_freq.back());
            }
            // accumulate Faraday rotation
            if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
                (*gint->tmp_fd_map)[ipix] += observables.fd;
            }
        }
    }// shell
}

void Integrator::radial_integration (struct_shell *shell_ref,
                                     pointing &ptg_in,
                                     struct_observables &pixobs,
                                     Breg *breg,
                                     Brnd *brnd,
                                     FEreg *fereg,
                                     FErnd *fernd,
                                     CRE *cre,
                                     Grid_breg *gbreg,
                                     Grid_brnd *gbrnd,
                                     Grid_fereg *gfereg,
                                     Grid_fernd *gfernd,
                                     Grid_cre *gcre,
                                     const Param *par) const{
    // pass in fd, zero others
    double inner_shells_fd {0.};
    if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {inner_shells_fd=pixobs.fd;}
    pixobs.dm=0.;pixobs.fd=0.;
    pixobs.Is=0.;pixobs.Us=0.;pixobs.Qs=0.;
    // angular position
    const double THE {ptg_in.theta};
    const double PHI {ptg_in.phi};
    double lambda_square=0.,i2bt_sync=0.,fd_forefactor=0.;
    if (par->grid_int.do_sync.back()){
        // for calculating synchrotron emission
        lambda_square = (CGS_U_C_light/par->grid_int.sim_sync_freq.back())*(CGS_U_C_light/par->grid_int.sim_sync_freq.back());
        // convert sync intensity(freq) to brightness temperature, Rayleigh-Jeans law
        i2bt_sync = CGS_U_C_light*CGS_U_C_light/(2.*CGS_U_kB*par->grid_int.sim_sync_freq.back()*par->grid_int.sim_sync_freq.back());
    }
    if (par->grid_int.do_fd){
        fd_forefactor = -(CGS_U_qe*CGS_U_qe*CGS_U_qe)/(2.*CGS_U_pi*CGS_U_MEC2*CGS_U_MEC2);
    }
    // radial accumulation
    for(decltype(shell_ref->step) looper=0;looper<shell_ref->step;++looper){
        // ec and gc position
        hvec<3,double> ec_pos {toolkit::los_versor(THE,PHI)*shell_ref->dist[looper]};
        hvec<3,double> pos {ec_pos + par->observer};
        // check LoS depth limit
        if (check_simulation_lower_limit (pos.length(),par->grid_int.gc_r_min)) {continue;}
        if (check_simulation_upper_limit (pos.length(),par->grid_int.gc_r_max)) {continue;}
        if (check_simulation_lower_limit (std::fabs(pos[2]),par->grid_int.gc_z_min)) {continue;}
        if (check_simulation_upper_limit (std::fabs(pos[2]),par->grid_int.gc_z_max)) {continue;}
        // B field
        hvec<3,double> B_vec {breg->get_breg(pos,par,gbreg)};
        // add random field
        B_vec += brnd->get_brnd (pos,par,gbrnd);
        const double B_par {toolkit::par2los(B_vec,THE,PHI)};
        // un-scaled B_per, if random B is active, we scale up this
        // in calculating emissivity, so it is not constant
        double B_per {toolkit::perp2los (B_vec,THE,PHI)};
        // FE field
        double te {fereg->get_density (pos,par,gfereg)};
        // add random field
        te += fernd->get_fernd (pos,par,gfernd);
        // to avoid negative value
        te *= double(te>0.);
        // DM
        if (par->grid_int.do_dm){
            pixobs.dm += te*shell_ref->delta_d;
        }
        // Faraday depth
        if (par->grid_int.do_fd or par->grid_int.do_sync.back()){
            pixobs.fd += te*B_par*fd_forefactor*shell_ref->delta_d;
        }
        // Synchrotron Emission
        if (par->grid_int.do_sync.back()){
            pixobs.Is += cre->get_emissivity_t(pos,par,gcre,B_per)*shell_ref->delta_d*i2bt_sync;
            // J_pol receives no contribution from unresolved random field
            const double Jpol {cre->get_emissivity_p(pos,par,gcre,B_per)*shell_ref->delta_d*i2bt_sync};
            // intrinsic polarization angle, following IAU definition
            const double qui {(inner_shells_fd+pixobs.fd)*lambda_square+toolkit::intr_pol_ang(B_vec,THE,PHI)};
            pixobs.Qs += cos(2.*qui)*Jpol;
            pixobs.Us += sin(2.*qui)*Jpol;
        }
    }// precalc
}//end of radial_integrate

// TOOLS
//---------------------------------------------------------

// assembling shell_ref structure
void Integrator::assemble_shell_ref (struct_shell *target,
                                     const Param *par,
                                     const std::size_t &shell_num) const{
    target->shell_num = shell_num;
    target->d_start = par->grid_int.radii_shell[shell_num-1];
    target->d_stop = par->grid_int.radii_shell[shell_num];
    target->delta_d = par->grid_int.radial_res;
    target->step = floor((target->d_stop/target->delta_d - target->d_start/target->delta_d))+1;
    // get rid of error in the previous step
    target->delta_d = (target->d_stop - target->d_start)/(target->step-1);
    for (std::size_t i=0;i<target->step;++i){
        target->dist.push_back(target->d_start+i*target->delta_d);
    }
}

//END ALL
