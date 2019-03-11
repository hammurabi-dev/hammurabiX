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
        gint->dm_map->fill (0.);
    }
    if (par->grid_int.do_sync.back()) {
        gint->Is_map->fill (0.);
        gint->Qs_map->fill (0.);
        gint->Us_map->fill (0.);
    }
    if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
        gint->fd_map->fill (0.);
    }
    auto shell_ref = std::make_unique<struct_shell>();
    for (decltype(par->grid_int.total_shell) current_shell=1;current_shell!=(par->grid_int.total_shell+1);++current_shell) {
        auto current_Is_map = std::make_unique<Healpix_Map<double>> ();
        auto current_Qs_map = std::make_unique<Healpix_Map<double>> ();
        auto current_Us_map = std::make_unique<Healpix_Map<double>> ();
        auto current_fd_map = std::make_unique<Healpix_Map<double>> ();
        auto current_dm_map = std::make_unique<Healpix_Map<double>> ();
        
        const std::size_t current_nside {par->grid_int.nside_shell[current_shell-1]};
        const std::size_t current_npix {12*current_nside*current_nside};
        if (par->grid_int.do_dm) {
            current_dm_map->SetNside (current_nside,
                                     RING);
            current_dm_map->fill (0.);
        }
        if (par->grid_int.do_sync.back()) {
            current_Is_map->SetNside (current_nside,
                                     RING);
            current_Qs_map->SetNside (current_nside,
                                     RING);
            current_Us_map->SetNside (current_nside,
                                     RING);
            current_Is_map->fill (0.);
            current_Qs_map->fill (0.);
            current_Us_map->fill (0.);
        }
        if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
            current_fd_map->SetNside (current_nside,
                                     RING);
            current_fd_map->fill (0.);
        }
        // setting for radial_integration
        // call auxiliary function assemble_shell_ref
        // use unique_ptr to avoid stack overflow
        assemble_shell_ref (shell_ref.get(),
                            par,
                            current_shell);
        
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::size_t ipix=0;ipix<current_npix;++ipix) {
            struct_observables observables;
            observables.Is = 0.;
            observables.Qs = 0.;
            observables.Us = 0.;
            observables.dm = 0.;
            // notice that either do_dm or do_fd should be true
            // if include dust emission, remember to complete logic for ptg assignment!
            pointing ptg;
            if (par->grid_int.do_dm) {
                ptg = current_dm_map->pix2ang (ipix);
            }
            // get fd out from inner shells
            if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
                ptg = current_fd_map->pix2ang (ipix);
                observables.fd = gint->fd_map->interpolated_value (ptg);
            }
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
                (*current_dm_map)[ipix] = observables.dm;
            }
            if (par->grid_int.do_sync.back()) {
                (*current_Is_map)[ipix] = toolkit::temp_convert (observables.Is,par->grid_int.sim_sync_freq.back());
                (*current_Qs_map)[ipix] = toolkit::temp_convert (observables.Qs,par->grid_int.sim_sync_freq.back());
                (*current_Us_map)[ipix] = toolkit::temp_convert (observables.Us,par->grid_int.sim_sync_freq.back());
            }
            if (par->grid_int.do_fd or par->grid_int.do_sync.back()) {
                (*current_fd_map)[ipix] = observables.fd;
            }
        }
        //adding up new shell map to sim map
        if (par->grid_int.do_dm) {
            std::size_t npix_dm = gint->dm_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (std::size_t ipix=0;ipix<npix_dm;++ipix) {
                pointing ptg {gint->dm_map->pix2ang (ipix)};
                (*gint->dm_map)[ipix] += current_dm_map->interpolated_value (ptg);
            }
        }
        if (par->grid_int.do_sync.back()) {
            std::size_t npix_sync = gint->Is_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (std::size_t ipix=0;ipix<npix_sync;++ipix) {
                pointing ptg = {gint->Is_map->pix2ang (ipix)};
                (*gint->Is_map)[ipix] += current_Is_map->interpolated_value (ptg);
                (*gint->Qs_map)[ipix] += current_Qs_map->interpolated_value (ptg);
                (*gint->Us_map)[ipix] += current_Us_map->interpolated_value (ptg);
            }
        }
        if (par->grid_int.do_fd) {
            std::size_t npix_fd = gint->fd_map->Npix();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
            for (std::size_t ipix=0;ipix<npix_fd;++ipix) {
                pointing ptg = {gint->fd_map->pix2ang (ipix)};
                (*gint->fd_map)[ipix] += current_fd_map->interpolated_value (ptg);
            }
        }
    }//end shell iteration
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
    if (check_simulation_lower_limit (0.5*CGS_U_pi-THE,par->grid_int.lat_min)) {return;}
    if (check_simulation_upper_limit (0.5*CGS_U_pi-THE,par->grid_int.lat_max)) {return;}
    if (check_simulation_lower_limit (PHI,par->grid_int.lon_min)) {return;}
    if (check_simulation_upper_limit (PHI,par->grid_int.lon_max)) {return;}
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
    // for Simpson's rule
    std::vector<double> F_dm,F_fd,F_Jtot,F_Jpol,intr_pol_ang;
    decltype(shell_ref->step) looper;
    for(looper=0;looper<shell_ref->step;++looper){
        // ec and gc position
        hvec<3,double> ec_pos {toolkit::los_versor(THE,PHI)*shell_ref->dist[looper]};
        hvec<3,double> pos {ec_pos + par->observer};
        // check LOS depth limit
        if (check_simulation_lower_limit (pos.length(),par->grid_int.gc_r_min)) {break;}
        if (check_simulation_upper_limit (pos.length(),par->grid_int.gc_r_max)) {break;}
        if (check_simulation_lower_limit (std::fabs(pos[2]),par->grid_int.gc_z_min)) {break;}
        if (check_simulation_upper_limit (std::fabs(pos[2]),par->grid_int.gc_z_max)) {break;}
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
        if (te<0){
            te = 0.;
        }
        // DM
        if (par->grid_int.do_dm){
            F_dm.push_back (te*shell_ref->delta_d);
        }
        // Faraday depth
        if (par->grid_int.do_fd or par->grid_int.do_sync.back()){
            F_fd.push_back (te*B_par*fd_forefactor*shell_ref->delta_d);
        }
        // Synchrotron Emission
        if (par->grid_int.do_sync.back()){
            F_Jtot.push_back (cre->get_emissivity_t(pos,par,gcre,B_per)*shell_ref->delta_d*i2bt_sync);
            // J_pol receive no contribution from missing random
            F_Jpol.push_back (cre->get_emissivity_p(pos,par,gcre,B_per)*shell_ref->delta_d*i2bt_sync);
            // intrinsic polarization angle, following IAU definition
            intr_pol_ang.push_back (toolkit::intr_pol_ang(B_vec,THE,PHI));
        }
    }// first iteration end
    // applying Simpson's rule
    for (decltype(shell_ref->step) i=1;i<looper-1;i+=2){
        // DM
        if (par->grid_int.do_dm){
            pixobs.dm += (F_dm[i-1]+4.*F_dm[i]+F_dm[i+1])*0.16666667;
        }
        // FD
        if (par->grid_int.do_fd or par->grid_int.do_sync.back()){
            pixobs.fd += (F_fd[i-1]+4.*F_fd[i]+F_fd[i+1])*0.16666667;
        }
        // Sync
        if (par->grid_int.do_sync.back()){
            // pol. angle after Faraday rotation
            double qui_base {(inner_shells_fd+pixobs.fd)*lambda_square};
            double qui_0 {qui_base+intr_pol_ang[i-1]};
            double qui_1 {qui_base+intr_pol_ang[i]};
            double qui_2 {qui_base+intr_pol_ang[i+1]};
            
            assert (abs(qui_0)+abs(qui_1)+abs(qui_2)<1e30);
            assert (F_Jtot[i-1]>=0 and F_Jtot[i]>=0 and F_Jtot[i+1]>=0);
            
            pixobs.Is += (F_Jtot[i-1]+4.*F_Jtot[i]+F_Jtot[i+1])*0.16666667;
            pixobs.Qs += (cos(2.*qui_0)*F_Jpol[i-1]+4.*cos(2.*qui_1)*F_Jpol[i]+cos(2.*qui_2)*F_Jpol[i+1])*0.16666667;
            pixobs.Us += (sin(2.*qui_0)*F_Jpol[i-1]+4.*sin(2.*qui_1)*F_Jpol[i]+sin(2.*qui_2)*F_Jpol[i+1])*0.16666667;
        }
    }
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
	// dragging target->delta_d outside will lose precision
    target->step = floor(2.*(target->d_stop/target->delta_d - target->d_start/target->delta_d))+1;
    for (std::size_t i=0;i<target->step;++i){
        target->dist.push_back(target->d_start+i*0.5*target->delta_d);
    }
}

//END ALL
