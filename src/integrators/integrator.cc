#include <iostream>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_errno.h>
#include <omp.h>
#include <cmath>
#include <cassert>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include <healpix_base.h>
#include <healpix_map.h>
#include <pointing.h>
#include <vec3.h>
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

using namespace std;

void Integrator::write_grid(Breg *breg,Brnd *brnd,FEreg *fereg,FErnd *fernd,CRE *cre,Grid_breg *gbreg,Grid_brnd *gbrnd,Grid_fereg *gfereg,Grid_fernd *gfernd,Grid_cre *gcre,Grid_int *gint,Param *par) {
    // unsigned int, pre-calculated in gint
    size_t npix_sim {gint->npix_sim};
    if (gint->do_dm) {
        gint->dm_map.SetNside(gint->nside_sim, RING);
        gint->dm_map.fill(0.);
    }
    if (gint->do_sync) {
        gint->Is_map.SetNside(gint->nside_sim, RING);
        gint->Qs_map.SetNside(gint->nside_sim, RING);
        gint->Us_map.SetNside(gint->nside_sim, RING);
        gint->Is_map.fill(0.);
        gint->Qs_map.fill(0.);
        gint->Us_map.fill(0.);
    }
    if (gint->do_fd or gint->do_sync) {
        gint->fd_map.SetNside(gint->nside_sim, RING);
        gint->fd_map.fill(0.);
    }
    unique_ptr<struct_shell> shell_ref = unique_ptr<struct_shell>(new struct_shell);
    for (decltype(gint->total_shell) current_shell=1;current_shell!=(gint->total_shell+1);++current_shell) {
        Healpix_Map<double> current_Is_map;
        Healpix_Map<double> current_Qs_map;
        Healpix_Map<double> current_Us_map;
        Healpix_Map<double> current_fd_map;
        Healpix_Map<double> current_dm_map;
        
        size_t current_nside {gint->nside_shell[current_shell-1]};
        size_t current_npix {12*current_nside*current_nside};
        if (gint->do_dm) {
            current_dm_map.SetNside(current_nside, RING);
            current_dm_map.fill(0.);
        }
        if (gint->do_sync) {
            current_Is_map.SetNside(current_nside, RING);
            current_Qs_map.SetNside(current_nside, RING);
            current_Us_map.SetNside(current_nside, RING);
            current_Is_map.fill(0.);
            current_Qs_map.fill(0.);
            current_Us_map.fill(0.);
        }
        if (gint->do_fd or gint->do_sync) {
            current_fd_map.SetNside(current_nside, RING);
            current_fd_map.fill(0.);
        }
        // setting for radial_integration
        // call auxiliary function assemble_shell_ref
        // use unique_ptr to avoid stack overflow
        assemble_shell_ref(shell_ref.get(),gint,current_shell);
        
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (decltype(current_npix) ipix=0;ipix<current_npix;++ipix) {
            struct_observables observables;
            observables.Is = observables.Qs = observables.Us = 0.;
            observables.dm = 0.;
            // notice that either do_dm or do_fd should be true
            // if include dust emission, remember to complete logic for ptg assignment!
            pointing ptg;
            if (gint->do_dm) {
                ptg = current_dm_map.pix2ang(ipix);
            }
            // get fd out from inner shells
            if (gint->do_fd or gint->do_sync) {
                ptg = current_fd_map.pix2ang(ipix);
                observables.fd = gint->fd_map.interpolated_value(ptg);
            }
            // core function!
            radial_integration(shell_ref.get(),ptg,observables,breg,brnd,fereg,fernd,cre,gbreg,gbrnd,gfereg,gfernd,gcre,gint,par);
            
            // assembling new shell
            if (gint->do_dm) {
                current_dm_map[ipix] = observables.dm;
            }
            if (gint->do_sync) {
                current_Is_map[ipix] = toolkit::temp_convert(observables.Is,par->sim_freq);
                current_Qs_map[ipix] = toolkit::temp_convert(observables.Qs,par->sim_freq);
                current_Us_map[ipix] = toolkit::temp_convert(observables.Us,par->sim_freq);
            }
            if (gint->do_fd or gint->do_sync) {
                current_fd_map[ipix] = observables.fd;
            }
        }
        //adding up new shell map to sim map
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (decltype(npix_sim) ipix=0;ipix<npix_sim;++ipix) {
            pointing ptg;
            if (gint->do_dm) {
                ptg = gint->dm_map.pix2ang(ipix);
                gint->dm_map[ipix] += current_dm_map.interpolated_value(ptg);
            }
            if (gint->do_sync) {
                ptg = gint->Is_map.pix2ang(ipix);
                gint->Is_map[ipix] += current_Is_map.interpolated_value(ptg);
                gint->Qs_map[ipix] += current_Qs_map.interpolated_value(ptg);
                gint->Us_map[ipix] += current_Us_map.interpolated_value(ptg);
            }
            if (gint->do_fd or gint->do_sync) {
                ptg = gint->fd_map.pix2ang(ipix);
                gint->fd_map[ipix] += current_fd_map.interpolated_value(ptg);
            }
        }
    }//end shell iteration
}

void Integrator::radial_integration(struct_shell *shell_ref,pointing &ptg_in, struct_observables &pixobs,Breg *breg,Brnd *brnd,FEreg *fereg,FErnd *fernd,CRE *cre,Grid_breg *gbreg,Grid_brnd *gbrnd,Grid_fereg *gfereg,Grid_fernd *gfernd,Grid_cre *gcre,Grid_int *gint,Param *par) {
    // pass in fd, zero others
    double inner_shells_fd {0.};
    if (gint->do_fd or gint->do_sync) {inner_shells_fd=pixobs.fd;}
    pixobs.dm=0.;pixobs.fd=0.;
    pixobs.Is=0.;pixobs.Us=0.;pixobs.Qs=0.;
    // angular position
    const double THE {ptg_in.theta};
    const double PHI {ptg_in.phi};
    if (check_simulation_lower_limit(fabs(0.5*CGS_U_pi-THE),gint->lat_lim)) {return;}
    // for calculating synchrotron emission
    const double lambda_square{(CGS_U_C_light/par->sim_freq)*(CGS_U_C_light/par->sim_freq)};
    // convert intensity(freq) to brightness temperature, Rayleigh-Jeans law
    const double i2bt {CGS_U_C_light*CGS_U_C_light/(2.*CGS_U_kB*par->sim_freq*par->sim_freq)};
    const double fd_forefactor {-(CGS_U_qe*CGS_U_qe*CGS_U_qe)/(2.*CGS_U_pi*CGS_U_MEC2*CGS_U_MEC2)};
    // for Simpson's rule
    vector<double> F_dm,F_fd,F_Jtot,F_Jpol,intr_pol_ang;
    decltype(shell_ref->step) looper;
    for(looper=0;looper<shell_ref->step;++looper){
        // ec and gc position
        vec3_t<double> ec_pos {toolkit::los_versor(THE,PHI)*shell_ref->dist[looper]};
        vec3_t<double> pos {ec_pos + par->SunPosition};
        // check LOS depth limit
        if (check_simulation_upper_limit(pos.Length(),gint->gc_r_max)) {break;}
        if (check_simulation_upper_limit(fabs(pos.z),gint->gc_z_max)) {break;}
        // B field
        vec3_t<double> B_vec {breg->get_breg(pos,par,gbreg)};
        // add random field
        B_vec += brnd->get_brnd(pos,gbrnd);
        const double B_par {toolkit::par2los(B_vec,THE,PHI)};
        // un-scaled B_per, if random B is active, we scale up this
        // in calculating emissivity, so it is not constant
        double B_per {toolkit::perp2los(B_vec,THE,PHI)};
        // FE field
        double te {fereg->get_density(pos,par,gfereg)};
        // add random field
        te += fernd->get_fernd(pos,gfernd);
        // to avoid negative value
        if(te<0){
            te = 0.;
        }
        // DM
        if(gint->do_dm){
            F_dm.push_back(te*shell_ref->delta_d);
        }
        // Faraday depth
        if(gint->do_fd or gint->do_sync){
            F_fd.push_back(te*B_par*fd_forefactor*shell_ref->delta_d);
        }
        // Synchrotron Emission
        if(gint->do_sync){
            F_Jtot.push_back(cre->get_emissivity_t(pos,par,gcre,B_per)*shell_ref->delta_d*i2bt);
            // J_pol receive no contribution from missing random
            F_Jpol.push_back(cre->get_emissivity_p(pos,par,gcre,B_per)*shell_ref->delta_d*i2bt);
            // intrinsic polarization angle, following IAU definition
            intr_pol_ang.push_back(toolkit::intr_pol_ang(B_vec,THE,PHI));
        }
    }// first iteration end
    // applying Simpson's rule
    for(decltype(shell_ref->step)i=1;i<looper-1;i+=2){
        // DM
        if(gint->do_dm){
            pixobs.dm += (F_dm[i-1]+4.*F_dm[i]+F_dm[i+1])*0.16666667;
        }
        // FD
        if(gint->do_fd or gint->do_sync){
            pixobs.fd += (F_fd[i-1]+4.*F_fd[i]+F_fd[i+1])*0.16666667;
        }
        // Sync
        if(gint->do_sync){
            // pol. angle after Faraday rotation
            double qui_base {(inner_shells_fd+pixobs.fd)*lambda_square};
            double qui_0 {qui_base+intr_pol_ang[i-1]};
            double qui_1 {qui_base+intr_pol_ang[i]};
            double qui_2 {qui_base+intr_pol_ang[i+1]};
            
            assert(abs(qui_0)+abs(qui_1)+abs(qui_2)<1e30);
            assert(F_Jtot[i-1]>=0 and F_Jtot[i]>=0 and F_Jtot[i+1]>=0);
            
            pixobs.Is += (F_Jtot[i-1]+4.*F_Jtot[i]+F_Jtot[i+1])*0.16666667;
            pixobs.Qs += (cos(2.*qui_0)*F_Jpol[i-1]+4.*cos(2.*qui_1)*F_Jpol[i]+cos(2.*qui_2)*F_Jpol[i+1])*0.16666667;
            pixobs.Us += (sin(2.*qui_0)*F_Jpol[i-1]+4.*sin(2.*qui_1)*F_Jpol[i]+sin(2.*qui_2)*F_Jpol[i+1])*0.16666667;
        }
    }
}//end of radial_integrate

// TOOLS
//---------------------------------------------------------

// compute upper radial boundary at given shell
double Integrator::get_max_shell_radius(const size_t &shell_numb,const size_t &total_shell,const double &radius) const {
    assert(shell_numb>=1 and shell_numb<=total_shell);
    double max_shell_radius {radius};
    for (size_t n=total_shell;n!=shell_numb;--n) {
        max_shell_radius *= 0.5;
    }
    return max_shell_radius;
}

// compute lower radial boundary at given shell
double Integrator::get_min_shell_radius(const size_t &shell_numb,const size_t &total_shell,const double &radius) const {
    assert(shell_numb>=1 and shell_numb<=total_shell);
    // min_shell_radius for the first shell is always zero.
    if(shell_numb==1){
        return 0.;
    }
    double min_shell_radius {radius};
    for (size_t n=total_shell;n!=(shell_numb-1);--n) {
        min_shell_radius *= 0.5;
    }
    return min_shell_radius;
}

// assembling shell_ref structure
void Integrator::assemble_shell_ref(struct_shell *target,const Grid_int *info,const size_t &shell_num){
    target->shell_num = shell_num;
    target->d_start = get_min_shell_radius(shell_num,info->total_shell,info->ec_r_max);
    target->d_stop = get_max_shell_radius(shell_num,info->total_shell,info->ec_r_max);
    target->delta_d = info->radial_res;
    target->step = floor((2.*target->d_stop - 2.*target->d_start)/(target->delta_d))+1;
    for(size_t i=0;i<target->step;++i){
        target->dist.push_back(target->d_start+i*0.5*target->delta_d);
    }
}

//END ALL
