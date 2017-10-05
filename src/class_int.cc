#include <iostream>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_synchrotron.h>
#include <gsl/gsl_errno.h>
#include <omp.h>
#include <cmath>
#include <fitshandle.h>
#include <healpix_map_fitsio.h>
#include <healpix_base.h>
#include <healpix_map.h>
#include <vec3.h>

#include "class_breg.h"
#include "class_brnd.h"
#include "class_cre.h"
#include "class_int.h"
#include "class_fe.h"
#include "class_fernd.h"
#include "class_pond.h"
#include "class_grid.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

void Integrator::write_grid(Breg *breg,Brnd *brnd,FE *fe,FErnd *fernd,CRE *cre,Grid_breg *gbreg,Grid_brnd *gbrnd,Grid_fe *gfe,Grid_fernd *gfernd,Grid_cre *gcre,Grid_int *gint,Pond *par) {
    // unsigned int, pre-calculated in gint
    unsigned int npix_sim {gint->npix_sim};
    if (gint->do_dm) {
        gint->dm_map.SetNside(gint->nside_sim, NEST);
        gint->dm_map.fill(0.);
    }
    if (gint->do_sync) {
        gint->Is_map.SetNside(gint->nside_sim, NEST);
        gint->Qs_map.SetNside(gint->nside_sim, NEST);
        gint->Us_map.SetNside(gint->nside_sim, NEST);
        gint->Is_map.fill(0.);
        gint->Qs_map.fill(0.);
        gint->Us_map.fill(0.);
    }
    if (gint->do_fd or gint->do_sync) {
        gint->fd_map.SetNside(gint->nside_sim, NEST);
        gint->fd_map.fill(0.);
    }
    
    for (decltype(gint->total_shell) current_shell=1;current_shell!=(gint->total_shell+1);++current_shell) {
        cout<<"INTEGRATOR: AT SHELL "<<current_shell<<endl;
        Healpix_Map<double> current_Is_map;
        Healpix_Map<double> current_Qs_map;
        Healpix_Map<double> current_Us_map;
        Healpix_Map<double> current_fd_map;
        Healpix_Map<double> current_dm_map;
        
        unsigned int current_nside {gint->nside_shell[current_shell-1]};
        unsigned int current_npix{12*current_nside*current_nside};
        if (gint->do_dm) {
            current_dm_map.SetNside(current_nside, NEST);
            current_dm_map.fill(0.);
        }
        if (gint->do_sync) {
            current_Is_map.SetNside(current_nside, NEST);
            current_Qs_map.SetNside(current_nside, NEST);
            current_Us_map.SetNside(current_nside, NEST);
            current_Is_map.fill(0.);
            current_Qs_map.fill(0.);
            current_Us_map.fill(0.);
        }
        if (gint->do_fd or gint->do_sync) {
            current_fd_map.SetNside(current_nside, NEST);
            current_fd_map.fill(0.);
        }
        
        struct_shell shell_ref;
        // setting for radial_integration
        shell_ref.shell_num = current_shell;
        shell_ref.d_start = get_min_shell_radius(current_shell,gint->total_shell,gint->ec_r_max);
        shell_ref.d_stop = get_max_shell_radius(current_shell,gint->total_shell,gint->ec_r_max);
        shell_ref.delta_d = (shell_ref.d_stop-shell_ref.d_start)/gint->bin_num;
        /*
         * logic for bin number calculation for each shell
         * if only one shell, bin number is set by bin_num
         * if more than one shells, first inner shell set by bin_num
         * outer shells are set by maintaining the same bin length
         * for example the 2nd inner shell set by bin_num, while the 3rd inner shell set by 2*bin_num
         */
        if (current_shell>1) {
            shell_ref.delta_d = (shell_ref.d_stop-shell_ref.d_start)/(gint->bin_num*(pow(2,current_shell-1)-pow(2,current_shell-2)));
        }
#pragma omp parallel for
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
            radial_integration(shell_ref,ptg,observables,breg,brnd,fe,fernd,cre,gbreg,gbrnd,gfe,gfernd,gcre,gint,par);
            
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
#pragma omp parallel for
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


void Integrator::radial_integration(struct_shell &shell_ref,pointing &ptg_in, struct_observables &pixobs,Breg *breg,Brnd *brnd,FE *fe,FErnd *fernd,CRE *cre,Grid_breg *gbreg,Grid_brnd *gbrnd,Grid_fe *gfe,Grid_fernd *gfernd,Grid_cre *gcre,Grid_int *gint,Pond *par) {
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
    // convert intensity to brightness temperature and include j(omega) to j(nu)=2pi*j(omega)
    const double i2bt {CGS_U_pi*CGS_U_C_light*CGS_U_C_light/(CGS_U_kB*par->sim_freq*par->sim_freq)};
    
    // for Simpson's rule
    vector<double> F_dm;
    vector<double> F_fd;
    vector<double> F_Jtot, F_Jpol, intr_pol_ang;
    // first iteration in distance
    for(double dist=shell_ref.d_start;dist<=shell_ref.d_stop;dist+=shell_ref.delta_d/2){
        // ec and gc position
        vec3 ec_pos {toolkit::get_LOS_unit_vec(THE,PHI)*dist};
        vec3 pos {ec_pos + par->SunPosition};
        // check LOS depth limit
        if (check_simulation_upper_limit(pos.Length(),gint->gc_r_max)) {break;}
        if (check_simulation_upper_limit(fabs(pos.z),gint->gc_z_max)) {break;}
        if (check_simulation_upper_limit(ec_pos.Length(),gint->ec_r_max)) {break;}
        
        // B field
        vec3 B_vec {breg->get_breg(pos,par,gbreg)};
        // add random field
        B_vec += brnd->get_brnd(pos,gbrnd);
        
        const double B_par {toolkit::get_par2LOS(B_vec,THE,PHI)};
        // un-scaled B_per, if random B is active, we scale up this
        // in calculating emissivity, so it is not constant
        double B_per {toolkit::get_perp2LOS(B_vec,THE,PHI)};
        
        // FE field
        double te {fe->get_density(pos,par,gfe)};
        // add random field
        te += fernd->get_fernd(pos,gfernd);
        // to avoid negative value
        if(te<0){
            /*
             cerr<<"ERR:"<<__FILE__
             <<" : in function "<<__func__<<endl
             <<" at line "<<__LINE__<<endl
             <<" NEGATIVE FE DENSITY"<<endl;
             exit(1);
             */
            te = 0.;
        }
        
        // DM
        if(gint->do_dm){
            F_dm.push_back(te*shell_ref.delta_d);
        }
        // Faraday depth
        if(gint->do_fd or gint->do_sync){
            const double fd_forefactor {-(CGS_U_qe*CGS_U_qe*CGS_U_qe)/(2.*CGS_U_pi*CGS_U_MEC2*CGS_U_MEC2)};
            F_fd.push_back(te*B_par*fd_forefactor*shell_ref.delta_d);
        }
        
        // Synchrotron Emission
        if(gint->do_sync){
            F_Jtot.push_back(cre->get_emissivity(pos,par,gcre,B_per,1)*shell_ref.delta_d*i2bt);
            // J_pol receive no contribution from missing random
            F_Jpol.push_back(cre->get_emissivity(pos,par,gcre,B_per,0)*shell_ref.delta_d*i2bt);
            // intrinsic polarization angle, following IAU definition
            intr_pol_ang.push_back(toolkit::get_intr_pol_ang(B_vec,THE,PHI));
        }
    }// first iteration end
    
    // second iteration in distance
    // applying Simpson's rule
    // be careful with the logic for assigining simpson_size!
    int simpson_size {0};
    if(gint->do_dm){
        simpson_size = F_dm.size();
    }
    if(gint->do_fd or gint->do_sync){
        simpson_size = F_fd.size();
    }
    // if nothing is calculated
    if(simpson_size==0){
        return;
    }
    
    for(decltype(simpson_size)i=1;i<simpson_size-1;i+=2){
        // DM
        if(gint->do_dm){
            pixobs.dm += (F_dm[i-1]+4.*F_dm[i]+F_dm[i+1])/6.;
        }
        // FD
        if(gint->do_fd or gint->do_sync){
            pixobs.fd += (F_fd[i-1]+4.*F_fd[i]+F_fd[i+1])/6.;
        }
        // Sync
        if(gint->do_sync){
            // pol. angle after Faraday rotation
            double qui {(inner_shells_fd+pixobs.fd)*lambda_square+intr_pol_ang[i]};
#ifndef NDEBUG
            if (abs(qui)>1e30) {
                cerr<<"ERR: ODD VALUE"<<endl;
                exit(1);
            }
            if (F_Jtot[i]<0) {
                cerr<<"ERR: J_tot NEGATIVE"<<endl;
                exit(1);
            }
#endif
            pixobs.Is += (F_Jtot[i-1]+4.*F_Jtot[i]+F_Jtot[i+1])/6.;
            pixobs.Qs += cos(2.*qui)*(F_Jpol[i-1]+4.*F_Jpol[i]+F_Jpol[i+1])/6.;
            pixobs.Us += sin(2.*qui)*(F_Jpol[i-1]+4.*F_Jpol[i]+F_Jpol[i+1])/6.;
        }
    }
}//end of radial_integrate

// TOOLS
//---------------------------------------------------------

// compute upper radial boundary at given shell
double Integrator::get_max_shell_radius(const unsigned int &shell_numb,const unsigned int &total_shell,const double &radius) const {
    if (shell_numb<1 or shell_numb>(total_shell)) {
        cerr<<"INVALID shell_numb"<<endl;
        exit(1);
    }
    double max_shell_radius {radius};
    for (unsigned int n=total_shell;n!=shell_numb;--n) {
        max_shell_radius /= 2.;
    }
    return max_shell_radius;
}

// compute lower radial boundary at given shell
double Integrator::get_min_shell_radius(const unsigned int &shell_numb,const unsigned int &total_shell,const double &radius) const {
    if (shell_numb<1 or shell_numb>total_shell) {
        cerr<<"INVALID shell_numb"<<endl;
        exit(1);
    }
    // min_shell_radius for the first shell is always zero.
    if(shell_numb==1){
        return 0.;
    }
    double min_shell_radius {radius};
    for (unsigned int n=total_shell;n!=(shell_numb-1);--n) {
        min_shell_radius /= 2.;
    }
    return min_shell_radius;
}

//END ALL
