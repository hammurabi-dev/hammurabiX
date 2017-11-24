#include <iostream>
#include <vec3.h>
#include <vector>
#include <array>
#include <cmath>

#include "pond.h"
#include "grid.h"
#include "breg.h"
#include "cgs_units_file.h"
#include "namespace_toolkit.h"

using namespace std;

vec3_t<double> Breg_jaffe::breg(const vec3_t<double> &pos,Pond *par){
    double inner_b {0};
    if(par->breg_jaffe.ring)
        inner_b = par->breg_jaffe.ring_amp;
    else if(par->breg_jaffe.bar)
        inner_b = par->breg_jaffe.bar_amp;
    
    vec3_t<double> bhat {versor(pos,par)};
    vec3_t<double> btot {0,0,0};
    btot = bhat*radial_scaling(pos,par)*(par->breg_jaffe.disk_amp*disk_scaling(pos,par) + par->breg_jaffe.halo_amp*halo_scaling(pos,par));
    // compress factor for each arm or for ring/bar
    vector<double> arm {arm_compress(pos,par)};
    // only inner region
    if(arm.size()==1){
        btot += bhat*arm[0]*inner_b;
    }
    // spiral arm region
    else{
        for(decltype(arm.size())i=0;i<arm.size();++i){
            btot += bhat*arm[i]*par->breg_jaffe.arm_amp[i];
        }
    }
    return btot;
}

vec3_t<double> Breg_jaffe::versor(const vec3_t<double> &pos,Pond *par){
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)}; // cylindrical frame
    const double r_lim {par->breg_jaffe.ring_r};
    const double bar_lim {par->breg_jaffe.bar_a + 0.5*par->breg_jaffe.comp_d};
    const double cos_p {cos(par->breg_jaffe.arm_pitch)}; const double sin_p {sin(par->breg_jaffe.arm_pitch)}; // pitch angle
    vec3_t<double> tmp {0,0,0}; double quadruple {1};
    if(r<0.5*CGS_U_kpc) // forbiden region
        return tmp;
    if(pos.z>par->breg_jaffe.disk_z0)
        quadruple = (1-2*par->breg_jaffe.quadruple);
    // molecular ring
    if(par->breg_jaffe.ring){
        // inside spiral arm
        if(r>r_lim){
            tmp.x = (cos_p*(pos.y/r)-sin_p*(pos.x/r))*quadruple; //sin(t-p)
            tmp.y = (-cos_p*(pos.x/r)-sin_p*(pos.y/r))*quadruple; //-cos(t-p)
        }
        // inside molecular ring
        else{
            tmp.x = (1-2*par->breg_jaffe.bss)*pos.y/r; //sin(phi)
            tmp.y = (2*par->breg_jaffe.bss-1)*pos.x/r; //-cos(phi)
        }
    }
    // elliptical bar (replace molecular ring)
    else if(par->breg_jaffe.bar){
        const double cos_phi {cos(par->breg_jaffe.bar_phi0)};
        const double sin_phi {sin(par->breg_jaffe.bar_phi0)};
        const double x {cos_phi*pos.x-sin_phi*pos.y};
        const double y {sin_phi*pos.x+cos_phi*pos.y};
        // inside spiral arm
        if(r>bar_lim){
            tmp.x = (cos_p*(pos.y/r)-sin_p*(pos.x/r))*quadruple; //sin(t-p)
            tmp.y = (-cos_p*(pos.x/r)-sin_p*(pos.y/r))*quadruple; //-cos(t-p)
        }
        // inside elliptical bar
        else{
            if(y!=0){
                const double new_x {copysign(1,y)};
                const double new_y {-copysign(1,y)*(x/y)*par->breg_jaffe.bar_b*par->breg_jaffe.bar_b/(par->breg_jaffe.bar_a*par->breg_jaffe.bar_a)};
                tmp.x = (cos_phi*new_x+sin_phi*new_y)*(1-2*par->breg_jaffe.bss);
                tmp.y = (-sin_phi*new_x+cos_phi*new_y)*(1-2*par->breg_jaffe.bss);
                tmp.Normalize();
            }
            else{
                tmp.x = (2*par->breg_jaffe.bss-1)*copysign(1,x)*sin_phi;
                tmp.y = (2*par->breg_jaffe.bss-1)*copysign(1,x)*cos_phi;
            }
        }
    }
    return tmp;
}

double Breg_jaffe::radial_scaling(const vec3_t<double> &pos,Pond *par){
    const double r2 {pos.x*pos.x+pos.y*pos.y};
    // separate into 3 parts for better view
    const double s1 {1.-exp(-r2/(par->breg_jaffe.r_inner*par->breg_jaffe.r_inner))};
    const double s2 {exp(-r2/(par->breg_jaffe.r_scale*par->breg_jaffe.r_scale))};
    const double s3 {exp(-r2*r2/(par->breg_jaffe.r_peak*par->breg_jaffe.r_peak*par->breg_jaffe.r_peak*par->breg_jaffe.r_peak))};
    return s1*( s2 + s3 );
}

vector<double> Breg_jaffe::arm_compress(const vec3_t<double> &pos,Pond *par){
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)/par->breg_jaffe.comp_r};
    const double c0 {1./par->breg_jaffe.comp_c -1.};
    vector<double> a0 = dist2arm(pos,par);
    const double r_scaling {radial_scaling(pos,par)};
    const double z_scaling {arm_scaling(pos,par)};
    // for saving computing time
    const double d0_inv {(r_scaling*z_scaling)/par->breg_jaffe.comp_d};
    double factor {c0*r_scaling*z_scaling};
    if(r>1){
        double cdrop {pow(r,-par->breg_jaffe.comp_p)};
        for(decltype(a0.size())i=0;i<a0.size();++i){
            a0[i] = factor*cdrop*exp(-a0[i]*a0[i]*cdrop*cdrop*d0_inv*d0_inv);
        }
    }
    else{
        for(decltype(a0.size())i=0;i<a0.size();++i){
            a0[i] = factor*exp(-a0[i]*a0[i]*d0_inv*d0_inv);
        }
    }
    return a0;
}

vector<double> Breg_jaffe::arm_compress_dust(const vec3_t<double> &pos,Pond *par){
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)/par->breg_jaffe.comp_r};
    const double c0 {1./par->breg_jaffe.comp_c -1.};
    vector<double> a0 = dist2arm(pos,par);
    const double r_scaling {radial_scaling(pos,par)};
    const double z_scaling {arm_scaling(pos,par)};
    // only difference from normal arm_compress
    const double d0_inv {(r_scaling)/par->breg_jaffe.comp_d};
    double factor {c0*r_scaling*z_scaling};
    if(r>1){
        double cdrop {pow(r,-par->breg_jaffe.comp_p)};
        for(decltype(a0.size())i=0;i<a0.size();++i){
            a0[i] = factor*cdrop*exp(-a0[i]*a0[i]*cdrop*cdrop*d0_inv*d0_inv);
        }
    }
    else{
        for(decltype(a0.size())i=0;i<a0.size();++i){
            a0[i] = factor*exp(-a0[i]*a0[i]*d0_inv*d0_inv);
        }
    }
    return a0;
}

vector<double> Breg_jaffe::dist2arm(const vec3_t<double> &pos,Pond *par){
    vector<double> d;
    const double r {sqrt(pos.x*pos.x+pos.y*pos.y)};
    const double r_lim {par->breg_jaffe.ring_r};
    const double bar_lim {par->breg_jaffe.bar_a + 0.5*par->breg_jaffe.comp_d};
    const double cos_p {cos(par->breg_jaffe.arm_pitch)}; const double sin_p {sin(par->breg_jaffe.arm_pitch)}; // pitch angle
    const double beta_inv {-sin_p/cos_p};
    double theta {atan2(pos.y,pos.x)};
    if(theta<0) theta += 2*CGS_U_pi;
    // if molecular ring
    if(par->breg_jaffe.ring){
        // in molecular ring, return single element vector
        if(r<r_lim){
            d.push_back(fabs(par->breg_jaffe.ring_r-r));
        }
        // in spiral arm, return vector with arm_num elements
        else{
            // loop through arms
            for(unsigned int i=0;i<par->breg_jaffe.arm_num;++i){
                double d_ang {par->breg_jaffe.arm_phi0[i]-theta};
                double d_rad {fabs(par->breg_jaffe.arm_r0*exp(d_ang*beta_inv)-r)};
                double d_rad_p {fabs(par->breg_jaffe.arm_r0*exp((d_ang+2*CGS_U_pi)*beta_inv)-r)};
                double d_rad_m {fabs(par->breg_jaffe.arm_r0*exp((d_ang-2*CGS_U_pi)*beta_inv)-r)};
                d.push_back(min(min(d_rad,d_rad_p),d_rad_m)*cos_p);
            }
        }
    }
    // if elliptical bar
    else if(par->breg_jaffe.bar){
        const double cos_tmp {cos(par->breg_jaffe.bar_phi0)*pos.x/r - sin(par->breg_jaffe.bar_phi0)*pos.y/r}; // cos(phi)cos(phi0) - sin(phi)sin(phi0)
        const double sin_tmp {cos(par->breg_jaffe.bar_phi0)*pos.y/r + sin(par->breg_jaffe.bar_phi0)*pos.x/r}; // sin(phi)cos(phi0) + cos(phi)sin(phi0)
        // in bar, return single element vector
        if(r<bar_lim){
            d.push_back(fabs(par->breg_jaffe.bar_a*par->breg_jaffe.bar_b/sqrt(par->breg_jaffe.bar_a*par->breg_jaffe.bar_a*sin_tmp*sin_tmp+par->breg_jaffe.bar_b*par->breg_jaffe.bar_b*cos_tmp*cos_tmp)-r));
        }
        // in spiral arm, return vector with arm_num elements
        else{
            // loop through arms
            for(unsigned int i=0;i<par->breg_jaffe.arm_num;++i){
                double d_ang {par->breg_jaffe.arm_phi0[i]-theta};
                double d_rad {fabs(par->breg_jaffe.arm_r0*exp(d_ang*beta_inv)-r)};
                double d_rad_p {fabs(par->breg_jaffe.arm_r0*exp((d_ang+2*CGS_U_pi)*beta_inv)-r)};
                double d_rad_m {fabs(par->breg_jaffe.arm_r0*exp((d_ang-2*CGS_U_pi)*beta_inv)-r)};
                d.push_back(min(min(d_rad,d_rad_p),d_rad_m)*cos_p);
            }
        }
    }
    return d;
}

// END
