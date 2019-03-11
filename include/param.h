// parameters (fixed or free) for physical models

#ifndef HAMMURABI_PARAM_H
#define HAMMURABI_PARAM_H

#include <string>
#include <vector>

#include <tinyxml2.h>
#include <hvec.h>

class Param {
public:
    Param (const std::string);
	Param () = default;
    Param (const Param &) = delete;
    Param (Param &&) = delete;
    Param& operator= (const Param &) = delete;
    Param& operator= (Param &&) = delete;
    virtual ~Param () = default;
//----------------------- GRID PARAMETERS --------------------------------------
    // observer
    hvec<3,double> observer;
    // regular GMF grid
    struct param_breg_grid{
        std::string filename;
        bool build_permission=false, read_permission=false, write_permission=false;
        double x_max, x_min, y_max, y_min, z_max, z_min;
        std::size_t nx, ny, nz, full_size;
    }grid_breg;
    // turbulent GMF grid
    struct param_brnd_grid{
        std::string filename;
        bool read_permission=false, write_permission=false, build_permission=false;
        double x_max, x_min, y_max, y_min, z_max, z_min;
        std::size_t nx, ny, nz, full_size;
    }grid_brnd;
    // regular FE grid
    struct param_fereg_grid{
        std::string filename;
        bool build_permission=false, read_permission=false, write_permission=false;
        double x_max, x_min, y_max, y_min, z_max, z_min;
        std::size_t nx, ny, nz, full_size;
    }grid_fereg;
    // turbulent FE grid
    struct param_fernd_grid{
        std::string filename;
        bool read_permission=false, write_permission=false, build_permission=false;
        double x_max, x_min, y_max, y_min, z_max, z_min;
        std::size_t nx, ny, nz, full_size;
    }grid_fernd;
    // cre grid
    struct param_cre_grid{
        std::string filename;
        bool build_permission=false, read_permission=false, write_permission=false;
        std::size_t nE, nz, nx, ny, cre_size;
        double x_max, x_min, y_max, y_min, z_max, z_min;
        double E_min, E_max, E_fact;
    }grid_cre;
    // observable grid
    struct param_int_grid{
        bool write_permission=false;
        // shell parameters
        std::size_t nside_dm, nside_fd, total_shell;
        std::vector<std::size_t> nside_sync, nside_shell;
        // storing upper limit of shell raidus ratio to max radius, except the last shell
        std::vector<double> cut_shell;
        std::vector<double> radii_shell;
        // shell boundary
        double gc_r_min, gc_r_max, ec_r_min, ec_r_max, gc_z_min, gc_z_max;
        double radial_res;
        double lat_min, lat_max, lon_min, lon_max;
        // switches
        bool do_dm=false, do_fd=false;
        std::vector<bool> do_sync;
        std::string sim_fd_name, sim_dm_name;
        std::vector<std::string> sim_sync_name;
        std::vector<double> sim_sync_freq;
    }grid_int;
//----------------------- FIELD PARAMETERS -------------------------------------
    // GMF
    std::string breg_type, brnd_type, brnd_method;
    // wmap lsa
    struct param_breg_wmap{
        double b0;
        double psi0;
        double psi1;
        double chi0;
    }breg_wmap;
    // test
#ifndef NDEBUG
    struct param_breg_unif{
        double b0;
        double l0;
        double r;
    }breg_unif;
#endif
    // jaffe
    struct param_breg_jaffe{
        bool quadruple,bss;
        double disk_amp,disk_z0;
        double halo_amp,halo_z0;
        double r_inner,r_scale,r_peak; // radial profile
        // ring/bar
        bool ring,bar;
        double ring_amp,bar_amp;
        double ring_r,bar_a,bar_b,bar_phi0;
        // spiral arms
        unsigned arm_num;
        std::vector<double> arm_amp,arm_phi0;
        double arm_pitch,arm_r0,arm_z0;
        // arm compress
        double comp_r,comp_c,comp_d,comp_p;
    }breg_jaffe;
    // random seed
    std::size_t brnd_seed;
    // global
    struct param_brnd_global_es{
        double rms;
        double k0;
        double a0;
        double rho;
        double r0,z0;
    }brnd_es;
    // local
    struct param_brnd_local_mhd{
        double pa0,pf0,ps0;
        double aa0,af0,as0;
        double k0;
        double ma,beta;
    }brnd_mhd;
    // FE
    std::string fereg_type, fernd_type, fernd_method;
    // ymw16
    struct param_fereg_ymw16{
        double r_warp, r0;
        double t0_gamma_w;
        double t1_ad, t1_bd, t1_n1, t1_h1;
        double t2_a2, t2_b2, t2_n2, t2_k2;
        double t3_b2s, t3_ka, t3_narm[5], t3_warm[5], t3_aa, t3_ncn, t3_wcn,
        t3_thetacn, t3_nsg, t3_wsg, t3_thetasg, t3_rmin[5], t3_phimin[5],
        t3_tpitch[5], t3_cpitch[5];
        double t4_ngc, t4_agc, t4_hgc;
        double t5_kgn, t5_ngn, t5_wgn, t5_agn;
        double t6_j_lb, t6_nlb1, t6_detlb1, t6_wlb1, t6_hlb1, t6_thetalb1,
        t6_nlb2, t6_detlb2, t6_wlb2, t6_hlb2, t6_thetalb2;
        double t7_nli, t7_rli, t7_wli, t7_detthetali, t7_thetali;
    }fereg_ymw16;
    // test
#ifndef NDEBUG
    struct param_fereg_unif{
        double n0;
        double r0;
    }fereg_unif;
#endif
    // random seed
    std::size_t fernd_seed;
    // isotropic
    struct param_fernd_global_dft{
        double rms;
        double k0;
        double a0;
        double r0;
        double z0;
    }fernd_dft;
    // CRE
    std::string cre_type;
    // analytical
    struct param_cre_ana{
        double alpha,beta,theta;
        double r0,z0;
        double E0,j0;
    }cre_ana;
    // test
#ifndef NDEBUG
    struct param_cre_unif{
        double alpha;
        double r0;
        double E0,j0;
    }cre_unif;
#endif
protected:
    // collect observable related parameters
    void obs_param (tinyxml2::XMLDocument *);
    // collect magnetic field related parameters
    void breg_param (tinyxml2::XMLDocument *);
    void brnd_param (tinyxml2::XMLDocument *);
    // collect free/thermal electron related parameters
    void fereg_param (tinyxml2::XMLDocument *);
    void fernd_param (tinyxml2::XMLDocument *);
    // collect CRE related parameters
    void cre_param (tinyxml2::XMLDocument *);
};

#endif

// END
