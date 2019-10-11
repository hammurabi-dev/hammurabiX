// parameters

#ifndef HAMMURABI_PARAM_H
#define HAMMURABI_PARAM_H

#include <string>
#include <vector>

#include <hamvec.h>
#include <tinyxml2.h>

class Param {
public:
  Param(const std::string);
  Param() = default;
  Param(const Param &) = delete;
  Param(Param &&) = delete;
  Param &operator=(const Param &) = delete;
  Param &operator=(Param &&) = delete;
  virtual ~Param() = default;
  // galactic centric Cartesian position of the observer
  hamvec<3, double> observer;
  // regular magnetic field grid
  struct param_breg_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool build_permission = false, read_permission = false,
         write_permission = false;
    // galactic centric Cartesian limit
    double x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    std::size_t nx, ny, nz, full_size;
  } grid_breg;
  // random magnetic field grid
  struct param_brnd_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool read_permission = false, write_permission = false,
         build_permission = false;
    // galactic centric Cartesian limit
    double x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    std::size_t nx, ny, nz, full_size;
  } grid_brnd;
  // regular thermal electron grid
  struct param_tereg_grid {
    std::string filename;
    // grid build/read/write controller
    bool build_permission = false, read_permission = false,
         write_permission = false;
    // galactic centric Cartesian limit
    double x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    std::size_t nx, ny, nz, full_size;
  } grid_tereg;
  // random thermal electron grid
  struct param_ternd_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool read_permission = false, write_permission = false,
         build_permission = false;
    // galactic centric Cartesian limit
    double x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    std::size_t nx, ny, nz, full_size;
  } grid_ternd;
  // cosmic ray electron grid
  struct param_cre_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool build_permission = false, read_permission = false,
         write_permission = false;
    // number Cartesian grid support points
    std::size_t nE, nz, nx, ny, cre_size;
    // galactic centric Cartesian limit
    double x_max, x_min, y_max, y_min, z_max, z_min;
    // cosmic ray electron logarithmic space energy limit
    // E_fact is used for assigning energy bin length
    double E_min, E_max, E_fact;
  } grid_cre;
  // observable grid
  struct param_obs_grid {
    // grid write controller
    bool write_permission = false;
    // shell resolution controllers
    std::size_t nside_dm, nside_fd, total_shell;
    // multi-frequency synchrotron shell resolution controllers
    std::vector<std::size_t> nside_sync, nside_shell;
    // nested shell radius controllers
    std::vector<double> cut_shell;
    std::vector<double> radii_shell;
    // simulation boundary
    // galactic centric cylindrical radius limit
    double gc_r_min, gc_r_max;
    // galactic centric cylindrical height limit
    double gc_z_min, gc_z_max;
    // observer centric spherical radius limit
    double oc_r_min, oc_r_max;
    // LoS integration radial resolution
    double oc_r_res;
    // simulation controllers
    bool do_dm = false, do_fd = false;
    std::vector<bool> do_sync;
    // in/output file name
    std::string sim_fd_name, sim_dm_name;
    std::vector<std::string> sim_sync_name;
    // synchrotron frequencies
    std::vector<double> sim_sync_freq;
    // mask controllers
    bool do_mask = false;
    std::string mask_name;
  } grid_obs;
  // magnetic field parameters
  std::string breg_type, brnd_type, brnd_method;
  // LSA model parameters
  struct param_breg_lsa {
    double b0;
    double psi0;
    double psi1;
    double chi0;
  } breg_lsa;
  // uniform model parameters
  struct param_breg_unif {
    double bp, bv;
    double l0;
  } breg_unif;
  // Jaffe model parameters
  struct param_breg_jaffe {
    bool quadruple, bss;
    double disk_amp, disk_z0;
    double halo_amp, halo_z0;
    double r_inner, r_scale, r_peak; // radial profile
    // ring/bar
    bool ring, bar;
    double ring_amp, bar_amp;
    double ring_r, bar_a, bar_b, bar_phi0;
    // spiral arms
    unsigned arm_num;
    std::vector<double> arm_amp, arm_phi0;
    double arm_pitch, arm_r0, arm_z0;
    // arm compress
    double comp_r, comp_c, comp_d, comp_p;
  } breg_jaffe;
  // random magnetic field generation seed
  std::size_t brnd_seed;
  // global ES model parameters
  struct param_brnd_global_es {
    double rms;
    double k0, k1;
    double a0, a1;
    double rho;
    double r0, z0;
  } brnd_es;
  // local MHD model parameters
  struct param_brnd_local_mhd {
    double pa0, pf0, ps0;
    double aa0, af0, as0, a1;
    double k0, k1;
    double ma, beta;
  } brnd_mhd;
  // thermal electron
  std::string tereg_type, ternd_type, ternd_method;
  // ymw16 model parameters
  struct param_tereg_ymw16 {
    double r_warp, r0;
    double t0_gamma_w;
    double t1_ad, t1_bd, t1_n1, t1_h1;
    double t2_a2, t2_b2, t2_n2, t2_k2;
    double t3_b2s, t3_ka, t3_narm[5], t3_warm[5], t3_aa, t3_ncn, t3_wcn,
        t3_thetacn, t3_nsg, t3_wsg, t3_thetasg, t3_rmin[5], t3_phimin[5],
        t3_tpitch[5], t3_cpitch[5];
    double t4_ngc, t4_agc, t4_hgc;
    double t5_kgn, t5_ngn, t5_wgn, t5_agn;
    double t6_j_lb, t6_nlb1, t6_detlb1, t6_wlb1, t6_hlb1, t6_thetalb1, t6_nlb2,
        t6_detlb2, t6_wlb2, t6_hlb2, t6_thetalb2;
    double t7_nli, t7_rli, t7_wli, t7_detthetali, t7_thetali;
  } tereg_ymw16;
  // uniform model parameters
  struct param_tereg_unif {
    double n0;
    double r0;
  } tereg_unif;
  // random thermal electron generation seed
  std::size_t ternd_seed;
  // isotropic default model parameters
  struct param_ternd_global_dft {
    double rms;
    double k0;
    double a0;
    double r0;
    double z0;
  } ternd_dft;
  // cosmic ray electron
  std::string cre_type;
  // analytical model parameters
  struct param_cre_ana {
    double alpha, beta, theta;
    double r0, z0;
    double E0, j0;
  } cre_ana;
  // uniform model parameters
  struct param_cre_unif {
    double alpha;
    double r0;
    double E0, j0;
  } cre_unif;

protected:
  // collect observable related parameters
  void obs_param(tinyxml2::XMLDocument *);
  // collect magnetic field related parameters
  void breg_param(tinyxml2::XMLDocument *);
  void brnd_param(tinyxml2::XMLDocument *);
  // collect thermal electron related parameters
  void tereg_param(tinyxml2::XMLDocument *);
  void ternd_param(tinyxml2::XMLDocument *);
  // collect cosmic ray electron related parameters
  void cre_param(tinyxml2::XMLDocument *);
};

#endif
