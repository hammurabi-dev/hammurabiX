// parameters

#ifndef HAMMURABI_PARAM_H
#define HAMMURABI_PARAM_H

#include <string>
#include <vector>

#include <hamtype.h>
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
  Hamvec<3, ham_float> observer;
  // regular magnetic field grid
  struct param_breg_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool build_permission = false, read_permission = false,
         write_permission = false;
    // galactic centric Cartesian limit
    ham_float x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    ham_uint nx, ny, nz, full_size;
  } grid_breg;
  // random magnetic field grid
  struct param_brnd_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool read_permission = false, write_permission = false,
         build_permission = false;
    // galactic centric Cartesian limit
    ham_float x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    ham_uint nx, ny, nz, full_size;
  } grid_brnd;
  // regular thermal electron grid
  struct param_tereg_grid {
    std::string filename;
    // grid build/read/write controller
    bool build_permission = false, read_permission = false,
         write_permission = false;
    // galactic centric Cartesian limit
    ham_float x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    ham_uint nx, ny, nz, full_size;
  } grid_tereg;
  // random thermal electron grid
  struct param_ternd_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool read_permission = false, write_permission = false,
         build_permission = false;
    // galactic centric Cartesian limit
    ham_float x_max, x_min, y_max, y_min, z_max, z_min;
    // number Cartesian grid support points
    ham_uint nx, ny, nz, full_size;
  } grid_ternd;
  // cosmic ray electron grid
  struct param_cre_grid {
    // in/output file name
    std::string filename;
    // grid build/read/write controller
    bool build_permission = false, read_permission = false,
         write_permission = false;
    // number Cartesian grid support points
    ham_uint nE, nz, nx, ny, cre_size;
    // galactic centric Cartesian limit
    ham_float x_max, x_min, y_max, y_min, z_max, z_min;
    // cosmic ray electron logarithmic space energy limit
    // E_fact is used for assigning energy bin length
    ham_float E_min, E_max, E_fact;
  } grid_cre;
  // observable grid
  struct param_obs_grid {
    // grid write controller
    bool write_permission = false;
    // shell resolution controllers
    ham_uint nside_dm, nside_fd, nside_mask, total_shell;
    // multi-frequency synchrotron shell resolution controllers
    std::vector<ham_uint> nside_sync, nside_shell;
    // nested shell radius controllers
    std::vector<ham_float> cut_shell;
    std::vector<ham_float> radii_shell;
    // simulation boundary
    // galactic centric cylindrical radius limit
    ham_float gc_r_min, gc_r_max;
    // galactic centric cylindrical height limit
    ham_float gc_z_min, gc_z_max;
    // observer centric spherical radius limit
    ham_float oc_r_min, oc_r_max;
    // LoS integration radial resolution
    ham_float oc_r_res;
    // simulation controllers
    bool do_dm = false, do_fd = false;
    std::vector<bool> do_sync;
    // in/output file name
    std::string sim_fd_name, sim_dm_name;
    std::vector<std::string> sim_sync_name;
    // synchrotron frequencies
    std::vector<ham_float> sim_sync_freq;
    // mask controllers
    bool do_mask = false;
    std::string mask_name;
  } grid_obs;
  // magnetic field parameters
  std::string breg_type, brnd_type, brnd_method;
  // LSA model parameters
  struct param_breg_lsa {
    ham_float b0;
    ham_float psi0;
    ham_float psi1;
    ham_float chi0;
  } breg_lsa;
  // uniform model parameters
  struct param_breg_unif {
    ham_float bp, bv;
    ham_float l0;
  } breg_unif;
  // Jaffe model parameters
  struct param_breg_jaffe {
    bool quadruple, bss;
    ham_float disk_amp, disk_z0;
    ham_float halo_amp, halo_z0;
    ham_float r_inner, r_scale, r_peak; // radial profile
    // ring/bar
    bool ring, bar;
    ham_float ring_amp, bar_amp;
    ham_float ring_r, bar_a, bar_b, bar_phi0;
    // spiral arms
    ham_uint arm_num;
    std::vector<ham_float> arm_amp, arm_phi0;
    ham_float arm_pitch, arm_r0, arm_z0;
    // arm compress
    ham_float comp_r, comp_c, comp_d, comp_p;
  } breg_jaffe;
  struct param_breg_cart {
    ham_float bx, by;
    ham_float bz;
  } breg_cart;
  struct param_breg_helix {
    ham_float bx, by, bz;
    ham_float r_min, r_max;
  } breg_helix;
  struct param_breg_jf12 {
    ham_float b_arm1, b_arm2, b_arm3, b_arm4, b_arm5, b_arm6, b_arm7;
    ham_float b_ring, h_disk, w_disk;
    ham_float Bn, Bs, rn, rs, wh, z0;
    ham_float B0_X, Xtheta, rpc_X, r0_X;
  } breg_jf12;

  // random magnetic field generation seed
  ham_uint brnd_seed;
  // global ES model parameters
  struct param_brnd_global_es {
    ham_float rms;
    ham_float k0, k1;
    ham_float a0, a1;
    ham_float rho;
    ham_float r0, z0;
  } brnd_es;
  struct param_brnd_global_jf12 {
    ham_float rms;
    ham_float k0, k1;
    ham_float a0, a1;
    ham_float rho;
    ham_float b0_1, b0_2, b0_3, b0_4, b0_5, b0_6, b0_7, b0_8;
    ham_float b0_int, z0_spiral;
    ham_float b0_halo, r0_halo, z0_halo;
  } brnd_jf12;
  // local MHD model parameters
  struct param_brnd_local_mhd {
    ham_float pa0, pf0, ps0;
    ham_float aa0, af0, as0, a1;
    ham_float k0, k1;
    ham_float ma, beta;
  } brnd_mhd;
  // thermal electron
  std::string tereg_type, ternd_type, ternd_method;
  // ymw16 model parameters
  struct param_tereg_ymw16 {
    ham_float r_warp, r0;
    ham_float t0_gamma_w;
    ham_float t1_ad, t1_bd, t1_n1, t1_h1;
    ham_float t2_a2, t2_b2, t2_n2, t2_k2;
    ham_float t3_b2s, t3_ka, t3_narm[5], t3_warm[5], t3_aa, t3_ncn, t3_wcn,
        t3_thetacn, t3_nsg, t3_wsg, t3_thetasg, t3_rmin[5], t3_phimin[5],
        t3_tpitch[5], t3_cpitch[5];
    ham_float t4_ngc, t4_agc, t4_hgc;
    ham_float t5_kgn, t5_ngn, t5_wgn, t5_agn;
    ham_float t6_j_lb, t6_nlb1, t6_detlb1, t6_wlb1, t6_hlb1, t6_thetalb1,
        t6_nlb2, t6_detlb2, t6_wlb2, t6_hlb2, t6_thetalb2;
    ham_float t7_nli, t7_rli, t7_wli, t7_detthetali, t7_thetali;
  } tereg_ymw16;
  // uniform model parameters
  struct param_tereg_unif {
    ham_float n0;
    ham_float r0;
  } tereg_unif;
  // random thermal electron generation seed
  ham_uint ternd_seed;
  // isotropic default model parameters
  struct param_ternd_global_dft {
    ham_float rms;
    ham_float k0;
    ham_float a0;
    ham_float r0;
    ham_float z0;
  } ternd_dft;
  // cosmic ray electron
  std::string cre_type;
  // analytical model parameters
  struct param_cre_ana {
    ham_float alpha, beta, theta;
    ham_float r0, z0;
    ham_float E0, j0;
  } cre_ana;
  // uniform model parameters
  struct param_cre_unif {
    ham_float alpha;
    ham_float r0;
    ham_float E0, j0;
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
