#include <cmath>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <memory>
#include <param.h>
#include <string>
#include <tinyxml2.h>
#include <toolkit.h>

Param::Param(const std::string file_name) {
  // load xml file
  std::unique_ptr<tinyxml2::XMLDocument> doc{toolkit::loadxml(file_name)};
  // parse parameters, function call order is critical
  parse_obs(doc.get());
  parse_field(doc.get());
  parse_bmodel(doc.get());
  parse_temodel(doc.get());
  parse_cremodel(doc.get());
  parse_grid(doc.get());
}

void Param::parse_obs(tinyxml2::XMLDocument *doc) {
  // observable base path
  const std::string binary_suffix(".bin");
  const std::string fits_suffix(".fits");
  tinyxml2::XMLElement *ele{toolkit::tracexml(doc, {"obsio"})};
  // controller for reading shell parameters
  grid_obs.write_permission = false;
  // if dispersion measure is required
  if (ele->FirstChildElement("dm") != nullptr) {
    grid_obs.write_permission = true;
    grid_obs.do_dm = true;
    grid_obs.sim_dm_name = toolkit::fetchstring(ele, "dm", "filename");
      if (grid_obs.sim_dm_name.find(binary_suffix) == std::string::npos and grid_obs.sim_dm_name.find(fits_suffix) == std::string::npos)
        throw std::runtime_error("wrong suffix");
    grid_obs.nside_dm = toolkit::fetchuint(ele, "dm", "nside");
  } else {
    grid_obs.do_dm = false;
  }
  // if faraday depth is required
  if (ele->FirstChildElement("faraday") != nullptr) {
    grid_obs.write_permission = true;
    grid_obs.do_fd = true;
    grid_obs.sim_fd_name = toolkit::fetchstring(ele, "faraday", "filename");
      if (grid_obs.sim_fd_name.find(binary_suffix) == std::string::npos and grid_obs.sim_fd_name.find(fits_suffix) == std::string::npos)
          throw std::runtime_error("wrong suffix");
    grid_obs.nside_fd = toolkit::fetchuint(ele, "faraday", "nside");
  } else {
    grid_obs.do_fd = false;
  }
  // if synchrotron emission is required
  if (ele->FirstChildElement("sync") != nullptr) {
    grid_obs.write_permission = true;
    tinyxml2::XMLElement *subele{ele->FirstChildElement("sync")};
    grid_obs.do_sync.push_back(true);
    grid_obs.sim_sync_freq.push_back(toolkit::fetchfloat(subele, "freq") *
                                     cgs::GHz);
    grid_obs.sim_sync_name.push_back(toolkit::fetchstring(subele, "filename"));
      if (grid_obs.sim_sync_name.back().find(binary_suffix) == std::string::npos and grid_obs.sim_sync_name.back().find(fits_suffix) == std::string::npos)
          throw std::runtime_error("wrong suffix");
    grid_obs.nside_sync.push_back(toolkit::fetchuint(subele, "nside"));
    for (auto e = subele->NextSiblingElement("sync"); e != nullptr;
         e = e->NextSiblingElement("sync")) {
      grid_obs.do_sync.push_back(true);
      grid_obs.sim_sync_freq.push_back(toolkit::fetchfloat(e, "freq") *
                                       cgs::GHz);
      grid_obs.sim_sync_name.push_back(toolkit::fetchstring(e, "filename"));
        if (grid_obs.sim_sync_name.back().find(binary_suffix) == std::string::npos and grid_obs.sim_sync_name.back().find(fits_suffix) == std::string::npos)
            throw std::runtime_error("wrong suffix");
      grid_obs.nside_sync.push_back(toolkit::fetchuint(e, "nside"));
    }
  } else {
    grid_obs.do_sync.push_back(false);
  }
  // if mask is required
  if (ele->FirstChildElement("mask") != nullptr) {
    grid_obs.do_mask = true;
    grid_obs.mask_name = toolkit::fetchstring(ele, "mask", "filename");
      if (grid_obs.mask_name.find(binary_suffix) == std::string::npos and grid_obs.mask_name.find(fits_suffix) == std::string::npos)
          throw std::runtime_error("wrong suffix");
    grid_obs.nside_mask = toolkit::fetchuint(ele, "mask", "nside");
  } else {
    grid_obs.do_mask = false;
  }
}

void Param::parse_field(tinyxml2::XMLDocument *doc) {
    const std::string binary_suffix(".bin");
  // bfield io
  tinyxml2::XMLElement *ele{toolkit::tracexml(doc, {"fieldio"})};
  if (ele->FirstChildElement("bfield") != nullptr) {
    grid_b.read_permission = toolkit::fetchbool(ele, "bfield", "read");
    grid_b.write_permission = toolkit::fetchbool(ele, "bfield", "write");
    grid_b.filename = toolkit::fetchstring(ele, "bfield", "filename");
      if (grid_b.filename.find(binary_suffix) == std::string::npos)
          throw std::runtime_error("wrong suffix");
  }
  if (ele->FirstChildElement("tefield") != nullptr) {
    grid_te.read_permission = toolkit::fetchbool(ele, "tefield", "read");
    grid_te.write_permission = toolkit::fetchbool(ele, "tefield", "write");
    grid_te.filename = toolkit::fetchstring(ele, "tefield", "filename");
      if (grid_te.filename.find(binary_suffix) == std::string::npos)
          throw std::runtime_error("wrong suffix");
  }
  if (ele->FirstChildElement("crefield") != nullptr) {
    grid_cre.read_permission = toolkit::fetchbool(ele, "crefield", "read");
    grid_cre.write_permission = toolkit::fetchbool(ele, "crefield", "write");
    grid_cre.filename = toolkit::fetchstring(ele, "crefield", "filename");
      if (grid_cre.filename.find(binary_suffix) == std::string::npos)
          throw std::runtime_error("wrong suffix");
  }
}

void Param::parse_bmodel(tinyxml2::XMLDocument *doc) {
  // bfield model selection
  tinyxml2::XMLElement *ele{toolkit::tracexml(doc, {"bmodel"})};
  // uniform model
  tinyxml2::XMLElement *subele{ele->FirstChildElement("unif")};
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    bmodel_list.push_back("unif");
    bmodel_unif.bp = toolkit::fetchfloat(subele, "bp", "value") * cgs::muGauss;
    bmodel_unif.bv = toolkit::fetchfloat(subele, "bv", "value") * cgs::muGauss;
    bmodel_unif.l0 = toolkit::fetchfloat(subele, "l0", "value") * cgs::rad;
  }
  // LSA model
  subele = ele->FirstChildElement("lsa");
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    bmodel_list.push_back("lsa");
    bmodel_lsa.b0 = toolkit::fetchfloat(subele, "b0", "value") * cgs::muGauss;
    bmodel_lsa.psi0 = toolkit::fetchfloat(subele, "psi0", "value") * cgs::rad;
    bmodel_lsa.psi1 = toolkit::fetchfloat(subele, "psi1", "value") * cgs::rad;
    bmodel_lsa.chi0 = toolkit::fetchfloat(subele, "chi0", "value") * cgs::rad;
  }
  // Jaffe model
  subele = ele->FirstChildElement("jaffe");
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    bmodel_list.push_back("jaffe");
    bmodel_jaffe.quadruple = toolkit::fetchbool(subele, "quadruple", "cue");
    bmodel_jaffe.bss = toolkit::fetchbool(subele, "bss", "cue");
    bmodel_jaffe.disk_amp =
        toolkit::fetchfloat(subele, "disk_amp", "value") * cgs::muGauss;
    bmodel_jaffe.disk_z0 =
        toolkit::fetchfloat(subele, "disk_z0", "value") * cgs::kpc;
    bmodel_jaffe.halo_amp =
        toolkit::fetchfloat(subele, "halo_amp", "value") * cgs::muGauss;
    bmodel_jaffe.halo_z0 =
        toolkit::fetchfloat(subele, "halo_z0", "value") * cgs::kpc;
    bmodel_jaffe.r_inner =
        toolkit::fetchfloat(subele, "r_inner", "value") * cgs::kpc;
    bmodel_jaffe.r_scale =
        toolkit::fetchfloat(subele, "r_scale", "value") * cgs::kpc;
    bmodel_jaffe.r_peak =
        toolkit::fetchfloat(subele, "r_peak", "value") * cgs::kpc;
    bmodel_jaffe.ring = toolkit::fetchbool(subele, "ring", "cue");
    bmodel_jaffe.ring_amp =
        toolkit::fetchfloat(subele, "ring_amp", "value") * cgs::muGauss;
    bmodel_jaffe.ring_r =
        toolkit::fetchfloat(subele, "ring_r", "value") * cgs::kpc;
    bmodel_jaffe.bar = toolkit::fetchbool(subele, "bar", "cue");
    bmodel_jaffe.bar_amp =
        toolkit::fetchfloat(subele, "bar_amp", "value") * cgs::muGauss;
    bmodel_jaffe.bar_a =
        toolkit::fetchfloat(subele, "bar_a", "value") * cgs::kpc;
    bmodel_jaffe.bar_b =
        toolkit::fetchfloat(subele, "bar_b", "value") * cgs::kpc;
    bmodel_jaffe.bar_phi0 =
        toolkit::fetchfloat(subele, "bar_phi0", "value") * cgs::rad;
    bmodel_jaffe.arm_num = toolkit::fetchuint(subele, "arm_num", "value");
    bmodel_jaffe.arm_r0 =
        toolkit::fetchfloat(subele, "arm_r0", "value") * cgs::kpc;
    bmodel_jaffe.arm_z0 =
        toolkit::fetchfloat(subele, "arm_z0", "value") * cgs::kpc;
    bmodel_jaffe.arm_phi0.push_back(
        toolkit::fetchfloat(subele, "arm_phi1", "value") * cgs::rad);
    bmodel_jaffe.arm_phi0.push_back(
        toolkit::fetchfloat(subele, "arm_phi2", "value") * cgs::rad);
    bmodel_jaffe.arm_phi0.push_back(
        toolkit::fetchfloat(subele, "arm_phi3", "value") * cgs::rad);
    bmodel_jaffe.arm_phi0.push_back(
        toolkit::fetchfloat(subele, "arm_phi4", "value") * cgs::rad);
    bmodel_jaffe.arm_amp.push_back(
        toolkit::fetchfloat(subele, "arm_amp1", "value") * cgs::muGauss);
    bmodel_jaffe.arm_amp.push_back(
        toolkit::fetchfloat(subele, "arm_amp2", "value") * cgs::muGauss);
    bmodel_jaffe.arm_amp.push_back(
        toolkit::fetchfloat(subele, "arm_amp3", "value") * cgs::muGauss);
    bmodel_jaffe.arm_amp.push_back(
        toolkit::fetchfloat(subele, "arm_amp4", "value") * cgs::muGauss);
    bmodel_jaffe.arm_pitch =
        toolkit::fetchfloat(subele, "arm_pitch", "value") * cgs::rad;
    bmodel_jaffe.comp_c = toolkit::fetchfloat(subele, "comp_c", "value");
    bmodel_jaffe.comp_d =
        toolkit::fetchfloat(subele, "comp_d", "value") * cgs::kpc;
    bmodel_jaffe.comp_r =
        toolkit::fetchfloat(subele, "comp_r", "value") * cgs::kpc;
    bmodel_jaffe.comp_p = toolkit::fetchfloat(subele, "comp_p", "value");
  }
  // ES model
  subele = ele->FirstChildElement("es");
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    bmodel_list.push_back("es");
    grid_b.build_permission = true;
    bmodel_es.seed = toolkit::fetchuint(subele, "seed");
    bmodel_es.rms = toolkit::fetchfloat(subele, "rms", "value") * cgs::muGauss;
    bmodel_es.k0 = toolkit::fetchfloat(subele, "k0", "value");
    bmodel_es.a0 = toolkit::fetchfloat(subele, "a0", "value");
    bmodel_es.k1 = toolkit::fetchfloat(subele, "k1", "value");
    bmodel_es.a1 = toolkit::fetchfloat(subele, "a1", "value");
    bmodel_es.rho = toolkit::fetchfloat(subele, "rho", "value");
    bmodel_es.r0 = toolkit::fetchfloat(subele, "r0", "value") * cgs::kpc;
    bmodel_es.z0 = toolkit::fetchfloat(subele, "z0", "value") * cgs::kpc;
  }
  // MHD model
  subele = ele->FirstChildElement("mhd");
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    bmodel_list.push_back("mhd");
    grid_b.build_permission = true;
    bmodel_mhd.seed = toolkit::fetchuint(subele, "seed");
    bmodel_mhd.pa0 = toolkit::fetchfloat(subele, "pa0", "value") *
                     cgs::muGauss * cgs::muGauss;
    bmodel_mhd.pf0 = toolkit::fetchfloat(subele, "pf0", "value") *
                     cgs::muGauss * cgs::muGauss;
    bmodel_mhd.ps0 = toolkit::fetchfloat(subele, "ps0", "value") *
                     cgs::muGauss * cgs::muGauss;
    bmodel_mhd.k0 = toolkit::fetchfloat(subele, "k0", "value");
    bmodel_mhd.aa0 = toolkit::fetchfloat(subele, "aa0", "value");
    bmodel_mhd.af0 = toolkit::fetchfloat(subele, "af0", "value");
    bmodel_mhd.as0 = toolkit::fetchfloat(subele, "as0", "value");
    bmodel_mhd.k1 = toolkit::fetchfloat(subele, "k1", "value");
    bmodel_mhd.a1 = toolkit::fetchfloat(subele, "a1", "value");
    bmodel_mhd.beta = toolkit::fetchfloat(subele, "beta", "value");
    bmodel_mhd.ma = toolkit::fetchfloat(subele, "ma", "value");
  }
}

void Param::parse_temodel(tinyxml2::XMLDocument *doc) {
  // te model selection
  tinyxml2::XMLElement *ele = toolkit::tracexml(doc, {"temodel"});
  // uniform model
  tinyxml2::XMLElement *subele{ele->FirstChildElement("unif")};
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    temodel_list.push_back("unif");
    temodel_unif.n0 = toolkit::fetchfloat(subele, "n0", "value");
    temodel_unif.r0 = toolkit::fetchfloat(subele, "r0", "value") * cgs::kpc;
  }
  // YMW16 model
  subele = ele->FirstChildElement("ymw16");
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    temodel_list.push_back("ymw16");
    // Warp_Sun
    tinyxml2::XMLElement *ssubele{subele->FirstChildElement("warp")};
    temodel_ymw16.r_warp =
        toolkit::fetchfloat(ssubele, "r_warp", "value") * cgs::kpc;
    temodel_ymw16.r0 = toolkit::fetchfloat(ssubele, "r0", "value") * cgs::kpc;
    temodel_ymw16.t0_gamma_w = toolkit::fetchfloat(ssubele, "gamma_w", "value");
    // thick disk
    ssubele = subele->FirstChildElement("thickdisk");
    temodel_ymw16.t1_ad = toolkit::fetchfloat(ssubele, "ad", "value") * cgs::pc;
    temodel_ymw16.t1_bd = toolkit::fetchfloat(ssubele, "bd", "value") * cgs::pc;
    temodel_ymw16.t1_n1 = toolkit::fetchfloat(ssubele, "n1", "value");
    temodel_ymw16.t1_h1 = toolkit::fetchfloat(ssubele, "h1", "value") * cgs::pc;
    // thin disk
    ssubele = subele->FirstChildElement("thindisk");
    temodel_ymw16.t2_a2 = toolkit::fetchfloat(ssubele, "a2", "value") * cgs::pc;
    temodel_ymw16.t2_b2 = toolkit::fetchfloat(ssubele, "b2", "value") * cgs::pc;
    temodel_ymw16.t2_n2 = toolkit::fetchfloat(ssubele, "n2", "value");
    temodel_ymw16.t2_k2 = toolkit::fetchfloat(ssubele, "k2", "value");
    // spiral arm
    ssubele = subele->FirstChildElement("spiralarm");
    temodel_ymw16.t3_b2s =
        toolkit::fetchfloat(ssubele, "b2s", "value") * cgs::pc;
    temodel_ymw16.t3_narm[0] =
        toolkit::fetchfloat(ssubele, "ele_arm_0", "value");
    temodel_ymw16.t3_narm[1] =
        toolkit::fetchfloat(ssubele, "ele_arm_1", "value");
    temodel_ymw16.t3_narm[2] =
        toolkit::fetchfloat(ssubele, "ele_arm_2", "value");
    temodel_ymw16.t3_narm[3] =
        toolkit::fetchfloat(ssubele, "ele_arm_3", "value");
    temodel_ymw16.t3_narm[4] =
        toolkit::fetchfloat(ssubele, "ele_arm_4", "value");
    temodel_ymw16.t3_warm[0] =
        toolkit::fetchfloat(ssubele, "wid_arm_0", "value") * cgs::pc;
    temodel_ymw16.t3_warm[1] =
        toolkit::fetchfloat(ssubele, "wid_arm_1", "value") * cgs::pc;
    temodel_ymw16.t3_warm[2] =
        toolkit::fetchfloat(ssubele, "wid_arm_2", "value") * cgs::pc;
    temodel_ymw16.t3_warm[3] =
        toolkit::fetchfloat(ssubele, "wid_arm_3", "value") * cgs::pc;
    temodel_ymw16.t3_warm[4] =
        toolkit::fetchfloat(ssubele, "wid_arm_4", "value") * cgs::pc;
    temodel_ymw16.t3_rmin[0] =
        toolkit::fetchfloat(ssubele, "rref_arm_0", "value") * cgs::kpc;
    temodel_ymw16.t3_rmin[1] =
        toolkit::fetchfloat(ssubele, "rref_arm_1", "value") * cgs::kpc;
    temodel_ymw16.t3_rmin[2] =
        toolkit::fetchfloat(ssubele, "rref_arm_2", "value") * cgs::kpc;
    temodel_ymw16.t3_rmin[3] =
        toolkit::fetchfloat(ssubele, "rref_arm_3", "value") * cgs::kpc;
    temodel_ymw16.t3_rmin[4] =
        toolkit::fetchfloat(ssubele, "rref_arm_4", "value") * cgs::kpc;
    temodel_ymw16.t3_phimin[0] =
        toolkit::fetchfloat(ssubele, "phiref_arm_0", "value") * cgs::rad;
    temodel_ymw16.t3_phimin[1] =
        toolkit::fetchfloat(ssubele, "phiref_arm_1", "value") * cgs::rad;
    temodel_ymw16.t3_phimin[2] =
        toolkit::fetchfloat(ssubele, "phiref_arm_2", "value") * cgs::rad;
    temodel_ymw16.t3_phimin[3] =
        toolkit::fetchfloat(ssubele, "phiref_arm_3", "value") * cgs::rad;
    temodel_ymw16.t3_phimin[4] =
        toolkit::fetchfloat(ssubele, "phiref_arm_4", "value") * cgs::rad;
    temodel_ymw16.t3_tpitch[0] =
        tan(toolkit::fetchfloat(ssubele, "pitch_arm_0", "value") * cgs::rad);
    temodel_ymw16.t3_tpitch[1] =
        tan(toolkit::fetchfloat(ssubele, "pitch_arm_1", "value") * cgs::rad);
    temodel_ymw16.t3_tpitch[2] =
        tan(toolkit::fetchfloat(ssubele, "pitch_arm_2", "value") * cgs::rad);
    temodel_ymw16.t3_tpitch[3] =
        tan(toolkit::fetchfloat(ssubele, "pitch_arm_3", "value") * cgs::rad);
    temodel_ymw16.t3_tpitch[4] =
        tan(toolkit::fetchfloat(ssubele, "pitch_arm_4", "value") * cgs::rad);
    temodel_ymw16.t3_cpitch[0] =
        cos(toolkit::fetchfloat(ssubele, "pitch_arm_0", "value") * cgs::rad);
    temodel_ymw16.t3_cpitch[1] =
        cos(toolkit::fetchfloat(ssubele, "pitch_arm_1", "value") * cgs::rad);
    temodel_ymw16.t3_cpitch[2] =
        cos(toolkit::fetchfloat(ssubele, "pitch_arm_2", "value") * cgs::rad);
    temodel_ymw16.t3_cpitch[3] =
        cos(toolkit::fetchfloat(ssubele, "pitch_arm_3", "value") * cgs::rad);
    temodel_ymw16.t3_cpitch[4] =
        cos(toolkit::fetchfloat(ssubele, "pitch_arm_4", "value") * cgs::rad);
    temodel_ymw16.t3_aa = toolkit::fetchfloat(ssubele, "aa", "value") * cgs::pc;
    temodel_ymw16.t3_ka = toolkit::fetchfloat(ssubele, "ka", "value");
    temodel_ymw16.t3_ncn = toolkit::fetchfloat(ssubele, "ncn", "value");
    temodel_ymw16.t3_thetacn = toolkit::fetchfloat(ssubele, "thetacn", "value");
    temodel_ymw16.t3_wcn = toolkit::fetchfloat(ssubele, "wcn", "value");
    temodel_ymw16.t3_nsg = toolkit::fetchfloat(ssubele, "nsg", "value");
    temodel_ymw16.t3_thetasg = toolkit::fetchfloat(ssubele, "thetasg", "value");
    temodel_ymw16.t3_wsg = toolkit::fetchfloat(ssubele, "wsg", "value");
    // gc
    ssubele = subele->FirstChildElement("galcenter");
    temodel_ymw16.t4_ngc = toolkit::fetchfloat(ssubele, "ngc", "value");
    temodel_ymw16.t4_agc =
        toolkit::fetchfloat(ssubele, "agc", "value") * cgs::pc;
    temodel_ymw16.t4_hgc =
        toolkit::fetchfloat(ssubele, "hgc", "value") * cgs::pc;
    // Gum
    ssubele = subele->FirstChildElement("gumnebula");
    temodel_ymw16.t5_ngn = toolkit::fetchfloat(ssubele, "ngn", "value");
    temodel_ymw16.t5_wgn =
        toolkit::fetchfloat(ssubele, "wgn", "value") * cgs::pc;
    temodel_ymw16.t5_agn =
        toolkit::fetchfloat(ssubele, "agn", "value") * cgs::pc;
    temodel_ymw16.t5_kgn = toolkit::fetchfloat(ssubele, "kgn", "value");
    // Local Bubble
    ssubele = subele->FirstChildElement("localbubble");
    temodel_ymw16.t6_j_lb = toolkit::fetchfloat(ssubele, "j_lb", "value");
    temodel_ymw16.t6_nlb1 = toolkit::fetchfloat(ssubele, "nlb1", "value");
    temodel_ymw16.t6_thetalb1 =
        toolkit::fetchfloat(ssubele, "thetalb1", "value");
    temodel_ymw16.t6_detlb1 = toolkit::fetchfloat(ssubele, "detlb1", "value");
    temodel_ymw16.t6_wlb1 =
        toolkit::fetchfloat(ssubele, "wlb1", "value") * cgs::pc;
    temodel_ymw16.t6_hlb1 =
        toolkit::fetchfloat(ssubele, "hlb1", "value") * cgs::pc;
    temodel_ymw16.t6_nlb2 = toolkit::fetchfloat(ssubele, "nlb2", "value");
    temodel_ymw16.t6_thetalb2 =
        toolkit::fetchfloat(ssubele, "thetalb2", "value");
    temodel_ymw16.t6_detlb2 = toolkit::fetchfloat(ssubele, "detlb2", "value");
    temodel_ymw16.t6_wlb2 =
        toolkit::fetchfloat(ssubele, "wlb2", "value") * cgs::pc;
    temodel_ymw16.t6_hlb2 =
        toolkit::fetchfloat(ssubele, "hlb2", "value") * cgs::pc;
    // Loop I
    ssubele = subele->FirstChildElement("loopi");
    temodel_ymw16.t7_nli = toolkit::fetchfloat(ssubele, "nli", "value");
    temodel_ymw16.t7_rli =
        toolkit::fetchfloat(ssubele, "rli", "value") * cgs::pc;
    temodel_ymw16.t7_wli =
        toolkit::fetchfloat(ssubele, "wli", "value") * cgs::pc;
    temodel_ymw16.t7_detthetali =
        toolkit::fetchfloat(ssubele, "detthetali", "value");
    temodel_ymw16.t7_thetali = toolkit::fetchfloat(ssubele, "thetali", "value");
  }
  // dft model
  subele = ele->FirstChildElement("dft");
  if (subele != nullptr and toolkit::fetchbool(subele, "cue")) {
    temodel_list.push_back("dft");
    grid_te.build_permission = true;
    temodel_dft.seed = toolkit::fetchuint(subele, "seed");
    temodel_dft.rms = toolkit::fetchfloat(subele, "rms", "value");
    temodel_dft.k0 = toolkit::fetchfloat(subele, "k0", "value");
    temodel_dft.a0 = toolkit::fetchfloat(subele, "a0", "value");
    temodel_dft.r0 = toolkit::fetchfloat(subele, "r0", "value") * cgs::kpc;
    temodel_dft.z0 = toolkit::fetchfloat(subele, "z0", "value") * cgs::kpc;
  }
}

void Param::parse_cremodel(tinyxml2::XMLDocument *doc) {
  // cre model
  tinyxml2::XMLElement *ele{toolkit::tracexml(doc, {"cremodel"})};
  // uniform
  tinyxml2::XMLElement *subele{ele->FirstChildElement("unif")};
  if (subele != nullptr and toolkit::fetchfloat(subele, "cue")) {
    cremodel_list.push_back("unif");
    cremodel_unif.alpha = toolkit::fetchfloat(subele, "alpha", "value");
    cremodel_unif.r0 =
        toolkit::fetchfloat(subele, "r0", "value") * cgs::kpc; // kpc
    cremodel_unif.e0 =
        toolkit::fetchfloat(subele, "e0", "value") * cgs::GeV; // GeV
    cremodel_unif.j0 = toolkit::fetchfloat(subele, "j0", "value");
  }
  // analytic
  subele = ele->FirstChildElement("pwrlaw");
  if (subele != nullptr and toolkit::fetchfloat(subele, "cue")) {
    cremodel_list.push_back("pwrlaw");
    cremodel_ana.alpha = toolkit::fetchfloat(subele, "alpha", "value");
    cremodel_ana.beta = toolkit::fetchfloat(subele, "beta", "value");
    cremodel_ana.theta = toolkit::fetchfloat(subele, "theta", "value");
    cremodel_ana.r0 =
        toolkit::fetchfloat(subele, "r0", "value") * cgs::kpc; // kpc
    cremodel_ana.z0 =
        toolkit::fetchfloat(subele, "z0", "value") * cgs::kpc; // kpc
    cremodel_ana.e0 =
        toolkit::fetchfloat(subele, "e0", "value") * cgs::GeV; // GeV
    cremodel_ana.j0 = toolkit::fetchfloat(subele, "j0", "value");
  }
}

void Param::parse_grid(tinyxml2::XMLDocument *doc) {
  // observer position
  tinyxml2::XMLElement *ele{toolkit::tracexml(doc, {"grid"})};
  tinyxml2::XMLElement *subele{ele->FirstChildElement("observer")};
  observer = Hamvec<3, ham_float>{
      cgs::kpc * toolkit::fetchfloat(subele, "x", "value"),
      cgs::kpc * toolkit::fetchfloat(subele, "y", "value"),
      cgs::kpc * toolkit::fetchfloat(subele, "z", "value")};
  if (grid_obs.write_permission) {
    // shell parameters
    subele = ele->FirstChildElement("shell");
    grid_obs.oc_r_min =
        toolkit::fetchfloat(subele, "oc_r_min", "value") * cgs::kpc;
    grid_obs.oc_r_max =
        toolkit::fetchfloat(subele, "oc_r_max", "value") * cgs::kpc;
    grid_obs.gc_r_min =
        toolkit::fetchfloat(subele, "gc_r_min", "value") * cgs::kpc;
    grid_obs.gc_r_max =
        toolkit::fetchfloat(subele, "gc_r_max", "value") * cgs::kpc;
    grid_obs.gc_z_min =
        toolkit::fetchfloat(subele, "gc_z_min", "value") * cgs::kpc;
    grid_obs.gc_z_max =
        toolkit::fetchfloat(subele, "gc_z_max", "value") * cgs::kpc;
    grid_obs.oc_r_res =
        toolkit::fetchfloat(subele, "oc_r_res", "value") * cgs::kpc;
    // auto dividing
    if (toolkit::fetchstring(subele, "type") == "auto") {
      tinyxml2::XMLElement *ssubele{subele->FirstChildElement("auto")};
      grid_obs.total_shell = toolkit::fetchuint(ssubele, "shell_num", "value");
      const ham_uint nside_sim{
          toolkit::fetchuint(ssubele, "nside_sim", "value")};
      // calculate shell arrangement
      for (ham_uint i = 0; i != grid_obs.total_shell; ++i) {
        grid_obs.nside_shell.push_back(pow(2, i) * nside_sim);
      }
      grid_obs.radii_shell.push_back(grid_obs.oc_r_min);
      for (ham_uint i = 0; i < grid_obs.total_shell; ++i) {
        grid_obs.radii_shell.push_back(
            (grid_obs.oc_r_max - grid_obs.oc_r_min) *
                std::pow(0.5, grid_obs.total_shell - i - 1) +
            grid_obs.oc_r_min);
      }
    }
    // manual dividing
    else if (toolkit::fetchstring(subele, "type") == "manual") {
      tinyxml2::XMLElement *ssubele{subele->FirstChildElement("manual")};
      grid_obs.total_shell = 0;
      for (auto e = ssubele->FirstChildElement("nside_sim"); e != nullptr;
           e = e->NextSiblingElement("nside_sim")) {
        grid_obs.total_shell++;
        grid_obs.nside_shell.push_back(toolkit::fetchuint(e, "value"));
      }
      for (auto e = ssubele->FirstChildElement("cut"); e != nullptr;
           e = e->NextSiblingElement("cut")) {
        grid_obs.cut_shell.push_back(toolkit::fetchfloat(e, "value"));
      }
      grid_obs.cut_shell.push_back(1.);
      assert(grid_obs.cut_shell.size() == grid_obs.total_shell);
      assert(grid_obs.nside_shell.size() == grid_obs.total_shell);
      grid_obs.radii_shell.push_back(grid_obs.oc_r_min);
      for (auto &i : grid_obs.cut_shell) {
        grid_obs.radii_shell.push_back(
            (grid_obs.oc_r_max - grid_obs.oc_r_min) * i + grid_obs.oc_r_min);
      }
    } else {
      throw std::runtime_error("unsupported shell option");
    }
  }
  // bfield box
  if (grid_b.write_permission or grid_b.read_permission or
      grid_b.build_permission) {
    tinyxml2::XMLElement *subele{ele->FirstChildElement("bfield")};
    grid_b.nx = toolkit::fetchuint(subele, "nx", "value");
    grid_b.ny = toolkit::fetchuint(subele, "ny", "value");
    grid_b.nz = toolkit::fetchuint(subele, "nz", "value");
    grid_b.full_size = grid_b.nx * grid_b.ny * grid_b.nz;
    grid_b.x_max = cgs::kpc * toolkit::fetchfloat(subele, "x_max", "value");
    grid_b.x_min = cgs::kpc * toolkit::fetchfloat(subele, "x_min", "value");
    grid_b.y_max = cgs::kpc * toolkit::fetchfloat(subele, "y_max", "value");
    grid_b.y_min = cgs::kpc * toolkit::fetchfloat(subele, "y_min", "value");
    grid_b.z_max = cgs::kpc * toolkit::fetchfloat(subele, "z_max", "value");
    grid_b.z_min = cgs::kpc * toolkit::fetchfloat(subele, "z_min", "value");
  }
  // tefield box
  if (grid_te.write_permission or grid_te.read_permission or
      grid_te.build_permission) {
    tinyxml2::XMLElement *subele{ele->FirstChildElement("tefield")};
    grid_te.nx = toolkit::fetchuint(subele, "nx", "value");
    grid_te.ny = toolkit::fetchuint(subele, "ny", "value");
    grid_te.nz = toolkit::fetchuint(subele, "nz", "value");
    grid_te.full_size = grid_te.nx * grid_te.ny * grid_te.nz;
    grid_te.x_max = cgs::kpc * toolkit::fetchfloat(subele, "x_max", "value");
    grid_te.x_min = cgs::kpc * toolkit::fetchfloat(subele, "x_min", "value");
    grid_te.y_max = cgs::kpc * toolkit::fetchfloat(subele, "y_max", "value");
    grid_te.y_min = cgs::kpc * toolkit::fetchfloat(subele, "y_min", "value");
    grid_te.z_max = cgs::kpc * toolkit::fetchfloat(subele, "z_max", "value");
    grid_te.z_min = cgs::kpc * toolkit::fetchfloat(subele, "z_min", "value");
  }
  // crefield box
  if (grid_cre.write_permission or grid_cre.read_permission or
      grid_cre.build_permission) {
    tinyxml2::XMLElement *subele{ele->FirstChildElement("crefield")};
    grid_cre.e_min = cgs::GeV * toolkit::fetchfloat(subele, "e_min", "value");
    grid_cre.e_max = cgs::GeV * toolkit::fetchfloat(subele, "e_max", "value");
    grid_cre.ne = toolkit::fetchuint(subele, "ne", "value");
    grid_cre.e_fact =
        std::log(grid_cre.e_max / grid_cre.e_min) / (grid_cre.ne - 1);
    grid_cre.nz = toolkit::fetchuint(subele, "nz", "value");
    grid_cre.nx = toolkit::fetchuint(subele, "nx", "value");
    grid_cre.ny = toolkit::fetchuint(subele, "ny", "value");
    grid_cre.x_max = cgs::kpc * toolkit::fetchfloat(subele, "x_max", "value");
    grid_cre.x_min = cgs::kpc * toolkit::fetchfloat(subele, "x_min", "value");
    grid_cre.y_max = cgs::kpc * toolkit::fetchfloat(subele, "y_max", "value");
    grid_cre.y_min = cgs::kpc * toolkit::fetchfloat(subele, "y_min", "value");
    grid_cre.z_max = cgs::kpc * toolkit::fetchfloat(subele, "z_max", "value");
    grid_cre.z_min = cgs::kpc * toolkit::fetchfloat(subele, "z_min", "value");
    grid_cre.cre_size = grid_cre.ne * grid_cre.nx * grid_cre.ny * grid_cre.nz;
  }
}
