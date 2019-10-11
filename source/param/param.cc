#include <cmath>
#include <memory>
#include <string>

#include <cgs_units.h>
#include <hamvec.h>
#include <param.h>
#include <tinyxml2.h>
#include <toolkit.h>

Param::Param(const std::string file_name) {
  // load xml file
  std::unique_ptr<tinyxml2::XMLDocument> doc{toolkit::loadxml(file_name)};
  // observer position
  tinyxml2::XMLElement *ptr{toolkit::tracexml(doc.get(), {"grid", "observer"})};
  observer =
      hamvec<3, double>{cgs_kpc * toolkit::fetchdouble(ptr, "value", "x"),
                        cgs_kpc * toolkit::fetchdouble(ptr, "value", "y"),
                        cgs_kpc * toolkit::fetchdouble(ptr, "value", "z")};
  // collect parameters
  obs_param(doc.get());
  breg_param(doc.get());
  brnd_param(doc.get());
  tereg_param(doc.get());
  ternd_param(doc.get());
  cre_param(doc.get());
}

void Param::obs_param(tinyxml2::XMLDocument *doc) {
  // observable base path
  tinyxml2::XMLElement *ptr{toolkit::tracexml(doc, {"observable"})};
  // controller for reading shell parameters
  grid_obs.write_permission = false;
  // if dispersion measure is required
  if (ptr->FirstChildElement("dm") != nullptr) {
    grid_obs.write_permission = true;
    grid_obs.do_dm = toolkit::fetchbool(ptr, "cue", "dm");
    grid_obs.sim_dm_name = toolkit::fetchstring(ptr, "filename", "dm");
    grid_obs.nside_dm = toolkit::fetchunsigned(ptr, "nside", "dm");
  } else {
    grid_obs.do_dm = false;
  }
  // if faraday depth is required
  if (ptr->FirstChildElement("faraday") != nullptr) {
    grid_obs.write_permission = true;
    grid_obs.do_fd = toolkit::fetchbool(ptr, "cue", "faraday");
    grid_obs.sim_fd_name = toolkit::fetchstring(ptr, "filename", "faraday");
    grid_obs.nside_fd = toolkit::fetchunsigned(ptr, "nside", "faraday");
  } else {
    grid_obs.do_fd = false;
  }
  // if synchrotron emission is required
  if (ptr->FirstChildElement("sync") != nullptr) {
    grid_obs.write_permission = true;
    tinyxml2::XMLElement *subptr{
        toolkit::tracexml(doc, {"observable", "sync"})};
    grid_obs.do_sync.push_back(toolkit::fetchbool(subptr, "cue"));
    grid_obs.sim_sync_freq.push_back(toolkit::fetchdouble(subptr, "freq") *
                                     cgs_GHz);
    grid_obs.sim_sync_name.push_back(toolkit::fetchstring(subptr, "filename"));
    grid_obs.nside_sync.push_back(toolkit::fetchunsigned(subptr, "nside"));
    for (auto e = subptr->NextSiblingElement("sync"); e != nullptr;
         e = e->NextSiblingElement("sync")) {
      grid_obs.do_sync.push_back(toolkit::fetchbool(e, "cue"));
      grid_obs.sim_sync_freq.push_back(toolkit::fetchdouble(e, "freq") *
                                       cgs_GHz);
      grid_obs.sim_sync_name.push_back(toolkit::fetchstring(e, "filename"));
      grid_obs.nside_sync.push_back(toolkit::fetchunsigned(e, "nside"));
    }
  } else {
    grid_obs.do_sync.push_back(false);
  }
  // if any observable is requried
  if (grid_obs.write_permission) {
    ptr = toolkit::tracexml(doc, {"grid", "shell"});
    grid_obs.oc_r_min =
        toolkit::fetchdouble(ptr, "value", "oc_r_min") * cgs_kpc;
    grid_obs.oc_r_max =
        toolkit::fetchdouble(ptr, "value", "oc_r_max") * cgs_kpc;
    grid_obs.gc_r_min =
        toolkit::fetchdouble(ptr, "value", "gc_r_min") * cgs_kpc;
    grid_obs.gc_r_max =
        toolkit::fetchdouble(ptr, "value", "gc_r_max") * cgs_kpc;
    grid_obs.gc_z_min =
        toolkit::fetchdouble(ptr, "value", "gc_z_min") * cgs_kpc;
    grid_obs.gc_z_max =
        toolkit::fetchdouble(ptr, "value", "gc_z_max") * cgs_kpc;
    grid_obs.oc_r_res =
        toolkit::fetchdouble(ptr, "value", "oc_r_res") * cgs_kpc;
    // auto shell cutting
    if (toolkit::fetchstring(ptr, "type", "layer") == "auto") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"grid", "shell", "layer", "auto"})};
      grid_obs.total_shell =
          toolkit::fetchunsigned(subptr, "value", "shell_num");
      const unsigned int nside_sim{
          toolkit::fetchunsigned(subptr, "value", "nside_sim")};
      for (std::size_t i = 0; i != grid_obs.total_shell; ++i) {
        grid_obs.nside_shell.push_back(pow(2, i) * nside_sim);
      }
      grid_obs.radii_shell.push_back(grid_obs.oc_r_min);
      for (std::size_t i = 0; i < grid_obs.total_shell; ++i) {
        grid_obs.radii_shell.push_back(
            (grid_obs.oc_r_max - grid_obs.oc_r_min) *
                std::pow(0.5, grid_obs.total_shell - i - 1) +
            grid_obs.oc_r_min);
      }
    }
    // manual shell cutting
    else if (toolkit::fetchstring(ptr, "type", "layer") == "manual") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"grid", "shell", "layer", "manual"})};
      grid_obs.total_shell = 0;
      for (auto e = subptr->FirstChildElement("nside_sim"); e != nullptr;
           e = e->NextSiblingElement("nside_sim")) {
        grid_obs.total_shell++;
        grid_obs.nside_shell.push_back(toolkit::fetchunsigned(e, "value"));
      }
      for (auto e = subptr->FirstChildElement("cut"); e != nullptr;
           e = e->NextSiblingElement("cut")) {
        grid_obs.cut_shell.push_back(toolkit::fetchdouble(e, "value"));
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
      throw std::runtime_error("unsupported layer option");
    }
    // get mask map upon request
    ptr = toolkit::tracexml(doc, {"mask"});
    if (toolkit::fetchbool(ptr, "cue")) {
      grid_obs.do_mask = true;
      grid_obs.mask_name = toolkit::fetchstring(ptr, "filename");
    } else {
      grid_obs.do_mask = false;
    }
  }
}

void Param::breg_param(tinyxml2::XMLDocument *doc) {
  // breg io
  tinyxml2::XMLElement *ptr{toolkit::tracexml(doc, {"fieldio"})};
  if (ptr->FirstChildElement("breg") != nullptr) {
    grid_breg.read_permission = toolkit::fetchbool(ptr, "read", "breg");
    grid_breg.write_permission = toolkit::fetchbool(ptr, "write", "breg");
    grid_breg.filename = toolkit::fetchstring(ptr, "filename", "breg");
  }
  // breg internal
  ptr = toolkit::tracexml(doc, {"magneticfield"});
  grid_breg.build_permission = toolkit::fetchbool(ptr, "cue", "regular");
  // if no external read and internal model is active
  if (grid_breg.build_permission and not grid_breg.read_permission) {
    breg_type = toolkit::fetchstring(ptr, "type", "regular");
    // blsa
    if (breg_type == "lsa") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"magneticfield", "regular", "lsa"})};
      breg_lsa.b0 = toolkit::fetchdouble(subptr, "value", "b0") *
                    cgs_muGauss; // microGauss
      breg_lsa.psi0 =
          toolkit::fetchdouble(subptr, "value", "psi0") * cgs_rad; // rad
      breg_lsa.psi1 =
          toolkit::fetchdouble(subptr, "value", "psi1") * cgs_rad; // rad
      breg_lsa.chi0 =
          toolkit::fetchdouble(subptr, "value", "chi0") * cgs_rad; // rad
    }
    // bjaffe
    else if (breg_type == "jaffe") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"magneticfield", "regular", "jaffe"})};
      breg_jaffe.quadruple = toolkit::fetchbool(subptr, "cue", "quadruple");
      breg_jaffe.bss = toolkit::fetchbool(subptr, "cue", "bss");
      breg_jaffe.disk_amp = toolkit::fetchdouble(subptr, "value", "disk_amp") *
                            cgs_muGauss; // microG
      breg_jaffe.disk_z0 =
          toolkit::fetchdouble(subptr, "value", "disk_z0") * cgs_kpc; // kpc
      breg_jaffe.halo_amp = toolkit::fetchdouble(subptr, "value", "halo_amp") *
                            cgs_muGauss; // microG
      breg_jaffe.halo_z0 =
          toolkit::fetchdouble(subptr, "value", "halo_z0") * cgs_kpc; // kpc
      breg_jaffe.r_inner =
          toolkit::fetchdouble(subptr, "value", "r_inner") * cgs_kpc; // kpc
      breg_jaffe.r_scale =
          toolkit::fetchdouble(subptr, "value", "r_scale") * cgs_kpc; // kpc
      breg_jaffe.r_peak =
          toolkit::fetchdouble(subptr, "value", "r_peak") * cgs_kpc; // kpc
      breg_jaffe.ring = toolkit::fetchbool(subptr, "cue", "ring");
      breg_jaffe.ring_amp = toolkit::fetchdouble(subptr, "value", "ring_amp") *
                            cgs_muGauss; // microG
      breg_jaffe.ring_r =
          toolkit::fetchdouble(subptr, "value", "ring_r") * cgs_kpc; // kpc
      breg_jaffe.bar = toolkit::fetchbool(subptr, "cue", "bar");
      breg_jaffe.bar_amp = toolkit::fetchdouble(subptr, "value", "bar_amp") *
                           cgs_muGauss; // microG
      breg_jaffe.bar_a =
          toolkit::fetchdouble(subptr, "value", "bar_a") * cgs_kpc; // kpc
      breg_jaffe.bar_b =
          toolkit::fetchdouble(subptr, "value", "bar_b") * cgs_kpc; // kpc
      breg_jaffe.bar_phi0 =
          toolkit::fetchdouble(subptr, "value", "bar_phi0") * cgs_rad; // rad
      breg_jaffe.arm_num = toolkit::fetchunsigned(subptr, "value", "arm_num");
      breg_jaffe.arm_r0 =
          toolkit::fetchdouble(subptr, "value", "arm_r0") * cgs_kpc; // kpc
      breg_jaffe.arm_z0 =
          toolkit::fetchdouble(subptr, "value", "arm_z0") * cgs_kpc; // kpc
      breg_jaffe.arm_phi0.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_phi1") * cgs_rad); // rad
      breg_jaffe.arm_phi0.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_phi2") * cgs_rad); // rad
      breg_jaffe.arm_phi0.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_phi3") * cgs_rad); // rad
      breg_jaffe.arm_phi0.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_phi4") * cgs_rad); // rad
      breg_jaffe.arm_amp.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_amp1") *
          cgs_muGauss); // microG
      breg_jaffe.arm_amp.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_amp2") *
          cgs_muGauss); // microG
      breg_jaffe.arm_amp.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_amp3") *
          cgs_muGauss); // microG
      breg_jaffe.arm_amp.push_back(
          toolkit::fetchdouble(subptr, "value", "arm_amp4") *
          cgs_muGauss); // microG
      breg_jaffe.arm_pitch =
          toolkit::fetchdouble(subptr, "value", "arm_pitch") * cgs_rad; // rad
      breg_jaffe.comp_c = toolkit::fetchdouble(subptr, "value", "comp_c");
      breg_jaffe.comp_d =
          toolkit::fetchdouble(subptr, "value", "comp_d") * cgs_kpc; // kpc
      breg_jaffe.comp_r =
          toolkit::fetchdouble(subptr, "value", "comp_r") * cgs_kpc; // kpc
      breg_jaffe.comp_p = toolkit::fetchdouble(subptr, "value", "comp_p");
    } else if (breg_type == "unif") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"magneticfield", "regular", "unif"})};
      breg_unif.bp = toolkit::fetchdouble(subptr, "value", "bp") *
                     cgs_muGauss; // microGauss
      breg_unif.bv = toolkit::fetchdouble(subptr, "value", "bv") *
                     cgs_muGauss; // microGauss
      breg_unif.l0 =
          toolkit::fetchdouble(subptr, "value", "l0") * cgs_rad; // rad
    } else {
      throw std::runtime_error("unsupported breg model");
    }
  }
  // breg io box
  if (grid_breg.read_permission or grid_breg.write_permission) {
    // breg box
    tinyxml2::XMLElement *subptr{toolkit::tracexml(doc, {"grid", "box_breg"})};
    grid_breg.nx = toolkit::fetchunsigned(subptr, "value", "nx");
    grid_breg.ny = toolkit::fetchunsigned(subptr, "value", "ny");
    grid_breg.nz = toolkit::fetchunsigned(subptr, "value", "nz");
    grid_breg.full_size = grid_breg.nx * grid_breg.ny * grid_breg.nz;
    grid_breg.x_max = cgs_kpc * toolkit::fetchdouble(subptr, "value", "x_max");
    grid_breg.x_min = cgs_kpc * toolkit::fetchdouble(subptr, "value", "x_min");
    grid_breg.y_max = cgs_kpc * toolkit::fetchdouble(subptr, "value", "y_max");
    grid_breg.y_min = cgs_kpc * toolkit::fetchdouble(subptr, "value", "y_min");
    grid_breg.z_max = cgs_kpc * toolkit::fetchdouble(subptr, "value", "z_max");
    grid_breg.z_min = cgs_kpc * toolkit::fetchdouble(subptr, "value", "z_min");
  }
}

void Param::brnd_param(tinyxml2::XMLDocument *doc) {
  // brnd io
  tinyxml2::XMLElement *ptr{toolkit::tracexml(doc, {"fieldio"})};
  if (ptr->FirstChildElement("brnd") != nullptr) {
    grid_brnd.read_permission = toolkit::fetchbool(ptr, "read", "brnd");
    grid_brnd.write_permission = toolkit::fetchbool(ptr, "write", "brnd");
    grid_brnd.filename = toolkit::fetchstring(ptr, "filename", "brnd");
  }
  // brnd internal
  ptr = toolkit::tracexml(doc, {"magneticfield"});
  grid_brnd.build_permission = toolkit::fetchbool(ptr, "cue", "random");
  // if no external read and internal model is active
  if (grid_brnd.build_permission and not grid_brnd.read_permission) {
    // random seed
    brnd_seed = toolkit::fetchunsigned(ptr, "seed", "random");
    brnd_type = toolkit::fetchstring(ptr, "type", "random");
    // brnd_global
    if (brnd_type == "global") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"magneticfield", "random", "global"})};
      brnd_method = toolkit::fetchstring(subptr, "type");
      if (brnd_method == "es") {
        subptr =
            toolkit::tracexml(doc, {"magneticfield", "random", "global", "es"});
        brnd_es.rms =
            toolkit::fetchdouble(subptr, "value", "rms") * cgs_muGauss;
        brnd_es.k0 = toolkit::fetchdouble(subptr, "value", "k0");
        brnd_es.a0 = toolkit::fetchdouble(subptr, "value", "a0");
        brnd_es.k1 = toolkit::fetchdouble(subptr, "value", "k1");
        brnd_es.a1 = toolkit::fetchdouble(subptr, "value", "a1");
        brnd_es.rho = toolkit::fetchdouble(subptr, "value", "rho");
        brnd_es.r0 = toolkit::fetchdouble(subptr, "value", "r0") * cgs_kpc;
        brnd_es.z0 = toolkit::fetchdouble(subptr, "value", "z0") * cgs_kpc;
      } else if (brnd_method == "jaffe") {
        subptr = toolkit::tracexml(
            doc, {"magneticfield", "random", "global", "jaffe"});
        // to be implemented
      } else {
        throw std::runtime_error("unsupported brnd model");
      }
    }
    // brnd_local
    else if (brnd_type == "local") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"magneticfield", "random", "local"})};
      brnd_method = toolkit::fetchstring(subptr, "type");
      if (brnd_method == "mhd") {
        subptr =
            toolkit::tracexml(doc, {"magneticfield", "random", "local", "mhd"});
        brnd_mhd.pa0 = toolkit::fetchdouble(subptr, "value", "pa0") *
                       cgs_muGauss * cgs_muGauss;
        brnd_mhd.pf0 = toolkit::fetchdouble(subptr, "value", "pf0") *
                       cgs_muGauss * cgs_muGauss;
        brnd_mhd.ps0 = toolkit::fetchdouble(subptr, "value", "ps0") *
                       cgs_muGauss * cgs_muGauss;
        brnd_mhd.k0 = toolkit::fetchdouble(subptr, "value", "k0");
        brnd_mhd.aa0 = toolkit::fetchdouble(subptr, "value", "aa0");
        brnd_mhd.af0 = toolkit::fetchdouble(subptr, "value", "af0");
        brnd_mhd.as0 = toolkit::fetchdouble(subptr, "value", "as0");
        brnd_mhd.k1 = toolkit::fetchdouble(subptr, "value", "k1");
        brnd_mhd.a1 = toolkit::fetchdouble(subptr, "value", "a1");
        brnd_mhd.beta = toolkit::fetchdouble(subptr, "value", "beta");
        brnd_mhd.ma = toolkit::fetchdouble(subptr, "value", "ma");
      } else {
        throw std::runtime_error("unsupported brnd model");
      }
    } else {
      throw std::runtime_error("unsupported brnd type");
    }
  }
  // brnd io box
  if (grid_brnd.read_permission or grid_brnd.write_permission or
      grid_brnd.build_permission) {
    // brnd box
    ptr = toolkit::tracexml(doc, {"grid", "box_brnd"});
    grid_brnd.nx = toolkit::fetchunsigned(ptr, "value", "nx");
    grid_brnd.ny = toolkit::fetchunsigned(ptr, "value", "ny");
    grid_brnd.nz = toolkit::fetchunsigned(ptr, "value", "nz");
    grid_brnd.full_size = grid_brnd.nx * grid_brnd.ny * grid_brnd.nz;
    grid_brnd.x_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_max");
    grid_brnd.x_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_min");
    grid_brnd.y_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_max");
    grid_brnd.y_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_min");
    grid_brnd.z_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_max");
    grid_brnd.z_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_min");
  }
}

void Param::tereg_param(tinyxml2::XMLDocument *doc) {
  // tereg io
  tinyxml2::XMLElement *ptr{toolkit::tracexml(doc, {"fieldio"})};
  if (ptr->FirstChildElement("tereg") != nullptr) {
    grid_tereg.read_permission = toolkit::fetchbool(ptr, "read", "tereg");
    grid_tereg.write_permission = toolkit::fetchbool(ptr, "write", "tereg");
    grid_tereg.filename = toolkit::fetchstring(ptr, "filename", "tereg");
  }
  // tereg internal
  ptr = toolkit::tracexml(doc, {"thermalelectron"});
  grid_tereg.build_permission = toolkit::fetchbool(ptr, "cue", "regular");
  if (grid_tereg.build_permission) {
    tereg_type = toolkit::fetchstring(ptr, "type", "regular");
    // YMW16
    if (tereg_type == "ymw16") {
      // Warp_Sun
      tinyxml2::XMLElement *subptr{toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "warp"})};
      tereg_ymw16.r_warp =
          toolkit::fetchdouble(subptr, "value", "r_warp") * cgs_kpc; // kpc
      tereg_ymw16.r0 =
          toolkit::fetchdouble(subptr, "value", "r0") * cgs_kpc; // kpc
      tereg_ymw16.t0_gamma_w = toolkit::fetchdouble(subptr, "value", "gamma_w");
      // thick disk
      subptr = toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "thickdisk"});
      tereg_ymw16.t1_ad =
          toolkit::fetchdouble(subptr, "value", "ad") * cgs_pc; // pc
      tereg_ymw16.t1_bd =
          toolkit::fetchdouble(subptr, "value", "bd") * cgs_pc;        // pc
      tereg_ymw16.t1_n1 = toolkit::fetchdouble(subptr, "value", "n1"); // pccm
      tereg_ymw16.t1_h1 =
          toolkit::fetchdouble(subptr, "value", "h1") * cgs_pc; // pc
      // thin disk
      subptr = toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "thindisk"});
      tereg_ymw16.t2_a2 =
          toolkit::fetchdouble(subptr, "value", "a2") * cgs_pc; // pc
      tereg_ymw16.t2_b2 =
          toolkit::fetchdouble(subptr, "value", "b2") * cgs_pc;        // pc
      tereg_ymw16.t2_n2 = toolkit::fetchdouble(subptr, "value", "n2"); // pccm
      tereg_ymw16.t2_k2 = toolkit::fetchdouble(subptr, "value", "k2");
      // spiral arm
      subptr = toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "spiralarm"});
      tereg_ymw16.t3_b2s =
          toolkit::fetchdouble(subptr, "value", "b2s") * cgs_pc; // pc
      tereg_ymw16.t3_narm[0] =
          toolkit::fetchdouble(subptr, "value", "ele_arm_0"); // pccm
      tereg_ymw16.t3_narm[1] =
          toolkit::fetchdouble(subptr, "value", "ele_arm_1");
      tereg_ymw16.t3_narm[2] =
          toolkit::fetchdouble(subptr, "value", "ele_arm_2");
      tereg_ymw16.t3_narm[3] =
          toolkit::fetchdouble(subptr, "value", "ele_arm_3");
      tereg_ymw16.t3_narm[4] =
          toolkit::fetchdouble(subptr, "value", "ele_arm_4");
      tereg_ymw16.t3_warm[0] =
          toolkit::fetchdouble(subptr, "value", "wid_arm_0") * cgs_pc; // pc
      tereg_ymw16.t3_warm[1] =
          toolkit::fetchdouble(subptr, "value", "wid_arm_1") * cgs_pc;
      tereg_ymw16.t3_warm[2] =
          toolkit::fetchdouble(subptr, "value", "wid_arm_2") * cgs_pc;
      tereg_ymw16.t3_warm[3] =
          toolkit::fetchdouble(subptr, "value", "wid_arm_3") * cgs_pc;
      tereg_ymw16.t3_warm[4] =
          toolkit::fetchdouble(subptr, "value", "wid_arm_4") * cgs_pc;
      tereg_ymw16.t3_rmin[0] =
          toolkit::fetchdouble(subptr, "value", "rref_arm_0") * cgs_kpc; // kpc
      tereg_ymw16.t3_rmin[1] =
          toolkit::fetchdouble(subptr, "value", "rref_arm_1") * cgs_kpc;
      tereg_ymw16.t3_rmin[2] =
          toolkit::fetchdouble(subptr, "value", "rref_arm_2") * cgs_kpc;
      tereg_ymw16.t3_rmin[3] =
          toolkit::fetchdouble(subptr, "value", "rref_arm_3") * cgs_kpc;
      tereg_ymw16.t3_rmin[4] =
          toolkit::fetchdouble(subptr, "value", "rref_arm_4") * cgs_kpc;
      tereg_ymw16.t3_phimin[0] =
          toolkit::fetchdouble(subptr, "value", "phiref_arm_0") *
          cgs_rad; // rad
      tereg_ymw16.t3_phimin[1] =
          toolkit::fetchdouble(subptr, "value", "phiref_arm_1") *
          cgs_rad; // rad
      tereg_ymw16.t3_phimin[2] =
          toolkit::fetchdouble(subptr, "value", "phiref_arm_2") *
          cgs_rad; // rad
      tereg_ymw16.t3_phimin[3] =
          toolkit::fetchdouble(subptr, "value", "phiref_arm_3") *
          cgs_rad; // rad
      tereg_ymw16.t3_phimin[4] =
          toolkit::fetchdouble(subptr, "value", "phiref_arm_4") *
          cgs_rad; // rad
      tereg_ymw16.t3_tpitch[0] =
          tan(toolkit::fetchdouble(subptr, "value", "pitch_arm_0") * cgs_rad);
      tereg_ymw16.t3_tpitch[1] =
          tan(toolkit::fetchdouble(subptr, "value", "pitch_arm_1") * cgs_rad);
      tereg_ymw16.t3_tpitch[2] =
          tan(toolkit::fetchdouble(subptr, "value", "pitch_arm_2") * cgs_rad);
      tereg_ymw16.t3_tpitch[3] =
          tan(toolkit::fetchdouble(subptr, "value", "pitch_arm_3") * cgs_rad);
      tereg_ymw16.t3_tpitch[4] =
          tan(toolkit::fetchdouble(subptr, "value", "pitch_arm_4") * cgs_rad);
      tereg_ymw16.t3_cpitch[0] =
          cos(toolkit::fetchdouble(subptr, "value", "pitch_arm_0") * cgs_rad);
      tereg_ymw16.t3_cpitch[1] =
          cos(toolkit::fetchdouble(subptr, "value", "pitch_arm_1") * cgs_rad);
      tereg_ymw16.t3_cpitch[2] =
          cos(toolkit::fetchdouble(subptr, "value", "pitch_arm_2") * cgs_rad);
      tereg_ymw16.t3_cpitch[3] =
          cos(toolkit::fetchdouble(subptr, "value", "pitch_arm_3") * cgs_rad);
      tereg_ymw16.t3_cpitch[4] =
          cos(toolkit::fetchdouble(subptr, "value", "pitch_arm_4") * cgs_rad);
      tereg_ymw16.t3_aa =
          toolkit::fetchdouble(subptr, "value", "aa") * cgs_pc; // pc
      tereg_ymw16.t3_ka = toolkit::fetchdouble(subptr, "value", "ka");
      tereg_ymw16.t3_ncn = toolkit::fetchdouble(subptr, "value", "ncn");
      tereg_ymw16.t3_thetacn =
          toolkit::fetchdouble(subptr, "value", "thetacn");              // deg
      tereg_ymw16.t3_wcn = toolkit::fetchdouble(subptr, "value", "wcn"); // deg
      tereg_ymw16.t3_nsg = toolkit::fetchdouble(subptr, "value", "nsg");
      tereg_ymw16.t3_thetasg =
          toolkit::fetchdouble(subptr, "value", "thetasg");              // deg
      tereg_ymw16.t3_wsg = toolkit::fetchdouble(subptr, "value", "wsg"); // deg
      // gc
      subptr = toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "galcenter"});
      tereg_ymw16.t4_ngc = toolkit::fetchdouble(subptr, "value", "ngc"); // pccm
      tereg_ymw16.t4_agc =
          toolkit::fetchdouble(subptr, "value", "agc") * cgs_pc; // pc
      tereg_ymw16.t4_hgc =
          toolkit::fetchdouble(subptr, "value", "hgc") * cgs_pc; // pc
      // Gum
      subptr = toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "gumnebula"});
      tereg_ymw16.t5_ngn = toolkit::fetchdouble(subptr, "value", "ngn"); // pccm
      tereg_ymw16.t5_wgn =
          toolkit::fetchdouble(subptr, "value", "wgn") * cgs_pc; // pc
      tereg_ymw16.t5_agn =
          toolkit::fetchdouble(subptr, "value", "agn") * cgs_pc; // pc
      tereg_ymw16.t5_kgn = toolkit::fetchdouble(subptr, "value", "kgn");
      // Local Bubble
      subptr = toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "localbubble"});
      tereg_ymw16.t6_j_lb = toolkit::fetchdouble(subptr, "value", "j_lb");
      tereg_ymw16.t6_nlb1 =
          toolkit::fetchdouble(subptr, "value", "nlb1"); // pccm
      tereg_ymw16.t6_thetalb1 =
          toolkit::fetchdouble(subptr, "value", "thetalb1"); // deg
      tereg_ymw16.t6_detlb1 =
          toolkit::fetchdouble(subptr, "value", "detlb1"); // deg
      tereg_ymw16.t6_wlb1 =
          toolkit::fetchdouble(subptr, "value", "wlb1") * cgs_pc; // pc
      tereg_ymw16.t6_hlb1 =
          toolkit::fetchdouble(subptr, "value", "hlb1") * cgs_pc; // pc
      tereg_ymw16.t6_nlb2 =
          toolkit::fetchdouble(subptr, "value", "nlb2"); // pccm
      tereg_ymw16.t6_thetalb2 =
          toolkit::fetchdouble(subptr, "value", "thetalb2"); // deg
      tereg_ymw16.t6_detlb2 =
          toolkit::fetchdouble(subptr, "value", "detlb2"); // deg
      tereg_ymw16.t6_wlb2 =
          toolkit::fetchdouble(subptr, "value", "wlb2") * cgs_pc; // pc
      tereg_ymw16.t6_hlb2 =
          toolkit::fetchdouble(subptr, "value", "hlb2") * cgs_pc; // pc
      // Loop I
      subptr = toolkit::tracexml(
          doc, {"thermalelectron", "regular", "ymw16", "loopi"});
      tereg_ymw16.t7_nli = toolkit::fetchdouble(subptr, "value", "nli"); // pccm
      tereg_ymw16.t7_rli =
          toolkit::fetchdouble(subptr, "value", "rli") * cgs_pc; // pc
      tereg_ymw16.t7_wli =
          toolkit::fetchdouble(subptr, "value", "wli") * cgs_pc; // pc
      tereg_ymw16.t7_detthetali =
          toolkit::fetchdouble(subptr, "value", "detthetali"); // deg
      tereg_ymw16.t7_thetali =
          toolkit::fetchdouble(subptr, "value", "thetali"); // deg
    }
    // uniform
    else if (tereg_type == "unif") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"thermalelectron", "regular", "unif"})};
      tereg_unif.n0 = toolkit::fetchdouble(subptr, "value", "n0");
      tereg_unif.r0 =
          toolkit::fetchdouble(subptr, "value", "r0") * cgs_kpc; // kpc
    } else {
      throw std::runtime_error("unsupported tereg model");
    }
  }
  // breg io box
  if (grid_tereg.read_permission or grid_tereg.write_permission) {
    ptr = toolkit::tracexml(doc, {"grid", "box_tereg"});
    grid_tereg.nx = toolkit::fetchunsigned(ptr, "value", "nx");
    grid_tereg.ny = toolkit::fetchunsigned(ptr, "value", "ny");
    grid_tereg.nz = toolkit::fetchunsigned(ptr, "value", "nz");
    grid_tereg.full_size = grid_tereg.nx * grid_tereg.ny * grid_tereg.nz;
    grid_tereg.x_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_max");
    grid_tereg.x_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_min");
    grid_tereg.y_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_max");
    grid_tereg.y_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_min");
    grid_tereg.z_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_max");
    grid_tereg.z_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_min");
  }
}

void Param::ternd_param(tinyxml2::XMLDocument *doc) {
  // ternd io
  tinyxml2::XMLElement *ptr{toolkit::tracexml(doc, {"fieldio"})};
  if (ptr->FirstChildElement("ternd") != nullptr) {
    grid_ternd.read_permission = toolkit::fetchbool(ptr, "read", "ternd");
    grid_ternd.write_permission = toolkit::fetchbool(ptr, "write", "ternd");
    grid_ternd.filename = toolkit::fetchstring(ptr, "filename", "ternd");
  }
  // ternd internal
  ptr = toolkit::tracexml(doc, {"thermalelectron"});
  grid_ternd.build_permission = toolkit::fetchbool(ptr, "cue", "random");
  if (grid_ternd.build_permission) {
    // random seed
    ternd_seed = toolkit::fetchunsigned(ptr, "seed", "random");
    ternd_type = toolkit::fetchstring(ptr, "type", "random");
    // global turbulent
    if (ternd_type == "global") {
      tinyxml2::XMLElement *subptr{
          toolkit::tracexml(doc, {"thermalelectron", "random", "global"})};
      ternd_method = toolkit::fetchstring(subptr, "type");
      if (ternd_method == "dft") {
        subptr = toolkit::tracexml(
            doc, {"thermalelectron", "random", "global", "dft"});
        ternd_dft.rms = toolkit::fetchdouble(subptr, "value", "rms");
        ternd_dft.k0 = toolkit::fetchdouble(subptr, "value", "k0");
        ternd_dft.a0 = toolkit::fetchdouble(subptr, "value", "a0");
        ternd_dft.r0 = toolkit::fetchdouble(subptr, "value", "r0") * cgs_kpc;
        ternd_dft.z0 = toolkit::fetchdouble(subptr, "value", "z0") * cgs_kpc;
      } else {
        throw std::runtime_error("unsupported ternd model");
      }
    } else {
      throw std::runtime_error("unsupported ternd type");
    }
  }
  // ternd io grid
  if (grid_ternd.read_permission or grid_ternd.write_permission or
      grid_ternd.build_permission) {
    ptr = toolkit::tracexml(doc, {"grid", "box_ternd"});
    grid_ternd.nx = toolkit::fetchunsigned(ptr, "value", "nx");
    grid_ternd.ny = toolkit::fetchunsigned(ptr, "value", "ny");
    grid_ternd.nz = toolkit::fetchunsigned(ptr, "value", "nz");
    grid_ternd.full_size = grid_ternd.nx * grid_ternd.ny * grid_ternd.nz;
    grid_ternd.x_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_max");
    grid_ternd.x_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_min");
    grid_ternd.y_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_max");
    grid_ternd.y_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_min");
    grid_ternd.z_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_max");
    grid_ternd.z_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_min");
  }
}

void Param::cre_param(tinyxml2::XMLDocument *doc) {
  // cre io
  tinyxml2::XMLElement *ptr{toolkit::tracexml(doc, {"fieldio"})};
  if (ptr->FirstChildElement("cre") != nullptr) {
    grid_cre.read_permission = toolkit::fetchbool(ptr, "read", "cre");
    grid_cre.write_permission = toolkit::fetchbool(ptr, "write", "cre");
    grid_cre.filename = toolkit::fetchstring(ptr, "filename", "cre");
  }
  // cre internal
  ptr = toolkit::tracexml(doc, {"cre"});
  grid_cre.build_permission = toolkit::fetchbool(ptr, "cue");
  if (grid_cre.build_permission) {
    cre_type = ptr->Attribute("type");
    // analytical
    if (cre_type == "analytic") {
      tinyxml2::XMLElement *subptr{toolkit::tracexml(doc, {"cre", "analytic"})};
      cre_ana.alpha = toolkit::fetchdouble(subptr, "value", "alpha");
      cre_ana.beta = toolkit::fetchdouble(subptr, "value", "beta");
      cre_ana.theta = toolkit::fetchdouble(subptr, "value", "theta");
      cre_ana.r0 = toolkit::fetchdouble(subptr, "value", "r0") * cgs_kpc; // kpc
      cre_ana.z0 = toolkit::fetchdouble(subptr, "value", "z0") * cgs_kpc; // kpc
      cre_ana.E0 = toolkit::fetchdouble(subptr, "value", "E0") * cgs_GeV; // GeV
      cre_ana.j0 = toolkit::fetchdouble(subptr, "value", "j0");
    }
    // uniform
    else if (cre_type == "unif") {
      tinyxml2::XMLElement *subptr{toolkit::tracexml(doc, {"cre", "unif"})};
      cre_unif.alpha = toolkit::fetchdouble(subptr, "value", "alpha");
      cre_unif.r0 =
          toolkit::fetchdouble(subptr, "value", "r0") * cgs_kpc; // kpc
      cre_unif.E0 =
          toolkit::fetchdouble(subptr, "value", "E0") * cgs_GeV; // GeV
      cre_unif.j0 = toolkit::fetchdouble(subptr, "value", "j0");
    } else {
      throw std::runtime_error("unsupported cre model");
    }
  }
  // cre io grid
  if (grid_cre.read_permission or grid_cre.write_permission) {
    ptr = toolkit::tracexml(doc, {"grid", "box_cre"});
    grid_cre.E_min = cgs_GeV * toolkit::fetchdouble(ptr, "value", "E_min");
    grid_cre.E_max = cgs_GeV * toolkit::fetchdouble(ptr, "value", "E_max");
    grid_cre.nE = toolkit::fetchunsigned(ptr, "value", "nE");
    grid_cre.E_fact =
        std::log(grid_cre.E_max / grid_cre.E_min) / (grid_cre.nE - 1);
    grid_cre.nz = toolkit::fetchunsigned(ptr, "value", "nz");
    grid_cre.nx = toolkit::fetchunsigned(ptr, "value", "nx");
    grid_cre.ny = toolkit::fetchunsigned(ptr, "value", "ny");
    grid_cre.x_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_max");
    grid_cre.x_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "x_min");
    grid_cre.y_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_max");
    grid_cre.y_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "y_min");
    grid_cre.z_max = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_max");
    grid_cre.z_min = cgs_kpc * toolkit::fetchdouble(ptr, "value", "z_min");
    grid_cre.cre_size = grid_cre.nE * grid_cre.nx * grid_cre.ny * grid_cre.nz;
  }
}
