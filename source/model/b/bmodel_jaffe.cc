#include <bfield.h>
#include <cmath>
#include <grid.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>
#include <vector>

Hamvec<3, ham_float> Bmodel_jaffe::read_model(const Hamvec<3, ham_float> &pos,
                                              const Param *par) const {
  ham_float inner_b{0};
  if (par->bmodel_jaffe.ring)
    inner_b = par->bmodel_jaffe.ring_amp;
  else if (par->bmodel_jaffe.bar)
    inner_b = par->bmodel_jaffe.bar_amp;

  Hamvec<3, ham_float> bhat{orientation(pos, par)};
  Hamvec<3, ham_float> btot;
  btot = bhat * radial_scaling(pos, par) *
         (par->bmodel_jaffe.disk_amp * disk_scaling(pos, par) +
          par->bmodel_jaffe.halo_amp * halo_scaling(pos, par));
  // compress factor for each arm or for ring/bar
  std::vector<ham_float> arm{arm_compress(pos, par)};
  // only inner region
  if (arm.size() == 1) {
    btot += bhat * arm[0] * inner_b;
  }
  // spiral arm region
  else {
    for (decltype(arm.size()) i = 0; i < arm.size(); ++i) {
      btot += bhat * arm[i] * par->bmodel_jaffe.arm_amp[i];
    }
  }
  return btot;
}

Hamvec<3, ham_float> Bmodel_jaffe::orientation(const Hamvec<3, ham_float> &pos,
                                               const Param *par) const {
  // cylindrical frame
  const ham_float r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  const ham_float r_lim{par->bmodel_jaffe.ring_r};
  const ham_float bar_lim{par->bmodel_jaffe.bar_a +
                          0.5 * par->bmodel_jaffe.comp_d};
  const ham_float cos_p{std::cos(par->bmodel_jaffe.arm_pitch)};
  // pitch angle
  const ham_float sin_p{std::sin(par->bmodel_jaffe.arm_pitch)};
  Hamvec<3, ham_float> tmp;
  ham_float quadruple{1};
  // forbiden region
  if (r < 0.5 * cgs::kpc)
    return tmp;
  if (pos[2] > par->bmodel_jaffe.disk_z0)
    quadruple = (1 - 2 * par->bmodel_jaffe.quadruple);
  // molecular ring
  if (par->bmodel_jaffe.ring) {
    // inside spiral arm
    if (r > r_lim) {
      // sin(t-p)
      tmp[0] = (cos_p * (pos[1] / r) - sin_p * (pos[0] / r)) * quadruple;
      //-cos(t-p)
      tmp[1] = (-cos_p * (pos[0] / r) - sin_p * (pos[1] / r)) * quadruple;
    }
    // inside molecular ring
    else {
      // sin(phi)
      tmp[0] = (1 - 2 * par->bmodel_jaffe.bss) * pos[1] / r;
      //-cos(phi)
      tmp[1] = (2 * par->bmodel_jaffe.bss - 1) * pos[0] / r;
    }
  }
  // elliptical bar (replace molecular ring)
  else if (par->bmodel_jaffe.bar) {
    const ham_float cos_phi{std::cos(par->bmodel_jaffe.bar_phi0)};
    const ham_float sin_phi{std::sin(par->bmodel_jaffe.bar_phi0)};
    const ham_float x{cos_phi * pos[0] - sin_phi * pos[1]};
    const ham_float y{sin_phi * pos[0] + cos_phi * pos[1]};
    // inside spiral arm
    if (r > bar_lim) {
      // sin(t-p)
      tmp[0] = (cos_p * (pos[1] / r) - sin_p * (pos[0] / r)) * quadruple;
      //-cos(t-p)
      tmp[1] = (-cos_p * (pos[0] / r) - sin_p * (pos[1] / r)) * quadruple;
    }
    // inside elliptical bar
    else {
      if (y != 0) {
        const ham_float new_x{copysign(1, y)};
        const ham_float new_y{
            -copysign(1, y) * (x / y) * par->bmodel_jaffe.bar_b *
            par->bmodel_jaffe.bar_b /
            (par->bmodel_jaffe.bar_a * par->bmodel_jaffe.bar_a)};
        tmp[0] = (cos_phi * new_x + sin_phi * new_y) *
                 (1 - 2 * par->bmodel_jaffe.bss);
        tmp[1] = (-sin_phi * new_x + cos_phi * new_y) *
                 (1 - 2 * par->bmodel_jaffe.bss);
        tmp = tmp.versor();
      } else {
        tmp[0] = (2 * par->bmodel_jaffe.bss - 1) * copysign(1, x) * sin_phi;
        tmp[1] = (2 * par->bmodel_jaffe.bss - 1) * copysign(1, x) * cos_phi;
      }
    }
  }
  return tmp;
}

ham_float Bmodel_jaffe::radial_scaling(const Hamvec<3, ham_float> &pos,
                                       const Param *par) const {
  const ham_float r2{pos[0] * pos[0] + pos[1] * pos[1]};
  // separate into 3 parts for better view
  const ham_float s1{1. - std::exp(-r2 / (par->bmodel_jaffe.r_inner *
                                          par->bmodel_jaffe.r_inner))};
  const ham_float s2{
      std::exp(-r2 / (par->bmodel_jaffe.r_scale * par->bmodel_jaffe.r_scale))};
  const ham_float s3{
      std::exp(-r2 * r2 /
               (par->bmodel_jaffe.r_peak * par->bmodel_jaffe.r_peak *
                par->bmodel_jaffe.r_peak * par->bmodel_jaffe.r_peak))};
  return s1 * (s2 + s3);
}

std::vector<ham_float>
Bmodel_jaffe::arm_compress(const Hamvec<3, ham_float> &pos,
                           const Param *par) const {
  const ham_float r{sqrt(pos[0] * pos[0] + pos[1] * pos[1]) /
                    par->bmodel_jaffe.comp_r};
  const ham_float c0{1. / par->bmodel_jaffe.comp_c - 1.};
  std::vector<ham_float> a0 = dist2arm(pos, par);
  const ham_float r_scaling{radial_scaling(pos, par)};
  const ham_float z_scaling{arm_scaling(pos, par)};
  // for saving computing time
  const ham_float d0_inv{(r_scaling * z_scaling) / par->bmodel_jaffe.comp_d};
  ham_float factor{c0 * r_scaling * z_scaling};
  if (r > 1) {
    ham_float cdrop{std::pow(r, -par->bmodel_jaffe.comp_p)};
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
      a0[i] = factor * cdrop *
              std::exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
    }
  } else {
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
      a0[i] = factor * std::exp(-a0[i] * a0[i] * d0_inv * d0_inv);
    }
  }
  return a0;
}

std::vector<ham_float>
Bmodel_jaffe::arm_compress_dust(const Hamvec<3, ham_float> &pos,
                                const Param *par) const {
  const ham_float r{sqrt(pos[0] * pos[0] + pos[1] * pos[1]) /
                    par->bmodel_jaffe.comp_r};
  const ham_float c0{1. / par->bmodel_jaffe.comp_c - 1.};
  std::vector<ham_float> a0 = dist2arm(pos, par);
  const ham_float r_scaling{radial_scaling(pos, par)};
  const ham_float z_scaling{arm_scaling(pos, par)};
  // only difference from normal arm_compress
  const ham_float d0_inv{(r_scaling) / par->bmodel_jaffe.comp_d};
  ham_float factor{c0 * r_scaling * z_scaling};
  if (r > 1) {
    ham_float cdrop{std::pow(r, -par->bmodel_jaffe.comp_p)};
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
      a0[i] = factor * cdrop *
              std::exp(-a0[i] * a0[i] * cdrop * cdrop * d0_inv * d0_inv);
    }
  } else {
    for (decltype(a0.size()) i = 0; i < a0.size(); ++i) {
      a0[i] = factor * std::exp(-a0[i] * a0[i] * d0_inv * d0_inv);
    }
  }
  return a0;
}

std::vector<ham_float> Bmodel_jaffe::dist2arm(const Hamvec<3, ham_float> &pos,
                                              const Param *par) const {
  std::vector<ham_float> d;
  const ham_float r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  const ham_float r_lim{par->bmodel_jaffe.ring_r};
  const ham_float bar_lim{par->bmodel_jaffe.bar_a +
                          0.5 * par->bmodel_jaffe.comp_d};
  const ham_float cos_p{std::cos(par->bmodel_jaffe.arm_pitch)};
  // pitch angle
  const ham_float sin_p{std::sin(par->bmodel_jaffe.arm_pitch)};
  const ham_float beta_inv{-sin_p / cos_p};
  ham_float theta{atan2(pos[1], pos[0])};
  if (theta < 0)
    theta += 2 * cgs::pi;
  // if molecular ring
  if (par->bmodel_jaffe.ring) {
    // in molecular ring, return single element vector
    if (r < r_lim) {
      d.push_back(std::fabs(par->bmodel_jaffe.ring_r - r));
    }
    // in spiral arm, return vector with arm_num elements
    else {
      // loop through arms
      for (ham_uint i = 0; i < par->bmodel_jaffe.arm_num; ++i) {
        ham_float d_ang{par->bmodel_jaffe.arm_phi0[i] - theta};
        ham_float d_rad{std::fabs(
            par->bmodel_jaffe.arm_r0 * std::exp(d_ang * beta_inv) - r)};
        ham_float d_rad_p{
            std::fabs(par->bmodel_jaffe.arm_r0 *
                          std::exp((d_ang + 2 * cgs::pi) * beta_inv) -
                      r)};
        ham_float d_rad_m{
            std::fabs(par->bmodel_jaffe.arm_r0 *
                          std::exp((d_ang - 2 * cgs::pi) * beta_inv) -
                      r)};
        d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
      }
    }
  }
  // if elliptical bar
  else if (par->bmodel_jaffe.bar) {
    // cos(phi)cos(phi0) - sin(phi)sin(phi0)
    const ham_float cos_tmp{std::cos(par->bmodel_jaffe.bar_phi0) * pos[0] / r -
                            std::sin(par->bmodel_jaffe.bar_phi0) * pos[1] / r};
    // sin(phi)cos(phi0) + cos(phi)sin(phi0)
    const ham_float sin_tmp{std::cos(par->bmodel_jaffe.bar_phi0) * pos[1] / r +
                            std::sin(par->bmodel_jaffe.bar_phi0) * pos[0] / r};
    // in bar, return single element vector
    if (r < bar_lim) {
      d.push_back(
          std::fabs(par->bmodel_jaffe.bar_a * par->bmodel_jaffe.bar_b /
                        sqrt(par->bmodel_jaffe.bar_a * par->bmodel_jaffe.bar_a *
                                 sin_tmp * sin_tmp +
                             par->bmodel_jaffe.bar_b * par->bmodel_jaffe.bar_b *
                                 cos_tmp * cos_tmp) -
                    r));
    }
    // in spiral arm, return vector with arm_num elements
    else {
      // loop through arms
      for (ham_uint i = 0; i < par->bmodel_jaffe.arm_num; ++i) {
        ham_float d_ang{par->bmodel_jaffe.arm_phi0[i] - theta};
        ham_float d_rad{std::fabs(
            par->bmodel_jaffe.arm_r0 * std::exp(d_ang * beta_inv) - r)};
        ham_float d_rad_p{
            std::fabs(par->bmodel_jaffe.arm_r0 *
                          std::exp((d_ang + 2 * cgs::pi) * beta_inv) -
                      r)};
        ham_float d_rad_m{
            std::fabs(par->bmodel_jaffe.arm_r0 *
                          std::exp((d_ang - 2 * cgs::pi) * beta_inv) -
                      r)};
        d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
      }
    }
  }
  return d;
}
