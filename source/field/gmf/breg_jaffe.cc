#include <cmath>
#include <vector>

#include <hvec.h>

#include <breg.h>
#include <cgs_units_file.h>
#include <grid.h>
#include <namespace_toolkit.h>
#include <param.h>

hvec<3, double> Breg_jaffe::breg(const hvec<3, double> &pos,
                                 const Param *par) const {
  double inner_b{0};
  if (par->breg_jaffe.ring)
    inner_b = par->breg_jaffe.ring_amp;
  else if (par->breg_jaffe.bar)
    inner_b = par->breg_jaffe.bar_amp;

  hvec<3, double> bhat{orientation(pos, par)};
  hvec<3, double> btot;
  btot = bhat * radial_scaling(pos, par) *
         (par->breg_jaffe.disk_amp * disk_scaling(pos, par) +
          par->breg_jaffe.halo_amp * halo_scaling(pos, par));
  // compress factor for each arm or for ring/bar
  std::vector<double> arm{arm_compress(pos, par)};
  // only inner region
  if (arm.size() == 1) {
    btot += bhat * arm[0] * inner_b;
  }
  // spiral arm region
  else {
    for (decltype(arm.size()) i = 0; i < arm.size(); ++i) {
      btot += bhat * arm[i] * par->breg_jaffe.arm_amp[i];
    }
  }
  return btot;
}

hvec<3, double> Breg_jaffe::orientation(const hvec<3, double> &pos,
                                        const Param *par) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])}; // cylindrical frame
  const double r_lim{par->breg_jaffe.ring_r};
  const double bar_lim{par->breg_jaffe.bar_a + 0.5 * par->breg_jaffe.comp_d};
  const double cos_p{std::cos(par->breg_jaffe.arm_pitch)};
  const double sin_p{std::sin(par->breg_jaffe.arm_pitch)}; // pitch angle
  hvec<3, double> tmp;
  double quadruple{1};
  if (r < 0.5 * CGS_U_kpc) // forbiden region
    return tmp;
  if (pos[2] > par->breg_jaffe.disk_z0)
    quadruple = (1 - 2 * par->breg_jaffe.quadruple);
  // molecular ring
  if (par->breg_jaffe.ring) {
    // inside spiral arm
    if (r > r_lim) {
      tmp[0] =
          (cos_p * (pos[1] / r) - sin_p * (pos[0] / r)) * quadruple; // sin(t-p)
      tmp[1] = (-cos_p * (pos[0] / r) - sin_p * (pos[1] / r)) *
               quadruple; //-cos(t-p)
    }
    // inside molecular ring
    else {
      tmp[0] = (1 - 2 * par->breg_jaffe.bss) * pos[1] / r; // sin(phi)
      tmp[1] = (2 * par->breg_jaffe.bss - 1) * pos[0] / r; //-cos(phi)
    }
  }
  // elliptical bar (replace molecular ring)
  else if (par->breg_jaffe.bar) {
    const double cos_phi{std::cos(par->breg_jaffe.bar_phi0)};
    const double sin_phi{std::sin(par->breg_jaffe.bar_phi0)};
    const double x{cos_phi * pos[0] - sin_phi * pos[1]};
    const double y{sin_phi * pos[0] + cos_phi * pos[1]};
    // inside spiral arm
    if (r > bar_lim) {
      tmp[0] =
          (cos_p * (pos[1] / r) - sin_p * (pos[0] / r)) * quadruple; // sin(t-p)
      tmp[1] = (-cos_p * (pos[0] / r) - sin_p * (pos[1] / r)) *
               quadruple; //-cos(t-p)
    }
    // inside elliptical bar
    else {
      if (y != 0) {
        const double new_x{copysign(1, y)};
        const double new_y{-copysign(1, y) * (x / y) * par->breg_jaffe.bar_b *
                           par->breg_jaffe.bar_b /
                           (par->breg_jaffe.bar_a * par->breg_jaffe.bar_a)};
        tmp[0] =
            (cos_phi * new_x + sin_phi * new_y) * (1 - 2 * par->breg_jaffe.bss);
        tmp[1] = (-sin_phi * new_x + cos_phi * new_y) *
                 (1 - 2 * par->breg_jaffe.bss);
        tmp = tmp.versor();
      } else {
        tmp[0] = (2 * par->breg_jaffe.bss - 1) * copysign(1, x) * sin_phi;
        tmp[1] = (2 * par->breg_jaffe.bss - 1) * copysign(1, x) * cos_phi;
      }
    }
  }
  return tmp;
}

double Breg_jaffe::radial_scaling(const hvec<3, double> &pos,
                                  const Param *par) const {
  const double r2{pos[0] * pos[0] + pos[1] * pos[1]};
  // separate into 3 parts for better view
  const double s1{
      1. - std::exp(-r2 / (par->breg_jaffe.r_inner * par->breg_jaffe.r_inner))};
  const double s2{
      std::exp(-r2 / (par->breg_jaffe.r_scale * par->breg_jaffe.r_scale))};
  const double s3{std::exp(-r2 * r2 /
                           (par->breg_jaffe.r_peak * par->breg_jaffe.r_peak *
                            par->breg_jaffe.r_peak * par->breg_jaffe.r_peak))};
  return s1 * (s2 + s3);
}

std::vector<double> Breg_jaffe::arm_compress(const hvec<3, double> &pos,
                                             const Param *par) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1]) /
                 par->breg_jaffe.comp_r};
  const double c0{1. / par->breg_jaffe.comp_c - 1.};
  std::vector<double> a0 = dist2arm(pos, par);
  const double r_scaling{radial_scaling(pos, par)};
  const double z_scaling{arm_scaling(pos, par)};
  // for saving computing time
  const double d0_inv{(r_scaling * z_scaling) / par->breg_jaffe.comp_d};
  double factor{c0 * r_scaling * z_scaling};
  if (r > 1) {
    double cdrop{std::pow(r, -par->breg_jaffe.comp_p)};
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

std::vector<double> Breg_jaffe::arm_compress_dust(const hvec<3, double> &pos,
                                                  const Param *par) const {
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1]) /
                 par->breg_jaffe.comp_r};
  const double c0{1. / par->breg_jaffe.comp_c - 1.};
  std::vector<double> a0 = dist2arm(pos, par);
  const double r_scaling{radial_scaling(pos, par)};
  const double z_scaling{arm_scaling(pos, par)};
  // only difference from normal arm_compress
  const double d0_inv{(r_scaling) / par->breg_jaffe.comp_d};
  double factor{c0 * r_scaling * z_scaling};
  if (r > 1) {
    double cdrop{std::pow(r, -par->breg_jaffe.comp_p)};
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

std::vector<double> Breg_jaffe::dist2arm(const hvec<3, double> &pos,
                                         const Param *par) const {
  std::vector<double> d;
  const double r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  const double r_lim{par->breg_jaffe.ring_r};
  const double bar_lim{par->breg_jaffe.bar_a + 0.5 * par->breg_jaffe.comp_d};
  const double cos_p{std::cos(par->breg_jaffe.arm_pitch)};
  const double sin_p{std::sin(par->breg_jaffe.arm_pitch)}; // pitch angle
  const double beta_inv{-sin_p / cos_p};
  double theta{atan2(pos[1], pos[0])};
  if (theta < 0)
    theta += 2 * CGS_U_pi;
  // if molecular ring
  if (par->breg_jaffe.ring) {
    // in molecular ring, return single element vector
    if (r < r_lim) {
      d.push_back(std::fabs(par->breg_jaffe.ring_r - r));
    }
    // in spiral arm, return vector with arm_num elements
    else {
      // loop through arms
      for (unsigned int i = 0; i < par->breg_jaffe.arm_num; ++i) {
        double d_ang{par->breg_jaffe.arm_phi0[i] - theta};
        double d_rad{
            std::fabs(par->breg_jaffe.arm_r0 * std::exp(d_ang * beta_inv) - r)};
        double d_rad_p{
            std::fabs(par->breg_jaffe.arm_r0 *
                          std::exp((d_ang + 2 * CGS_U_pi) * beta_inv) -
                      r)};
        double d_rad_m{
            std::fabs(par->breg_jaffe.arm_r0 *
                          std::exp((d_ang - 2 * CGS_U_pi) * beta_inv) -
                      r)};
        d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
      }
    }
  }
  // if elliptical bar
  else if (par->breg_jaffe.bar) {
    const double cos_tmp{std::cos(par->breg_jaffe.bar_phi0) * pos[0] / r -
                         std::sin(par->breg_jaffe.bar_phi0) * pos[1] /
                             r}; // cos(phi)cos(phi0) - sin(phi)sin(phi0)
    const double sin_tmp{std::cos(par->breg_jaffe.bar_phi0) * pos[1] / r +
                         std::sin(par->breg_jaffe.bar_phi0) * pos[0] /
                             r}; // sin(phi)cos(phi0) + cos(phi)sin(phi0)
    // in bar, return single element vector
    if (r < bar_lim) {
      d.push_back(
          std::fabs(par->breg_jaffe.bar_a * par->breg_jaffe.bar_b /
                        sqrt(par->breg_jaffe.bar_a * par->breg_jaffe.bar_a *
                                 sin_tmp * sin_tmp +
                             par->breg_jaffe.bar_b * par->breg_jaffe.bar_b *
                                 cos_tmp * cos_tmp) -
                    r));
    }
    // in spiral arm, return vector with arm_num elements
    else {
      // loop through arms
      for (unsigned int i = 0; i < par->breg_jaffe.arm_num; ++i) {
        double d_ang{par->breg_jaffe.arm_phi0[i] - theta};
        double d_rad{
            std::fabs(par->breg_jaffe.arm_r0 * std::exp(d_ang * beta_inv) - r)};
        double d_rad_p{
            std::fabs(par->breg_jaffe.arm_r0 *
                          std::exp((d_ang + 2 * CGS_U_pi) * beta_inv) -
                      r)};
        double d_rad_m{
            std::fabs(par->breg_jaffe.arm_r0 *
                          std::exp((d_ang - 2 * CGS_U_pi) * beta_inv) -
                      r)};
        d.push_back(std::min(std::min(d_rad, d_rad_p), d_rad_m) * cos_p);
      }
    }
  }
  return d;
}

// END
