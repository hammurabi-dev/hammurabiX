#include <cmath>
#include <grid.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>
#include <tefield.h>

ham_float TEmodel_ymw16::read_model(const Hamvec<3, ham_float> &pos,
                                    const Param *par) const {
  // YMW16 using a different Cartesian frame from our default one
  Hamvec<3, ham_float> gc_pos{pos[1], -pos[0], pos[2]};
  // sylindrical r
  ham_float r_cyl{std::sqrt(gc_pos[0] * gc_pos[0] + gc_pos[1] * gc_pos[1])};
  // warp
  if (r_cyl >= par->temodel_ymw16.r_warp) {
    ham_float theta_warp{std::atan2(gc_pos[1], gc_pos[0])};
    gc_pos[2] -= par->temodel_ymw16.t0_gamma_w *
                 (r_cyl - par->temodel_ymw16.r_warp) * std::cos(theta_warp);
  }
  if (gc_pos.length() > 25 * cgs::kpc) {
    return 0.;
  } else {
    ham_float ne{0.};
    ham_float ne_comp[8]{0.};
    ham_float weight_localbubble{0.};
    ham_float weight_gum{0.};
    ham_float weight_loop{0.};
    // longitude, in deg
    const ham_float ec_l{
        std::atan2(gc_pos[0], par->temodel_ymw16.r0 - gc_pos[1]) / cgs::rad};
    // call structure functions
    // since in YMW16, Fermi Bubble is not actually contributing, we ignore FB
    // for thick disk
    ne_comp[1] = thick(gc_pos[2], r_cyl, par);
    ne_comp[2] = thin(gc_pos[2], r_cyl, par);
    ne_comp[3] = spiral(gc_pos[0], gc_pos[1], gc_pos[2], r_cyl, par);
    ne_comp[4] = galcen(gc_pos[0], gc_pos[1], gc_pos[2], par);
    ne_comp[5] = gum(gc_pos[0], gc_pos[1], gc_pos[2], par);
    // localbubble boundary
    const ham_float localbubble_boundary{110. * cgs::pc};
    ne_comp[6] = localbubble(gc_pos[0], gc_pos[1], gc_pos[2], ec_l,
                             localbubble_boundary, par);
    ne_comp[7] = nps(gc_pos[0], gc_pos[1], gc_pos[2], par);
    // adding up rules
    ne_comp[0] = ne_comp[1] + std::max(ne_comp[2], ne_comp[3]);
    // distance to local bubble
    const ham_float rlb{std::sqrt(
        std::pow(((gc_pos[1] - 8.34 * cgs::kpc) * 0.94 - 0.34 * gc_pos[2]), 2) +
        gc_pos[0] * gc_pos[0])};
    if (rlb < localbubble_boundary) { // inside local bubble
      ne_comp[0] = par->temodel_ymw16.t6_j_lb * ne_comp[1] +
                   std::max(ne_comp[2], ne_comp[3]);
      if (ne_comp[6] > ne_comp[0]) {
        weight_localbubble = 1;
      }
    } else { // outside local bubble
      if (ne_comp[6] > ne_comp[0] and ne_comp[6] > ne_comp[5]) {
        weight_localbubble = 1;
      }
    }
    if (ne_comp[7] > ne_comp[0]) {
      weight_loop = 1;
    }
    if (ne_comp[5] > ne_comp[0]) {
      weight_gum = 1;
    }
    // final density
    ne =
        (1 - weight_localbubble) *
            ((1 - weight_gum) * ((1 - weight_loop) * (ne_comp[0] + ne_comp[4]) +
                                 weight_loop * ne_comp[7]) +
             weight_gum * ne_comp[5]) +
        (weight_localbubble) * (ne_comp[6]);
    assert(std::isfinite(ne));
    return ne;
  }
}

// thick disk
ham_float TEmodel_ymw16::thick(const ham_float &zz, const ham_float &rr,
                               const Param *par) const {
  if (zz > 10. * par->temodel_ymw16.t1_h1)
    return 0.; // timesaving
  ham_float gd{1.};
  if (rr > par->temodel_ymw16.t1_bd) {
    gd = std::pow(1. / std::cosh((rr - par->temodel_ymw16.t1_bd) /
                                 par->temodel_ymw16.t1_ad),
                  2);
  }
  return par->temodel_ymw16.t1_n1 * gd *
         std::pow(1. / std::cosh(zz / par->temodel_ymw16.t1_h1), 2);
}

// thin disk
ham_float TEmodel_ymw16::thin(const ham_float &zz, const ham_float &rr,
                              const Param *par) const {
  // z scaling, K_2*h0 in ref
  ham_float h0{par->temodel_ymw16.t2_k2 *
               (32 * cgs::pc + 1.6e-3 * rr + (4.e-7 / cgs::pc) * rr * rr)};
  if (zz > 10. * h0)
    return 0.; // timesaving
  ham_float gd{1.};
  if (rr > par->temodel_ymw16.t1_bd) {
    gd = std::pow(1. / std::cosh((rr - par->temodel_ymw16.t1_bd) /
                                 par->temodel_ymw16.t1_ad),
                  2);
  }
  return par->temodel_ymw16.t2_n2 * gd *
         std::pow(1. / std::cosh((rr - par->temodel_ymw16.t2_b2) /
                                 par->temodel_ymw16.t2_a2),
                  2) *
         std::pow(1. / std::cosh(zz / h0), 2);
}

// spiral arms
ham_float TEmodel_ymw16::spiral(const ham_float &xx, const ham_float &yy,
                                const ham_float &zz, const ham_float &rr,
                                const Param *par) const {
  // structure scaling
  ham_float scaling{1.};
  if (rr > par->temodel_ymw16.t1_bd) {
    if ((rr - par->temodel_ymw16.t1_bd) > 10. * par->temodel_ymw16.t1_ad)
      return 0.;
    scaling = std::pow(1. / std::cosh((rr - par->temodel_ymw16.t1_bd) /
                                      par->temodel_ymw16.t1_ad),
                       2);
  }
  // z scaling, K_a*h0 in ref
  const ham_float h0{
      par->temodel_ymw16.t3_ka *
      (32 * cgs::pc + 1.6e-3 * rr + (4.e-7 / cgs::pc) * pow(rr, 2))};
  if (zz > 10. * h0)
    return 0.; // timesaving
  scaling *= std::pow(1. / std::cosh(zz / h0), 2);
  if ((rr - par->temodel_ymw16.t3_b2s) > 10. * par->temodel_ymw16.t3_aa)
    return 0.; // timesaving
  // 2nd raidus scaling
  scaling *= std::pow(1. / std::cosh((rr - par->temodel_ymw16.t3_b2s) /
                                     par->temodel_ymw16.t3_aa),
                      2);
  ham_float smin;
  ham_float theta{std::atan2(yy, xx)};
  if (theta < 0)
    theta += 2 * cgs::pi;
  ham_float ne3s{0.};
  // looping through arms
  for (ham_uint i = 0; i < 4; ++i) {
    // get distance to arm center
    if (i != 4) {
      ham_float d_phi = theta - par->temodel_ymw16.t3_phimin[i];
      if (d_phi < 0)
        d_phi += 2. * cgs::pi;
      ham_float d =
          std::fabs(par->temodel_ymw16.t3_rmin[i] *
                        std::exp(d_phi * par->temodel_ymw16.t3_tpitch[i]) -
                    rr);
      ham_float d_p = std::fabs(par->temodel_ymw16.t3_rmin[i] *
                                    std::exp((d_phi + 2. * cgs::pi) *
                                             par->temodel_ymw16.t3_tpitch[i]) -
                                rr);
      smin = std::min(d, d_p) * par->temodel_ymw16.t3_cpitch[i];
    } else if (i == 4 and theta >= par->temodel_ymw16.t3_phimin[i] and
               theta < 2) { // Local arm
      smin = std::fabs(par->temodel_ymw16.t3_rmin[i] *
                           std::exp((theta + 2 * cgs::pi -
                                     par->temodel_ymw16.t3_phimin[i]) *
                                    par->temodel_ymw16.t3_tpitch[i]) -
                       rr) *
             par->temodel_ymw16.t3_cpitch[i];
    } else {
      continue;
    }
    if (smin > 10. * par->temodel_ymw16.t3_warm[i])
      continue; // timesaving
    // accumulate density
    if (i != 2) {
      ne3s += par->temodel_ymw16.t3_narm[i] * scaling *
              std::pow(1. / std::cosh(smin / par->temodel_ymw16.t3_warm[i]), 2);
    } else if (rr > 6 * cgs::kpc and
               theta * cgs::rad >
                   par->temodel_ymw16
                       .t3_thetacn) { // correction for Carina-Sagittarius
      const ham_float ga =
          (1. - (par->temodel_ymw16.t3_nsg) *
                    (std::exp(-std::pow(
                        (theta * cgs::rad - par->temodel_ymw16.t3_thetasg) /
                            par->temodel_ymw16.t3_wsg,
                        2)))) *
          (1. + par->temodel_ymw16.t3_ncn) *
          std::pow(1. / std::cosh(smin / par->temodel_ymw16.t3_warm[i]), 2);
      ne3s += par->temodel_ymw16.t3_narm[i] * scaling * ga;
    } else {
      const ham_float ga =
          (1. - (par->temodel_ymw16.t3_nsg) *
                    (std::exp(-std::pow(
                        (theta * cgs::rad - par->temodel_ymw16.t3_thetasg) /
                            par->temodel_ymw16.t3_wsg,
                        2)))) *
          (1. + par->temodel_ymw16.t3_ncn *
                    std::exp(-std::pow(
                        (theta * cgs::rad - par->temodel_ymw16.t3_thetacn) /
                            par->temodel_ymw16.t3_wcn,
                        2))) *
          std::pow(1. / std::cosh(smin / par->temodel_ymw16.t3_warm[i]), 2);
      ne3s += par->temodel_ymw16.t3_narm[i] * scaling * ga;
    }
  } // end of looping through arms
  return ne3s;
}

// galactic center
ham_float TEmodel_ymw16::galcen(const ham_float &xx, const ham_float &yy,
                                const ham_float &zz, const Param *par) const {
  // pos of center
  const ham_float Xgc{50. * cgs::pc};
  const ham_float Ygc{0.};
  const ham_float Zgc{-7. * cgs::pc};
  const ham_float R2gc{(xx - Xgc) * (xx - Xgc) + (yy - Ygc) * (yy - Ygc)};
  if (R2gc > 10. * par->temodel_ymw16.t4_agc * par->temodel_ymw16.t4_agc)
    return 0.; // timesaving
  const ham_float Ar{std::exp(
      -R2gc / (par->temodel_ymw16.t4_agc * par->temodel_ymw16.t4_agc))};
  if (std::fabs(zz - Zgc) > 10. * par->temodel_ymw16.t4_hgc)
    return 0.; // timesaving
  const ham_float Az{
      std::pow(1. / std::cosh((zz - Zgc) / par->temodel_ymw16.t4_hgc), 2)};
  return par->temodel_ymw16.t4_ngc * Ar * Az;
}

// gum nebula
ham_float TEmodel_ymw16::gum(const ham_float &xx, const ham_float &yy,
                             const ham_float &zz, const Param *par) const {
  if (yy < 0 or xx > 0)
    return 0.; // timesaving
  // center of Gum Nebula
  const ham_float lc{264. * cgs::rad};
  const ham_float bc{-4. * cgs::rad};
  const ham_float dc{450. * cgs::pc};
  const ham_float xc{dc * std::cos(bc) * std::sin(lc)};
  const ham_float yc{par->temodel_ymw16.r0 - dc * std::cos(bc) * std::cos(lc)};
  const ham_float zc{dc * std::sin(bc)};
  // theta is limited in I quadrant
  const ham_float theta{
      std::atan2(std::fabs(zz - zc),
                 std::sqrt((xx - xc) * (xx - xc) + (yy - yc) * (yy - yc)))};
  const ham_float tantheta = std::tan(theta);
  // zp is positive
  ham_float zp{(par->temodel_ymw16.t5_agn * par->temodel_ymw16.t5_agn *
                par->temodel_ymw16.t5_kgn * tantheta) /
               std::sqrt(par->temodel_ymw16.t5_agn * par->temodel_ymw16.t5_agn +
                         par->temodel_ymw16.t5_agn * par->temodel_ymw16.t5_agn *
                             par->temodel_ymw16.t5_kgn *
                             par->temodel_ymw16.t5_kgn)};
  // xyp is positive
  const ham_float xyp{zp / tantheta};
  // alpha is positive
  const ham_float xy_dist = {
      std::sqrt(par->temodel_ymw16.t5_agn * par->temodel_ymw16.t5_agn -
                xyp * xyp) *
      ham_float(par->temodel_ymw16.t5_agn > xyp)};
  const ham_float alpha{std::atan2(par->temodel_ymw16.t5_kgn * xyp, xy_dist) +
                        theta}; // add theta, timesaving
  const ham_float R2{(xx - xc) * (xx - xc) + (yy - yc) * (yy - yc) +
                     (zz - zc) * (zz - zc)};
  const ham_float r2{zp * zp + xyp * xyp};
  const ham_float D2min{(R2 + r2 - 2. * std::sqrt(R2 * r2)) * std::sin(alpha) *
                        std::sin(alpha)};
  if (D2min > 10. * par->temodel_ymw16.t5_wgn * par->temodel_ymw16.t5_wgn)
    return 0.;
  return par->temodel_ymw16.t5_ngn *
         std::exp(-D2min /
                  (par->temodel_ymw16.t5_wgn * par->temodel_ymw16.t5_wgn));
}

// local bubble
ham_float TEmodel_ymw16::localbubble(const ham_float &xx, const ham_float &yy,
                                     const ham_float &zz, const ham_float &ll,
                                     const ham_float &Rlb,
                                     const Param *par) const {
  if (yy < 0)
    return 0.; // timesaving
  ham_float nel{0.};
  // r_LB in ref
  const ham_float rLB{
      std::sqrt(std::pow(((yy - 8.34 * cgs::kpc) * 0.94 - 0.34 * zz), 2) +
                std::pow(xx, 2))};
  // l-l_LB1 in ref
  const ham_float dl1{
      std::min(std::fabs(ll + 360. - par->temodel_ymw16.t6_thetalb1),
               std::fabs(par->temodel_ymw16.t6_thetalb1 - (ll)))};
  if (dl1 < 10. * par->temodel_ymw16.t6_detlb1 or
      (rLB - Rlb) < 10. * par->temodel_ymw16.t6_wlb1 or
      zz < 10. * par->temodel_ymw16.t6_hlb1) // timesaving
    nel +=
        par->temodel_ymw16.t6_nlb1 *
        std::pow(1. / std::cosh(dl1 / par->temodel_ymw16.t6_detlb1), 2) *
        std::pow(1. / std::cosh((rLB - Rlb) / par->temodel_ymw16.t6_wlb1), 2) *
        std::pow(1. / std::cosh(zz / par->temodel_ymw16.t6_hlb1), 2);
  // l-l_LB2 in ref
  const ham_float dl2{
      std::min(std::fabs(ll + 360 - par->temodel_ymw16.t6_thetalb2),
               std::fabs(par->temodel_ymw16.t6_thetalb2 - (ll)))};
  if (dl2 < 10. * par->temodel_ymw16.t6_detlb2 or
      (rLB - Rlb) < 10. * par->temodel_ymw16.t6_wlb2 or
      zz < 10. * par->temodel_ymw16.t6_hlb2) // timesaving
    nel +=
        par->temodel_ymw16.t6_nlb2 *
        std::pow(1. / std::cosh(dl2 / par->temodel_ymw16.t6_detlb2), 2) *
        std::pow(1. / std::cosh((rLB - Rlb) / par->temodel_ymw16.t6_wlb2), 2) *
        std::pow(1. / std::cosh(zz / par->temodel_ymw16.t6_hlb2), 2);
  return nel;
}

// north spur
ham_float TEmodel_ymw16::nps(const ham_float &xx, const ham_float &yy,
                             const ham_float &zz, const Param *par) const {
  if (yy < 0)
    return 0.; // timesaving
  const ham_float theta_LI{(par->temodel_ymw16.t7_thetali) * cgs::rad};
  const ham_float x_c{-10.156 * cgs::pc};
  const ham_float y_c{8106.207 * cgs::pc};
  const ham_float z_c{10.467 * cgs::pc};
  // r_LI in ref
  const ham_float rLI{std::sqrt((xx - x_c) * (xx - x_c) +
                                (yy - y_c) * (yy - y_c) +
                                (zz - z_c) * (zz - z_c))};
  const ham_float theta{std::acos(((xx - x_c) * (std::cos(theta_LI)) +
                                   (zz - z_c) * (std::sin(theta_LI))) /
                                  rLI) /
                        cgs::rad};
  if (theta > 10. * par->temodel_ymw16.t7_detthetali or
      (rLI - par->temodel_ymw16.t7_rli) > 10. * par->temodel_ymw16.t7_wli)
    return 0.; // timesaving
  return (par->temodel_ymw16.t7_nli) *
         std::exp(-std::pow((rLI - par->temodel_ymw16.t7_rli) /
                                par->temodel_ymw16.t7_wli,
                            2)) *
         std::exp(-std::pow(theta / par->temodel_ymw16.t7_detthetali, 2));
}
