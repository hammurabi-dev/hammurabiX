#include <bfield.h>
#include <cmath>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>
#include <param.h>

Hamvec<3, ham_float> Breg_jf12::write_field(const Hamvec<3, ham_float> &pos,
                                            const Param *par) const {

  // input parameters
  //----------------------------------------------
  // disk parameters
  const ham_float b_arm_1{par->breg_jf12.b_arm1};
  const ham_float b_arm_2{par->breg_jf12.b_arm2};
  const ham_float b_arm_3{par->breg_jf12.b_arm3};
  const ham_float b_arm_4{par->breg_jf12.b_arm4};
  const ham_float b_arm_5{par->breg_jf12.b_arm5};
  const ham_float b_arm_6{par->breg_jf12.b_arm6};
  const ham_float b_arm_7{par->breg_jf12.b_arm7};
  const ham_float b_ring{par->breg_jf12.b_ring};
  const ham_float h_disk{par->breg_jf12.h_disk};
  const ham_float w_disk{par->breg_jf12.w_disk};

  // toroidal halo parameters
  const ham_float Bn{par->breg_jf12.Bn};
  const ham_float Bs{par->breg_jf12.Bs};
  const ham_float rn{par->breg_jf12.rn};
  const ham_float rs{par->breg_jf12.rs};
  const ham_float wh{par->breg_jf12.wh};
  const ham_float z0{par->breg_jf12.z0};

  // X-field parameters
  const ham_float B0_X{par->breg_jf12.B0_X};
  const ham_float Xtheta_const{(par->breg_jf12.Xtheta)};
  const ham_float rpc_X{par->breg_jf12.rpc_X};
  const ham_float r0_X{par->breg_jf12.r0_X};
  //---------------------------------------------

  // define fixed parameters
  const ham_float Rmax = 20 * cgs::kpc;   // outer boundary of GMF
  const ham_float rho_GC = 1. * cgs::kpc; // interior boundary of GMF

  // fixed disk parameters
  const ham_float inc = 11.5; // inclination, in degrees
  const ham_float rmin =
      5. * cgs::kpc; // outer boundary of the molecular ring region
  const ham_float rcent =
      3. * cgs::kpc; // inner boundary of the molecular ring region (field is
                     // zero within this region)
  const ham_float f[8] = {
      0.130, 0.165, 0.094, 0.122,
      0.13,  0.118, 0.084, 0.156}; // fractions of circumference spanned by each
                                   // spiral, sums to unity
  const ham_float rc_B[8] = {
      5.1,  6.3,  7.1, 8.3, 9.8,
      11.4, 12.7, 15.5}; // the radii where the spiral arm boundaries cross the
                         // negative x-axis

  const ham_float r{sqrt(pos[0] * pos[0] + pos[1] * pos[1])};
  const ham_float rho{
      sqrt(pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2])};
  const ham_float phi{atan2(pos[1], pos[0])};
  const ham_float z{pos[2]};

  // define boundaries for where magnetic field is zero (outside of galaxy)
  if (r > Rmax || rho < rho_GC) {
    return Hamvec<3, ham_float>{1.e-12, 1.e-12, 1.e-12};
  }

  //------------------------------------------------------------------------------
  // DISK COMPONENT (8 spiral regions, 7 free parameters with 8th set to
  // conserve flux)
  // B0 set to 1 at r=5kpc
  const ham_float B0 = (rmin / r); //
  // the logistic equation, to be multiplied to the toroidal halo field and
  // (1-zprofile) multiplied to the disk:
  const ham_float zprofile{1. /
                           (1 + exp(-2. / w_disk * (std::fabs(z) - h_disk)))};

  // printf("%g, %g \n", z, zprofile);
  ham_float B_cyl_disk[3] = {0, 0,
                             0}; // the disk field in cylindrical coordinates

  if ((r > rcent)) // disk field zero elsewhere
  {
    if (r < rmin) { // circular field in molecular ring
      B_cyl_disk[1] = B0 * b_ring * (1 - zprofile);
    } else {
      // use flux conservation to calculate the field strength in the 8th spiral
      // arm
      ham_float bv_B[8] = {b_arm_1, b_arm_2, b_arm_3, b_arm_4,
                           b_arm_5, b_arm_6, b_arm_7, 0.};
      ham_float b8 = 0.;

      for (int i = 0; i < 7; i++) {
        b8 -= f[i] * bv_B[i] / f[7];
      }
      bv_B[7] = b8;

      // iteratively figure out which spiral arm the current coordinates (r.phi)
      // correspond to
      ham_float b_disk = 0.;
      ham_float r_negx =
          r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi - M_PI));

      if (r_negx > rc_B[7] * cgs::kpc) {
        r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + M_PI));
      }
      if (r_negx > rc_B[7] * cgs::kpc) {
        r_negx = r * exp(-1 / tan(M_PI / 180. * (90 - inc)) * (phi + 3 * M_PI));
      }
      for (int i = 7; i >= 0; i--) {
        if (r_negx < rc_B[i] * cgs::kpc) {
          b_disk = bv_B[i];
        }
      } // "region 8,7,6,..,2"

      B_cyl_disk[0] = b_disk * B0 * sin(M_PI / 180. * inc) * (1 - zprofile);
      B_cyl_disk[1] = b_disk * B0 * cos(M_PI / 180. * inc) * (1 - zprofile);
    }
  }

  //-------------------------------------------------------------------------
  ////TOROIDAL HALO COMPONENT

  ham_float b1, rh;
  ham_float B_h = 0.;

  if (z >= 0) { // North
    b1 = Bn;
    rh = rn; // transition radius between inner-outer region
  } else {   // South
    b1 = Bs;
    rh = rs;
  }

  B_h = b1 * (1. - 1. / (1. + exp(-2. / wh * (r - rh)))) *
        exp(-(std::fabs(z)) / (z0)); // vertical exponential fall-off
  const ham_float B_cyl_h[3] = {0., B_h * zprofile, 0.};

  //------------------------------------------------------------------------
  // X- FIELD

  ham_float Xtheta = 0.;
  ham_float rp_X =
      0.; // the mid-plane radius for the field line that pass through r
  ham_float B_X = 0.;
  ham_float r_sign = 1.; // +1 for north, -1 for south
  if (z < 0) {
    r_sign = -1.;
  }

  // dividing line between region with constat elevation angle, and the
  // interior:
  ham_float rc_X = rpc_X + std::fabs(z) / tan(Xtheta_const);

  if (r < rc_X) { // interior region, with varying elevation angle
    rp_X = r * rpc_X / rc_X;
    B_X = B0_X * pow(rpc_X / rc_X, 2.) * exp(-rp_X / r0_X);
    Xtheta = atan(std::abs(z) /
                  (r - rp_X)); // modified elevation angle in interior region
    // printf("Xtheta %g at z %g , r %g , rc_X %g \n",Xtheta, z,r,rc_X);
    if (z == 0.) {
      Xtheta = M_PI / 2.;
    }      // to avoid some NaN
  } else { // exterior region with constant elevation angle
    Xtheta = Xtheta_const;
    rp_X = r - std::abs(z) / tan(Xtheta);
    B_X = B0_X * rp_X / r * exp(-rp_X / r0_X);
  }

  // X-field in cylindrical coordinates
  ham_float B_cyl_X[3] = {B_X * cos(Xtheta) * r_sign, 0., B_X * sin(Xtheta)};

  // add fields together
  ham_float B_cyl[3] = {0.0, 0.0, 0.0};
  B_cyl[0] = B_cyl_disk[0] + B_cyl_h[0] + B_cyl_X[0];
  B_cyl[1] = B_cyl_disk[1] + B_cyl_h[1] + B_cyl_X[1];
  B_cyl[2] = B_cyl_disk[2] + B_cyl_h[2] + B_cyl_X[2];

  // convert field to cartesian coordinates
  ham_float B_cart[3] = {0.0, 0.0, 0.0};
  B_cart[0] = B_cyl[0] * cos(phi) - B_cyl[1] * sin(phi);
  B_cart[1] = B_cyl[0] * sin(phi) + B_cyl[1] * cos(phi);
  B_cart[2] = B_cyl[2];

  return Hamvec<3, ham_float>{B_cart[0], B_cart[1], B_cart[2]};
}
