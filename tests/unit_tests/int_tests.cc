// unit tests for Integrator class member functions
// feel free to add more rational testing blocks

#include <gtest/gtest.h>

#include <cgs_units_file.h>
#include <cmath>
#include <grid.h>
#include <integrator.h>
#include <memory>

// testing:
// Integrator::check_simulation_upper_limit
// Integrator::check_simulation_lower_limit
TEST(integrator, boundary_check) {
  Integrator pipe;
  double R0 = 10.;
  double R_lim = (1.0 + 1.0e-5) * R0;
  EXPECT_FALSE(pipe.check_simulation_upper_limit(R0, R_lim));
  EXPECT_TRUE(pipe.check_simulation_lower_limit(R0, R_lim));
}

// testing:
// Integrator::assemble_shell_ref
TEST(integrator, shell_info_assembling) {
  unsigned int test_unsigned;
  double test_double;
  auto pipe = std::make_unique<Integrator>();
  auto ref = std::make_unique<Integrator::struct_shell>();
  //
  auto par = std::make_unique<Param>("reference/int_tests_01.xml");
  //
  test_unsigned = 2;
  EXPECT_EQ(par->grid_obs.sim_sync_freq.size(), test_unsigned);
  test_double = 23.0 * CGS_U_GHz;
  EXPECT_EQ(par->grid_obs.sim_sync_freq[0], test_double);
  test_double = 1.4 * CGS_U_GHz;
  EXPECT_EQ(par->grid_obs.sim_sync_freq[1], test_double);
  test_unsigned = 32;
  EXPECT_EQ(par->grid_obs.nside_shell[0], test_unsigned);
  test_unsigned = 64;
  EXPECT_EQ(par->grid_obs.nside_shell[1], test_unsigned);
  test_unsigned = 128;
  EXPECT_EQ(par->grid_obs.nside_shell[2], test_unsigned);
  //
  pipe->assemble_shell_ref(ref.get(), par.get(), 1);
  test_unsigned = 333;
  EXPECT_EQ(ref->step, test_unsigned);
  test_double = 5 * CGS_U_kpc;
  EXPECT_EQ(ref->d_start, test_double);
  test_double = 40 * CGS_U_kpc / 4 + 5 * CGS_U_kpc;
  EXPECT_EQ(ref->d_stop, test_double);
  //
  pipe->assemble_shell_ref(ref.get(), par.get(), 3);
  test_unsigned = 666;
  EXPECT_EQ(ref->step, test_unsigned);
  test_double = 40 * CGS_U_kpc / 2 + 5 * CGS_U_kpc;
  EXPECT_EQ(ref->d_start, test_double);
  test_double = 45 * CGS_U_kpc;
  EXPECT_EQ(ref->d_stop, test_double);
  //
  par = std::make_unique<Param>("reference/int_tests_02.xml");
  //
  test_unsigned = 32;
  EXPECT_EQ(par->grid_obs.nside_shell[0], test_unsigned);
  test_unsigned = 16;
  EXPECT_EQ(par->grid_obs.nside_shell[1], test_unsigned);
  test_unsigned = 8;
  EXPECT_EQ(par->grid_obs.nside_shell[2], test_unsigned);
  //
  pipe->assemble_shell_ref(ref.get(), par.get(), 1);
  test_unsigned = 100;
  EXPECT_EQ(ref->step, test_unsigned);
  test_double = 0.;
  EXPECT_EQ(ref->d_start, test_double);
  test_double = par->grid_obs.oc_r_max * 0.1;
  EXPECT_EQ(ref->d_stop, test_double);
  //
  pipe->assemble_shell_ref(ref.get(), par.get(), 2);
  test_unsigned = 600;
  EXPECT_EQ(ref->step, test_unsigned);
  test_double = par->grid_obs.oc_r_max * 0.1;
  EXPECT_EQ(ref->d_start, test_double);
  test_double = par->grid_obs.oc_r_max * 0.7;
  EXPECT_EQ(ref->d_stop, test_double);
  //
  pipe->assemble_shell_ref(ref.get(), par.get(), 3);
  test_unsigned = 300;
  EXPECT_EQ(ref->step, test_unsigned);
  test_double = par->grid_obs.oc_r_max * 0.7;
  EXPECT_EQ(ref->d_start, test_double);
  test_double = par->grid_obs.oc_r_max * 1.0;
  EXPECT_EQ(ref->d_stop, test_double);
}

// testing:
// los_versor
TEST(toolkit, los_versor) {
  auto pipe = std::make_unique<Integrator>();
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};

  EXPECT_LT(std::fabs(pipe->los_versor(theta[0], phi[0])[0] - 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[0], phi[0])[1] - 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[0], phi[0])[2] - 1.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[1], phi[1])[0] + 1.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[1], phi[1])[1] - 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[1], phi[1])[2] - 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[2], phi[2])[0] - 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[2], phi[2])[1] + 1.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_versor(theta[2], phi[2])[2] - 0.), 1e-10);
}

// testing:
// los_parproj
TEST(toolkit, los_parproj) {
  auto pipe = std::make_unique<Integrator>();
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};
  const hamvec<3, double> A{1., 0., 0.};

  EXPECT_LT(std::fabs(pipe->los_parproj(A, theta[0], phi[0]) - 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_parproj(A, theta[1], phi[1]) + 1.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_parproj(A, theta[2], phi[2]) - 0.), 1e-10);
}

// testing:
// los_perproj
TEST(toolkit, los_perproj) {
  auto pipe = std::make_unique<Integrator>();
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};
  const hamvec<3, double> A{1., 0., 0.};

  EXPECT_LT(std::fabs(pipe->los_perproj(A, theta[0], phi[0]) - 1.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_perproj(A, theta[1], phi[1]) + 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_perproj(A, theta[2], phi[2]) - 1.), 1e-10);
}

// testing:
// sync_ipa
TEST(toolkit, sync_ipa) {
  auto pipe = std::make_unique<Integrator>();
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};
  const hamvec<3, double> A{1., 0., 0.};

  EXPECT_LT(std::fabs(pipe->sync_ipa(A, theta[0], phi[0]) + 90. * CGS_U_rad),
            1e-10);
  EXPECT_LT(std::fabs(pipe->sync_ipa(A, theta[2], phi[2]) - 180. * CGS_U_rad),
            1e-10);
}
