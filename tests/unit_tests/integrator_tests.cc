// unit tests for Integrator class

#include <gtest/gtest.h>

#include <cmath>
#include <grid.h>
#include <hamtype.h>
#include <hamunits.h>
#include <integrator.h>
#include <memory>
#include <random>

// testing:
// Integrator::check_simulation_upper_limit
// Integrator::check_simulation_lower_limit
TEST(integrator, boundary_check) {
  Integrator pipe;
  std::default_random_engine rng;
  std::uniform_real_distribution<ham_float> smp(0.0, 10.0);
  ham_float R0 = smp(rng);
  ham_float R_lim = (1.0 + 1.0e-5) * R0;
  EXPECT_FALSE(pipe.check_simulation_upper_limit(R0, R_lim));
  EXPECT_TRUE(pipe.check_simulation_lower_limit(R0, R_lim));
}

// testing:
// Integrator::assemble_shell_ref
TEST(integrator, shell_info_assembling) {
  ham_uint test_unsigned;
  ham_float test_ham_float;
  auto pipe = std::make_unique<Integrator>();
  auto ref = std::make_unique<Integrator::struct_shell>();
  // auto shell
  auto par = std::make_unique<Param>("reference/int_tests_01.xml");
  // frequency
  test_unsigned = 2;
  EXPECT_EQ(par->grid_obs.sim_sync_freq.size(), test_unsigned);
  test_ham_float = 23.0 * cgs::GHz;
  EXPECT_EQ(par->grid_obs.sim_sync_freq[0], test_ham_float);
  test_ham_float = 1.4 * cgs::GHz;
  EXPECT_EQ(par->grid_obs.sim_sync_freq[1], test_ham_float);
  // shell resolution
  test_unsigned = 32;
  EXPECT_EQ(par->grid_obs.nside_shell[0], test_unsigned);
  test_unsigned = 64;
  EXPECT_EQ(par->grid_obs.nside_shell[1], test_unsigned);
  test_unsigned = 128;
  EXPECT_EQ(par->grid_obs.nside_shell[2], test_unsigned);
  // shell info
  pipe->assemble_shell_ref(ref.get(), par.get(), 0); // 1st shell
  test_unsigned = 333;
  EXPECT_EQ(ref->step, test_unsigned);
  test_ham_float = 5 * cgs::kpc;
  EXPECT_EQ(ref->d_start, test_ham_float);
  test_ham_float = 40 * cgs::kpc / 4 + 5 * cgs::kpc;
  EXPECT_EQ(ref->d_stop, test_ham_float);
  test_ham_float = (10 * cgs::kpc) / 333;
  EXPECT_NEAR(ref->delta_d, test_ham_float, 1.0e-10 * cgs::kpc);
  // radial quadrature point
  std::default_random_engine rng;
  std::uniform_int_distribution<ham_int> smp(1, 333);
  ham_int idx = smp(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  idx = smp(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  // shell info
  pipe->assemble_shell_ref(ref.get(), par.get(), 2); // 3rd shell
  test_unsigned = 666;
  EXPECT_EQ(ref->step, test_unsigned);
  test_ham_float = 40 * cgs::kpc / 2 + 5 * cgs::kpc;
  EXPECT_EQ(ref->d_start, test_ham_float);
  test_ham_float = 45 * cgs::kpc;
  EXPECT_EQ(ref->d_stop, test_ham_float);
  test_ham_float = (20 * cgs::kpc) / 666;
  EXPECT_NEAR(ref->delta_d, test_ham_float, 1.0e-10 * cgs::kpc);
  // radial quadrature point
  std::uniform_int_distribution<ham_int> smp2(1, 666);
  idx = smp2(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  idx = smp2(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  // manual shell
  par = std::make_unique<Param>("reference/int_tests_02.xml");
  // shell resolution
  test_unsigned = 32;
  EXPECT_EQ(par->grid_obs.nside_shell[0], test_unsigned);
  test_unsigned = 16;
  EXPECT_EQ(par->grid_obs.nside_shell[1], test_unsigned);
  test_unsigned = 8;
  EXPECT_EQ(par->grid_obs.nside_shell[2], test_unsigned);
  // 1st shell
  pipe->assemble_shell_ref(ref.get(), par.get(), 0);
  test_unsigned = 100;
  EXPECT_EQ(ref->step, test_unsigned);
  test_ham_float = 0.;
  EXPECT_EQ(ref->d_start, test_ham_float);
  test_ham_float = par->grid_obs.oc_r_max * 0.1;
  EXPECT_EQ(ref->d_stop, test_ham_float);
  test_ham_float = (3 * cgs::kpc) / 100;
  EXPECT_NEAR(ref->delta_d, test_ham_float, 1.0e-10 * cgs::kpc);
  // radial quadrature point
  std::uniform_int_distribution<ham_int> smp3(1, 100);
  idx = smp3(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  idx = smp3(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  // 2nd shell
  pipe->assemble_shell_ref(ref.get(), par.get(), 1);
  test_unsigned = 600;
  EXPECT_EQ(ref->step, test_unsigned);
  test_ham_float = par->grid_obs.oc_r_max * 0.1;
  EXPECT_EQ(ref->d_start, test_ham_float);
  test_ham_float = par->grid_obs.oc_r_max * 0.7;
  EXPECT_EQ(ref->d_stop, test_ham_float);
  test_ham_float = (18 * cgs::kpc) / 600;
  EXPECT_NEAR(ref->delta_d, test_ham_float, 1.0e-10 * cgs::kpc);
  // radial quadrature point
  std::uniform_int_distribution<ham_int> smp4(1, 600);
  idx = smp4(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  idx = smp4(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  // 3rd shell
  pipe->assemble_shell_ref(ref.get(), par.get(), 2);
  test_unsigned = 300;
  EXPECT_EQ(ref->step, test_unsigned);
  test_ham_float = par->grid_obs.oc_r_max * 0.7;
  EXPECT_EQ(ref->d_start, test_ham_float);
  test_ham_float = par->grid_obs.oc_r_max * 1.0;
  EXPECT_EQ(ref->d_stop, test_ham_float);
  test_ham_float = (9 * cgs::kpc) / 300;
  EXPECT_NEAR(ref->delta_d, test_ham_float, 1.0e-10 * cgs::kpc);
  // radial quadrature point
  std::uniform_int_distribution<ham_int> smp5(1, 300);
  idx = smp5(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
  idx = smp5(rng);
  EXPECT_EQ(ref->dist[idx], ref->d_start + (idx + 0.5) * ref->delta_d);
}

// testing:
// los_versor
TEST(toolkit, los_versor) {
  auto pipe = std::make_unique<Integrator>();
  const ham_float theta[3] = {0., 90. * cgs::rad, 90. * cgs::rad};
  const ham_float phi[3] = {0., 180. * cgs::rad, 270. * cgs::rad};

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
  const ham_float theta[3] = {0., 90. * cgs::rad, 90. * cgs::rad};
  const ham_float phi[3] = {0., 180. * cgs::rad, 270. * cgs::rad};
  const Hamvec<3, ham_float> A{1., 0., 0.};

  EXPECT_LT(std::fabs(pipe->los_parproj(A, theta[0], phi[0]) - 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_parproj(A, theta[1], phi[1]) + 1.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_parproj(A, theta[2], phi[2]) - 0.), 1e-10);
}

// testing:
// los_perproj
TEST(toolkit, los_perproj) {
  auto pipe = std::make_unique<Integrator>();
  const ham_float theta[3] = {0., 90. * cgs::rad, 90. * cgs::rad};
  const ham_float phi[3] = {0., 180. * cgs::rad, 270. * cgs::rad};
  const Hamvec<3, ham_float> A{1., 0., 0.};

  EXPECT_LT(std::fabs(pipe->los_perproj(A, theta[0], phi[0]) - 1.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_perproj(A, theta[1], phi[1]) + 0.), 1e-10);
  EXPECT_LT(std::fabs(pipe->los_perproj(A, theta[2], phi[2]) - 1.), 1e-10);
}

// testing:
// sync_ipa
TEST(toolkit, sync_ipa) {
  auto pipe = std::make_unique<Integrator>();
  const ham_float theta[3] = {0., 90. * cgs::rad, 90. * cgs::rad};
  const ham_float phi[3] = {0., 180. * cgs::rad, 270. * cgs::rad};
  const Hamvec<3, ham_float> A{1., 0., 0.};

  EXPECT_LT(std::fabs(pipe->sync_ipa(A, theta[0], phi[0]) + 90. * cgs::rad),
            1e-10);
  EXPECT_LT(std::fabs(pipe->sync_ipa(A, theta[2], phi[2]) - 180. * cgs::rad),
            1e-10);
}
