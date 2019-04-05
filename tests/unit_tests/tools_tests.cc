// unit tests for namespace toolkit
// feel free to add more rational testing blocks

#include <gtest/gtest.h>

#include <cgs_units_file.h>
#include <cmath>
#include <fftw3.h>
#include <hvec.h>
#include <memory>
#include <namespace_toolkit.h>
#include <tinyxml2.h>
#include <vector>

TEST(toolkit, los_versor) {
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};

  EXPECT_LT(std::fabs(toolkit::los_versor(theta[0], phi[0])[0] - 0.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[0], phi[0])[1] - 0.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[0], phi[0])[2] - 1.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[1], phi[1])[0] + 1.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[1], phi[1])[1] - 0.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[1], phi[1])[2] - 0.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[2], phi[2])[0] - 0.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[2], phi[2])[1] + 1.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::los_versor(theta[2], phi[2])[2] - 0.), 1e-10);
}

TEST(toolkit, par2los) {
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};
  const hvec<3, double> A{1., 0., 0.};

  EXPECT_LT(std::fabs(toolkit::par2los(A, theta[0], phi[0]) - 0.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::par2los(A, theta[1], phi[1]) + 1.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::par2los(A, theta[2], phi[2]) - 0.), 1e-10);
}

TEST(toolkit, perp2los) {
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};
  const hvec<3, double> A{1., 0., 0.};

  EXPECT_LT(std::fabs(toolkit::perp2los(A, theta[0], phi[0]) - 1.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::perp2los(A, theta[1], phi[1]) + 0.), 1e-10);
  EXPECT_LT(std::fabs(toolkit::perp2los(A, theta[2], phi[2]) - 1.), 1e-10);
}

TEST(toolkit, intr_pol_ang) {
  const double theta[3] = {0., 90. * CGS_U_rad, 90. * CGS_U_rad};
  const double phi[3] = {0., 180. * CGS_U_rad, 270. * CGS_U_rad};
  const hvec<3, double> A{1., 0., 0.};

  EXPECT_LT(
      std::fabs(toolkit::intr_pol_ang(A, theta[0], phi[0]) + 90. * CGS_U_rad),
      1e-10);
  EXPECT_LT(
      std::fabs(toolkit::intr_pol_ang(A, theta[2], phi[2]) - 180. * CGS_U_rad),
      1e-10);
}

TEST(toolkit, cart_coord2cyl_coord) {
  const hvec<3, double> A{1., 0., 0.};
  const hvec<3, double> B{0., 0., 1.};
  const hvec<3, double> C{0., 1., 0.};
  hvec<3, double> tmp;

  toolkit::cart_coord2cyl_coord(A, tmp);

  EXPECT_LT(std::fabs(tmp[0] - 1.), 1e-10);
  EXPECT_LT(std::fabs(tmp[1] - 0.), 1e-10);
  EXPECT_LT(std::fabs(tmp[2] - 0.), 1e-10);

  toolkit::cart_coord2cyl_coord(B, tmp);

  EXPECT_LT(std::fabs(tmp[0] - 0.), 1e-10);
  EXPECT_LT(std::fabs(tmp[2] - 1.), 1e-10);

  toolkit::cart_coord2cyl_coord(C, tmp);

  EXPECT_LT(std::fabs(tmp[0] - 1.), 1e-10);
  EXPECT_LT(std::fabs(tmp[1] - 90. * CGS_U_rad), 1e-10);
  EXPECT_LT(std::fabs(tmp[2] - 0.), 1e-10);

  double tmp_r, tmp_phi, tmp_z;

  toolkit::cart_coord2cyl_coord(A, tmp_r, tmp_phi, tmp_z);

  EXPECT_LT(std::fabs(tmp_r - 1.), 1e-10);
  EXPECT_LT(std::fabs(tmp_phi - 0.), 1e-10);
  EXPECT_LT(std::fabs(tmp_z - 0.), 1e-10);

  toolkit::cart_coord2cyl_coord(B, tmp_r, tmp_phi, tmp_z);

  EXPECT_LT(std::fabs(tmp_r - 0.), 1e-10);
  EXPECT_LT(std::fabs(tmp_z - 1.), 1e-10);

  toolkit::cart_coord2cyl_coord(C, tmp_r, tmp_phi, tmp_z);

  EXPECT_LT(std::fabs(tmp_r - 1.), 1e-10);
  EXPECT_LT(std::fabs(tmp_phi - 90. * CGS_U_rad), 1e-10);
  EXPECT_LT(std::fabs(tmp_z - 0.), 1e-10);
}

TEST(toolkit, index3d) {
  std::size_t test_idx{53};
  EXPECT_EQ(toolkit::index3d(0, 5, 4, 2, 3, 1), test_idx);
}

TEST(toolkit, index4d) {
  std::size_t test_idx{461};
  EXPECT_EQ(toolkit::index4d(0, 4, 7, 3, 5, 1, 6, 2), test_idx);
}

TEST(toolkit, mean) {
  const double test_array[3] = {3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::mean(test_array, 3), 4.);

  const std::vector<double> test_vector{3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::mean(test_vector), 4.);
}

TEST(toolkit, variance) {
  const double test_array[3] = {3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::variance(test_array, 3), 2. / 3.);

  const std::vector<double> test_vector{3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::variance(test_vector), 2. / 3.);
}

TEST(toolkit, covariance) {
  const double test_array1[3] = {1, 2, 3};
  const double test_array2[3] = {3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::covariance(test_array1, test_array2, 3), 2. / 3.);

  const std::vector<double> test_vector1{1, 2, 3};
  const std::vector<double> test_vector2{3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::covariance(test_vector1, test_vector2), 2. / 3.);
}

TEST(toolkit, complex_stripping) {
  const auto test_complex = fftw_alloc_complex(2);
  test_complex[0][0] = 1.;
  test_complex[0][1] = 2.;
  test_complex[1][0] = 3.;
  test_complex[1][1] = 4.;
  double test_real[2];
  double test_imag[2];

  toolkit::complex2real(test_complex, test_real, 2);
  EXPECT_DOUBLE_EQ(test_real[0], test_complex[0][0]);
  EXPECT_DOUBLE_EQ(test_real[1], test_complex[1][0]);

  toolkit::complex2imag(test_complex, test_imag, 2);
  EXPECT_DOUBLE_EQ(test_imag[0], test_complex[0][1]);
  EXPECT_DOUBLE_EQ(test_imag[1], test_complex[1][1]);

  toolkit::complex2rni(test_complex, test_real, test_imag, 2);
  EXPECT_DOUBLE_EQ(test_real[0], test_complex[0][0]);
  EXPECT_DOUBLE_EQ(test_real[1], test_complex[1][0]);
  EXPECT_DOUBLE_EQ(test_imag[0], test_complex[0][1]);
  EXPECT_DOUBLE_EQ(test_imag[1], test_complex[1][1]);

  fftw_free(test_complex);
}

TEST(toolkit, xml_parser) {
  auto doc = toolkit::loadxml("reference/tools_tests.xml");

  tinyxml2::XMLElement *el = toolkit::tracexml(doc.get(), {"double"});

  EXPECT_EQ(toolkit::fetchdouble(el, "value"), double(3.14));

  EXPECT_EQ(toolkit::fetchdouble(el, "value", 2.0), double(3.14));

  EXPECT_EQ(toolkit::fetchdouble(el, "dft_value", 2.0), double(2.0));

  EXPECT_EQ(toolkit::fetchdouble(el, "dft_value"), double(0));

  el = toolkit::tracexml(doc.get(), {"empty_double"});

  EXPECT_EQ(toolkit::fetchdouble(el, "value", 5.0), double(5.0));

  el = toolkit::tracexml(doc.get(), {"integer"});

  EXPECT_EQ(toolkit::fetchunsigned(el, "value"), unsigned(23));

  EXPECT_EQ(toolkit::fetchunsigned(el, "value", 30), unsigned(23));

  EXPECT_EQ(toolkit::fetchunsigned(el, "dft_value", 30), unsigned(30));

  EXPECT_EQ(toolkit::fetchunsigned(el, "dft_value"), unsigned(0));

  el = toolkit::tracexml(doc.get(), {});

  EXPECT_EQ(toolkit::fetchstring(el, "value", "string"), std::string("string"));

  el = toolkit::tracexml(doc.get(), {});

  EXPECT_EQ(toolkit::fetchbool(el, "value", "false"), false);

  EXPECT_EQ(toolkit::fetchbool(el, "value", "true"), true);

  EXPECT_EQ(toolkit::fetchbool(el, "dft_value", "true"), false);
}
