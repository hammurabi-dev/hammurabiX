// unit tests for namespace toolkit
// feel free to add more rational testing blocks

#include <gtest/gtest.h>

#include <memory>
#include <random>
#include <tinyxml2.h>
#include <toolkit.h>
#include <vector>

// testing:
// toolkit::index3d
TEST(toolkit, index3d) {
  std::default_random_engine rng;
  std::uniform_int_distribution<int> smp(1, 10);

  const int nx = smp(rng);
  const int ny = smp(rng);
  const int nz = smp(rng);
  const int ix = smp(rng) % nx;
  const int iy = smp(rng) % ny;
  const int iz = smp(rng) % nz;
  EXPECT_EQ(toolkit::index3d(nx, ny, nz, ix, iy, iz),
            ix * ny * nz + iy * nz + iz);
}

// testing:
// toolkit::index4d
TEST(toolkit, index4d) {
  std::default_random_engine rng;
  std::uniform_int_distribution<int> smp(1, 10);

  const int ne = smp(rng);
  const int nx = smp(rng);
  const int ny = smp(rng);
  const int nz = smp(rng);
  const int ie = smp(rng) % ne;
  const int ix = smp(rng) % nx;
  const int iy = smp(rng) % ny;
  const int iz = smp(rng) % nz;
  EXPECT_EQ(toolkit::index4d(ne, nx, ny, nz, ie, ix, iy, iz),
            ie * nx * ny * nz + ix * ny * nz + iy * nz + iz);
}

// testing:
// toolkit::mean
TEST(toolkit, mean) {
  std::default_random_engine rng;
  std::uniform_real_distribution<double> smp(0.0, 1.0);

  const double arr[3] = {smp(rng), smp(rng), smp(rng)};
  EXPECT_DOUBLE_EQ(toolkit::mean(arr, 3), (arr[0] + arr[1] + arr[2]) / 3.);

  const std::vector<double> vec{smp(rng), smp(rng), smp(rng)};
  EXPECT_DOUBLE_EQ(toolkit::mean(vec), (vec[0] + vec[1] + vec[2]) / 3.);
}

// testing:
// toolkit::variance
TEST(toolkit, variance) {
  const double test_array[3] = {3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::variance(test_array, 3), 2. / 3.);

  const std::vector<double> test_vector{3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::variance(test_vector), 2. / 3.);
}

// testing:
// toolkit::covariance
TEST(toolkit, covariance) {
  const double test_array1[3] = {1, 2, 3};
  const double test_array2[3] = {3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::covariance(test_array1, test_array2, 3), 2. / 3.);

  const std::vector<double> test_vector1{1, 2, 3};
  const std::vector<double> test_vector2{3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::covariance(test_vector1, test_vector2), 2. / 3.);
}

// testing:
// toolkit::xml_parser
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
