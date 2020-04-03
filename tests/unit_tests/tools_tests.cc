// unit tests for namespace toolkit

#include <gtest/gtest.h>
#include <hamtype.h>
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
  std::uniform_real_distribution<ham_float> smp(0.0, 1.0);

  const ham_float arr[3] = {smp(rng), smp(rng), smp(rng)};
  EXPECT_DOUBLE_EQ(toolkit::mean(arr, 3), (arr[0] + arr[1] + arr[2]) / 3.);

  const std::vector<ham_float> vec{smp(rng), smp(rng), smp(rng)};
  EXPECT_DOUBLE_EQ(toolkit::mean(vec), (vec[0] + vec[1] + vec[2]) / 3.);
}

// testing:
// toolkit::variance
TEST(toolkit, variance) {
  const ham_float test_array[3] = {3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::variance(test_array, 3), 2. / 3.);

  const std::vector<ham_float> test_vector{3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::variance(test_vector), 2. / 3.);
}

// testing:
// toolkit::covariance
TEST(toolkit, covariance) {
  const ham_float test_array1[3] = {1, 2, 3};
  const ham_float test_array2[3] = {3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::covariance(test_array1, test_array2, 3), 2. / 3.);

  const std::vector<ham_float> test_vector1{1, 2, 3};
  const std::vector<ham_float> test_vector2{3, 4, 5};
  EXPECT_DOUBLE_EQ(toolkit::covariance(test_vector1, test_vector2), 2. / 3.);
}

// testing:
// toolkit::xml_parser
TEST(toolkit, xml_parser) {
  auto doc = toolkit::loadxml("reference/tools_tests.xml");

  tinyxml2::XMLElement *el = toolkit::tracexml(doc.get(), {"float"});

  EXPECT_EQ(toolkit::fetchfloat(el, "value"), ham_float(3.14));

  el = toolkit::tracexml(doc.get(), {"integer"});

  EXPECT_EQ(toolkit::fetchuint(el, "value"), ham_uint(23));

  el = toolkit::tracexml(doc.get(), {});

  EXPECT_EQ(toolkit::fetchstring(el, "string", "value"), std::string("string"));

  el = toolkit::tracexml(doc.get(), {});

  EXPECT_EQ(toolkit::fetchbool(el, "false", "value"), false);

  EXPECT_EQ(toolkit::fetchbool(el, "true", "value"), true);
}
