// unit tests for Hamvec class

#include <cmath>
#include <gtest/gtest.h>
#include <hamtype.h>
#include <hamvec.h>

TEST(vector, basic) {
  // default ctor
  Hamvec<1, float> vec_dft1;
  EXPECT_EQ(vec_dft1[0], float(0));

  Hamvec<2, int> vec_dft2;
  EXPECT_EQ(vec_dft2[0], int(0));
  EXPECT_EQ(vec_dft2[1], int(0));

  Hamvec<3, double> vec_dft3;
  EXPECT_EQ(vec_dft3[0], double(0));
  EXPECT_EQ(vec_dft3[1], double(0));
  EXPECT_EQ(vec_dft3[2], double(0));

  // argument ctor
  Hamvec<1, float> vec1(0.0);
  EXPECT_EQ(vec1[0], 0.);

  // implicit ctor
  vec1 = 0.1;
  EXPECT_EQ(vec1[0], float(0.1));

  // mv assign, operator=
  Hamvec<1, float> vec1_mva(0.2);
  vec1 = std::move(vec1_mva);
  EXPECT_EQ(vec1[0], float(0.2));

  // cp assign, operator=
  Hamvec<1, float> vec1_cpa(0.2);
  vec1 = vec1_cpa;
  EXPECT_EQ(vec1[0], float(0.2));

  // cp ctor
  Hamvec<1, float> vec1_cpc(vec1);
  EXPECT_EQ(vec1_cpc[0], vec1[0]);

  // mv ctor
  Hamvec<1, float> vec1_mvc = std::move(vec1_cpc);
  EXPECT_EQ(vec1_mvc[0], vec1[0]);

  // list ctor
  Hamvec<2, double> vec2{0.3, 0.4};
  EXPECT_EQ(vec2[0], double(0.3));
  EXPECT_EQ(vec2[1], double(0.4));

  // operator +
  auto vecp = vec2 + Hamvec<2, double>{0.4, 0.8};
  EXPECT_EQ(vecp[0], vec2[0] + double(0.4));
  EXPECT_EQ(vecp[1], vec2[1] + double(0.8));

  // operator +=
  vec2 += vec2;
  EXPECT_EQ(vec2[0], double(0.3) + double(0.3));
  EXPECT_EQ(vec2[1], double(0.4) + double(0.4));

  // resign
  vec2 = {0.3, 0.4};
  EXPECT_EQ(vec2[0], double(0.3));
  EXPECT_EQ(vec2[1], double(0.4));

  // operator -
  auto vecm = vec2 - Hamvec<2, double>{0.3, 0.8};
  EXPECT_EQ(vecm[0], double(0.3) - double(0.3));
  EXPECT_EQ(vecm[1], double(0.4) - double(0.8));

  // operator -=
  vec2 = {0.3, 0.4};
  vec2 -= vec2;
  EXPECT_EQ(vec2[0], double(0.0));
  EXPECT_EQ(vec2[1], double(0.0));

  // operator *
  vec2 = {0.3, 0.4};
  vec2 = vec2 * 3;
  EXPECT_EQ(vec2[0], double(0.3) * double(3));
  EXPECT_EQ(vec2[1], double(0.4) * double(3));

  // operator *=
  vec2 = {0.3, 0.4};
  vec2 *= 0.3;
  EXPECT_EQ(vec2[0], double(0.3) * double(0.3));
  EXPECT_EQ(vec2[1], double(0.4) * double(0.3));

  // operator /
  vec2 = {0.27, 0.36};
  vec2 = vec2 / 0.3;
  EXPECT_EQ(vec2[0], double(0.27) / double(0.3));
  EXPECT_EQ(vec2[1], double(0.36) / double(0.3));

  // operator /=
  vec2 = {0.9, 1.2};
  vec2 /= 3;
  EXPECT_EQ(vec2[0], double(0.9) / double(3));
  EXPECT_EQ(vec2[1], double(1.2) / double(3));

  // operator !=
  auto test = Hamvec<2, double>(1., 1.);
  EXPECT_TRUE(vec2 != test);

  Hamvec<3, int> vec3{1, 2, 3};

  // function lengthsq
  EXPECT_EQ(vec3.lengthsq(), ham_float(14.0));

  // function length
  EXPECT_EQ(vec3.length(), std::sqrt(ham_float(14.0)));

  // function versor, non-zero
  auto norm = vec3.versor();
  EXPECT_EQ(norm[0], ham_float(1. / std::sqrt(14.)));
  EXPECT_EQ(norm[1], ham_float(2. / std::sqrt(14.)));
  EXPECT_EQ(norm[2], ham_float(3. / std::sqrt(14.)));

  // function versor, zero
  vec3 = {0, 0, 0};
  norm = vec3.versor();
  EXPECT_EQ(norm[0], ham_float(0));
  EXPECT_EQ(norm[1], ham_float(0));
  EXPECT_EQ(norm[2], ham_float(0));

  // function flip
  vec3 = {1, 2, 3};
  auto flipt = Hamvec<3, int>(-1, -2, -3);
  vec3.flip();
  EXPECT_TRUE(vec3 == flipt);

  Hamvec<3, double> vec3a{0.1, 0.2, 0.3};
  Hamvec<3, double> vec3b{0.4, 0.5, 0.6};

  // function dotprod
  auto rslt = vec3a.dotprod(vec3b);
  EXPECT_EQ(rslt, ham_float(0.32));

  // function crossprod
  Hamvec<3, double> prodt(0, 0, 0);
  prodt[0] = double(0.2) * double(0.6) - double(0.3) * double(0.5);
  prodt[1] = double(0.3) * double(0.4) - double(0.1) * double(0.6);
  prodt[2] = double(0.1) * double(0.5) - double(0.2) * double(0.4);
  auto vec3c = vec3a.crossprod(vec3b);
  EXPECT_TRUE(vec3c == prodt);
}
