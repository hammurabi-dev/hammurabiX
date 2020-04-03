// unit tests for Hamp class

#include <cmath>
#include <gtest/gtest.h>
#include <hamp.h>
#include <hamtype.h>
#include <hamunits.h>
#include <hamvec.h>

TEST(Hamp, basic) {
  // init
  Hamp ptr1;
  EXPECT_EQ(ptr1.theta(), ham_float(0));
  EXPECT_EQ(ptr1.phi(), ham_float(0));

  // return func
  EXPECT_EQ(ptr1.theta(), ptr1.theta());
  EXPECT_EQ(ptr1.phi(), ptr1.phi());

  // explicit init
  Hamp ptr2(0.1, 0.2);
  EXPECT_EQ(ptr2.theta(), ham_float(0.1));
  EXPECT_EQ(ptr2.phi(), ham_float(0.2));

  // list init
  Hamp ptr3 = {0.3, 0.4};
  EXPECT_EQ(ptr3.theta(), ham_float(0.3));
  EXPECT_EQ(ptr3.phi(), ham_float(0.4));

  // copy
  Hamp ptr4(ptr2);
  EXPECT_EQ(ptr4.theta(), ham_float(0.1));
  EXPECT_EQ(ptr4.phi(), ham_float(0.2));
  ptr4 = ptr3;
  EXPECT_EQ(ptr4.theta(), ham_float(0.3));
  EXPECT_EQ(ptr4.phi(), ham_float(0.4));

  // move
  Hamp ptr5(std::move(ptr2));
  EXPECT_EQ(ptr5.theta(), ham_float(0.1));
  EXPECT_EQ(ptr5.phi(), ham_float(0.2));
  ptr5 = std::move(ptr4);
  EXPECT_EQ(ptr5.theta(), ham_float(0.3));
  EXPECT_EQ(ptr5.phi(), ham_float(0.4));

  // norm
  Hamp ptr6(1.2 * cgs::pi, 2.1 * cgs::pi);
  EXPECT_NEAR(ptr6.theta(), ham_float(0.2 * cgs::pi), 1.0e-10);
  EXPECT_NEAR(ptr6.phi(), ham_float(0.1 * cgs::pi), 1.0e-10);
  ptr6.theta(-1.2 * cgs::pi);
  EXPECT_NEAR(ptr6.theta(), ham_float(0.8 * cgs::pi), 1.0e-10);
  ptr6.theta(-0.2 * cgs::pi);
  EXPECT_NEAR(ptr6.theta(), ham_float(0.8 * cgs::pi), 1.0e-10);
  ptr6.theta(-0.9 * cgs::pi);
  EXPECT_NEAR(ptr6.theta(), ham_float(0.1 * cgs::pi), 1.0e-10);
  ptr6.phi(-2.3 * cgs::pi);
  EXPECT_NEAR(ptr6.phi(), ham_float(1.7 * cgs::pi), 1.0e-10);
  ptr6.phi(-0.8 * cgs::pi);
  EXPECT_NEAR(ptr6.phi(), ham_float(1.2 * cgs::pi), 1.0e-10);
  ptr6.phi(-1.4 * cgs::pi);
  EXPECT_NEAR(ptr6.phi(), ham_float(0.6 * cgs::pi), 1.0e-10);
}

TEST(Hamp, vector) {
  Hamp ptr1(0.3 * cgs::pi, 1.7 * cgs::pi);
  const auto v1 = ptr1.versor();
  EXPECT_NEAR(v1[0], std::sin(0.3 * cgs::pi) * std::cos(1.7 * cgs::pi),
              1.0e-10);
  EXPECT_NEAR(v1[1], std::sin(0.3 * cgs::pi) * std::sin(1.7 * cgs::pi),
              1.0e-10);
  EXPECT_NEAR(v1[2], std::cos(0.3 * cgs::pi), 1.0e-10);

  Hamp ptr2(v1);
  EXPECT_NEAR(ptr2.theta(), ham_float(0.3 * cgs::pi), 1.0e-10);
  EXPECT_NEAR(ptr2.phi(), ham_float(1.7 * cgs::pi), 1.0e-10);
}
