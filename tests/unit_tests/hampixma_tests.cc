// unit tests for hampixma class

#include <gtest/gtest.h>
#include <iostream>
#include <memory>

#include <hampixma.h>
#include <param.h>

TEST(hampixma, basic) {
  // init
  auto testma = std::make_unique<hampixma>();
  // check pivot resolution
  EXPECT_EQ(testma->pivot_nside, 0);
  // read HEALPix file
  auto testpar = std::make_unique<Param>();
  testpar->grid_obs.do_mask = true;
  // the reference map is prepared in Nside 8
  testpar->grid_obs.mask_name = "reference/mask_stripe_60deg.fits";
  testma->import(testpar.get());
  // check pivot resolution
  EXPECT_EQ(testma->pivot_nside, 8);
  // check actual Nside
  EXPECT_EQ((testma->maps->at(testma->pivot_nside)).Nside(), 8);
  // get a high-res mask copy
  testma->duplicate(16);
  EXPECT_EQ((testma->maps->at(16)).Nside(), 16);
  // get a low-res mask copy
  testma->duplicate(4);
  EXPECT_EQ((testma->maps->at(4)).Nside(), 4);
  // check the pessimistic degrading
  EXPECT_TRUE( (testma->maps->at(4)).average() > (testma->maps->at(8)).average() );
  // check the upgrading
  EXPECT_TRUE( (testma->maps->at(8)).average() == (testma->maps->at(16)).average() );
  // a more precise check
  for (int i=0;i<768;++i) {
    const double pixval = testma->info(8,i);
    if (pixval == 1.0) {
      EXPECT_TRUE( pixval == testma->info(4,i%4) );
    }
  }
}
