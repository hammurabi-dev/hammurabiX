// unit tests for Hamsk class

#include <gtest/gtest.h>
#include <hamdis.h>
#include <hamio.h>
#include <hamsk.h>
#include <hamtype.h>
#include <string>

TEST(Hampisk, basic) {
  // init
  Hampisk<ham_float> mask_dft;
  EXPECT_EQ(mask_dft.pivot(), 0);

  Hampix<ham_float> map_n(4);
  Hampisk<ham_float> mask_n(map_n);
  EXPECT_EQ(mask_n.pivot(), 4);
  for (ham_uint i = 0; i < 192; ++i) {
    EXPECT_EQ(mask_n.data(4, i), ham_float(0));
  }

  // duplicate
  mask_n.duplicate(2);
  for (ham_uint i = 0; i < 48; ++i) {
    EXPECT_EQ(mask_n.data(2, i), ham_float(0));
  }
  mask_n.duplicate(8);
  for (ham_uint i = 0; i < 768; ++i) {
    EXPECT_EQ(mask_n.data(8, i), ham_float(0));
  }

  // import external map
  Hamio<ham_float> io_dft("reference/mask_map_nside8.bin");
  map_n.reset(8);
  io_dft.load(map_n);
  Hampisk<ham_float> mask_ext(map_n);
  mask_ext.duplicate(2);
  mask_ext.duplicate(64);

  // test mask downgrade
  map_n.reset(2);
  io_dft.filename("reference/mask_map_nside2.bin");
  io_dft.load(map_n);
  for (ham_uint i = 0; i < 48; ++i) {
    EXPECT_EQ(mask_ext.data(2, i), map_n.data(i));
  }

  // test mask upgrade
  map_n.reset(64);
  io_dft.filename("reference/mask_map_nside64.bin");
  io_dft.load(map_n);
  for (ham_uint i = 0; i < 49152; ++i) {
    EXPECT_EQ(mask_ext.data(64, i), map_n.data(i));
  }
}
