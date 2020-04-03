// unit tests for Hamio class

#include <gtest/gtest.h>
#include <hamdis.h>
#include <hamio.h>
#include <hamtype.h>
#include <string>

TEST(Hamio, basic) {
  ham_uint base_nside = 4;
  ham_uint base_npix = 192;
  // init
  Hamio<ham_float> io_dft("reference/random_map_nside4.bin");
  EXPECT_EQ(io_dft.filename(), std::string("reference/random_map_nside4.bin"));

  // manual read
  std::vector<ham_float> base_input;
  std::fstream base_file("reference/random_map_nside4.bin",
                         std::ios::in | std::ios::binary);
  if (base_file.is_open()) {
    ham_float tmp;
    for (ham_uint i = 0; i != base_npix; ++i) {
      base_file.read(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
      base_input.push_back(tmp);
    }
    base_file.close();
  }
  Hampix<ham_float> base_map(base_nside, base_input);
  // hamio load
  Hampix<ham_float> test_map(base_nside);
  io_dft.load(test_map);
  for (ham_uint i = 0; i < base_npix; ++i) {
    EXPECT_EQ(test_map.data(i), base_map.data(i));
    EXPECT_EQ(test_map.index(i), base_map.index(i));
    EXPECT_EQ(test_map.pointing(i).theta(), base_map.pointing(i).theta());
    EXPECT_EQ(test_map.pointing(i).phi(), base_map.pointing(i).phi());
  }
  // hamio dump
  io_dft.filename("reference/test_map_nside4.bin");
  io_dft.dump(test_map);
  // reload
  io_dft.load(test_map);
  for (ham_uint i = 0; i < base_npix; ++i) {
    EXPECT_EQ(test_map.data(i), base_map.data(i));
    EXPECT_EQ(test_map.index(i), base_map.index(i));
    EXPECT_EQ(test_map.pointing(i).theta(), base_map.pointing(i).theta());
    EXPECT_EQ(test_map.pointing(i).phi(), base_map.pointing(i).phi());
  }
  // dynamic bind
  Hamdis<ham_float> *bind_map = &test_map;
  io_dft.dump(*bind_map);
  io_dft.load(*bind_map);
  for (ham_uint i = 0; i < base_npix; ++i) {
    EXPECT_EQ(bind_map->data(i), base_map.data(i));
    EXPECT_EQ(bind_map->index(i), base_map.index(i));
    EXPECT_EQ(bind_map->pointing(i).theta(), base_map.pointing(i).theta());
    EXPECT_EQ(bind_map->pointing(i).phi(), base_map.pointing(i).phi());
  }
}
