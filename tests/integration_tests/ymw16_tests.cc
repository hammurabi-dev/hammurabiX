// integrated test for thermal electron I/O

#include <gtest/gtest.h>

#include <cgs_units_file.h>
#include <cmath>
#include <param.h>
#include <tefield.h>

TEST(YMW16, precision) {
  // load template params
  std::unique_ptr<Param> par =
      std::make_unique<Param>("reference/ymw16_unitest.xml");

  std::unique_ptr<Grid_tereg> test_grid =
      std::make_unique<Grid_tereg>(par.get());
  std::unique_ptr<TEreg> test_field = std::make_unique<TEreg_ymw16>();
  test_field->write_grid(par.get(), test_grid.get());
  // grid->export_grid();

  // read reference file
  std::unique_ptr<Grid_tereg> ref_grid =
      std::make_unique<Grid_tereg>(par.get());
  std::unique_ptr<TEreg> ref_field = std::make_unique<TEreg_ymw16>();
  ref_grid->import_grid(par.get());

  for (decltype(par->grid_tereg.full_size) i = 0;
       i != par->grid_tereg.full_size; ++i) {
    EXPECT_LT(fabs(ref_grid->te[i] - test_grid->te[i]), 1e-5);
  }
}
