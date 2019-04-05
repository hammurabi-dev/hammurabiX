/**
 * unit-test for namespace toolkit
 * feel free to add more rational testing blocks
 */
#include <gtest/gtest.h>

#include <cgs_units_file.h>
#include <cmath>
#include <fereg.h>
#include <param.h>

TEST(YMW16, precision) {
  // load template params
  std::unique_ptr<Param> par =
      std::make_unique<Param>("reference/ymw16_unitest.xml");

  std::unique_ptr<Grid_fereg> test_grid =
      std::make_unique<Grid_fereg>(par.get());
  std::unique_ptr<FEreg> test_field = std::make_unique<FEreg_ymw16>();
  test_field->write_grid(par.get(), test_grid.get());
  // grid->export_grid();

  // read reference file
  std::unique_ptr<Grid_fereg> ref_grid =
      std::make_unique<Grid_fereg>(par.get());
  std::unique_ptr<FEreg> ref_field = std::make_unique<FEreg_ymw16>();
  ref_grid->import_grid(par.get());

  for (decltype(par->grid_fereg.full_size) i = 0;
       i != par->grid_fereg.full_size; ++i) {
    EXPECT_LT(fabs(ref_grid->fe[i] - test_grid->fe[i]), 1e-5);
  }
}
