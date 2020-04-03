// unit tests for Grid class
// for testing the ''read_grid'' function
// we setup field with linear relation to coordinates
// so that an arbitrary position can be precisely interpoalted and tested

#include <bfield.h>
#include <cmath>
#include <crefield.h>
#include <grid.h>
#include <gtest/gtest.h>
#include <hamtype.h>
#include <hamvec.h>
#include <memory>
#include <random>
#include <tefield.h>
#include <toolkit.h>

void fill_grid_b(const Param *par, Grid_b *grid);
void fill_grid_te(const Param *par, Grid_te *grid);
void fill_grid_cre(const Param *par, Grid_cre *grid);

// assign field vector by position
void fill_grid_b(const Param *par, Grid_b *grid) {
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_b.x_max - par->grid_b.x_min};
  ham_float ly{par->grid_b.y_max - par->grid_b.y_min};
  ham_float lz{par->grid_b.z_max - par->grid_b.z_min};
  for (decltype(par->grid_b.nx) i = 0; i != par->grid_b.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_b.nx - 1) + par->grid_b.x_min;
    for (decltype(par->grid_b.ny) j = 0; j != par->grid_b.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_b.ny - 1) + par->grid_b.y_min;
      for (decltype(par->grid_b.nz) k = 0; k != par->grid_b.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_b.nz - 1) + par->grid_b.z_min;
        ham_uint idx{toolkit::index3d(par->grid_b.nx, par->grid_b.ny,
                                      par->grid_b.nz, i, j, k)};
        grid->bx[idx] = gc_pos[0];
        grid->by[idx] = gc_pos[1];
        grid->bz[idx] = gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_grid_te(const Param *par, Grid_te *grid) {
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_te.x_max - par->grid_te.x_min};
  ham_float ly{par->grid_te.y_max - par->grid_te.y_min};
  ham_float lz{par->grid_te.z_max - par->grid_te.z_min};
  for (decltype(par->grid_te.nx) i = 0; i != par->grid_te.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_te.nx - 1) + par->grid_te.x_min;
    for (decltype(par->grid_te.ny) j = 0; j != par->grid_te.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_te.ny - 1) + par->grid_te.y_min;
      for (decltype(par->grid_te.nz) k = 0; k != par->grid_te.nz; ++k) {
        ham_uint idx{toolkit::index3d(par->grid_te.nx, par->grid_te.ny,
                                      par->grid_te.nz, i, j, k)};
        gc_pos[2] = lz * k / (par->grid_te.nz - 1) + par->grid_te.z_min;
        grid->te[idx] = gc_pos[0] + gc_pos[1] + gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_grid_cre(const Param *par, Grid_cre *grid) {
  Hamvec<3, ham_float> gc_pos;
  ham_float E;
  ham_float lx{par->grid_cre.x_max - par->grid_cre.x_min};
  ham_float ly{par->grid_cre.y_max - par->grid_cre.y_min};
  ham_float lz{par->grid_cre.z_max - par->grid_cre.z_min};
  for (decltype(par->grid_cre.nx) i = 0; i != par->grid_cre.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_cre.nx - 1) + par->grid_cre.x_min;
    for (decltype(par->grid_cre.ny) j = 0; j != par->grid_cre.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_cre.ny - 1) + par->grid_cre.y_min;
      for (decltype(par->grid_cre.nz) k = 0; k != par->grid_cre.nz; ++k) {
        gc_pos[2] = lz * k / (par->grid_cre.nz - 1) + par->grid_cre.z_min;
        for (decltype(par->grid_cre.ne) m = 0; m != par->grid_cre.ne; ++m) {
          E = par->grid_cre.e_min * std::exp(m * par->grid_cre.e_fact);
          ham_uint idx{toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny,
                                        par->grid_cre.nz, par->grid_cre.ne, i,
                                        j, k, m)};
          grid->cre_flux[idx] = gc_pos[0] + gc_pos[1] + gc_pos[2] + E;
        }
      }
    }
  }
}

// testing:
// Grid_b::read_field
TEST(grid, grid_b) {
  // initialize parameter set manually
  auto test_par = std::make_unique<Param>();
  test_par->grid_b.nx = 10;
  test_par->grid_b.ny = 8;
  test_par->grid_b.nz = 29;
  test_par->grid_b.x_max = 1;
  test_par->grid_b.x_min = 0;
  test_par->grid_b.y_max = 1;
  test_par->grid_b.y_min = 0;
  test_par->grid_b.z_max = 1;
  test_par->grid_b.z_min = 0;
  test_par->grid_b.full_size = 2320;
  test_par->grid_b.read_permission = true;
  // initialize testing grid manually
  auto test_grid = std::make_unique<Grid_b>(test_par.get());
  // fill testing grid
  fill_grid_b(test_par.get(), test_grid.get());
  // initialize empty field
  auto test_bfield = std::make_unique<Bfield>(test_par.get());
  // choose an arbitrary position and test read_grid function
  // will be used to obtain a seed for the random number engine
  std::random_device rd;
  // standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> position(dis(gen), dis(gen), dis(gen));
  auto test_b =
      test_bfield->read_field(position, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_b[0], position[0], 1.0e-10);
  EXPECT_NEAR(test_b[1], position[1], 1.0e-10);
  EXPECT_NEAR(test_b[2], position[2], 1.0e-10);
}

// testing:
// Grid_te::read_field
TEST(grid, grid_te) {
  // initialize parameter set manually
  auto test_par = std::make_unique<Param>();
  test_par->grid_te.nx = 10;
  test_par->grid_te.ny = 8;
  test_par->grid_te.nz = 29;
  test_par->grid_te.x_max = 1;
  test_par->grid_te.x_min = 0;
  test_par->grid_te.y_max = 1;
  test_par->grid_te.y_min = 0;
  test_par->grid_te.z_max = 1;
  test_par->grid_te.z_min = 0;
  test_par->grid_te.full_size = 2320;
  test_par->grid_te.read_permission = true;
  // initialize testing grid manually
  auto test_grid = std::make_unique<Grid_te>(test_par.get());
  // fill testing grid
  fill_grid_te(test_par.get(), test_grid.get());
  // initialize empty field
  auto test_tefield = std::make_unique<TEfield>(test_par.get());
  // choose an arbitrary position and test read_grid function
  // will be used to obtain a seed for the random number engine
  std::random_device rd;
  // standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> position(dis(gen), dis(gen), dis(gen));
  auto test_te =
      test_tefield->read_field(position, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_te, position[0] + position[1] + position[2], 1.0e-10);
}

// testing:
// Grid_cre::read_field
TEST(grid, grid_cre) {
  // initialize parameter set manually
  auto test_par = std::make_unique<Param>();
  test_par->grid_cre.nx = 10;
  test_par->grid_cre.ny = 8;
  test_par->grid_cre.nz = 29;
  test_par->grid_cre.ne = 19;
  test_par->grid_cre.x_max = 1;
  test_par->grid_cre.x_min = 0;
  test_par->grid_cre.y_max = 1;
  test_par->grid_cre.y_min = 0;
  test_par->grid_cre.z_max = 1;
  test_par->grid_cre.z_min = 0;
  test_par->grid_cre.e_min = 0.01;
  test_par->grid_cre.e_max = 1.;
  test_par->grid_cre.e_fact =
      std::log(test_par->grid_cre.e_max / test_par->grid_cre.e_min) /
      (test_par->grid_cre.ne - 1);
  test_par->grid_cre.cre_size = 44080;
  test_par->grid_cre.read_permission = true;
  auto test_grid = std::make_unique<Grid_cre>(test_par.get());
  // fill test_grid
  fill_grid_cre(test_par.get(), test_grid.get());
  // initialize empty field
  auto test_crefield = std::make_unique<CREfield>(test_par.get());
  // choose an arbitrary position and test read_grid function
  // will be used to obtain a seed for the random number engine
  std::random_device rd;
  // standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> position(dis(gen), dis(gen), dis(gen));
  std::uniform_int_distribution<> disi(0, test_par->grid_cre.ne - 1);
  auto idxe = disi(gen);
  const ham_float e =
      test_par->grid_cre.e_min * std::exp(idxe * test_par->grid_cre.e_fact);
  auto test_c = test_crefield->read_field(position, idxe, test_par.get(),
                                          test_grid.get());
  EXPECT_NEAR(test_c, position[0] + position[1] + position[2] + e, 1.0e-10);
}
