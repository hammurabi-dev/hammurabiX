// unit tests for Grid class
// for testing the ''read_grid'' function
// we setup field with linear relation to coordinates
// so that an arbitrary position can be precisely interpoalted and tested

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <random>

#include <bfield.h>
#include <crefield.h>
#include <grid.h>
#include <hamtype.h>
#include <hamvec.h>
#include <tefield.h>
#include <toolkit.h>

void fill_breg_grid(const Param *par, Grid_breg *grid);
void fill_brnd_grid(const Param *par, Grid_brnd *grid);
void fill_tereg_grid(const Param *par, Grid_tereg *grid);
void fill_ternd_grid(const Param *par, Grid_tereg *grid);
void fill_cre_grid(const Param *par, Grid_cre *grid);

// assign field vector by position
void fill_breg_grid(const Param *par, Grid_breg *grid) {
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_breg.x_max - par->grid_breg.x_min};
  ham_float ly{par->grid_breg.y_max - par->grid_breg.y_min};
  ham_float lz{par->grid_breg.z_max - par->grid_breg.z_min};
  for (decltype(par->grid_breg.nx) i = 0; i != par->grid_breg.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_breg.nx - 1) + par->grid_breg.x_min;
    for (decltype(par->grid_breg.ny) j = 0; j != par->grid_breg.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_breg.ny - 1) + par->grid_breg.y_min;
      for (decltype(par->grid_breg.nz) k = 0; k != par->grid_breg.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_breg.nz - 1) + par->grid_breg.z_min;
        ham_uint idx{toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
                                      par->grid_breg.nz, i, j, k)};
        grid->bx[idx] = gc_pos[0];
        grid->by[idx] = gc_pos[1];
        grid->bz[idx] = gc_pos[2];
      }
    }
  }
}

// assign field vector by position
void fill_brnd_grid(const Param *par, Grid_brnd *grid) {
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_brnd.x_max - par->grid_brnd.x_min};
  ham_float ly{par->grid_brnd.y_max - par->grid_brnd.y_min};
  ham_float lz{par->grid_brnd.z_max - par->grid_brnd.z_min};
  for (decltype(par->grid_brnd.nx) i = 0; i != par->grid_brnd.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_brnd.nx - 1) + par->grid_brnd.x_min;
    for (decltype(par->grid_brnd.ny) j = 0; j != par->grid_brnd.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_brnd.ny - 1) + par->grid_brnd.y_min;
      for (decltype(par->grid_brnd.nz) k = 0; k != par->grid_brnd.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_brnd.nz - 1) + par->grid_brnd.z_min;
        ham_uint idx{toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                                      par->grid_brnd.nz, i, j, k)};
        grid->bx[idx] = gc_pos[0];
        grid->by[idx] = gc_pos[1];
        grid->bz[idx] = gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_ternd_grid(const Param *par, Grid_ternd *grid) {
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_ternd.x_max - par->grid_ternd.x_min};
  ham_float ly{par->grid_ternd.y_max - par->grid_ternd.y_min};
  ham_float lz{par->grid_ternd.z_max - par->grid_ternd.z_min};
  for (decltype(par->grid_ternd.nx) i = 0; i != par->grid_ternd.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_ternd.nx - 1) + par->grid_ternd.x_min;
    for (decltype(par->grid_ternd.ny) j = 0; j != par->grid_ternd.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_ternd.ny - 1) + par->grid_ternd.y_min;
      for (decltype(par->grid_ternd.nz) k = 0; k != par->grid_ternd.nz; ++k) {
        ham_uint idx{toolkit::index3d(par->grid_ternd.nx, par->grid_ternd.ny,
                                      par->grid_ternd.nz, i, j, k)};
        gc_pos[2] = lz * k / (par->grid_ternd.nz - 1) + par->grid_ternd.z_min;
        grid->te[idx] = gc_pos[0] + gc_pos[1] + gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_tereg_grid(const Param *par, Grid_tereg *grid) {
  Hamvec<3, ham_float> gc_pos;
  ham_float lx{par->grid_tereg.x_max - par->grid_tereg.x_min};
  ham_float ly{par->grid_tereg.y_max - par->grid_tereg.y_min};
  ham_float lz{par->grid_tereg.z_max - par->grid_tereg.z_min};
  for (decltype(par->grid_tereg.nx) i = 0; i != par->grid_tereg.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_tereg.nx - 1) + par->grid_tereg.x_min;
    for (decltype(par->grid_tereg.ny) j = 0; j != par->grid_tereg.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_tereg.ny - 1) + par->grid_tereg.y_min;
      for (decltype(par->grid_tereg.nz) k = 0; k != par->grid_tereg.nz; ++k) {
        ham_uint idx{toolkit::index3d(par->grid_tereg.nx, par->grid_tereg.ny,
                                      par->grid_tereg.nz, i, j, k)};
        gc_pos[2] = lz * k / (par->grid_tereg.nz - 1) + par->grid_tereg.z_min;
        grid->te[idx] = gc_pos[0] + gc_pos[1] + gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_cre_grid(const Param *par, Grid_cre *grid) {
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
        for (decltype(par->grid_cre.nE) m = 0; m != par->grid_cre.nE; ++m) {
          E = par->grid_cre.E_min * std::exp(m * par->grid_cre.E_fact);
          ham_uint idx{toolkit::index4d(par->grid_cre.nx, par->grid_cre.ny,
                                        par->grid_cre.nz, par->grid_cre.nE, i,
                                        j, k, m)};
          grid->cre_flux[idx] = gc_pos[0] + gc_pos[1] + gc_pos[2] + E;
        }
      }
    }
  }
}

// testing:
// Breg::read_grid
TEST(grid, breg_grid) {
  auto test_par = std::make_unique<Param>();
  test_par->grid_breg.nx = 10;
  test_par->grid_breg.ny = 8;
  test_par->grid_breg.nz = 29;
  test_par->grid_breg.x_max = 1;
  test_par->grid_breg.x_min = 0;
  test_par->grid_breg.y_max = 1;
  test_par->grid_breg.y_min = 0;
  test_par->grid_breg.z_max = 1;
  test_par->grid_breg.z_min = 0;
  test_par->grid_breg.full_size = 2320;
  test_par->grid_breg.read_permission = true;
  auto test_grid = std::make_unique<Grid_breg>(test_par.get());
  // fill test_grid
  fill_breg_grid(test_par.get(), test_grid.get());
  //
  auto test_breg = std::make_unique<Breg>();
  // choose an arbitrary position and test read_grid function
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> baseline(dis(gen), dis(gen), dis(gen));
  auto test_b = test_breg->read_grid(baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_b[0], baseline[0], 1.0e-10);
  EXPECT_NEAR(test_b[1], baseline[1], 1.0e-10);
  EXPECT_NEAR(test_b[2], baseline[2], 1.0e-10);
}

// testing:
// Brnd::read_grid
TEST(grid, brnd_grid) {
  auto test_par = std::make_unique<Param>();
  test_par->grid_brnd.nx = 10;
  test_par->grid_brnd.ny = 8;
  test_par->grid_brnd.nz = 29;
  test_par->grid_brnd.x_max = 1;
  test_par->grid_brnd.x_min = 0;
  test_par->grid_brnd.y_max = 1;
  test_par->grid_brnd.y_min = 0;
  test_par->grid_brnd.z_max = 1;
  test_par->grid_brnd.z_min = 0;
  test_par->grid_brnd.full_size = 2320;
  test_par->grid_brnd.read_permission = true;
  auto test_grid = std::make_unique<Grid_brnd>(test_par.get());
  // fill test_grid
  fill_brnd_grid(test_par.get(), test_grid.get());
  //
  auto test_brnd = std::make_unique<Brnd>();
  // choose an arbitrary position and test read_grid function
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> baseline(dis(gen), dis(gen), dis(gen));
  auto test_b = test_brnd->read_grid(baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_b[0], baseline[0], 1.0e-10);
  EXPECT_NEAR(test_b[1], baseline[1], 1.0e-10);
  EXPECT_NEAR(test_b[2], baseline[2], 1.0e-10);
}

// testing:
// TEreg::read_grid
TEST(grid, tereg_grid) {
  auto test_par = std::make_unique<Param>();
  test_par->grid_tereg.nx = 10;
  test_par->grid_tereg.ny = 8;
  test_par->grid_tereg.nz = 29;
  test_par->grid_tereg.x_max = 1;
  test_par->grid_tereg.x_min = 0;
  test_par->grid_tereg.y_max = 1;
  test_par->grid_tereg.y_min = 0;
  test_par->grid_tereg.z_max = 1;
  test_par->grid_tereg.z_min = 0;
  test_par->grid_tereg.full_size = 2320;
  test_par->grid_tereg.read_permission = true;
  auto test_grid = std::make_unique<Grid_tereg>(test_par.get());
  // fill test_grid
  fill_tereg_grid(test_par.get(), test_grid.get());
  //
  auto test_tereg = std::make_unique<TEreg>();
  // choose an arbitrary position and test read_grid function
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> baseline(dis(gen), dis(gen), dis(gen));
  auto test_te =
      test_tereg->read_grid(baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_te, baseline[0] + baseline[1] + baseline[2], 1.0e-10);
}

// testing:
// TErnd::read_grid
TEST(grid, ternd_grid) {
  auto test_par = std::make_unique<Param>();
  test_par->grid_ternd.nx = 10;
  test_par->grid_ternd.ny = 8;
  test_par->grid_ternd.nz = 29;
  test_par->grid_ternd.x_max = 1;
  test_par->grid_ternd.x_min = 0;
  test_par->grid_ternd.y_max = 1;
  test_par->grid_ternd.y_min = 0;
  test_par->grid_ternd.z_max = 1;
  test_par->grid_ternd.z_min = 0;
  test_par->grid_ternd.full_size = 2320;
  test_par->grid_ternd.read_permission = true;
  auto test_grid = std::make_unique<Grid_ternd>(test_par.get());
  // fill test_grid
  fill_ternd_grid(test_par.get(), test_grid.get());
  //
  auto test_ternd = std::make_unique<TErnd>();
  // choose an arbitrary position and test read_grid function
  // will be used to obtain a seed for the random number engine
  std::random_device rd;
  // standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> baseline(dis(gen), dis(gen), dis(gen));
  auto test_te =
      test_ternd->read_grid(baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_te, baseline[0] + baseline[1] + baseline[2], 1.0e-10);
}

// testing:
// CRE::read_grid
TEST(grid, cre_grid) {
  auto test_par = std::make_unique<Param>();
  test_par->grid_cre.nx = 10;
  test_par->grid_cre.ny = 8;
  test_par->grid_cre.nz = 29;
  test_par->grid_cre.nE = 19;
  test_par->grid_cre.x_max = 1;
  test_par->grid_cre.x_min = 0;
  test_par->grid_cre.y_max = 1;
  test_par->grid_cre.y_min = 0;
  test_par->grid_cre.z_max = 1;
  test_par->grid_cre.z_min = 0;
  test_par->grid_cre.E_min = 0.01;
  test_par->grid_cre.E_max = 1;
  test_par->grid_cre.E_fact =
      std::log(test_par->grid_cre.E_max / test_par->grid_cre.E_min) /
      (test_par->grid_cre.nE - 1);
  test_par->grid_cre.cre_size = 44080;
  test_par->grid_cre.read_permission = true;
  auto test_grid = std::make_unique<Grid_cre>(test_par.get());
  // fill test_grid
  fill_cre_grid(test_par.get(), test_grid.get());
  //
  auto test_cre = std::make_unique<CRE_num>();
  // choose an arbitrary position and test read_grid function
  // will be used to obtain a seed for the random number engine
  std::random_device rd;
  // standard mersenne_twister_engine seeded with rd()
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  Hamvec<3, ham_float> baseline(dis(gen), dis(gen), dis(gen));
  std::uniform_int_distribution<> disi(0, test_par->grid_cre.nE - 1);
  auto idxE = disi(gen);
  const ham_float E =
      test_par->grid_cre.E_min * std::exp(idxE * test_par->grid_cre.E_fact);
  auto test_c =
      test_cre->read_grid_num(baseline, idxE, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_c, baseline[0] + baseline[1] + baseline[2] + E, 1.0e-10);
}
