// unit tests for grid interpolation
// for testing read_grid, we define field with linear relation to coordinates
// so that an arbitrary position get precisely interpoalted and tested
// feel free to add more rational testing blocks

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <random>

#include <breg.h>
#include <brnd.h>
#include <cre.h>
#include <fereg.h>
#include <fernd.h>
#include <grid.h>
#include <hvec.h>
#include <namespace_toolkit.h>

void fill_breg_grid(const Param *par, Grid_breg *grid);
void fill_brnd_grid(const Param *par, Grid_brnd *grid);
void fill_fereg_grid(const Param *par, Grid_fereg *grid);
void fill_fernd_grid(const Param *par, Grid_fereg *grid);
void fill_cre_grid(const Param *par, Grid_cre *grid);

// assign field vector by position
void fill_breg_grid(const Param *par, Grid_breg *grid) {
  hvec<3, double> gc_pos;
  double lx{par->grid_breg.x_max - par->grid_breg.x_min};
  double ly{par->grid_breg.y_max - par->grid_breg.y_min};
  double lz{par->grid_breg.z_max - par->grid_breg.z_min};
  for (decltype(par->grid_breg.nx) i = 0; i != par->grid_breg.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_breg.nx - 1) + par->grid_breg.x_min;
    for (decltype(par->grid_breg.ny) j = 0; j != par->grid_breg.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_breg.ny - 1) + par->grid_breg.y_min;
      for (decltype(par->grid_breg.nz) k = 0; k != par->grid_breg.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_breg.nz - 1) + par->grid_breg.z_min;
        std::size_t idx{toolkit::index3d(par->grid_breg.nx, par->grid_breg.ny,
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
  hvec<3, double> gc_pos;
  double lx{par->grid_brnd.x_max - par->grid_brnd.x_min};
  double ly{par->grid_brnd.y_max - par->grid_brnd.y_min};
  double lz{par->grid_brnd.z_max - par->grid_brnd.z_min};
  for (decltype(par->grid_brnd.nx) i = 0; i != par->grid_brnd.nx; ++i) {
    gc_pos[0] = i * lx / (par->grid_brnd.nx - 1) + par->grid_brnd.x_min;
    for (decltype(par->grid_brnd.ny) j = 0; j != par->grid_brnd.ny; ++j) {
      gc_pos[1] = j * ly / (par->grid_brnd.ny - 1) + par->grid_brnd.y_min;
      for (decltype(par->grid_brnd.nz) k = 0; k != par->grid_brnd.nz; ++k) {
        gc_pos[2] = k * lz / (par->grid_brnd.nz - 1) + par->grid_brnd.z_min;
        std::size_t idx{toolkit::index3d(par->grid_brnd.nx, par->grid_brnd.ny,
                                         par->grid_brnd.nz, i, j, k)};
        grid->bx[idx] = gc_pos[0];
        grid->by[idx] = gc_pos[1];
        grid->bz[idx] = gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_fernd_grid(const Param *par, Grid_fernd *grid) {
  hvec<3, double> gc_pos;
  double lx{par->grid_fernd.x_max - par->grid_fernd.x_min};
  double ly{par->grid_fernd.y_max - par->grid_fernd.y_min};
  double lz{par->grid_fernd.z_max - par->grid_fernd.z_min};
  for (decltype(par->grid_fernd.nx) i = 0; i != par->grid_fernd.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_fernd.nx - 1) + par->grid_fernd.x_min;
    for (decltype(par->grid_fernd.ny) j = 0; j != par->grid_fernd.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_fernd.ny - 1) + par->grid_fernd.y_min;
      for (decltype(par->grid_fernd.nz) k = 0; k != par->grid_fernd.nz; ++k) {
        std::size_t idx{toolkit::index3d(par->grid_fernd.nx, par->grid_fernd.ny,
                                         par->grid_fernd.nz, i, j, k)};
        gc_pos[2] = lz * k / (par->grid_fernd.nz - 1) + par->grid_fernd.z_min;
        grid->fe[idx] = gc_pos[0] + gc_pos[1] + gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_fereg_grid(const Param *par, Grid_fereg *grid) {
  hvec<3, double> gc_pos;
  double lx{par->grid_fereg.x_max - par->grid_fereg.x_min};
  double ly{par->grid_fereg.y_max - par->grid_fereg.y_min};
  double lz{par->grid_fereg.z_max - par->grid_fereg.z_min};
  for (decltype(par->grid_fereg.nx) i = 0; i != par->grid_fereg.nx; ++i) {
    gc_pos[0] = lx * i / (par->grid_fereg.nx - 1) + par->grid_fereg.x_min;
    for (decltype(par->grid_fereg.ny) j = 0; j != par->grid_fereg.ny; ++j) {
      gc_pos[1] = ly * j / (par->grid_fereg.ny - 1) + par->grid_fereg.y_min;
      for (decltype(par->grid_fereg.nz) k = 0; k != par->grid_fereg.nz; ++k) {
        std::size_t idx{toolkit::index3d(par->grid_fereg.nx, par->grid_fereg.ny,
                                         par->grid_fereg.nz, i, j, k)};
        gc_pos[2] = lz * k / (par->grid_fereg.nz - 1) + par->grid_fereg.z_min;
        grid->fe[idx] = gc_pos[0] + gc_pos[1] + gc_pos[2];
      }
    }
  }
}

// fill field scalar by sum of position coordinates
void fill_cre_grid(const Param *par, Grid_cre *grid) {
  hvec<3, double> gc_pos;
  double E;
  double lx{par->grid_cre.x_max - par->grid_cre.x_min};
  double ly{par->grid_cre.y_max - par->grid_cre.y_min};
  double lz{par->grid_cre.z_max - par->grid_cre.z_min};
  for (decltype(par->grid_cre.nE) i = 0; i != par->grid_cre.nE; ++i) {
    E = par->grid_cre.E_min * std::exp(i * par->grid_cre.E_fact);
    for (decltype(par->grid_cre.nx) j = 0; j != par->grid_cre.nx; ++j) {
      gc_pos[0] = lx * j / (par->grid_cre.nx - 1) + par->grid_cre.x_min;
      for (decltype(par->grid_cre.ny) k = 0; k != par->grid_cre.ny; ++k) {
        gc_pos[1] = ly * k / (par->grid_cre.ny - 1) + par->grid_cre.y_min;
        for (decltype(par->grid_cre.nz) m = 0; m != par->grid_cre.nz; ++m) {
          gc_pos[2] = lz * m / (par->grid_cre.nz - 1) + par->grid_cre.z_min;
          std::size_t idx{toolkit::index4d(par->grid_cre.nE, par->grid_cre.nx,
                                           par->grid_cre.ny, par->grid_cre.nz,
                                           i, j, k, m)};
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
  hvec<3, double> baseline(dis(gen), dis(gen), dis(gen));
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
  hvec<3, double> baseline(dis(gen), dis(gen), dis(gen));
  auto test_b = test_brnd->read_grid(baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_b[0], baseline[0], 1.0e-10);
  EXPECT_NEAR(test_b[1], baseline[1], 1.0e-10);
  EXPECT_NEAR(test_b[2], baseline[2], 1.0e-10);
}

// testing:
// FEreg::read_grid
TEST(grid, fereg_grid) {
  auto test_par = std::make_unique<Param>();
  test_par->grid_fereg.nx = 10;
  test_par->grid_fereg.ny = 8;
  test_par->grid_fereg.nz = 29;
  test_par->grid_fereg.x_max = 1;
  test_par->grid_fereg.x_min = 0;
  test_par->grid_fereg.y_max = 1;
  test_par->grid_fereg.y_min = 0;
  test_par->grid_fereg.z_max = 1;
  test_par->grid_fereg.z_min = 0;
  test_par->grid_fereg.full_size = 2320;
  test_par->grid_fereg.read_permission = true;
  auto test_grid = std::make_unique<Grid_fereg>(test_par.get());
  // fill test_grid
  fill_fereg_grid(test_par.get(), test_grid.get());
  //
  auto test_fereg = std::make_unique<FEreg>();
  // choose an arbitrary position and test read_grid function
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  hvec<3, double> baseline(dis(gen), dis(gen), dis(gen));
  auto test_fe =
      test_fereg->read_grid(baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_fe, baseline[0] + baseline[1] + baseline[2], 1.0e-10);
}

// testing:
// FErnd::read_grid
TEST(grid, fernd_grid) {
  auto test_par = std::make_unique<Param>();
  test_par->grid_fernd.nx = 10;
  test_par->grid_fernd.ny = 8;
  test_par->grid_fernd.nz = 29;
  test_par->grid_fernd.x_max = 1;
  test_par->grid_fernd.x_min = 0;
  test_par->grid_fernd.y_max = 1;
  test_par->grid_fernd.y_min = 0;
  test_par->grid_fernd.z_max = 1;
  test_par->grid_fernd.z_min = 0;
  test_par->grid_fernd.full_size = 2320;
  test_par->grid_fernd.read_permission = true;
  auto test_grid = std::make_unique<Grid_fernd>(test_par.get());
  // fill test_grid
  fill_fernd_grid(test_par.get(), test_grid.get());
  //
  auto test_fernd = std::make_unique<FErnd>();
  // choose an arbitrary position and test read_grid function
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  hvec<3, double> baseline(dis(gen), dis(gen), dis(gen));
  auto test_fe =
      test_fernd->read_grid(baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_fe, baseline[0] + baseline[1] + baseline[2], 1.0e-10);
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
  auto test_cre = std::make_unique<CRE>();
  // choose an arbitrary position and test read_grid function
  std::random_device
      rd; // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0.0, 1.0);
  hvec<3, double> baseline(dis(gen), dis(gen), dis(gen));
  std::uniform_int_distribution<> disi(0, test_par->grid_cre.nE - 1);
  auto idxE = disi(gen);
  const double E =
      test_par->grid_cre.E_min * std::exp(idxE * test_par->grid_cre.E_fact);
  auto test_c =
      test_cre->read_grid(idxE, baseline, test_par.get(), test_grid.get());
  EXPECT_NEAR(test_c, baseline[0] + baseline[1] + baseline[2] + E, 1.0e-10);
}
