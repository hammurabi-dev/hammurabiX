// random magnetic vector field grid

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <omp.h>
#include <sstream>
#include <string>
#include <vector>

#include <fftw3.h>

#include <grid.h>
#include <hamtype.h>
#include <param.h>

Grid_brnd::Grid_brnd(const Param *par) {
  if (par->grid_brnd.build_permission or par->grid_brnd.read_permission) {
    build_grid(par);
    clean_switch = true;
  }
}

void Grid_brnd::build_grid(const Param *par) {
  // allocate spatial domian magnetic field
  bx = std::make_unique<ham_float[]>(par->grid_brnd.full_size);
  by = std::make_unique<ham_float[]>(par->grid_brnd.full_size);
  bz = std::make_unique<ham_float[]>(par->grid_brnd.full_size);
  // Fourier domain complex field
  c0 = fftw_alloc_complex(par->grid_brnd.full_size);
  c0[0][0] = 0;
  c0[0][1] = 0; // 0th term should be zero
  c1 = fftw_alloc_complex(par->grid_brnd.full_size);
  c1[0][0] = 0;
  c1[0][1] = 0; // 0th term should be zero
#ifdef _OPENMP
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  // backword in-place plans
  plan_c0_bw =
      fftw_plan_dft_3d(par->grid_brnd.nx, par->grid_brnd.ny, par->grid_brnd.nz,
                       c0, c0, FFTW_BACKWARD, FFTW_ESTIMATE);
  plan_c1_bw =
      fftw_plan_dft_3d(par->grid_brnd.nx, par->grid_brnd.ny, par->grid_brnd.nz,
                       c1, c1, FFTW_BACKWARD, FFTW_ESTIMATE);
  // forward in-place plans
  plan_c0_fw =
      fftw_plan_dft_3d(par->grid_brnd.nx, par->grid_brnd.ny, par->grid_brnd.nz,
                       c0, c0, FFTW_FORWARD, FFTW_ESTIMATE);
  plan_c1_fw =
      fftw_plan_dft_3d(par->grid_brnd.nx, par->grid_brnd.ny, par->grid_brnd.nz,
                       c1, c1, FFTW_FORWARD, FFTW_ESTIMATE);
}

void Grid_brnd::export_grid(const Param *par) {
  assert(!par->grid_brnd.filename.empty());
  std::ofstream output(par->grid_brnd.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  ham_float tmp;
  for (decltype(par->grid_brnd.full_size) i = 0; i != par->grid_brnd.full_size;
       ++i) {
    assert(!output.eof());
    tmp = bx[i];
    output.write(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
    tmp = by[i];
    output.write(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
    tmp = bz[i];
    output.write(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
  }
  output.close();
}

void Grid_brnd::import_grid(const Param *par) {
  assert(!par->grid_brnd.filename.empty());
  std::ifstream input(par->grid_brnd.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  ham_float tmp;
  for (decltype(par->grid_brnd.full_size) i = 0; i != par->grid_brnd.full_size;
       ++i) {
    assert(!input.eof());
    input.read(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
    bx[i] = tmp;
    input.read(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
    by[i] = tmp;
    input.read(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
    bz[i] = tmp;
  }
#ifndef NDEBUG
  auto eof = input.tellg();
  input.seekg(0, input.end);
#endif
  assert(eof == input.tellg());
  input.close();
}
