// random thermal electron density field grid

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <fftw3.h>
#include <omp.h>

#include <grid.h>
#include <param.h>

// turbulent free electron density field
Grid_ternd::Grid_ternd(const Param *par) {
  if (par->grid_ternd.build_permission or par->grid_ternd.read_permission) {
    build_grid(par);
    clean_switch = true;
  }
}

void Grid_ternd::build_grid(const Param *par) {
  // allocate spatial domain thermal electron field
  te = std::make_unique<double[]>(par->grid_ternd.full_size);
  // allocate Fourier domain thermal electron field
  te_k = fftw_alloc_complex(par->grid_ternd.full_size);
  te_k[0][0] = 0;
  te_k[0][1] = 0; // 0th term should be zero
#ifdef _OPENMP
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  plan_te_bw = fftw_plan_dft_3d(par->grid_ternd.nx, par->grid_ternd.ny,
                                par->grid_ternd.nz, te_k, te_k, FFTW_BACKWARD,
                                FFTW_ESTIMATE);
}

void Grid_ternd::export_grid(const Param *par) {
  assert(!par->grid_ternd.filename.empty());
  std::ofstream output(par->grid_ternd.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  double tmp;
  for (decltype(par->grid_ternd.full_size) i = 0;
       i != par->grid_ternd.full_size; ++i) {
    assert(!output.eof());
    tmp = te[i];
    output.write(reinterpret_cast<char *>(&tmp), sizeof(double));
  }
  output.close();
}

void Grid_ternd::import_grid(const Param *par) {
  assert(!par->grid_ternd.filename.empty());
  std::ifstream input(par->grid_ternd.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  double tmp;
  for (decltype(par->grid_ternd.full_size) i = 0;
       i != par->grid_ternd.full_size; ++i) {
    assert(!input.eof());
    input.read(reinterpret_cast<char *>(&tmp), sizeof(double));
    te[i] = tmp;
  }
#ifndef NDEBUG
  auto eof = input.tellg();
  input.seekg(0, input.end);
#endif
  assert(eof == input.tellg());
  input.close();
}
