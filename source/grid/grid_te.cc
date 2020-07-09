#include <array>
#include <cassert>
#include <fftw3.h>
#include <fstream>
#include <grid.h>
#include <hamtype.h>
#include <memory>
#include <omp.h>
#include <param.h>
#include <sstream>
#include <string>
#include <vector>

// turbulent free electron density field
Grid_te::Grid_te(const Param *par) {
  if (par->grid_te.build_permission or par->grid_te.read_permission)
    build_grid(par);
}

void Grid_te::build_grid(const Param *par) {
  // allocate spatial domain thermal electron field
  te = std::make_unique<ham_float[]>(par->grid_te.full_size);
  // activate FFT allocation upon build_permission
  if (par->grid_te.build_permission) {
    clean_fft = true;
    // allocate Fourier domain thermal electron field
    te_k = fftw_alloc_complex(par->grid_te.full_size);
    te_k[0][0] = 0;
    te_k[0][1] = 0; // 0th term should be zero
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    plan_te_bw =
        fftw_plan_dft_3d(par->grid_te.nx, par->grid_te.ny, par->grid_te.nz,
                         te_k, te_k, FFTW_BACKWARD, FFTW_ESTIMATE);
  }
}

void Grid_te::export_grid(const Param *par) {
  assert(!par->grid_te.filename.empty());
  std::ofstream output(par->grid_te.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  hamio_float tmp;
  for (decltype(par->grid_te.full_size) i = 0; i != par->grid_te.full_size;
       ++i) {
    assert(!output.eof());
    tmp = static_cast<hamio_float>(te[i]);
    output.write(reinterpret_cast<char *>(&tmp), sizeof(hamio_float));
  }
  output.close();
}

void Grid_te::import_grid(const Param *par) {
  assert(!par->grid_te.filename.empty());
  std::ifstream input(par->grid_te.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  hamio_float tmp;
  for (decltype(par->grid_te.full_size) i = 0; i != par->grid_te.full_size;
       ++i) {
    assert(!input.eof());
    input.read(reinterpret_cast<char *>(&tmp), sizeof(hamio_float));
    te[i] = static_cast<ham_float>(tmp);
  }
#ifndef NDEBUG
  auto eof = input.tellg();
  input.seekg(0, input.end);
#endif
  assert(eof == input.tellg());
  input.close();
}
