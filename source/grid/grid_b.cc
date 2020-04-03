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

Grid_b::Grid_b(const Param *par) {
  if (par->grid_b.build_permission or par->grid_b.read_permission or
      par->grid_b.write_permission)
    build_grid(par);
}

void Grid_b::build_grid(const Param *par) {
  // allocate spatial domian magnetic field
  bx = std::make_unique<ham_float[]>(par->grid_b.full_size);
  by = std::make_unique<ham_float[]>(par->grid_b.full_size);
  bz = std::make_unique<ham_float[]>(par->grid_b.full_size);
  // activate FFT allocation upon build_permission
  if (par->grid_b.build_permission) {
    clean_fft = true;
    // Fourier domain complex field
    c0 = fftw_alloc_complex(par->grid_b.full_size);
    c0[0][0] = 0;
    c0[0][1] = 0; // 0th term should be zero
    c1 = fftw_alloc_complex(par->grid_b.full_size);
    c1[0][0] = 0;
    c1[0][1] = 0; // 0th term should be zero
#ifdef _OPENMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif
    // backword in-place plans
    plan_c0_bw =
        fftw_plan_dft_3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz, c0, c0,
                         FFTW_BACKWARD, FFTW_ESTIMATE);
    plan_c1_bw =
        fftw_plan_dft_3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz, c1, c1,
                         FFTW_BACKWARD, FFTW_ESTIMATE);
    // forward in-place plans
    plan_c0_fw =
        fftw_plan_dft_3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz, c0, c0,
                         FFTW_FORWARD, FFTW_ESTIMATE);
    plan_c1_fw =
        fftw_plan_dft_3d(par->grid_b.nx, par->grid_b.ny, par->grid_b.nz, c1, c1,
                         FFTW_FORWARD, FFTW_ESTIMATE);
  }
}

void Grid_b::export_grid(const Param *par) {
  assert(!par->grid_b.filename.empty());
  std::ofstream output(par->grid_b.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  ham_float tmp;
  for (ham_uint i = 0; i != par->grid_b.full_size; ++i) {
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

void Grid_b::import_grid(const Param *par) {
  assert(!par->grid_b.filename.empty());
  std::ifstream input(par->grid_b.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  ham_float tmp;
  for (ham_uint i = 0; i != par->grid_b.full_size; ++i) {
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
