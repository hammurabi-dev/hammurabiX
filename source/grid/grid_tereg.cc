// regular thermal electron density field grid

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <grid.h>
#include <hamtype.h>
#include <param.h>

Grid_tereg::Grid_tereg(const Param *par) {
  if (par->grid_tereg.read_permission or par->grid_tereg.write_permission) {
    build_grid(par);
  }
}

void Grid_tereg::build_grid(const Param *par) {
  // allocate spatial domain thermal electron field
  te = std::make_unique<ham_float[]>(par->grid_tereg.full_size);
}

void Grid_tereg::export_grid(const Param *par) {
  assert(!par->grid_tereg.filename.empty());
  std::ofstream output(par->grid_tereg.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  ham_float tmp;
  for (decltype(par->grid_tereg.full_size) i = 0;
       i != par->grid_tereg.full_size; ++i) {
    assert(!output.eof());
    tmp = te[i];
    assert(tmp >= 0);
    output.write(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
  }
  output.close();
}

void Grid_tereg::import_grid(const Param *par) {
  assert(!par->grid_tereg.filename.empty());
  std::ifstream input(par->grid_tereg.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  ham_float tmp;
  for (decltype(par->grid_tereg.full_size) i = 0;
       i != par->grid_tereg.full_size; ++i) {
    assert(!input.eof());
    input.read(reinterpret_cast<char *>(&tmp), sizeof(ham_float));
    te[i] = tmp;
  }
#ifndef NDEBUG
  auto eof = input.tellg();
  input.seekg(0, input.end);
#endif
  assert(eof == input.tellg());
  input.close();
}
