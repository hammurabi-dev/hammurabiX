// regular magnetic vector field grid

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

Grid_breg::Grid_breg(const Param *par) {
  // build up grid when have read or write permission
  if (par->grid_breg.read_permission or par->grid_breg.write_permission) {
    build_grid(par);
  }
}

void Grid_breg::build_grid(const Param *par) {
  // allocate spatial domain regular magnetic field
  bx = std::make_unique<ham_float[]>(par->grid_breg.full_size);
  by = std::make_unique<ham_float[]>(par->grid_breg.full_size);
  bz = std::make_unique<ham_float[]>(par->grid_breg.full_size);
}

void Grid_breg::export_grid(const Param *par) {
  assert(!par->grid_breg.filename.empty());
  std::ofstream output(par->grid_breg.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  ham_float tmp;
  for (decltype(par->grid_breg.full_size) i = 0; i != par->grid_breg.full_size;
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

void Grid_breg::import_grid(const Param *par) {
  assert(!par->grid_breg.filename.empty());
  std::ifstream input(par->grid_breg.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  ham_float tmp;
  for (decltype(par->grid_breg.full_size) i = 0; i != par->grid_breg.full_size;
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
