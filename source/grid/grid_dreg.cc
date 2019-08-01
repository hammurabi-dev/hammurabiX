// regular dust density field grid

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <grid.h>
#include <param.h>

Grid_dreg::Grid_dreg(const Param *par) {
  if (par->grid_dreg.read_permission or par->grid_dreg.write_permission) {
    build_grid(par);
  }
}

void Grid_dreg::build_grid(const Param *par) {
  // allocate spatial domain dust field
  d = std::make_unique<double[]>(par->grid_dreg.full_size);
}

void Grid_dreg::export_grid(const Param *par) {
  assert(!par->grid_dreg.filename.empty());
  std::ofstream output(par->grid_dreg.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  double tmp;
  for (decltype(par->grid_dreg.full_size) i = 0; i != par->grid_dreg.full_size;
       ++i) {
    assert(!output.eof());
    tmp = d[i];
    assert(tmp >= 0);
    output.write(reinterpret_cast<char *>(&tmp), sizeof(double));
  }
  output.close();
}

void Grid_dreg::import_grid(const Param *par) {
  assert(!par->grid_dreg.filename.empty());
  std::ifstream input(par->grid_dreg.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  double tmp;
  for (decltype(par->grid_dreg.full_size) i = 0; i != par->grid_dreg.full_size;
       ++i) {
    assert(!input.eof());
    input.read(reinterpret_cast<char *>(&tmp), sizeof(double));
    d[i] = tmp;
  }
#ifndef NDEBUG
  auto eof = input.tellg();
  input.seekg(0, input.end);
#endif
  assert(eof == input.tellg());
  input.close();
}
