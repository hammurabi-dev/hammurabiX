// cosmic ray electron flux phase-space density field grid

#include <array>
#include <cassert>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <grid.h>
#include <param.h>

Grid_cre::Grid_cre(const Param *par) {
  // build up grid when have read or write permission
  if (par->grid_cre.read_permission or par->grid_cre.write_permission) {
    build_grid(par);
  }
}

void Grid_cre::build_grid(const Param *par) {
  // allocate phase-space CRE flux
  cre_flux = std::make_unique<double[]>(par->grid_cre.cre_size);
}

void Grid_cre::export_grid(const Param *par) {
  assert(!par->grid_cre.filename.empty());
  std::ofstream output(par->grid_cre.filename.c_str(),
                       std::ios::out | std::ios::binary);
  assert(output.is_open());
  double tmp;
  for (decltype(par->grid_cre.cre_size) i = 0; i != par->grid_cre.cre_size;
       ++i) {
    assert(!output.eof());
    tmp = cre_flux[i];
    output.write(reinterpret_cast<char *>(&tmp), sizeof(double));
  }
  output.close();
}

void Grid_cre::import_grid(const Param *par) {
  assert(!par->grid_cre.filename.empty());
  std::ifstream input(par->grid_cre.filename.c_str(),
                      std::ios::in | std::ios::binary);
  assert(input.is_open());
  double tmp;
  for (decltype(par->grid_cre.cre_size) i = 0; i != par->grid_cre.cre_size;
       ++i) {
    assert(!input.eof());
    input.read(reinterpret_cast<char *>(&tmp), sizeof(double));
    cre_flux[i] = tmp;
  }
#ifndef NDEBUG
  auto eof = input.tellg();
  input.seekg(0, input.end);
#endif
  assert(eof == input.tellg());
  input.close();
}
