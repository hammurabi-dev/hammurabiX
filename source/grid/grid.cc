#include <cassert>
#include <stdexcept>

#include <grid.h>
#include <param.h>

void Grid::build_grid(const Param *) {
  throw std::runtime_error("wrong inheritance");
}

void Grid::export_grid(const Param *) {
  throw std::runtime_error("wrong inheritance");
}

void Grid::import_grid(const Param *) {
  throw std::runtime_error("wrong inheritance");
}

// END
