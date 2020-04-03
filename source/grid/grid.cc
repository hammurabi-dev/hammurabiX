#include <grid.h>
#include <param.h>
#include <stdexcept>

void Grid::build_grid(const Param *) {
  throw std::runtime_error("wrong inheritance");
}

void Grid::export_grid(const Param *) {
  throw std::runtime_error("wrong inheritance");
}

void Grid::import_grid(const Param *) {
  throw std::runtime_error("wrong inheritance");
}
