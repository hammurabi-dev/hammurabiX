#include <cstdlib>
#include <grid.h>
#include <iostream>
#include <memory>
#include <param.h>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  // helping
  if (argc != 2) {
    std::cout << "wrong input(s)!" << std::endl
              << "hammurabi X requires the path to the XML parameter file"
              << std::endl
              << "try hamx -h for more details." << std::endl;
    throw std::runtime_error("exit before execution");
  }
  const std::string input(argv[1]);
  if (input == "-h") {
    std::cout << "to execute hammurabi X you need to use" << std::endl
              << "hamx [XML parameter file path]" << std::endl
              << "an XML template file can be found in the templates directory"
              << std::endl;
    return EXIT_SUCCESS;
  }

  auto par = std::make_unique<Param>(input);
  auto grid_b = std::make_unique<Grid_b>(par.get());

  // to avoid seg-falt, turn on the "read" and "write" flags in the XML
  // parameter file.
  auto mock_input =
      std::make_unique<std::vector<ham_float>>(3 * par->grid_b.full_size);

  std::cout << "import grid" << std::endl;

  grid_b->import_grid(mock_input.get());

  std::cout << "export grid" << std::endl;

  mock_input.reset(grid_b->export_grid());

  return EXIT_SUCCESS;
}
