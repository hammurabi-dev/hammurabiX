#include <cstdlib>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <pipeline.h>
#include <timer.h>

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
#ifndef NTIMING
  auto tmr = std::make_unique<Timer>();
  tmr->start("main");
#endif
  auto run = std::make_unique<Pipeline>(input);
  run->assemble_grid();
  run->assemble_tereg();
  run->assemble_breg();
  run->assemble_ternd();
  run->assemble_brnd();
  run->assemble_cre();
  run->assemble_obs();
#ifndef NTIMING
  tmr->stop("main");
  tmr->print();
#endif
  return EXIT_SUCCESS;
}
