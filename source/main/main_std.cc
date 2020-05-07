#include <cstdlib>
#include <iostream>
#include <memory>
#include <pipeline.h>
#include <stdexcept>
#include <string>
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
  tmr->start("mainroutine");
  tmr->start("prepfield");
#endif
  auto pipe = std::make_unique<Pipeline>(input);
  pipe->assemble_grid();
  pipe->assemble_b();
  pipe->assemble_te();
  pipe->assemble_cre();
#ifndef NTIMING
  tmr->stop("prepfield");
  tmr->start("LoSintegral");
#endif
  pipe->assemble_obs();
#ifndef NTIMING
  tmr->stop("LoSintegral");
  tmr->stop("mainroutine");
  tmr->print();
#endif
  return EXIT_SUCCESS;
}
