// sanity tests for parameters

#include <gtest/gtest.h>
#include <memory>
#include <param.h>

// testing:
// template XML parameters
TEST(parameter, parse) {
  auto par = std::make_unique<Param>("templates/params_template.xml");
}
