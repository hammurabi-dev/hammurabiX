// unit tests for timer module
// timer record time elapse from Timer::start to Timer::stop
// timer support nested/overlaped recording

#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <thread>
#include <timer.h>

// testing:
// Timer::start
// Timer::stop
TEST(timing, precision) {
  Timer t;
  t.start("main");

  // sleep 5s
  for (int i = 0; i < 2; ++i) {
    t.start("sub");
    std::this_thread::sleep_for(std::chrono::milliseconds(2500));
    t.stop("sub");
  }
  // sleep another 5s
  std::this_thread::sleep_for(std::chrono::milliseconds(5000));
  t.stop("main");

  // docker is not precise
  EXPECT_LT(fabs(t.record.find("sub")->second.second - 5000), 50);
  EXPECT_LT(fabs(t.record.find("main")->second.second - 10000), 100);
}
