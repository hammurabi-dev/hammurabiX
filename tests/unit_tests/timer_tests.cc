/**
 * unit tests for timer module
 * timer record time elapse from Timer::start to Timer::stop
 * timer support nested/overlaped recording
 * feel free to add more rational testing blocks
 */
#include <gtest/gtest.h>

#include <timer.h>
#include <chrono>
#include <thread>
#include <cmath>

TEST(timing, precision){
    Timer t;
    t.start("main");
    
    t.start("sub");
    // sleep 5s
    std::this_thread::sleep_for(std::chrono::milliseconds(5000));
    t.stop("sub");
    
    // sleep another 5s
    std::this_thread::sleep_for(std::chrono::milliseconds(5000));
    t.stop("main");
    
    EXPECT_LT(fabs(t.record.find("sub")->second.second - 5000),1);
    EXPECT_LT(fabs(t.record.find("main")->second.second - 10000),1);
}
