// Timer class
//
// Timer records time elapse from Timer::start to Timer::stop
// Timer supports nested/overlaped recording

#ifndef HAMMURABI_TIMR_H
#define HAMMURABI_TIMR_H

#include <cassert>
#include <chrono>
#include <ctime>
#include <hamtype.h>
#include <iostream>
#include <map>
#include <string>

class Timer {
  typedef std::chrono::time_point<std::chrono::high_resolution_clock> tick;
  typedef std::chrono::duration<ham_float, std::milli> duration;
  typedef std::map<std::string, std::pair<tick, ham_float>> timecache;
#ifndef NDEBUG
public:
#endif
  timecache record;
  bool usato = false; // marker for multiple time recording
public:
  Timer() = default;
  virtual ~Timer() = default;
  Timer(const Timer &) = delete;
  Timer(Timer &&) = delete;
  Timer &operator=(const Timer &) = delete;
  Timer &operator=(Timer &&) = delete;
  // record start point of timer
  // 1st argument: name of timing record
  inline void start(std::string flag) {
    record[flag].first = std::chrono::high_resolution_clock::now();
  }
  // record end point of timer
  // then calculate corresponding elapsed time
  // 1st argument: name of timing record
  inline void stop(std::string flag) {
    duration diff =
        (std::chrono::high_resolution_clock::now() - record[flag].first);
    if (usato) { // multiple stop
      record[flag].second += diff.count();
    } else { // 1st stop
      record[flag].second = diff.count();
      usato = true;
    }
  }
  // print to stdout the elapsed time in ms resolution
  // 1st argument: (optional) name of timing record
  // if no name provided, print all records
  inline void print(std::string flag = std::string()) {
    if (flag.empty()) {
      for (auto &i : record)
        std::cout << i.first << " " << i.second.second << " ms" << std::endl;
    } else {
      auto target = record.find(flag);
      assert(target != record.end());
      std::cout << target->first << " " << target->second.second << " ms"
                << std::endl;
    }
  }
};

#endif
