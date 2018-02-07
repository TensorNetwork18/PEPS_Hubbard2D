#ifndef _TIMER_H_
#define _TIMER_H_

#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>

class Timer {
  protected:
    std::chrono::high_resolution_clock::time_point tp_;

  public:
    Timer() : tp_(std::chrono::high_resolution_clock::now()) {};

    // return duration in seconds
    double tick() {
      auto now = std::chrono::high_resolution_clock::now();
      double out = std::chrono::duration_cast<std::chrono::nanoseconds>(now - tp_).count()*1.0e-9;
      tp_ = now;
      return out;
    }
};
#endif /*_TIMER_H_*/
