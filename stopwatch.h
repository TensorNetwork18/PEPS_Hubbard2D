#ifndef _STOPWATCH_H
#define _STOPWATCH_H

#include <chrono>


template <typename TimeUnit = std::chrono::nanoseconds, typename Clock = std::chrono::high_resolution_clock>
class Stopwatch {
  public:
    Stopwatch();

    void start();
    void stop();
    void reset();
    TimeUnit elapsed() const;

  private:
    bool switch_;
    TimeUnit duration_;
    std::chrono::time_point<Clock> start_;
    std::chrono::time_point<Clock> stop_;
};

template <typename TimeUnit, typename Clock>
Stopwatch<TimeUnit, Clock>::Stopwatch() : switch_(false), duration_(TimeUnit(0)), start_(TimeUnit(0)), stop_(TimeUnit(0)) {
}

template <typename TimeUnit, typename Clock>
void Stopwatch<TimeUnit, Clock>::start() {
  if ( switch_ )
    return;
  switch_ = true;
  start_ = Clock::now();
}

template <typename TimeUnit, typename Clock>
void Stopwatch<TimeUnit, Clock>::stop() {
  if ( !switch_ )
    return;
  stop_ = Clock::now();
  duration_ += (stop_ - start_);
  switch_ = false;
}

template <typename TimeUnit, typename Clock>
void Stopwatch<TimeUnit, Clock>::reset() {
  switch_ = false;
  duration_ = TimeUnit(0);
}

template <typename TimeUnit, typename Clock>
TimeUnit Stopwatch<TimeUnit, Clock>::elapsed() const {
  return (switch_) ? duration_ + std::chrono::duration_cast<TimeUnit>( Clock::now() - start_ ) : duration_;
}
#endif
