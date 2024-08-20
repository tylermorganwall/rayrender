#ifndef RAYLOG_H
#define RAYLOG_H

#include <chrono>
#include <string>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <atomic>
#include <memory>
#include <iostream>
#include <iomanip>
#include <stack>

// #define ENABLE_RAYLOG

#include <cstddef>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#elif defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <unistd.h>
#include <sys/resource.h>
#else
#error "Unsupported platform for memory usage tracking"
#endif

size_t getCurrentMemoryUsage();

class RayLog {
public:
  struct TimerData {
    std::atomic<uint64_t> totalNanoseconds;
    std::atomic<size_t> count;
    std::string context;
    
    TimerData() : totalNanoseconds(0), count(0) {}
  };
  
  struct CounterData {
    std::atomic<size_t> value;
    std::string context;
    
    CounterData() : value(0) {}
  };
  
  static RayLog& getInstance() {
    static RayLog instance;
    return instance;
  }
  
  void startTimer(const std::string& name);
  
  void stopTimer(const std::string& name);
  
  void incrementCounter(const std::string& name);
  
  void setMemoryUsage(size_t bytes);
  
  void printReport(size_t n_cores);
  
  void reset();
  void printCurrentMemoryUsage(const std::string& point);
  
  void pushContext(const std::string& context);
  void popContext();
  
  std::string getCurrentContext() const;
  std::string getFullContext() const;
  
private:
  RayLog() = default;
  ~RayLog() = default;
  RayLog(const RayLog&) = delete;
  RayLog& operator=(const RayLog&) = delete;
  
  std::unordered_map<std::string, TimerData> timers;
  std::unordered_map<std::string, CounterData> counters;  
  std::atomic<size_t> peakMemoryUsage{0};
  
  // Thread-local storage for start times
  thread_local static std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> threadLocalStartTimes;
  thread_local static std::stack<std::string> contextStack;
  
  // We need a mutex for creating timers, but this should be infrequent
  std::mutex timerCreationMutex;
  std::mutex counterCreationMutex;
  
  TimerData& getOrCreateTimer(const std::string& name) {
    auto it = timers.find(name);
    if (it == timers.end()) {
      std::lock_guard<std::mutex> lock(timerCreationMutex);
      auto& timer = timers[name];
      timer.context = getCurrentContext();
      return timer;
    }
    return it->second;
  }
  
  CounterData& getOrCreateCounter(const std::string& name) {
    auto it = counters.find(name);
    if (it == counters.end()) {
      std::lock_guard<std::mutex> lock(timerCreationMutex);
      auto& counter = counters[name];
      counter.context = getCurrentContext();
      return counter;
    }
    return it->second;
  }
};

// New RAII class for automatic timer management
class ScopedTimer {
public:
  ScopedTimer(const std::string& name) : name_(name) {
    RayLog::getInstance().startTimer(name_);
  }
  
  ~ScopedTimer() {
    RayLog::getInstance().stopTimer(name_);
  }
  
private:
  std::string name_;
};

// New RAII class for automatic timer management
class ScopedTimerCounter {
public:
  ScopedTimerCounter(const std::string& name) : name_(name) {
    RayLog::getInstance().startTimer(name_);
    RayLog::getInstance().incrementCounter(name_);
  }
  
  ~ScopedTimerCounter() {
    RayLog::getInstance().stopTimer(name_);
  }
  
private:
  std::string name_;
};

// RAII class for automatic context management
class ScopedContext {
public:
  ScopedContext(const std::string& context) {
    RayLog::getInstance().pushContext(context);
  }
  
  ~ScopedContext() {
    RayLog::getInstance().popContext();
  }
};



#ifdef ENABLE_RAYLOG
  #define START_TIMER(name) RayLog::getInstance().startTimer(name)
  #define STOP_TIMER(name) RayLog::getInstance().stopTimer(name)
  #define INCREMENT_COUNTER(name) RayLog::getInstance().incrementCounter(name)
  #define QUERY_MEMORY_USAGE() RayLog::getInstance().setMemoryUsage(getCurrentMemoryUsage())
  #define PRINT_LOG_REPORT(numbercores) RayLog::getInstance().printReport(numbercores)
  #define RESET_RAYLOG() RayLog::getInstance().reset()
  #define PRINT_CURRENT_MEMORY(point) RayLog::getInstance().printCurrentMemoryUsage(point)
  #define SCOPED_TIMER(name) ScopedTimer scopedTimer##__LINE__(name)
  #define SCOPED_TIMER_COUNTER(name) ScopedTimerCounter scopedTimer##__LINE__(name)
  #define PUSH_CONTEXT(context) RayLog::getInstance().pushContext(context)
  #define POP_CONTEXT() RayLog::getInstance().popContext()
  #define SCOPED_CONTEXT(context) ScopedContext scopedContext##__LINE__(context)
#else
  #define START_TIMER(name)
  #define STOP_TIMER(name)
  #define INCREMENT_COUNTER(name)
  #define QUERY_MEMORY_USAGE(bytes)
  #define RESET_RAYLOG()
  #define PRINT_LOG_REPORT(numbercores)
  #define PRINT_CURRENT_MEMORY(point)
  #define SCOPED_TIMER(name)
  #define SCOPED_TIMER_COUNTER(name) 
  #define PUSH_CONTEXT(context)
  #define POP_CONTEXT()
  #define SCOPED_CONTEXT(context)
#endif
#endif // RAYLOG_H