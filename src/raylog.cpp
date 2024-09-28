#include "raylog.h"
#include "Rcpp.h"

size_t getCurrentMemoryUsage() {
#ifdef _WIN32
  PROCESS_MEMORY_COUNTERS_EX pmc;
  if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) {
    return pmc.WorkingSetSize;
  }
  return 0;
  
#elif defined(__APPLE__)
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  
  if (KERN_SUCCESS == task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count)) {
    return t_info.resident_size;
  }
  return 0;
  
#elif defined(__linux__)
  long rss = 0L;
  FILE* fp = NULL;
  if ((fp = fopen("/proc/self/statm", "r")) == NULL) {
    return 0;
  }
  if (fscanf(fp, "%*s%ld", &rss) != 1) {
    fclose(fp);
    return 0;
  }
  fclose(fp);
  return rss * sysconf(_SC_PAGESIZE);
  
#else
  return 0;
#endif
}

void RayLog::startTimer(const std::string& name) {
#ifdef ENABLE_RAYLOG
  auto now = std::chrono::high_resolution_clock::now();
  threadLocalStartTimes[name] = now;
  getOrCreateTimer(name); 
#endif
}

void RayLog::stopTimer(const std::string& name) {
#ifdef ENABLE_RAYLOG
  auto stop = std::chrono::high_resolution_clock::now();
  auto it = threadLocalStartTimes.find(name);
  if (it != threadLocalStartTimes.end()) {
    auto start = it->second;
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
    
    auto& timer = getOrCreateTimer(name);
    timer.totalNanoseconds.fetch_add(duration.count(), std::memory_order_relaxed);
    timer.count.fetch_add(1, std::memory_order_relaxed);
    
    threadLocalStartTimes.erase(it);
  }
#endif
}

void RayLog::incrementCounter(const std::string& name) {
#ifdef ENABLE_RAYLOG
  auto& counter = getOrCreateCounter(name);
  counter.value.fetch_add(1, std::memory_order_relaxed);
#endif
}

void RayLog::setMemoryUsage(size_t bytes) {
#ifdef ENABLE_RAYLOG
  size_t currentMax = peakMemoryUsage.load(std::memory_order_relaxed);
  while(bytes > currentMax && 
        !peakMemoryUsage.compare_exchange_weak(currentMax, bytes, 
                                               std::memory_order_relaxed, 
                                               std::memory_order_relaxed));
#endif
}

void RayLog::printReport(size_t n_cores) {
#ifdef ENABLE_RAYLOG
  Rcpp::Rcout << "=== Logging Report ===\n";
  Rcpp::Rcout << "CPU Cores Used: " << n_cores << std::endl;
  
  auto printTreeItem = [](const std::string& prefix, const std::string& name, const std::string& value, bool last) {
    Rcpp::Rcout << prefix;
    Rcpp::Rcout << (last ? "└─ " : "├─ ");
    Rcpp::Rcout << std::left << std::setw(20) << name << ": " << value << "\n";
  };
  
  
  // Group timers by context
  std::map<std::string, std::vector<std::pair<std::string, const TimerData&>>> groupedTimers;
  for (const auto& timer : timers) {
    groupedTimers[timer.second.context].emplace_back(timer.first, timer.second);
  }
  
  // Group counters by context
  std::map<std::string, std::vector<std::pair<std::string, const CounterData&>>> groupedCounters;
  for (const auto& counter : counters) {
    groupedCounters[counter.second.context].emplace_back(counter.first, counter.second);
  }
  
  Rcpp::Rcout << "CPU Time:\n";
  for (auto it = groupedTimers.begin(); it != groupedTimers.end(); ++it) {
    const auto& context = it->first;
    const auto& timers = it->second;
    
    bool lastContext = (std::next(it) == groupedTimers.end());
    Rcpp::Rcout << (lastContext ? "└─ " : "├─ ") << context << ":\n";
    
    std::string prefix = lastContext ? "    " : "│   ";
    
    for (size_t i = 0; i < timers.size(); ++i) {
      const auto& timer = timers[i];
      double seconds = timer.second.totalNanoseconds.load(std::memory_order_relaxed) / 1e9;
      if(seconds > 0 && timer.second.count > 0) {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(6)
           << seconds << "s [~" 
           << seconds/n_cores << "s/core + "  
           << int(seconds/timer.second.count * 1e9)  << std::setprecision(0) << "ns/iter]";
        
        printTreeItem(prefix, timer.first, ss.str(), i == timers.size() - 1);
      }
    }
  }
  
  Rcpp::Rcout << "Counters:\n";
  for (auto it = groupedCounters.begin(); it != groupedCounters.end(); ++it) {
    const auto& context = it->first;
    const auto& counters = it->second;
    
    bool lastContext = (std::next(it) == groupedCounters.end());
    Rcpp::Rcout << (lastContext ? "└─ " : "├─ ") << context << ":\n";
    
    std::string prefix = lastContext ? "    " : "│   ";
    
    for (size_t i = 0; i < counters.size(); ++i) {
      const auto& counter = counters[i];
      if(counter.second.value.load(std::memory_order_relaxed) > 0) {
        printTreeItem(prefix, counter.first, 
                      std::to_string(counter.second.value.load(std::memory_order_relaxed)),
                      i == counters.size() - 1);
      }
    }
  }
  
  Rcpp::Rcout << "Peak Memory Usage: " 
              << (peakMemoryUsage.load(std::memory_order_relaxed) / (1024.0 * 1024.0)) 
              << " MB\n";
#endif
}

void RayLog::reset() {
#ifdef ENABLE_RAYLOG
  std::lock_guard<std::mutex> timerLock(timerCreationMutex);
  std::lock_guard<std::mutex> counterLock(counterCreationMutex);
  for (auto& timer : timers) {
    timer.second.totalNanoseconds.store(0, std::memory_order_relaxed);
    timer.second.count.store(0, std::memory_order_relaxed);
  }
  for (auto& counter : counters) {
    counter.second.value.store(0, std::memory_order_relaxed);
  }
  peakMemoryUsage.store(0, std::memory_order_relaxed);
  threadLocalStartTimes.clear();
  while (!contextStack.empty()) {
    contextStack.pop();
  }
#endif
}

void RayLog::printCurrentMemoryUsage(const std::string& point) {
#ifdef ENABLE_RAYLOG
  size_t currentMemory = getCurrentMemoryUsage();
  Rcpp::Rcout << "Memory usage at " << point << ": " 
              << (currentMemory / (1024.0 * 1024.0)) 
              << " MB\n";
#endif
}

std::string RayLog::getCurrentContext() const {
  if (contextStack.empty()) {
    return "Global";
  }
  return contextStack.top();
}

void RayLog::pushContext(const std::string& context) {
#ifdef ENABLE_RAYLOG
  contextStack.push(context);
#endif
}

void RayLog::popContext() {
#ifdef ENABLE_RAYLOG
  if (!contextStack.empty()) {
    contextStack.pop();
  }
#endif
}

thread_local std::unordered_map<std::string, std::chrono::high_resolution_clock::time_point> RayLog::threadLocalStartTimes;
thread_local std::stack<std::string> RayLog::contextStack;
