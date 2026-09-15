#ifndef RAYRENDER_RENDER_JOBS_H
#define RAYRENDER_RENDER_JOBS_H

#include <atomic>
#include <chrono>
#include <exception>
#include <future>
#include <vector>
#include "RcppThread.h"

// Wait on completion rather than sleeping between readiness checks. Short
// sample passes wake immediately; long tiles still allow the main R thread to
// service preview cancellation and interrupts at least once every 10 ms.
template <typename PollCancel>
bool wait_for_render_jobs(std::vector<std::future<void>>& futures,
                          std::atomic<bool>& cancelled, PollCancel poll_cancel) {
  std::exception_ptr error;
  auto next_poll = std::chrono::steady_clock::now() + std::chrono::milliseconds(10);
  try {
    for (auto& future : futures) {
      // Share the deadline across tiles: a sequence of individually short
      // waits must not postpone interrupt processing for the whole pass.
      while (future.wait_until(next_poll) == std::future_status::timeout) {
        if (poll_cancel()) {
          cancelled.store(true, std::memory_order_relaxed);
        }
        RcppThread::checkUserInterrupt();
        next_poll = std::chrono::steady_clock::now() + std::chrono::milliseconds(10);
      }
      future.get();
    }
  } catch (...) {
    cancelled.store(true, std::memory_order_relaxed);
    error = std::current_exception();
  }

  // Workers capture scene and sampler state by reference. Drain every task
  // before propagating a worker or main-thread error, so unwinding cannot
  // destroy that state while a persistent pool is still using it.
  for (auto& future : futures) {
    if (!future.valid()) continue;
    try {
      future.get();
    } catch (...) {
      cancelled.store(true, std::memory_order_relaxed);
      if (!error) error = std::current_exception();
    }
  }
  if (error) std::rethrow_exception(error);
  return !cancelled.load(std::memory_order_relaxed);
}

#endif
