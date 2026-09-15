#ifdef NOT_CRAN
#include "render_jobs.h"
#include <stdexcept>
#include <string>
#include <thread>
#include <testthat.h>

namespace {
void finish_when_cancelled(const std::atomic<bool>& cancelled, std::atomic<int>& finished) {
  // Bound the worker even if a regression prevents cancellation being delivered.
  auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
  while (!cancelled.load() && std::chrono::steady_clock::now() < deadline)
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
  ++finished;
}
}

context("Render sample scheduling") {
  test_that("a persistent pool completes independent batches and empty passes") {
    std::atomic<bool> cancelled(false);
    RcppThread::ThreadPool pool(2);
    std::vector<int> values(8, 0);
    for (int sample = 0; sample < 3; ++sample) {
      std::vector<std::future<void>> futures;
      for (size_t i = 0; i < values.size(); ++i)
        futures.push_back(pool.pushReturn([&, i] { ++values[i]; }));
      expect_true(wait_for_render_jobs(futures, cancelled, [] { return false; }));
      for (int value : values) expect_true(value == sample + 1);
    }
    std::vector<std::future<void>> empty;
    expect_true(wait_for_render_jobs(empty, cancelled, [] { return false; }));
  }

  test_that("preview cancellation is polled on the caller and drains the batch") {
    std::atomic<bool> cancelled(false);
    std::atomic<int> finished(0);
    const auto caller = std::this_thread::get_id();
    RcppThread::ThreadPool pool(2);
    std::vector<std::future<void>> futures;
    for (int i = 0; i < 4; ++i)
      futures.push_back(pool.pushReturn([&] { finish_when_cancelled(cancelled, finished); }));
    bool on_caller = true;
    int polls = 0;
    bool completed = wait_for_render_jobs(futures, cancelled, [&] {
      on_caller = on_caller && std::this_thread::get_id() == caller;
      ++polls;
      return true;
    });
    expect_false(completed);
    expect_true((on_caller && polls > 0));
    expect_true(finished == 4);
  }

  test_that("worker and polling exceptions cancel and drain all captured work") {
    for (bool worker_error : {false, true}) {
      std::atomic<bool> cancelled(false);
      std::atomic<int> finished(0);
      RcppThread::ThreadPool pool(2);
      std::vector<std::future<void>> futures;
      if (worker_error)
        futures.push_back(pool.pushReturn([] { throw std::runtime_error("worker failure"); }));
      for (int i = 0; i < 3; ++i)
        futures.push_back(pool.pushReturn([&] { finish_when_cancelled(cancelled, finished); }));
      std::string error;
      try {
        wait_for_render_jobs(futures, cancelled, [&]() -> bool {
          if (!worker_error) throw std::runtime_error("polling failure");
          return false;
        });
      } catch (const std::runtime_error& caught) {
        error = caught.what();
      }
      expect_true(error == (worker_error ? "worker failure" : "polling failure"));
      expect_true((cancelled && finished == 3));
      // An error in one completed batch must not poison the persistent pool.
      cancelled = false;
      futures.clear();
      futures.push_back(pool.pushReturn([&] { ++finished; }));
      expect_true(wait_for_render_jobs(futures, cancelled, [] { return false; }));
      expect_true(finished == 4);
    }
  }
}
#endif
