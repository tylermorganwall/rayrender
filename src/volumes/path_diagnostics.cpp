#include "path_diagnostics.h"
#include <iomanip>
#include <sstream>

namespace {
const char *names[] = {"repeated_entry", "invalid_exit", "invalid_density",
                      "invalid_majorant", "invalid_roulette_weight",
                      "position_precision", "invalid_radiance"};
}

void PathDiagnostics::Record(const PathFailure &failure, const Ray &camera,
                             const Ray &current, size_t depth, uint64_t events,
                             uint64_t boundary, const char *stage) {
  const size_t kind = size_t(failure.kind);
  std::lock_guard<std::mutex> lock(mutex);
  ++counts[kind];
  if (examples[kind].size() >= examples_per_kind) return;
  std::ostringstream message;
  message << std::setprecision(17) << failure.what() << "\n  stage=" << stage
          << " depth=" << depth << " subsurface_events=" << events
          << " active_boundary=" << boundary << " time=" << current.time();
  auto describe = [&](const char *label, const Ray &ray) {
    message << "\n  " << label << " origin=(" << ray.o[0] << ", " << ray.o[1]
            << ", " << ray.o[2] << ") direction=(" << ray.d[0] << ", "
            << ray.d[1] << ", " << ray.d[2] << ")";
  };
  describe("camera", camera);
  describe("path", current);
  examples[kind].push_back(message.str());
}

Rcpp::List PathDiagnostics::Take() {
  std::lock_guard<std::mutex> lock(mutex);
  Rcpp::NumericVector totals(kinds);
  Rcpp::CharacterVector labels(kinds);
  Rcpp::List records(kinds);
  uint64_t total = 0;
  for (size_t k = 0; k < kinds; ++k) {
    total += counts[k];
    totals[k] = double(counts[k]);
    labels[k] = names[k];
    records[k] = Rcpp::wrap(examples[k]);
    counts[k] = 0;
    examples[k].clear();
  }
  totals.attr("names") = labels;
  records.attr("names") = labels;
  return Rcpp::List::create(Rcpp::Named("terminated_paths") = double(total),
                            Rcpp::Named("counts") = totals,
                            Rcpp::Named("examples") = records,
                            Rcpp::Named("examples_per_kind") = examples_per_kind);
}
