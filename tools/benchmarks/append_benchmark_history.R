script = sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
source(file.path(dirname(script), "benchmark_data.R"))
args = benchmark_cli(c("history_dir", "new_results"))
data_dir = file.path(args$history_dir, "data")
history_csv = file.path(data_dir, "render_benchmarks.csv")
existing = if (file.exists(history_csv)) {
  benchmark_read(history_csv)
} else {
  data.frame()
}
incoming = benchmark_read(args$new_results)
if (!nrow(incoming)) {
  stop("No attempt to publish")
}
combined = benchmark_merge(existing, incoming)
benchmark_write(combined, history_csv)
saveRDS(combined, file.path(data_dir, "render_benchmarks.rds"), version = 2)
dashboard = benchmark_dashboard(combined)
attempts = dashboard$attempts
latest = function(values) {
  if (!length(values)) {
    return(NULL)
  }
  values[[order(
    vapply(values, `[[`, character(1), "timestamp_utc"),
    decreasing = TRUE
  )[1]]]
}
writeLines(
  benchmark_json(list(
    row_count = nrow(combined),
    latest_attempt = latest(attempts),
    latest_complete = latest(Filter(
      function(x) x$status == "complete",
      attempts
    ))
  )),
  file.path(data_dir, "latest.json")
)
message("Benchmark history rows: ", nrow(combined))
