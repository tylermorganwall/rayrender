script = sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
benchmark_repo = normalizePath(file.path(dirname(script), "../../.."))
source(file.path(dirname(script), "helper-fixtures.R"))
destination = commandArgs(trailingOnly = TRUE)[1]
f = benchmark_fixture()
good = f$rows
good$attempt_status = "complete"
failed = good
failed$run_attempt = "2"
failed$timestamp_utc = "2026-09-18T01:00:00Z"
failed$status = "build_failed"
failed$error = "compiler failed: <test> & diagnostics </script>"
failed$attempt_status = "failed"
failed$validation_errors = "Both builds failed"
partial = failed
partial$status[partial$effective_backend == "sse"] = "ok"
partial$error[partial$effective_backend == "sse"] = NA_character_
partial$validation_errors = "Scalar build failed; comparison incomplete"
legacy = good
legacy$bvh_build_seconds = NULL
legacy$attempt_status = NULL
legacy$config_revision = NULL
fixtures = list(
  success = good,
  failed = failed,
  partial = partial,
  mixed = benchmark_bind(list(good, failed)),
  legacy = legacy,
  empty = good[FALSE, ]
)
for (name in names(fixtures)) {
  csv = file.path(destination, paste0(name, ".csv"))
  benchmark_write(fixtures[[name]], csv)
  result = benchmark_script(
    "build_benchmark_site.R",
    c("--history-csv", csv, "--site-dir", file.path(destination, name))
  )
  if (result$status != 0) stop(paste(result$output, collapse = "\n"))
}
