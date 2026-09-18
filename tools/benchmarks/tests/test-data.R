testthat::test_that("GNU time parses full labels and rejects malformed durations", {
  path = tempfile()
  on.exit(unlink(path))
  for (item in list(
    c("0:04.39", 4.39),
    c("12:34.50", 754.5),
    c("1:02:03.45", 3723.45)
  )) {
    writeLines(
      paste("Elapsed (wall clock) time (h:mm:ss or m:ss):", item[1]),
      path
    )
    testthat::expect_equal(
      as.numeric(parse_time_log(path)$process_elapsed_seconds),
      as.numeric(item[2])
    )
  }
  for (value in c(
    "",
    NA,
    "NA",
    "Inf",
    "NaN",
    "-1",
    "1:2:3:4",
    "1:60",
    "2:61:02",
    "1:x",
    "1:",
    "a:4"
  )) {
    testthat::expect_true(is.na(parse_elapsed_time(value)))
  }
  writeLines("unrelated text", path)
  testthat::expect_equal(parse_time_log(path)$process_elapsed_seconds, "NA")
  testthat::expect_equal(
    parse_time_log(paste0(path, "-absent"))$max_rss_mb,
    "NA"
  )
  writeLines(
    c(
      "Elapsed (wall clock) time (h:mm:ss or m:ss): bad",
      "Maximum resident set size (kbytes): 2048",
      "User time (seconds): 1.5"
    ),
    path
  )
  testthat::expect_equal(parse_time_log(path)$process_elapsed_seconds, "NA")
  testthat::expect_equal(as.numeric(parse_time_log(path)$max_rss_mb), 2)
})

testthat::test_that("only complete successful requested measurements validate", {
  f = benchmark_fixture()
  testthat::expect_length(benchmark_validate(f$rows, f$expected), 0)
  for (status in c("build_failed", "render_failed", "artifact_missing")) {
    rows = f$rows
    rows$status[1] = status
    testthat::expect_match(
      paste(benchmark_validate(rows, f$expected), collapse = " "),
      "Unsuccessful"
    )
  }
  testthat::expect_true(
    length(benchmark_validate(f$rows[-1, ], f$expected)) > 0
  )
  testthat::expect_true(
    length(benchmark_validate(rbind(f$rows, f$rows[1, ]), f$expected)) > 0
  )
  for (field in c(
    "run_id",
    "run_attempt",
    "commit_sha",
    "config_revision",
    "width",
    "samples",
    "seed",
    "time_build",
    "effective_backend",
    "compile_commands_valid",
    "benchmark_settings_json"
  )) {
    rows = f$rows
    rows[[field]][1] = "wrong"
    testthat::expect_true(
      length(benchmark_validate(rows, f$expected)) > 0,
      info = field
    )
  }
  for (metric in f$expected$required_metrics) {
    for (value in c(NA, "Inf", "NaN", "-1")) {
      rows = f$rows
      rows[[metric]][1] = value
      testthat::expect_true(
        length(benchmark_validate(rows, f$expected)) > 0,
        info = paste(metric, value)
      )
    }
  }
  rows = f$rows
  rows$render_seconds[1] = "0"
  testthat::expect_true(length(benchmark_validate(rows, f$expected)) > 0)
  rows = f$rows
  rows$bvh_build_seconds = "0"
  testthat::expect_length(benchmark_validate(rows, f$expected), 0)
})

testthat::test_that("collection fails with explicit diagnostics for missing configuration artifacts", {
  f = benchmark_fixture()
  input = tempfile()
  dir.create(input)
  output = tempfile(fileext = ".csv")
  on.exit(unlink(c(input, output), recursive = TRUE))
  result = collect_fixture(input, output, f)
  testthat::expect_gt(result$status, 0)
  rows = benchmark_read(output)
  testthat::expect_setequal(
    rows$build_config_name,
    c("o3_native_no_simd", "o3_native_simd_sse")
  )
  testthat::expect_true(all(rows$status == "artifact_missing"))
  for (name in unique(f$rows$build_config_name)) {
    benchmark_write(
      f$rows[f$rows$build_config_name == name, ],
      file.path(input, paste0(name, ".csv"))
    )
  }
  result = collect_fixture(input, output, f)
  testthat::expect_equal(
    result$status,
    0,
    info = paste(result$output, collapse = "\n")
  )
  testthat::expect_true(all(
    benchmark_read(output)$attempt_status == "complete"
  ))
  unlink(file.path(input, "o3_native_no_simd.csv"))
  testthat::expect_gt(collect_fixture(input, output, f)$status, 0)
  testthat::expect_true(all(benchmark_read(output)$attempt_status == "failed"))
})

testthat::test_that("history merging is idempotent and preserves rerun attempts and legacy rows", {
  f = benchmark_fixture()
  old = f$rows[
    1,
    c("timestamp_utc", "benchmark_name", "build_config_name", "status")
  ]
  old$status = "build_failed"
  first = benchmark_merge(rbind(old, old), f$rows)
  second = benchmark_merge(first, f$rows)
  testthat::expect_equal(second, first)
  rerun = f$rows
  rerun$run_attempt = "2"
  testthat::expect_equal(nrow(benchmark_merge(first, rerun)), 18)
  testthat::expect_equal(sum(is.na(first$run_id)), 2)
  # Current validation does not scan historical failures or missing new fields.
  testthat::expect_length(benchmark_validate(f$rows, f$expected), 0)
  testthat::expect_silent(benchmark_dashboard(first))
})

testthat::test_that("R summaries group iterations, isolate settings and count builds once", {
  f = benchmark_fixture()
  f$rows$attempt_status = "complete"
  data = benchmark_dashboard(f$rows)
  testthat::expect_length(data$summaries, 4)
  testthat::expect_true(all(vapply(
    data$summaries,
    function(x) x$render_seconds_n == 2 && x$render_seconds == 3,
    logical(1)
  )))
  testthat::expect_length(data$builds, 2)
  testthat::expect_true(all(vapply(
    data$builds,
    function(x) x$n == 1,
    logical(1)
  )))
  rows = f$rows
  rows$timestamp_utc[1] = "2026-09-17T05:00:00Z"
  testthat::expect_length(benchmark_dashboard(rows)$summaries, 4)
  rows$benchmark_settings_json[1] = '{"width":999}'
  testthat::expect_length(benchmark_dashboard(rows)$summaries, 5)
  rows = f$rows
  rows$status[1] = "render_failed"
  rows$render_seconds[1] = "9999"
  rows$bvh_build_seconds = NA
  data = benchmark_dashboard(rows)
  testthat::expect_true(all(vapply(
    data$summaries,
    function(x) is.na(x$bvh_build_seconds),
    logical(1)
  )))
  testthat::expect_true(all(vapply(
    data$summaries,
    function(x) x$render_seconds <= 4,
    logical(1)
  )))
  rows$render_seconds[2] = "Inf"
  testthat::expect_true(all(vapply(
    benchmark_dashboard(rows)$summaries,
    function(x) x$render_seconds <= 4,
    logical(1)
  )))
  testthat::expect_equal(benchmark_dashboard(rows[FALSE, ])$row_count, 0)
})

testthat::test_that("successful comparison is retained when a later attempt fails", {
  rows = benchmark_fixture()$rows
  rows$attempt_status = "complete"
  failed = rows
  failed$run_attempt = "2"
  failed$attempt_status = "failed"
  failed$status = "build_failed"
  data = benchmark_dashboard(benchmark_merge(rows, failed))
  testthat::expect_setequal(
    vapply(data$attempts, `[[`, character(1), "status"),
    c("complete", "failed")
  )
})

testthat::test_that("configuration flags are applied to the selected C++ standard", {
  config = load_config(file.path(
    benchmark_repo,
    "tools/benchmarks/configs/default.json"
  ))
  for (build in config$build_configs[2:3]) {
    testthat::expect_match(
      build$makevars$CXX20FLAGS,
      "-O3 -march=native",
      fixed = TRUE
    )
    testthat::expect_equal(build$revision, "2")
  }
  testthat::expect_equal(
    config$build_configs[[2]]$env$RAYRENDER_DISABLE_SIMD,
    "true"
  )
  testthat::expect_match(
    config$build_configs[[3]]$makevars$CXX20FLAGS,
    "-DRAYSIMD -DHAS_SSE",
    fixed = TRUE
  )
  path = tempfile()
  write_makevars(
    path,
    list(
      CXX20FLAGS = "-O3 -march=native",
      PKG_LIBS = "-ltest",
      PKG_CXXFLAGS = "-DTEST"
    )
  )
  testthat::expect_equal(
    readLines(path),
    c(
      "CXX20FLAGS = -O3 -march=native",
      "PKG_LIBS = -ltest",
      "PKG_CXXFLAGS = -DTEST"
    )
  )
  unlink(path)
})

testthat::test_that("collector isolates invalid identities and duplicate artifacts", {
  f = benchmark_fixture()
  input = tempfile()
  dir.create(input)
  output = tempfile(fileext = ".csv")
  on.exit(unlink(c(input, output), recursive = TRUE))
  for (name in unique(f$rows$build_config_name)) {
    benchmark_write(
      f$rows[f$rows$build_config_name == name, ],
      file.path(input, paste0(name, ".csv"))
    )
  }
  rows = f$rows[f$rows$effective_backend == "scalar", ]
  rows$run_id = "old-run"
  benchmark_write(rows, file.path(input, "o3_native_no_simd.csv"))
  testthat::expect_gt(collect_fixture(input, output, f)$status, 0)
  result = benchmark_read(output)
  testthat::expect_true(all(result$run_id == "123"))
  testthat::expect_true("artifact_invalid" %in% result$status)
  dir.create(file.path(input, "duplicate"))
  file.copy(
    file.path(input, "o3_native_simd_sse.csv"),
    file.path(input, "duplicate", "o3_native_simd_sse.csv")
  )
  testthat::expect_gt(collect_fixture(input, output, f)$status, 0)
  testthat::expect_true(all(
    benchmark_read(output)$status == "artifact_invalid"
  ))
})

testthat::test_that("the step summary reports failures and measured units", {
  f = benchmark_fixture()
  path = tempfile()
  on.exit(unlink(path))
  f$rows$status[1] = "render_failed"
  suppressMessages(benchmark_step_summary(f$rows, "render failed", path))
  text = paste(readLines(path), collapse = "\n")
  testthat::expect_match(text, "FAILED", fixed = TRUE)
  testthat::expect_match(text, "Render median (s)", fixed = TRUE)
  testthat::expect_match(text, "render_failed", fixed = TRUE)
})

testthat::test_that("publication and site generation retain the last valid comparison", {
  rows = benchmark_fixture()$rows
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  incoming = file.path(root, "incoming.csv")
  rows$attempt_status = "complete"
  benchmark_write(rows, incoming)
  args = c("--history-dir", root, "--new-results", incoming)
  testthat::expect_equal(
    benchmark_script("append_benchmark_history.R", args)$status,
    0
  )
  history = file.path(root, "data/render_benchmarks.csv")
  first_hash = tools::md5sum(history)
  testthat::expect_equal(
    benchmark_script("append_benchmark_history.R", args)$status,
    0
  )
  testthat::expect_equal(tools::md5sum(history), first_hash)
  rows$run_attempt = "2"
  rows$attempt_status = "failed"
  rows$status = "render_failed"
  rows$timestamp_utc = "2026-09-18T01:00:00Z"
  benchmark_write(rows, incoming)
  testthat::expect_equal(
    benchmark_script("append_benchmark_history.R", args)$status,
    0
  )
  latest = jsonlite::fromJSON(file.path(root, "data/latest.json"))
  testthat::expect_equal(latest$latest_attempt$run_attempt, "2")
  testthat::expect_equal(latest$latest_complete$run_attempt, "1")
  site_args = c("--history-csv", history, "--site-dir", file.path(root, "site"))
  testthat::expect_equal(
    benchmark_script("build_benchmark_site.R", site_args)$status,
    0
  )
  unlink(history)
  testthat::expect_gt(
    benchmark_script("build_benchmark_site.R", site_args)$status,
    0
  )
  writeLines("wrong,header\n1,2", history)
  testthat::expect_gt(
    benchmark_script("build_benchmark_site.R", site_args)$status,
    0
  )
})

testthat::test_that("effective compiler metadata detects optimization and backend mismatches", {
  compiler = Sys.which("clang++")
  if (!nzchar(compiler)) {
    compiler = Sys.which("g++")
  }
  testthat::skip_if(!nzchar(compiler), "No standalone compiler")
  path = tempfile()
  on.exit(unlink(path))
  config = list(name = "o3_native_no_simd", backend = "scalar")
  command = paste(
    shQuote(compiler),
    "-std=c++20 -O3 -march=native -c PrintClassSizes.cpp -o PrintClassSizes.o"
  )
  writeLines(command, path)
  metadata = effective_build_metadata(path, config, benchmark_repo, character())
  testthat::expect_true(metadata$compile_commands_valid)
  testthat::expect_equal(metadata$effective_backend, "scalar")
  testthat::expect_match(metadata$cxx_standard, "c++20", fixed = TRUE)
  testthat::expect_true(nzchar(metadata$compiler_id))
  writeLines(sub("-O3", "-O2", command, fixed = TRUE), path)
  testthat::expect_false(
    effective_build_metadata(
      path,
      config,
      benchmark_repo,
      character()
    )$compile_commands_valid
  )
  writeLines(command, path)
  config$backend = "sse"
  testthat::expect_false(
    effective_build_metadata(
      path,
      config,
      benchmark_repo,
      character()
    )$compile_commands_valid
  )
})

testthat::test_that("a worker render error exits nonzero and retains structured diagnostics", {
  testthat::skip_if_not_installed("rayrender")
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  scene = file.path(root, "failure.R")
  writeLines(
    'run_benchmark = function(settings) stop("deliberate render failure")',
    scene
  )
  output = file.path(root, "result.json")
  result = benchmark_script(
    "run_one_render_benchmark.R",
    c(
      "--benchmark-name",
      "failure",
      "--benchmark-source",
      scene,
      "--lib-path",
      dirname(find.package("rayrender")),
      "--iteration",
      "1",
      "--seed",
      "1",
      "--width",
      "8",
      "--height",
      "8",
      "--samples",
      "1",
      "--threads",
      "1",
      "--settings-json",
      "{}",
      "--output-json",
      output,
      "--artifact-dir",
      root
    )
  )
  testthat::expect_gt(result$status, 0)
  diagnostic = jsonlite::fromJSON(output)
  testthat::expect_identical(diagnostic$status, "render_failed")
  testthat::expect_identical(diagnostic$error, "deliberate render failure")
})
