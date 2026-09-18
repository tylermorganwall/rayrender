if (!exists("benchmark_repo")) {
  benchmark_repo = normalizePath(file.path(testthat::test_path(), "../../.."))
}
source(file.path(benchmark_repo, "tools/benchmarks/run_render_benchmarks.R"))

benchmark_fixture = function() {
  config = load_config(file.path(
    benchmark_repo,
    "tools/benchmarks/configs/default.json"
  ))
  configs = config$build_configs[2:3]
  settings = base_settings(config)
  settings$time_build = TRUE
  expected = benchmark_expectation(
    configs,
    c("bvh_many_spheres", "bvh_mixed_primitives"),
    2,
    1,
    settings,
    "123",
    "1",
    "abc123",
    "master",
    require_process = TRUE
  )
  rows = expand.grid(
    build_config_name = vapply(configs, `[[`, character(1), "name"),
    benchmark_name = expected$scenes,
    iteration_index = as.character(1:2),
    stringsAsFactors = FALSE
  )
  for (field in c("run_id", "run_attempt", "commit_sha", "branch_name")) {
    rows[[field]] = expected[[field]]
  }
  for (field in c("width", "height", "samples", "seed", "threads")) {
    rows[[field]] = as.character(settings[[field]])
  }
  rows$benchmark_settings_json = benchmark_json(settings)
  rows$timestamp_utc = sprintf("2026-09-17T01:00:%02dZ", seq_len(nrow(rows)))
  rows$iterations = "2"
  rows$warmup_iterations = "1"
  rows$time_build = "TRUE"
  rows$config_revision = "2"
  rows$status = "ok"
  rows$error = NA_character_
  rows$compiler_id = "test compiler"
  rows$cxx_standard = "-std=gnu++20"
  rows$effective_compile_command = "c++ -std=gnu++20 -O3 -march=native -c file.cpp -o file.o"
  rows$effective_flags = "-std=gnu++20 -O3 -march=native"
  rows$effective_backend = ifelse(
    rows$build_config_name == configs[[1]]$name,
    "scalar",
    "sse"
  )
  rows$compile_commands_valid = "TRUE"
  rows$build_id = paste(
    rows$run_id,
    rows$run_attempt,
    rows$build_config_name,
    sep = "/"
  )
  for (metric in expected$required_metrics) {
    rows[[metric]] = "1"
  }
  rows$render_seconds = rep(c("2", "4"), each = 4)
  rows$total_seconds = "5"
  rows$package_build_seconds = "20"
  list(rows = rows, expected = expected)
}

benchmark_script = function(name, args) {
  output = suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    shQuote(c(file.path(benchmark_repo, "tools/benchmarks", name), args)),
    stdout = TRUE,
    stderr = TRUE
  ))
  list(status = attr(output, "status") %||% 0L, output = output)
}

collect_fixture = function(input, output, fixture) {
  benchmark_script(
    "combine_benchmark_csvs.R",
    c(
      "--input-dir",
      input,
      "--output",
      output,
      "--config",
      file.path(benchmark_repo, "tools/benchmarks/configs/default.json"),
      "--configs",
      "o3_native_no_simd,o3_native_simd_sse",
      "--benchmarks",
      "bvh_many_spheres,bvh_mixed_primitives",
      "--iterations",
      "2",
      "--warmup",
      "1",
      "--extra-r-arg",
      "--time-build true",
      "--run-id",
      "123",
      "--run-attempt",
      "1",
      "--commit-sha",
      "abc123",
      "--branch-name",
      "master"
    )
  )
}
