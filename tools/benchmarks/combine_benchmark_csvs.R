script = sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
)
source(file.path(dirname(script), "run_render_benchmarks.R"))
args = benchmark_cli(c(
  "input_dir",
  "output",
  "config",
  "configs",
  "benchmarks",
  "iterations",
  "warmup"
))
config = load_config(args$config)
configs = select_configs(config, args$configs)
settings = base_settings(config)
for (name in names(extra_settings(extra_r_args(
  args$extra_r_arg %||% character()
)))) {
  settings[[name]] = extra_settings(extra_r_args(args$extra_r_arg))[[name]]
}
ci = github_context(".")
expected = benchmark_expectation(
  configs,
  benchmark_names(args$benchmarks),
  as.integer(args$iterations),
  as.integer(args$warmup),
  settings,
  args$run_id %||% ci$run_id,
  args$run_attempt %||% ci$run_attempt,
  args$commit_sha %||% git_metadata(".", "current")$commit_sha,
  args$branch_name %||% ci$branch_name
)

# Collect each expected artifact by name. A download/setup failure is represented
# explicitly even if that matrix job never reached the R harness.
frames = list()
for (build_config in configs) {
  paths = list.files(
    args$input_dir,
    pattern = "[.]csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  paths = paths[basename(paths) == paste0(build_config$name, ".csv")]
  error = NULL
  frame = tryCatch(
    {
      if (length(paths) != 1) {
        stop(
          "Expected exactly one configuration artifact; found ",
          length(paths)
        )
      }
      value = benchmark_read(paths)
      if (!nrow(value)) {
        stop("Configuration artifact is empty")
      }
      for (field in c("run_id", "run_attempt", "commit_sha", "branch_name")) {
        actual = benchmark_field(value, field)
        if (anyNA(actual) || any(actual != expected[[field]])) {
          stop(
            "Artifact identity does not match requested ",
            field,
            "; original rows retained in the configuration artifact"
          )
        }
      }
      if (
        anyNA(value$build_config_name) ||
          !"build_config_name" %in% names(value) ||
          any(value$build_config_name != build_config$name)
      ) {
        stop("Artifact configuration does not match filename")
      }
      value
    },
    error = function(e) {
      error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(frame)) {
    frame = data.frame(
      timestamp_utc = timestamp_utc(),
      run_id = expected$run_id,
      run_attempt = expected$run_attempt,
      commit_sha = expected$commit_sha,
      branch_name = expected$branch_name,
      build_config_name = build_config$name,
      config_revision = build_config$revision,
      benchmark_name = expected$scenes,
      benchmark_settings_json = benchmark_json(settings),
      iteration_index = NA_character_,
      status = if (length(paths)) "artifact_invalid" else "artifact_missing",
      error = error
    )
  }
  frames[[build_config$name]] = frame
}
combined = benchmark_bind(frames)
errors = benchmark_validate(combined, expected)
if (!is.null(args$matrix_result) && args$matrix_result != "success") {
  errors = c(errors, paste("Matrix jobs:", args$matrix_result))
}
combined$attempt_status = if (length(errors)) "failed" else "complete"
combined$validation_errors = if (length(errors)) {
  paste(errors, collapse = "; ")
} else {
  NA_character_
}
benchmark_write(combined, args$output)
benchmark_step_summary(combined, errors)
if (length(errors)) {
  stop(
    "Incomplete/invalid comparison; diagnostic CSV written for publication",
    call. = FALSE
  )
}
