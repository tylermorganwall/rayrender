# Shared R validation, history and dashboard summaries. No historical row is
# subjected to the current measurement contract.
benchmark_cli = function(required = character()) {
  argv = commandArgs(trailingOnly = TRUE)
  if (
    length(argv) %% 2 || any(!startsWith(argv[seq(1, length(argv), 2)], "--"))
  ) {
    stop("Arguments must be --name value pairs")
  }
  args = as.list(argv[seq(2, length(argv), 2)])
  names(args) = gsub("-", "_", sub("^--", "", argv[seq(1, length(argv), 2)]))
  if (length(setdiff(required, names(args)))) {
    stop("Missing: ", paste(setdiff(required, names(args)), collapse = ", "))
  }
  args
}

benchmark_read = function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    colClasses = "character",
    na.strings = c("NA", "")
  )
}

benchmark_write = function(frame, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(
    frame,
    path,
    row.names = FALSE,
    na = "NA",
    fileEncoding = "UTF-8"
  )
}

benchmark_bind = function(frames) {
  columns = unique(unlist(lapply(frames, names), use.names = FALSE))
  frames = lapply(frames, function(x) {
    for (name in setdiff(columns, names(x))) {
      x[[name]] = rep(NA_character_, nrow(x))
    }
    x[columns]
  })
  do.call(rbind, frames)
}

benchmark_field = function(frame, name, default = NA_character_) {
  if (name %in% names(frame)) {
    as.character(frame[[name]])
  } else {
    rep(default, nrow(frame))
  }
}

benchmark_json = function(value) {
  as.character(jsonlite::toJSON(
    value,
    auto_unbox = TRUE,
    na = "null",
    null = "null",
    digits = NA
  ))
}

benchmark_key = function(frame, fields) {
  columns = lapply(fields, function(field) benchmark_field(frame, field))
  vapply(
    seq_len(nrow(frame)),
    function(i) benchmark_json(lapply(columns, `[`, i)),
    character(1)
  )
}

benchmark_settings_key = function(value) {
  tryCatch(
    {
      x = jsonlite::fromJSON(value, simplifyVector = FALSE)
      benchmark_json(x[sort(names(x))])
    },
    error = function(e) NA_character_
  )
}

benchmark_expectation = function(
  configs,
  scenes,
  iterations,
  warmup,
  settings,
  run_id,
  run_attempt,
  commit_sha,
  branch_name,
  require_process = Sys.info()[["sysname"]] == "Linux"
) {
  list(
    run_id = as.character(run_id),
    run_attempt = as.character(run_attempt),
    commit_sha = commit_sha,
    branch_name = branch_name,
    configs = lapply(configs, function(x) {
      list(
        name = x$name,
        revision = x$revision %||% "unversioned",
        backend = x$backend %||% "auto"
      )
    }),
    scenes = scenes,
    iterations = iterations,
    warmup = warmup,
    settings = settings,
    required_metrics = c(
      "render_seconds",
      "total_seconds",
      "package_build_seconds",
      "bvh_build_seconds",
      "bvh_build_count",
      if (isTRUE(settings$time_build)) "scene_build_seconds",
      if (require_process) c("process_elapsed_seconds", "max_rss_mb")
    )
  )
}

benchmark_validate = function(rows, expected) {
  errors = character()
  add = function(message) errors <<- c(errors, message)
  eq = function(field, value) {
    actual = benchmark_field(rows, field)
    if (anyNA(actual) || any(actual != as.character(value))) {
      add(paste("Incorrect", field))
    }
  }
  if (!nrow(rows)) {
    add("No measurements")
  }
  for (field in c("run_id", "run_attempt", "commit_sha", "branch_name")) {
    eq(field, expected[[field]])
  }
  eq("iterations", expected$iterations)
  eq("warmup_iterations", expected$warmup)
  for (field in c("width", "height", "samples", "seed", "threads")) {
    eq(field, expected$settings[[field]])
  }
  eq("time_build", as.character(isTRUE(expected$settings$time_build)))
  settings = vapply(
    benchmark_field(rows, "benchmark_settings_json"),
    benchmark_settings_key,
    character(1)
  )
  if (
    anyNA(settings) ||
      any(settings != benchmark_settings_key(benchmark_json(expected$settings)))
  ) {
    add("Incorrect benchmark settings")
  }
  config_names = vapply(expected$configs, `[[`, character(1), "name")
  if (any(!benchmark_field(rows, "build_config_name") %in% config_names)) {
    add("Unexpected configuration")
  }
  if (any(!benchmark_field(rows, "benchmark_name") %in% expected$scenes)) {
    add("Unexpected scene")
  }
  for (config in expected$configs) {
    part = rows[
      which(benchmark_field(rows, "build_config_name") == config$name),
      ,
      drop = FALSE
    ]
    for (scene in expected$scenes) {
      indices = benchmark_field(part, "iteration_index")[which(
        benchmark_field(part, "benchmark_name") == scene
      )]
      if (
        !identical(
          sort(indices),
          sort(as.character(seq_len(expected$iterations)))
        )
      ) {
        add(paste(
          config$name,
          scene,
          "missing, duplicate or unexpected iterations"
        ))
      }
    }
    if (
      anyNA(benchmark_field(part, "config_revision")) ||
        any(benchmark_field(part, "config_revision") != config$revision)
    ) {
      add(paste(config$name, "incorrect configuration revision"))
    }
    if (
      config$backend != "auto" &&
        any(
          benchmark_field(part, "effective_backend", "unknown") !=
            config$backend,
          na.rm = TRUE
        )
    ) {
      add(paste(config$name, "incorrect backend"))
    }
  }
  statuses = benchmark_field(rows, "status")
  if (anyNA(statuses) || any(statuses != "ok", na.rm = TRUE)) {
    add(paste(
      "Unsuccessful measurements:",
      paste(
        unique(statuses[is.na(statuses) | statuses != "ok"]),
        collapse = ", "
      )
    ))
  }
  for (field in c(
    "compiler_id",
    "cxx_standard",
    "effective_compile_command",
    "effective_flags",
    "effective_backend",
    "build_id",
    "config_revision"
  )) {
    x = benchmark_field(rows, field)
    if (anyNA(x) || any(!nzchar(x))) add(paste("Missing", field))
  }
  eq("compile_commands_valid", "TRUE")
  for (metric in expected$required_metrics) {
    values = suppressWarnings(as.numeric(benchmark_field(rows, metric)))
    nonnegative = metric %in% c("bvh_build_seconds", "scene_build_seconds")
    if (
      any(!is.finite(values) | if (nonnegative) values < 0 else values <= 0)
    ) {
      add(paste("Missing, nonfinite or out-of-range", metric))
    }
  }
  counts = suppressWarnings(as.numeric(benchmark_field(
    rows,
    "bvh_build_count"
  )))
  if (any(is.finite(counts) & counts != floor(counts))) {
    add("Noninteger bvh_build_count")
  }
  # Related intervals have deliberately different scopes, but these bounds must hold.
  render = suppressWarnings(as.numeric(benchmark_field(rows, "render_seconds")))
  bvh = suppressWarnings(as.numeric(benchmark_field(rows, "bvh_build_seconds")))
  total = suppressWarnings(as.numeric(benchmark_field(rows, "total_seconds")))
  if (any(bvh > render + 0.01 | render > total + 0.01, na.rm = TRUE)) {
    add("Inconsistent timing intervals")
  }
  unique(errors)
}

benchmark_step_summary = function(
  rows,
  errors,
  path = Sys.getenv("GITHUB_STEP_SUMMARY")
) {
  lines = c(
    "## Render benchmark validation",
    if (length(errors)) "**FAILED**" else "**PASS**",
    "",
    "| Configuration | Scene | Status | Render median (s) | BVH median (s) | Successful samples |",
    "|---|---|---|---:|---:|---:|"
  )
  groups = split(
    seq_len(nrow(rows)),
    benchmark_key(rows, c("build_config_name", "benchmark_name"))
  )
  for (indices in groups) {
    part = rows[indices, , drop = FALSE]
    ok = benchmark_field(part, "status") %in% "ok"
    med = function(name) {
      values = suppressWarnings(as.numeric(benchmark_field(part, name)))
      values = values[ok & is.finite(values) & values >= 0]
      if (length(values)) {
        format(stats::median(values), digits = 6)
      } else {
        "unavailable"
      }
    }
    lines = c(
      lines,
      paste0(
        "| ",
        paste(
          c(
            benchmark_field(part, "build_config_name")[1],
            benchmark_field(part, "benchmark_name")[1],
            paste(unique(benchmark_field(part, "status")), collapse = ", "),
            med("render_seconds"),
            med("bvh_build_seconds"),
            sum(ok)
          ),
          collapse = " | "
        ),
        " |"
      )
    )
  }
  lines = c(
    lines,
    "",
    if (length(errors)) paste0("- ", errors),
    "",
    "CSV rows and diagnostic logs are retained in the workflow artifacts."
  )
  if (nzchar(path)) {
    cat(paste(lines, collapse = "\n"), "\n", file = path, append = TRUE)
  }
  message(paste(lines, collapse = "\n"))
}

benchmark_merge = function(existing, incoming) {
  keys = c(
    "run_id",
    "run_attempt",
    "commit_sha",
    "branch_name",
    "build_config_name",
    "config_revision",
    "benchmark_name",
    "benchmark_settings_json",
    "iteration_index"
  )
  # Preserve all established rows, including legacy duplicates; replace only exact
  # incoming keys. A rerun attempt has its own identity.
  frames = benchmark_bind(list(existing, incoming))
  incoming = tail(frames, nrow(incoming))
  existing = head(frames, nrow(existing))
  existing = existing[
    !benchmark_key(existing, keys) %in% benchmark_key(incoming, keys),
    ,
    drop = FALSE
  ]
  incoming = incoming[
    !duplicated(benchmark_key(incoming, keys), fromLast = TRUE),
    ,
    drop = FALSE
  ]
  result = rbind(existing, incoming)
  result = result[
    order(benchmark_field(result, "timestamp_utc"), na.last = TRUE),
    ,
    drop = FALSE
  ]
  rownames(result) = NULL
  result
}

benchmark_attempt_key = function(rows) {
  benchmark_key(rows, c("run_id", "run_attempt", "commit_sha", "branch_name"))
}

benchmark_dashboard = function(rows) {
  if (
    !all(
      c("timestamp_utc", "benchmark_name", "build_config_name", "status") %in%
        names(rows)
    )
  ) {
    stop("Invalid benchmark history schema")
  }
  if (!nrow(rows)) {
    return(list(
      row_count = 0,
      attempts = list(),
      summaries = list(),
      builds = list()
    ))
  }
  rows$settings_key = vapply(
    benchmark_field(rows, "benchmark_settings_json"),
    benchmark_settings_key,
    character(1)
  )
  rows$attempt_key = benchmark_attempt_key(rows)
  grouping = c(
    "attempt_key",
    "benchmark_name",
    "build_config_name",
    "config_revision",
    "settings_key",
    "width",
    "height",
    "samples",
    "seed",
    "threads",
    "time_build",
    "warmup_iterations",
    "iterations",
    "cxx_standard",
    "compiler_id",
    "effective_flags",
    "effective_backend",
    "platform_system",
    "platform_machine"
  )
  identity = function(part, fields) {
    setNames(
      lapply(fields, function(field) benchmark_field(part, field)[1]),
      fields
    )
  }
  timestamp = function(part) {
    max(benchmark_field(part, "timestamp_utc"), na.rm = TRUE)
  }
  summaries = lapply(
    split(seq_len(nrow(rows)), benchmark_key(rows, grouping)),
    function(i) {
      part = rows[i, , drop = FALSE]
      result = identity(
        part,
        unique(c(
          grouping,
          "run_id",
          "run_attempt",
          "commit_sha",
          "branch_name",
          "benchmark_settings_json"
        ))
      )
      result$timestamp_utc = timestamp(part)
      result$statuses = as.list(table(
        benchmark_field(part, "status"),
        useNA = "ifany"
      ))
      result$failures = unique(na.omit(benchmark_field(
        part[!benchmark_field(part, "status") %in% "ok", , drop = FALSE],
        "error"
      )))
      result$n = sum(benchmark_field(part, "status") %in% "ok")
      result$row_count = nrow(part)
      for (metric in c(
        "render_seconds",
        "bvh_build_seconds",
        "max_rss_mb",
        "total_seconds"
      )) {
        value = suppressWarnings(as.numeric(benchmark_field(part, metric)))
        value = value[
          benchmark_field(part, "status") %in%
            "ok" &
            is.finite(value) &
            if (metric == "bvh_build_seconds") value >= 0 else value > 0
        ]
        result[[paste0(metric, "_n")]] = length(value)
        result[[metric]] = if (length(value)) stats::median(value) else NA_real_
      }
      result
    }
  )
  attempts = lapply(split(seq_len(nrow(rows)), rows$attempt_key), function(i) {
    part = rows[i, , drop = FALSE]
    result = identity(
      part,
      c("attempt_key", "run_id", "run_attempt", "commit_sha", "branch_name")
    )
    result$timestamp_utc = timestamp(part)
    statuses = benchmark_field(part, "attempt_status")
    result$status = if (all(statuses %in% "complete")) {
      "complete"
    } else if (any(statuses %in% "failed")) {
      "failed"
    } else {
      "legacy / unvalidated"
    }
    result$validation_errors = unique(na.omit(benchmark_field(
      part,
      "validation_errors"
    )))
    result
  })
  attempts = attempts[order(
    vapply(attempts, `[[`, character(1), "timestamp_utc"),
    vapply(attempts, `[[`, character(1), "run_id"),
    suppressWarnings(as.numeric(vapply(
      attempts,
      `[[`,
      character(1),
      "run_attempt"
    ))),
    decreasing = TRUE,
    na.last = TRUE
  )]
  # One R CMD INSTALL per build_id, irrespective of scene/iteration replication.
  builds = lapply(
    split(
      seq_len(nrow(rows)),
      benchmark_key(
        rows,
        c("attempt_key", "build_config_name", "config_revision", "build_id")
      )
    ),
    function(i) {
      part = rows[i, , drop = FALSE]
      result = identity(
        part,
        c("attempt_key", "build_config_name", "config_revision", "build_id")
      )
      values = unique(suppressWarnings(as.numeric(benchmark_field(
        part,
        "package_build_seconds"
      ))))
      values = values[is.finite(values) & values > 0]
      result$seconds = if (
        length(values) == 1 && any(benchmark_field(part, "status") %in% "ok")
      ) {
        values
      } else {
        NA_real_
      }
      result$n = as.integer(is.finite(result$seconds))
      result
    }
  )
  list(
    row_count = nrow(rows),
    attempts = unname(attempts),
    summaries = unname(summaries),
    builds = unname(builds)
  )
}

`%||%` = function(x, y) if (is.null(x) || !length(x) || anyNA(x)) y else x
