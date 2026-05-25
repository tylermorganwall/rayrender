usage = function() {
  cat(
    paste(
      "Usage:",
      "Rscript tools/benchmarks/append_benchmark_history.R \\",
      "  --history-dir benchmark-history-worktree \\",
      "  --new-results combined-new-results.csv",
      sep = "\n"
    )
  )
}

stop_usage = function(message) {
  usage()
  stop(message, call. = FALSE)
}

parse_args = function(argv) {
  values = list()
  i = 1
  while (i <= length(argv)) {
    key = argv[[i]]
    if (key %in% c("--help", "-h")) {
      usage()
      quit(status = 0)
    }
    if (!startsWith(key, "--") || i == length(argv)) {
      stop_usage(paste0("Invalid argument: ", key))
    }
    name = gsub("-", "_", sub("^--", "", key))
    values[[name]] = argv[[i + 1]]
    i = i + 2
  }
  missing = setdiff(c("history_dir", "new_results"), names(values))
  if (length(missing) > 0) {
    stop_usage(paste0("Missing arguments: ", paste(missing, collapse = ", ")))
  }
  values
}

require_jsonlite = function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("append_benchmark_history.R requires jsonlite", call. = FALSE)
  }
}

read_csv_or_empty = function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE, check.names = FALSE))
  }
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
}

normalize_columns = function(frames) {
  columns = unique(unlist(lapply(frames, names), use.names = FALSE))
  if (length(columns) == 0) {
    return(frames)
  }
  lapply(frames, function(frame) {
    missing = setdiff(columns, names(frame))
    for (name in missing) {
      frame[[name]] = rep("NA", nrow(frame))
    }
    frame[columns]
  })
}

row_key = function(frame, keys) {
  for (key in keys) {
    if (!key %in% names(frame)) {
      frame[[key]] = "NA"
    }
  }
  do.call(paste, c(frame[keys], sep = "\r"))
}

sort_history = function(frame) {
  if (!"timestamp_utc" %in% names(frame)) {
    return(frame)
  }
  timestamps = suppressWarnings(as.POSIXct(
    frame$timestamp_utc,
    tz = "UTC",
    format = "%Y-%m-%dT%H:%M:%SZ"
  ))
  frame[order(timestamps, frame$timestamp_utc, na.last = TRUE), , drop = FALSE]
}

latest_summary = function(frame) {
  if (nrow(frame) == 0) {
    return(list(
      row_count = 0,
      latest_timestamp_utc = NA_character_,
      latest = list()
    ))
  }
  field = function(data, name) {
    if (name %in% names(data) && nrow(data) > 0) {
      return(data[[name]][[1]])
    }
    NA_character_
  }
  if ("timestamp_utc" %in% names(frame)) {
    frame = sort_history(frame)
  }
  latest_row = frame[nrow(frame), , drop = FALSE]
  run_id = field(latest_row, "run_id")
  run_attempt = field(latest_row, "run_attempt")
  same_run = rep(TRUE, nrow(frame))
  if (!is.na(run_id) && nzchar(run_id) && "run_id" %in% names(frame)) {
    same_run = same_run & frame$run_id == run_id
  }
  if (
    !is.na(run_attempt) &&
      nzchar(run_attempt) &&
      "run_attempt" %in% names(frame)
  ) {
    same_run = same_run & frame$run_attempt == run_attempt
  }
  latest = frame[same_run, , drop = FALSE]
  list(
    row_count = nrow(frame),
    latest_row_count = nrow(latest),
    latest_timestamp_utc = field(latest_row, "timestamp_utc"),
    latest_run_id = run_id,
    latest_run_attempt = run_attempt,
    latest_commit_sha = field(latest_row, "commit_sha"),
    latest_branch_name = field(latest_row, "branch_name"),
    latest_status_counts = as.list(table(latest$status %||% "NA"))
  )
}

write_json = function(path, value) {
  require_jsonlite()
  writeLines(
    jsonlite::toJSON(value, auto_unbox = TRUE, pretty = TRUE, null = "null"),
    path
  )
}

`%||%` = function(x, y) {
  if (is.null(x) || length(x) == 0 || any(is.na(x))) {
    y
  } else {
    x
  }
}

args = parse_args(commandArgs(trailingOnly = TRUE))
history_dir = normalizePath(args$history_dir, mustWork = FALSE)
data_dir = file.path(history_dir, "data")
history_csv = file.path(data_dir, "render_benchmarks.csv")
history_rds = file.path(data_dir, "render_benchmarks.rds")
latest_json = file.path(data_dir, "latest.json")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

existing = read_csv_or_empty(history_csv)
new_results = read_csv_or_empty(args$new_results)

if (nrow(new_results) == 0 && nrow(existing) == 0) {
  message("No benchmark rows to append.")
  quit(status = 0)
}

frames = normalize_columns(list(existing, new_results))
combined = do.call(rbind, frames)

dedupe_keys = c(
  "commit_sha",
  "branch_name",
  "build_config_name",
  "benchmark_name",
  "benchmark_settings_json",
  "iteration_index",
  "run_id",
  "run_attempt"
)
keys = row_key(combined, dedupe_keys)
combined = combined[!duplicated(keys, fromLast = TRUE), , drop = FALSE]
combined = sort_history(combined)

utils::write.table(
  combined,
  file = history_csv,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  qmethod = "double",
  fileEncoding = "UTF-8"
)
saveRDS(combined, history_rds, version = 2)
write_json(latest_json, latest_summary(combined))

message("Benchmark history rows: ", nrow(combined))
