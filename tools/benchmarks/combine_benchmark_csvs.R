usage = function() {
  cat(
    paste(
      "Usage:",
      "Rscript tools/benchmarks/combine_benchmark_csvs.R \\",
      "  --input-dir downloaded-benchmark-results \\",
      "  --output combined-new-results.csv",
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
  missing = setdiff(c("input_dir", "output"), names(values))
  if (length(missing) > 0) {
    stop_usage(paste0("Missing arguments: ", paste(missing, collapse = ", ")))
  }
  values
}

read_csv = function(path) {
  utils::read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
}

normalize_columns = function(frames) {
  columns = unique(unlist(lapply(frames, names), use.names = FALSE))
  lapply(frames, function(frame) {
    missing = setdiff(columns, names(frame))
    for (name in missing) {
      frame[[name]] = "NA"
    }
    frame[columns]
  })
}

args = parse_args(commandArgs(trailingOnly = TRUE))
input_dir = normalizePath(args$input_dir, mustWork = TRUE)
csv_paths = list.files(
  input_dir,
  pattern = "[.]csv$",
  recursive = TRUE,
  full.names = TRUE
)

dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)

if (length(csv_paths) == 0) {
  file.create(args$output)
  quit(status = 0)
}

frames = lapply(csv_paths, read_csv)
frames = frames[vapply(frames, nrow, integer(1)) > 0]

if (length(frames) == 0) {
  file.create(args$output)
  quit(status = 0)
}

combined = do.call(rbind, normalize_columns(frames))
utils::write.table(
  combined,
  file = args$output,
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  qmethod = "double",
  fileEncoding = "UTF-8"
)
