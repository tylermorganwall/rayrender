csv_columns = c(
  "timestamp_utc",
  "benchmark_name",
  "benchmark_source",
  "benchmark_settings_json",
  "iteration_index",
  "iterations",
  "warmup_iterations",
  "seed",
  "width",
  "height",
  "samples",
  "threads",
  "platform_system",
  "platform_release",
  "platform_machine",
  "platform_processor",
  "cpu_count_logical",
  "python_version",
  "r_version",
  "package_name",
  "package_version",
  "source_ref",
  "commit_label",
  "commit_sha",
  "git_dirty",
  "build_config_name",
  "compiler_id",
  "cc",
  "cxx",
  "cxxflags",
  "cxx17flags",
  "pkg_cxxflags",
  "makeflags",
  "extra_env_json",
  "build_seconds",
  "install_seconds",
  "render_seconds",
  "status",
  "error",
  "output_hash",
  "artifact_path",
  "build_log_path",
  "benchmark_log_path",
  "workdir_path"
)

supported_makevars = c(
  "CXX",
  "CXX11",
  "CXX14",
  "CXX17",
  "CC",
  "CFLAGS",
  "CXXFLAGS",
  "CXX11FLAGS",
  "CXX14FLAGS",
  "CXX17FLAGS",
  "PKG_CXXFLAGS",
  "PKG_LIBS",
  "MAKEFLAGS"
)

copy_exclude_names = c(
  ".git",
  ".Rproj.user",
  ".cache",
  "00LOCK-rayrender"
)

copy_exclude_suffixes = c(
  ".o",
  ".so",
  ".dll",
  ".dylib",
  ".a",
  ".rds",
  ".RData",
  ".Rhistory",
  ".exr"
)

usage = function() {
  cat(
    paste(
      "Usage:",
      "Rscript tools/benchmarks/run_render_benchmarks.R \\",
      "  --repo . \\",
      "  --ref current \\",
      "  --config tools/benchmarks/configs/default.json \\",
      "  --benchmarks bvh_many_spheres,bvh_mixed_primitives \\",
      "  --iterations 5 \\",
      "  --warmup 1 \\",
      "  --output tools/benchmarks/results/render_benchmark_times.csv",
      "",
      "Required arguments:",
      "  --repo PATH",
      "  --ref REF_OR_current",
      "  --config PATH",
      "  --benchmarks NAME[,NAME...]",
      "  --iterations N",
      "  --warmup N",
      "  --output PATH",
      "",
      "Optional arguments:",
      "  --keep-workdirs",
      "  --workdir-root PATH",
      "  --only-config NAME[,NAME...]",
      "  --timeout-seconds SECONDS",
      "  --extra-r-arg ARG",
      "  --append",
      "  --no-append",
      "  --dry-run",
      "  --help",
      sep = "\n"
    )
  )
}

stop_usage = function(message) {
  usage()
  stop(message, call. = FALSE)
}

parse_args = function(argv) {
  values = list(
    append = TRUE,
    keep_workdirs = FALSE,
    dry_run = FALSE,
    extra_r_arg = character()
  )
  i = 1
  while (i <= length(argv)) {
    key = argv[[i]]
    if (key %in% c("--help", "-h")) {
      usage()
      quit(status = 0)
    }
    if (!startsWith(key, "--")) {
      stop_usage(paste0("Unexpected positional argument: ", key))
    }
    name = gsub("-", "_", sub("^--", "", key))
    if (name %in% c("keep_workdirs", "dry_run", "append", "no_append")) {
      if (name == "no_append") {
        values$append = FALSE
      } else if (name == "append") {
        values$append = TRUE
      } else {
        values[[name]] = TRUE
      }
      i = i + 1
      next
    }
    if (i == length(argv)) {
      stop_usage(paste0("Missing value for ", key))
    }
    value = argv[[i + 1]]
    if (name == "extra_r_arg") {
      values$extra_r_arg = c(values$extra_r_arg, value)
    } else {
      values[[name]] = value
    }
    i = i + 2
  }

  required = c(
    "repo",
    "ref",
    "config",
    "benchmarks",
    "iterations",
    "warmup",
    "output"
  )
  missing = required[!required %in% names(values)]
  if (length(missing) > 0) {
    stop_usage(paste0(
      "Missing required arguments: ",
      paste(missing, collapse = ", ")
    ))
  }

  values$iterations = as.integer(values$iterations)
  values$warmup = as.integer(values$warmup)
  if (is.na(values$iterations) || values$iterations < 1) {
    stop_usage("--iterations must be at least 1")
  }
  if (is.na(values$warmup) || values$warmup < 0) {
    stop_usage("--warmup must be non-negative")
  }
  if ("timeout_seconds" %in% names(values)) {
    values$timeout_seconds = as.numeric(values$timeout_seconds)
    if (is.na(values$timeout_seconds) || values$timeout_seconds <= 0) {
      stop_usage("--timeout-seconds must be positive")
    }
  } else {
    values$timeout_seconds = NA_real_
  }
  values
}

require_jsonlite = function() {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "The benchmark harness requires the jsonlite R package.",
      call. = FALSE
    )
  }
}

read_json = function(path) {
  require_jsonlite()
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

to_json = function(value) {
  require_jsonlite()
  as.character(jsonlite::toJSON(value, auto_unbox = TRUE, null = "null"))
}

quote_command = function(parts) {
  paste(vapply(parts, shQuote, character(1)), collapse = " ")
}

absolute_path = function(path, base = getwd()) {
  if (grepl("^/", path) || grepl("^[A-Za-z]:[/\\\\]", path)) {
    return(path)
  }
  file.path(base, path)
}

safe_name = function(value) {
  safe = gsub("[^A-Za-z0-9_.-]+", "-", value)
  gsub("^-+|-+$", "", safe)
}

run_text = function(
  command,
  args = character(),
  cwd = ".",
  env = character(),
  check = TRUE
) {
  oldwd = setwd(cwd)
  on.exit(setwd(oldwd), add = TRUE)
  shell_args = if (length(args) > 0) {
    vapply(args, shQuote, character(1))
  } else {
    character()
  }
  output = tryCatch(
    suppressWarnings(system2(
      command,
      shell_args,
      stdout = TRUE,
      stderr = TRUE,
      env = env
    )),
    error = function(e) {
      structure(conditionMessage(e), status = 127)
    }
  )
  status = attr(output, "status")
  if (is.null(status)) {
    status = 0
  }
  text = paste(as.character(output), collapse = "\n")
  if (check && status != 0) {
    stop(
      paste0(
        "Command failed (",
        status,
        "): ",
        quote_command(c(command, args)),
        "\n",
        text
      ),
      call. = FALSE
    )
  }
  trimws(text)
}

run_logged = function(
  command,
  args,
  cwd,
  env,
  log_path,
  timeout_seconds = NA_real_
) {
  start = proc.time()[["elapsed"]]
  oldwd = setwd(cwd)
  on.exit(setwd(oldwd), add = TRUE)
  timeout = if (is.na(timeout_seconds)) 0 else timeout_seconds
  shell_args = if (length(args) > 0) {
    vapply(args, shQuote, character(1))
  } else {
    character()
  }
  output = tryCatch(
    suppressWarnings(
      system2(
        command,
        shell_args,
        stdout = TRUE,
        stderr = TRUE,
        env = env,
        timeout = timeout
      )
    ),
    error = function(e) {
      structure(conditionMessage(e), status = 127)
    }
  )
  elapsed = proc.time()[["elapsed"]] - start
  status = attr(output, "status")
  if (is.null(status)) {
    status = 0
  }
  log_dir = dirname(log_path)
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  }
  lines = c(
    paste0("$ ", quote_command(c(command, args))),
    as.character(output)
  )
  if (status == 124 && !is.na(timeout_seconds)) {
    lines = c(
      lines,
      paste0("Command timed out after ", timeout_seconds, " seconds")
    )
  }
  writeLines(lines, log_path)
  list(
    status = status,
    elapsed = elapsed,
    error = if (status == 124) tail(lines, 1) else ""
  )
}

read_tail = function(path, max_chars = 4000) {
  if (!file.exists(path)) {
    return("Log file does not exist")
  }
  text = paste(readLines(path, warn = FALSE), collapse = "\n")
  if (nchar(text) <= max_chars) {
    return(text)
  }
  substring(text, nchar(text) - max_chars + 1, nchar(text))
}

load_config = function(path) {
  config = read_json(path)
  if (!"build_configs" %in% names(config) || !is.list(config$build_configs)) {
    stop("Config JSON must contain a build_configs list", call. = FALSE)
  }
  if (length(config$build_configs) == 0) {
    stop("Config JSON must contain at least one build config", call. = FALSE)
  }
  if (
    !"benchmark_defaults" %in% names(config) ||
      is.null(config$benchmark_defaults)
  ) {
    config$benchmark_defaults = list()
  }
  for (i in seq_along(config$build_configs)) {
    build_config = config$build_configs[[i]]
    if (!"name" %in% names(build_config)) {
      stop("Each build config must have a name", call. = FALSE)
    }
    if (!"env" %in% names(build_config) || is.null(build_config$env)) {
      build_config$env = list()
    }
    if (
      !"makevars" %in% names(build_config) || is.null(build_config$makevars)
    ) {
      build_config$makevars = list()
    }
    unsupported = setdiff(names(build_config$makevars), supported_makevars)
    if (length(unsupported) > 0) {
      stop(
        paste0(
          "Unsupported makevars keys in ",
          build_config$name,
          ": ",
          paste(unsupported, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    config$build_configs[[i]] = build_config
  }
  config
}

select_configs = function(config, only_config) {
  configs = config$build_configs
  if (is.null(only_config) || !nzchar(only_config)) {
    return(configs)
  }
  requested = trimws(strsplit(only_config, ",", fixed = TRUE)[[1]])
  requested = requested[nzchar(requested)]
  selected = configs[vapply(
    configs,
    function(item) item$name %in% requested,
    logical(1)
  )]
  found = vapply(selected, function(item) item$name, character(1))
  missing = setdiff(requested, found)
  if (length(missing) > 0) {
    stop(
      paste0("Unknown build configs: ", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  selected
}

parse_description = function(path) {
  lines = readLines(path, warn = FALSE)
  fields = list()
  current_key = NULL
  for (line in lines) {
    if (!nzchar(line)) {
      next
    }
    if (grepl("^\\s", line) && !is.null(current_key)) {
      fields[[current_key]] = paste(fields[[current_key]], trimws(line))
      next
    }
    split_at = regexpr(":", line, fixed = TRUE)[[1]]
    if (split_at > 0) {
      key = substr(line, 1, split_at - 1)
      value = trimws(substr(line, split_at + 1, nchar(line)))
      fields[[key]] = value
      current_key = key
    }
  }
  fields
}

git_metadata = function(repo, ref) {
  sha = run_text("git", c("rev-parse", "HEAD"), cwd = repo)
  dirty = run_text("git", c("status", "--porcelain"), cwd = repo, check = FALSE)
  list(
    source_ref = ref,
    commit_label = ref,
    commit_sha = sha,
    git_dirty = if (nzchar(dirty)) "true" else "false"
  )
}

is_excluded_path = function(rel_path) {
  parts = strsplit(rel_path, "/", fixed = TRUE)[[1]]
  name = tail(parts, 1)
  if (any(parts %in% copy_exclude_names)) {
    return(TRUE)
  }
  if (
    length(parts) >= 3 &&
      identical(parts[1:3], c("tools", "benchmarks", "results"))
  ) {
    return(TRUE)
  }
  if (grepl("^rayrender_.*[.]tar[.]gz$", name)) {
    return(TRUE)
  }
  if (any(endsWith(name, copy_exclude_suffixes))) {
    return(TRUE)
  }
  if (endsWith(name, ".Rcheck")) {
    return(TRUE)
  }
  FALSE
}

copy_current_tree = function(repo, source) {
  all_paths = list.files(
    repo,
    all.files = TRUE,
    recursive = TRUE,
    no.. = TRUE,
    include.dirs = TRUE
  )
  if (!dir.exists(source)) {
    dir.create(source, recursive = TRUE, showWarnings = FALSE)
  }
  keep = !vapply(all_paths, is_excluded_path, logical(1))
  all_paths = all_paths[keep]
  full_paths = file.path(repo, all_paths)
  info = file.info(full_paths)
  dirs = all_paths[info$isdir %in% TRUE]
  dirs = dirs[order(nchar(dirs))]
  for (rel in dirs) {
    dir.create(file.path(source, rel), recursive = TRUE, showWarnings = FALSE)
  }
  files = all_paths[!(info$isdir %in% TRUE)]
  for (rel in files) {
    from = file.path(repo, rel)
    to = file.path(source, rel)
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(from, to, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
}

prepare_source_tree = function(repo, ref, workdir) {
  source = file.path(workdir, "source")
  if (identical(ref, "current")) {
    metadata = git_metadata(repo, "current")
    copy_current_tree(repo, source)
    return(list(source = source, metadata = metadata, used_worktree = FALSE))
  }

  dir.create(source, recursive = TRUE, showWarnings = FALSE)
  archive = file.path(workdir, "source.tar")
  run_text(
    "git",
    c("archive", "--format=tar", "--output", archive, ref),
    cwd = repo
  )
  utils::untar(archive, exdir = source)
  sha = run_text("git", c("rev-parse", paste0(ref, "^{commit}")), cwd = repo)
  metadata = list(
    source_ref = ref,
    commit_label = ref,
    commit_sha = sha,
    git_dirty = "false"
  )
  list(source = source, metadata = metadata, used_worktree = FALSE)
}

write_makevars = function(path, values) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  keys = names(values)
  lines = vapply(
    keys,
    function(key) paste0(key, " = ", values[[key]]),
    character(1)
  )
  writeLines(lines, path)
}

env_assignments = function(values) {
  if (length(values) == 0) {
    return(character())
  }
  keys = names(values)
  vapply(keys, function(key) paste0(key, "=", values[[key]]), character(1))
}

command_env = function(build_config, makevars_path) {
  values = build_config$env
  values$R_MAKEVARS_USER = makevars_path
  ccache_dir = Sys.getenv("CCACHE_DIR", unset = NA_character_)
  if (is.na(ccache_dir) || !nzchar(ccache_dir)) {
    values$CCACHE_DIR = file.path(dirname(makevars_path), "ccache")
  }
  env_assignments(values)
}

r_cmd_config = function(variable, cwd, env) {
  output = run_text(
    "R",
    c("CMD", "config", variable),
    cwd = cwd,
    env = env,
    check = FALSE
  )
  first = strsplit(output, "\n", fixed = TRUE)[[1]][[1]]
  first = trimws(first)
  if (!nzchar(first) || startsWith(first, "ERROR:")) {
    return("")
  }
  first
}

configured_value = function(build_config, env, name, source) {
  if (name %in% names(build_config$makevars)) {
    return(as.character(build_config$makevars[[name]]))
  }
  if (name %in% names(build_config$env)) {
    return(as.character(build_config$env[[name]]))
  }
  sys_value = Sys.getenv(name, unset = NA_character_)
  if (!is.na(sys_value) && nzchar(sys_value)) {
    return(sys_value)
  }
  if (
    name %in% c("CC", "CXX", "CXX17", "CXXFLAGS", "CXX17FLAGS", "PKG_CXXFLAGS")
  ) {
    return(r_cmd_config(name, source, env))
  }
  "NA"
}

split_command_line = function(value) {
  if (!nzchar(value)) {
    return(character())
  }
  scan(text = value, what = character(), quiet = TRUE)
}

compiler_id = function(build_config, env, source) {
  candidates = c(
    as.character(build_config$makevars$CXX17 %||% ""),
    as.character(build_config$makevars$CXX %||% ""),
    as.character(build_config$env$CXX17 %||% ""),
    as.character(build_config$env$CXX %||% ""),
    r_cmd_config("CXX17", source, env),
    r_cmd_config("CXX", source, env)
  )
  for (candidate in candidates) {
    if (!nzchar(candidate)) {
      next
    }
    parts = split_command_line(candidate)
    if (length(parts) == 0) {
      next
    }
    output = run_text(
      parts[[1]],
      "--version",
      cwd = source,
      env = env,
      check = FALSE
    )
    if (nzchar(output)) {
      return(strsplit(output, "\n", fixed = TRUE)[[1]][[1]])
    }
  }
  "NA"
}

benchmark_names = function(value) {
  names = trimws(strsplit(value, ",", fixed = TRUE)[[1]])
  names = names[nzchar(names)]
  if (length(names) == 0) {
    stop("--benchmarks must name at least one benchmark", call. = FALSE)
  }
  names
}

benchmark_source = function(repo, name) {
  source = file.path(repo, "tools", "benchmarks", "scenes", paste0(name, ".R"))
  if (!file.exists(source)) {
    stop(paste0("Benchmark scene not found: ", source), call. = FALSE)
  }
  source
}

parse_scalar = function(value) {
  value_trim = trimws(value)
  value_lower = tolower(value_trim)
  if (value_lower %in% c("true", "false")) {
    return(value_lower == "true")
  }
  if (value_lower %in% c("null", "na")) {
    return(NA)
  }
  numeric_value = suppressWarnings(as.numeric(value_trim))
  if (
    !is.na(numeric_value) && grepl("^-?[0-9.]+([eE][-+]?[0-9]+)?$", value_trim)
  ) {
    return(numeric_value)
  }
  value_trim
}

extra_r_args = function(values) {
  output = character()
  for (value in values) {
    output = c(output, split_command_line(value))
  }
  output
}

extra_settings = function(extra_args) {
  settings = list()
  i = 1
  while (i <= length(extra_args)) {
    key = extra_args[[i]]
    if (!startsWith(key, "--")) {
      i = i + 1
      next
    }
    name = gsub("-", "_", sub("^--", "", key))
    if (i == length(extra_args) || startsWith(extra_args[[i + 1]], "--")) {
      settings[[name]] = TRUE
      i = i + 1
    } else {
      settings[[name]] = parse_scalar(extra_args[[i + 1]])
      i = i + 2
    }
  }
  settings
}

base_settings = function(config) {
  defaults = config$benchmark_defaults
  defaults$width = as.integer(defaults$width %||% 200)
  defaults$height = as.integer(defaults$height %||% 120)
  defaults$samples = as.integer(defaults$samples %||% 8)
  defaults$seed = as.integer(defaults$seed %||% 1)
  defaults$threads = as.integer(defaults$threads %||% 1)
  defaults
}

timestamp_utc = function() {
  format(
    as.POSIXct(Sys.time(), tz = "UTC"),
    "%Y-%m-%dT%H:%M:%SZ",
    tz = "UTC",
    usetz = FALSE
  )
}

na = function(value) {
  if (is.null(value) || length(value) == 0) {
    return("NA")
  }
  if (is.numeric(value) && any(is.na(value))) {
    return("NA")
  }
  value = as.character(value[[1]])
  if (!nzchar(value)) {
    return("NA")
  }
  value
}

base_row = function(
  benchmark_name,
  benchmark_source_path,
  settings_json,
  args,
  settings,
  metadata,
  build_config,
  env,
  source,
  package_name,
  package_version,
  build_seconds,
  install_seconds,
  build_log_path,
  benchmark_log_path,
  workdir
) {
  sys = Sys.info()
  row = list(
    timestamp_utc = timestamp_utc(),
    benchmark_name = benchmark_name,
    benchmark_source = benchmark_source_path,
    benchmark_settings_json = settings_json,
    iteration_index = "NA",
    iterations = as.character(args$iterations),
    warmup_iterations = as.character(args$warmup),
    seed = as.character(settings$seed),
    width = as.character(settings$width),
    height = as.character(settings$height),
    samples = as.character(settings$samples),
    threads = as.character(settings$threads),
    platform_system = na(sys[["sysname"]]),
    platform_release = na(sys[["release"]]),
    platform_machine = na(sys[["machine"]]),
    platform_processor = na(Sys.getenv("PROCESSOR_IDENTIFIER", unset = "NA")),
    cpu_count_logical = as.character(
      parallel::detectCores(logical = TRUE) %||% "NA"
    ),
    python_version = "NA",
    r_version = R.version.string,
    package_name = package_name,
    package_version = package_version,
    source_ref = metadata$source_ref,
    commit_label = metadata$commit_label,
    commit_sha = metadata$commit_sha,
    git_dirty = metadata$git_dirty,
    build_config_name = build_config$name,
    compiler_id = compiler_id(build_config, env, source),
    cc = na(configured_value(build_config, env, "CC", source)),
    cxx = na(configured_value(build_config, env, "CXX", source)),
    cxxflags = na(configured_value(build_config, env, "CXXFLAGS", source)),
    cxx17flags = na(configured_value(build_config, env, "CXX17FLAGS", source)),
    pkg_cxxflags = na(configured_value(
      build_config,
      env,
      "PKG_CXXFLAGS",
      source
    )),
    makeflags = na(configured_value(build_config, env, "MAKEFLAGS", source)),
    extra_env_json = to_json(build_config$env),
    build_seconds = build_seconds,
    install_seconds = install_seconds,
    render_seconds = "NA",
    status = "NA",
    error = "NA",
    output_hash = "NA",
    artifact_path = "NA",
    build_log_path = build_log_path,
    benchmark_log_path = benchmark_log_path %||% "NA",
    workdir_path = workdir
  )
  row[csv_columns]
}

write_csv_row = function(output, append, row) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  write_header = !append || !file.exists(output) || file.info(output)$size == 0
  row = row[csv_columns]
  row = lapply(row, na)
  frame = as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
  utils::write.table(
    frame,
    file = output,
    sep = ",",
    row.names = FALSE,
    col.names = write_header,
    append = append && file.exists(output) && !write_header,
    qmethod = "double",
    fileEncoding = "UTF-8"
  )
}

worker_command = function(
  repo,
  benchmark_name,
  benchmark_source_path,
  lib_path,
  iteration,
  settings,
  settings_json,
  output_json,
  artifact_dir,
  extra_args
) {
  worker = file.path(repo, "tools", "benchmarks", "run_one_render_benchmark.R")
  c(
    worker,
    "--benchmark-name",
    benchmark_name,
    "--benchmark-source",
    benchmark_source_path,
    "--lib-path",
    lib_path,
    "--iteration",
    as.character(iteration),
    "--seed",
    as.character(settings$seed),
    "--width",
    as.character(settings$width),
    "--height",
    as.character(settings$height),
    "--samples",
    as.character(settings$samples),
    "--threads",
    as.character(settings$threads),
    "--settings-json",
    settings_json,
    "--output-json",
    output_json,
    "--artifact-dir",
    artifact_dir,
    extra_args
  )
}

dry_run = function(args, repo, config, configs, benchmarks, settings_json) {
  cat("Benchmark plan\n")
  cat("  repo: ", repo, "\n", sep = "")
  cat("  ref: ", args$ref, "\n", sep = "")
  cat("  iterations: ", args$iterations, "\n", sep = "")
  cat("  warmup: ", args$warmup, "\n", sep = "")
  cat("  output: ", absolute_path(args$output), "\n", sep = "")
  cat("  settings: ", settings_json, "\n", sep = "")
  for (build_config in configs) {
    cat("\nBuild config: ", build_config$name, "\n", sep = "")
    cat(
      "  install: ",
      quote_command(c(
        "R",
        "CMD",
        "INSTALL",
        "-l",
        "<temp-r-lib>",
        "<source>"
      )),
      "\n",
      sep = ""
    )
    if (length(build_config$makevars) > 0) {
      cat("  makevars:\n")
      for (name in names(build_config$makevars)) {
        cat("    ", name, " = ", build_config$makevars[[name]], "\n", sep = "")
      }
    }
    if (length(build_config$env) > 0) {
      cat("  env: ", to_json(build_config$env), "\n", sep = "")
    }
  }
  cat("\nBenchmarks:\n")
  for (benchmark in benchmarks) {
    cat("  - ", benchmark, "\n", sep = "")
  }
}

cleanup_workdir = function(repo, source_info, workdir, keep_workdirs) {
  if (keep_workdirs) {
    message("Keeping benchmark workdir: ", workdir)
    return()
  }
  if (isTRUE(source_info$used_worktree) && dir.exists(source_info$source)) {
    run_text(
      "git",
      c("worktree", "remove", "--force", source_info$source),
      cwd = repo,
      check = FALSE
    )
  }
  if (dir.exists(workdir)) {
    unlink(workdir, recursive = TRUE, force = TRUE)
  }
}

run_benchmarks = function() {
  args = parse_args(commandArgs(trailingOnly = TRUE))
  launch_cwd = getwd()
  repo = normalizePath(absolute_path(args$repo, launch_cwd), mustWork = TRUE)
  config_path = normalizePath(
    absolute_path(args$config, launch_cwd),
    mustWork = TRUE
  )
  output = absolute_path(args$output, launch_cwd)
  config = load_config(config_path)
  configs = select_configs(config, args$only_config %||% NULL)
  benchmarks = benchmark_names(args$benchmarks)
  invisible(lapply(benchmarks, function(name) benchmark_source(repo, name)))

  settings = base_settings(config)
  flattened_extra_args = extra_r_args(args$extra_r_arg)
  parsed_extra_settings = extra_settings(flattened_extra_args)
  for (name in names(parsed_extra_settings)) {
    settings[[name]] = parsed_extra_settings[[name]]
  }
  settings_json = to_json(settings)

  if (isTRUE(args$dry_run)) {
    dry_run(args, repo, config, configs, benchmarks, settings_json)
    return(invisible(TRUE))
  }

  if (!isTRUE(args$append) && file.exists(output)) {
    unlink(output)
  }
  args$append = TRUE

  root_parent = if (!is.null(args$workdir_root)) {
    normalizePath(
      absolute_path(args$workdir_root, launch_cwd),
      mustWork = FALSE
    )
  } else {
    os_tmp = Sys.getenv("TMPDIR", unset = "")
    if (!nzchar(os_tmp)) {
      os_tmp = "/tmp"
    }
    normalizePath(os_tmp, mustWork = FALSE)
  }
  dir.create(root_parent, recursive = TRUE, showWarnings = FALSE)
  root = tempfile("rayrender-bench-root-", tmpdir = root_parent)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  for (build_config in configs) {
    workdir = tempfile(
      paste0(
        "rayrender-",
        safe_name(args$ref),
        "-",
        safe_name(build_config$name),
        "-"
      ),
      tmpdir = root
    )
    dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
    source_info = NULL
    tryCatch(
      {
        source_info = prepare_source_tree(repo, args$ref, workdir)
        source = source_info$source
        metadata = source_info$metadata
        description = parse_description(file.path(source, "DESCRIPTION"))
        package_name = description$Package %||% "NA"
        package_version = description$Version %||% "NA"
        lib_path = file.path(workdir, "r-lib")
        makevars_path = file.path(workdir, "Makevars")
        log_dir = file.path(workdir, "logs")
        result_dir = file.path(workdir, "worker-results")
        artifact_dir = file.path(workdir, "artifacts")
        dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
        dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
        write_makevars(makevars_path, build_config$makevars)
        env = command_env(build_config, makevars_path)
        runtime_env = env_assignments(build_config$env)
        build_log_path = file.path(log_dir, "R-CMD-INSTALL.log")
        build = run_logged(
          "R",
          c("CMD", "INSTALL", "-l", lib_path, source),
          cwd = source,
          env = env,
          log_path = build_log_path
        )
        build_seconds = sprintf("%.6f", build$elapsed)
        install_seconds = "NA"

        if (build$status != 0) {
          for (benchmark in benchmarks) {
            row = base_row(
              benchmark,
              benchmark_source(repo, benchmark),
              settings_json,
              args,
              settings,
              metadata,
              build_config,
              env,
              source,
              package_name,
              package_version,
              build_seconds,
              install_seconds,
              build_log_path,
              NULL,
              workdir
            )
            row$status = "build_failed"
            row$error = read_tail(build_log_path)
            write_csv_row(output, args$append, row)
          }
        } else {
          for (benchmark in benchmarks) {
            scene_source = benchmark_source(repo, benchmark)
            if (args$warmup > 0) {
              for (warmup_index in seq_len(args$warmup)) {
                warmup_json = file.path(
                  result_dir,
                  paste0(benchmark, "-warmup-", warmup_index, ".json")
                )
                warmup_log = file.path(
                  log_dir,
                  paste0(benchmark, "-warmup-", warmup_index, ".log")
                )
                command_args = worker_command(
                  repo,
                  benchmark,
                  scene_source,
                  lib_path,
                  -warmup_index,
                  settings,
                  settings_json,
                  warmup_json,
                  artifact_dir,
                  flattened_extra_args
                )
                warmup = run_logged(
                  "Rscript",
                  command_args,
                  cwd = repo,
                  env = runtime_env,
                  log_path = warmup_log,
                  timeout_seconds = args$timeout_seconds
                )
                if (warmup$status != 0) {
                  row = base_row(
                    benchmark,
                    scene_source,
                    settings_json,
                    args,
                    settings,
                    metadata,
                    build_config,
                    env,
                    source,
                    package_name,
                    package_version,
                    build_seconds,
                    install_seconds,
                    build_log_path,
                    warmup_log,
                    workdir
                  )
                  row$iteration_index = paste0("warmup_", warmup_index)
                  row$status = "render_failed"
                  row$error = read_tail(warmup_log)
                  write_csv_row(output, args$append, row)
                  next
                }
              }
            }

            for (iteration in seq_len(args$iterations)) {
              output_json = file.path(
                result_dir,
                paste0(benchmark, "-", iteration, ".json")
              )
              benchmark_log = file.path(
                log_dir,
                paste0(benchmark, "-", iteration, ".log")
              )
              command_args = worker_command(
                repo,
                benchmark,
                scene_source,
                lib_path,
                iteration,
                settings,
                settings_json,
                output_json,
                artifact_dir,
                flattened_extra_args
              )
              render = run_logged(
                "Rscript",
                command_args,
                cwd = repo,
                env = runtime_env,
                log_path = benchmark_log,
                timeout_seconds = args$timeout_seconds
              )
              row = base_row(
                benchmark,
                scene_source,
                settings_json,
                args,
                settings,
                metadata,
                build_config,
                env,
                source,
                package_name,
                package_version,
                build_seconds,
                install_seconds,
                build_log_path,
                benchmark_log,
                workdir
              )
              row$iteration_index = as.character(iteration)
              if (render$status != 0 || !file.exists(output_json)) {
                row$status = "render_failed"
                row$error = if (nzchar(render$error)) render$error else
                  read_tail(benchmark_log)
              } else {
                result = read_json(output_json)
                row$status = result$status %||% "NA"
                row$error = na(result$error)
                row$render_seconds = na(result$render_seconds)
                row$output_hash = na(result$output_hash)
                row$artifact_path = na(result$image_path)
              }
              write_csv_row(output, args$append, row)
            }
          }
        }
      },
      finally = {
        if (!is.null(source_info)) {
          cleanup_workdir(repo, source_info, workdir, args$keep_workdirs)
        } else if (!isTRUE(args$keep_workdirs) && dir.exists(workdir)) {
          unlink(workdir, recursive = TRUE, force = TRUE)
        }
      }
    )
  }

  if (!isTRUE(args$keep_workdirs) && dir.exists(root)) {
    unlink(root, recursive = TRUE, force = TRUE)
  } else if (isTRUE(args$keep_workdirs)) {
    message("Benchmark workdir root: ", root)
  }
  invisible(TRUE)
}

`%||%` = function(x, y) {
  if (is.null(x) || length(x) == 0 || any(is.na(x))) {
    y
  } else {
    x
  }
}

run_benchmarks()
