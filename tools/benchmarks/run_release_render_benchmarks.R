usage = function() {
  cat(
    paste(
      "Usage:",
      "Rscript tools/benchmarks/run_release_render_benchmarks.R \\",
      "  --repo . \\",
      "  --config tools/benchmarks/configs/default.json \\",
      "  --benchmarks bvh_many_spheres \\",
      "  --iterations 1 \\",
      "  --warmup 0 \\",
      "  --output tools/benchmarks/results/release_render_benchmark_times.csv",
      "",
      "Runs render benchmarks for commits whose commit message contains a",
      "release version matching v0.x.y.",
      "",
      "Options:",
      "  --repo PATH                    Default: .",
      "  --config PATH                  Default: tools/benchmarks/configs/default.json",
      "  --benchmarks NAME[,NAME...]     Default: bvh_many_spheres",
      "  --iterations N                 Default: 1",
      "  --warmup N                     Default: 0",
      "  --output PATH                  Default: tools/benchmarks/results/release_render_benchmark_times.csv",
      "  --manifest-output PATH         Default: output path with _commits.csv suffix",
      "  --release-regex REGEX          Default: v0\\.[0-9]+\\.[0-9]+",
      "  --cxx-std STD                  Default: auto. One of auto, none, CXX11, CXX14, CXX17, CXX20",
      "  --max-releases N               Limit how many matching release commits run",
      "  --newest-first                 Run newest matching commits first",
      "  --all-refs                     Search branches, tags, and remotes instead of current first-parent history",
      "  --keep-workdirs                Forwarded to run_render_benchmarks.R",
      "  --workdir-root PATH            Forwarded to run_render_benchmarks.R",
      "  --only-config NAME[,NAME...]    Forwarded to run_render_benchmarks.R",
      "  --timeout-seconds SECONDS      Forwarded to run_render_benchmarks.R",
      "  --extra-r-arg ARG              Forwarded to run_render_benchmarks.R; may be repeated",
      "  --append                       Append to output CSV. Default.",
      "  --no-append                    Remove output CSV before starting this release sweep",
      "  --dry-run                      Print the planned release commits and benchmark commands",
      "  --stop-on-error                Stop after the first benchmark subprocess error",
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
    repo = ".",
    config = file.path("tools", "benchmarks", "configs", "default.json"),
    benchmarks = "bvh_many_spheres",
    iterations = "1",
    warmup = "0",
    output = file.path(
      "tools",
      "benchmarks",
      "results",
      "release_render_benchmark_times.csv"
    ),
    release_regex = "v0\\.[0-9]+\\.[0-9]+",
    cxx_std = "auto",
    append = TRUE,
    all_refs = FALSE,
    newest_first = FALSE,
    keep_workdirs = FALSE,
    dry_run = FALSE,
    stop_on_error = FALSE,
    extra_r_arg = character()
  )

  value_options = c(
    "repo",
    "config",
    "benchmarks",
    "iterations",
    "warmup",
    "output",
    "manifest_output",
    "release_regex",
    "cxx_std",
    "max_releases",
    "workdir_root",
    "only_config",
    "timeout_seconds",
    "extra_r_arg"
  )
  flag_options = c(
    "append",
    "no_append",
    "all_refs",
    "newest_first",
    "keep_workdirs",
    "dry_run",
    "stop_on_error"
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
    if (name %in% flag_options) {
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

    if (!name %in% value_options) {
      stop_usage(paste0("Unknown option: ", key))
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

  values$iterations = parse_count(values$iterations, "iterations")
  values$warmup = parse_count(values$warmup, "warmup")
  if (!is.null(values$max_releases)) {
    values$max_releases = parse_count(values$max_releases, "max-releases")
  }
  values$cxx_std = normalize_cxx_std(values$cxx_std)
  values
}

parse_count = function(value, name) {
  number = suppressWarnings(as.integer(value))
  if (
    is.na(number) || number < 0 || as.character(number) != as.character(value)
  ) {
    stop_usage(paste0("--", name, " must be a non-negative integer"))
  }
  number
}

normalize_cxx_std = function(value) {
  normalized = toupper(trimws(value))
  normalized = gsub("\\+", "X", normalized)
  normalized = gsub("[^A-Z0-9]", "", normalized)

  if (normalized %in% c("AUTO", "")) {
    return("auto")
  }
  if (normalized %in% c("NONE", "FALSE", "OFF")) {
    return("none")
  }
  if (normalized %in% c("11", "CXX11", "CXXSTD11", "CXX11STD", "CXX11X")) {
    return("CXX11")
  }
  if (normalized %in% c("14", "CXX14", "CXXSTD14", "CXX14STD", "CXX14X")) {
    return("CXX14")
  }
  if (normalized %in% c("17", "CXX17", "CXXSTD17", "CXX17STD", "CXX17X")) {
    return("CXX17")
  }
  if (normalized %in% c("20", "CXX20", "CXXSTD20", "CXX20STD", "CXX20X")) {
    return("CXX20")
  }

  stop_usage("--cxx-std must be one of auto, none, CXX11, CXX14, CXX17, CXX20")
}

absolute_path = function(path, base = getwd()) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path
  } else {
    file.path(base, path)
  }
}

derive_manifest_output = function(output) {
  if (grepl("\\.csv$", output, ignore.case = TRUE)) {
    sub("\\.csv$", "_commits.csv", output, ignore.case = TRUE)
  } else {
    paste0(output, "_commits.csv")
  }
}

current_script_path = function(base) {
  file_arg = commandArgs(FALSE)
  file_arg = file_arg[startsWith(file_arg, "--file=")]
  if (length(file_arg)) {
    return(normalizePath(
      absolute_path(sub("^--file=", "", file_arg[[1]]), base),
      mustWork = TRUE
    ))
  }

  normalizePath(
    absolute_path(
      file.path("tools", "benchmarks", "run_release_render_benchmarks.R"),
      base
    ),
    mustWork = TRUE
  )
}

run_git_text = function(args, allow_failure = FALSE) {
  stdout = tempfile("rayrender-release-git-stdout-")
  stderr = tempfile("rayrender-release-git-stderr-")
  on.exit(unlink(c(stdout, stderr)), add = TRUE)

  status = tryCatch(
    system2("git", args, stdout = stdout, stderr = stderr),
    error = function(error) {
      message(conditionMessage(error))
      1
    }
  )
  if (is.null(status)) {
    status = 0
  }
  if (status != 0) {
    if (isTRUE(allow_failure)) {
      return("")
    }
    error_size = file.info(stderr)$size
    error_text = if (!is.na(error_size) && error_size > 0) {
      readChar(stderr, error_size, useBytes = TRUE)
    } else {
      "git log failed"
    }
    stop(error_text, call. = FALSE)
  }

  output_size = file.info(stdout)$size
  if (is.na(output_size) || output_size == 0) {
    return("")
  }
  readChar(stdout, output_size, useBytes = TRUE)
}

discover_release_commits = function(
  repo,
  release_regex,
  all_refs,
  newest_first
) {
  args = c("-C", repo, "log")
  if (isTRUE(all_refs)) {
    args = c(args, "--branches", "--tags", "--remotes")
  } else {
    args = c(args, "--first-parent")
  }
  if (!isTRUE(newest_first)) {
    args = c(args, "--reverse")
  }
  args = c(args, "--format=%H%x09%B%x1e")

  output = run_git_text(args)
  if (!nzchar(output)) {
    return(data.frame(
      release_version = character(),
      commit_sha = character(),
      subject = character()
    ))
  }

  entries = strsplit(output, "\036", fixed = TRUE)[[1]]
  entries = entries[nzchar(trimws(entries))]
  entries = trimws(entries, which = "left")
  tab = regexpr("\t", entries, fixed = TRUE)
  entries = entries[tab > 0]
  tab = tab[tab > 0]

  sha = substr(entries, 1, tab - 1)
  message = substring(entries, tab + 1)
  message = sub("[\r\n]+$", "", message)
  subject = vapply(
    strsplit(message, "\n", fixed = TRUE),
    function(lines) {
      lines[[1]]
    },
    character(1)
  )

  version_match = regexpr(release_regex, message, perl = TRUE)
  matched = version_match > 0
  releases = data.frame(
    release_version = regmatches(message, version_match),
    commit_sha = sha[matched],
    subject = subject[matched],
    stringsAsFactors = FALSE
  )
  releases[!duplicated(releases$commit_sha), , drop = FALSE]
}

git_show_file = function(repo, commit_sha, path) {
  run_git_text(
    c("-C", repo, "show", paste0(commit_sha, ":", path)),
    allow_failure = TRUE
  )
}

detect_cxx_std = function(repo, commit_sha, mode) {
  if (identical(mode, "none")) {
    return(NA_character_)
  }
  if (!identical(mode, "auto")) {
    return(mode)
  }

  description = git_show_file(repo, commit_sha, "DESCRIPTION")
  if (grepl("C\\+\\+20", description)) {
    return("CXX20")
  }
  if (
    grepl("C\\+\\+17", description) ||
      grepl("C\\+\\+14", description) ||
      grepl("C\\+\\+11", description)
  ) {
    return("CXX17")
  }

  "CXX17"
}

write_release_config = function(config_path, output_dir, release) {
  cxx_std = release$cxx_std[[1]]
  if (is.na(cxx_std) || !nzchar(cxx_std)) {
    return(config_path)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The release benchmark wrapper requires the jsonlite R package.")
  }

  config = jsonlite::read_json(config_path, simplifyVector = FALSE)
  for (i in seq_along(config$build_configs)) {
    build_config = config$build_configs[[i]]
    if (
      !"makevars" %in% names(build_config) || is.null(build_config$makevars)
    ) {
      build_config$makevars = list()
    }
    build_config$makevars$CXX_STD = cxx_std
    config$build_configs[[i]] = build_config
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  release_config = file.path(
    output_dir,
    paste0(
      "config-",
      substr(release$commit_sha[[1]], 1, 12),
      "-",
      tolower(cxx_std),
      ".json"
    )
  )
  jsonlite::write_json(config, release_config, auto_unbox = TRUE, pretty = TRUE)
  release_config
}

write_manifest = function(releases, output, args, started_at) {
  manifest = data.frame(
    run_started_utc = started_at,
    repo = args$repo,
    release_regex = args$release_regex,
    cxx_std_mode = args$cxx_std,
    release_version = releases$release_version,
    commit_sha = releases$commit_sha,
    cxx_std = releases$cxx_std,
    subject = releases$subject,
    stringsAsFactors = FALSE
  )

  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(manifest, output, row.names = FALSE, na = "")
}

benchmark_args_for_release = function(
  args,
  repo,
  script,
  release,
  config_path
) {
  command_args = c(
    script,
    "--repo",
    repo,
    "--ref",
    release$commit_sha,
    "--config",
    config_path,
    "--benchmarks",
    args$benchmarks,
    "--iterations",
    as.character(args$iterations),
    "--warmup",
    as.character(args$warmup),
    "--output",
    args$output,
    "--append"
  )

  if (isTRUE(args$keep_workdirs)) {
    command_args = c(command_args, "--keep-workdirs")
  }
  if (isTRUE(args$dry_run)) {
    command_args = c(command_args, "--dry-run")
  }
  if (!is.null(args$workdir_root)) {
    command_args = c(command_args, "--workdir-root", args$workdir_root)
  }
  if (!is.null(args$only_config)) {
    command_args = c(command_args, "--only-config", args$only_config)
  }
  if (!is.null(args$timeout_seconds)) {
    command_args = c(command_args, "--timeout-seconds", args$timeout_seconds)
  }
  for (extra_arg in args$extra_r_arg) {
    command_args = c(command_args, "--extra-r-arg", shQuote(extra_arg))
  }

  command_args
}

shell_quote = function(value) {
  paste(shQuote(value), collapse = " ")
}

run_release_benchmarks = function() {
  launch_cwd = getwd()
  args = parse_args(commandArgs(trailingOnly = TRUE))
  args$repo = normalizePath(
    absolute_path(args$repo, launch_cwd),
    mustWork = TRUE
  )
  args$config = normalizePath(
    absolute_path(args$config, launch_cwd),
    mustWork = TRUE
  )
  args$output = absolute_path(args$output, launch_cwd)
  args$manifest_output = absolute_path(
    args$manifest_output %||% derive_manifest_output(args$output),
    launch_cwd
  )

  script = normalizePath(
    file.path(
      dirname(current_script_path(launch_cwd)),
      "run_render_benchmarks.R"
    ),
    mustWork = TRUE
  )

  releases = discover_release_commits(
    args$repo,
    args$release_regex,
    args$all_refs,
    args$newest_first
  )
  if (!nrow(releases)) {
    stop("No release commits matched ", args$release_regex, call. = FALSE)
  }
  if (!is.null(args$max_releases)) {
    releases = utils::head(releases, args$max_releases)
  }
  releases$cxx_std = vapply(
    releases$commit_sha,
    function(commit_sha) detect_cxx_std(args$repo, commit_sha, args$cxx_std),
    character(1)
  )

  cat("Release benchmark sweep\n")
  cat("  repo: ", args$repo, "\n", sep = "")
  cat("  releases: ", nrow(releases), "\n", sep = "")
  cat("  cxx-std: ", args$cxx_std, "\n", sep = "")
  cat("  output: ", args$output, "\n", sep = "")
  cat("  manifest: ", args$manifest_output, "\n", sep = "")

  for (i in seq_len(nrow(releases))) {
    cat(
      sprintf(
        "  %3d. %s %s [%s] %s\n",
        i,
        releases$release_version[[i]],
        substr(releases$commit_sha[[i]], 1, 12),
        releases$cxx_std[[i]],
        releases$subject[[i]]
      )
    )
  }

  if (!isTRUE(args$dry_run)) {
    if (!isTRUE(args$append) && file.exists(args$output)) {
      unlink(args$output)
    }
    write_manifest(
      releases,
      args$manifest_output,
      args,
      format(Sys.time(), tz = "UTC")
    )
  }

  failures = data.frame(
    release_version = character(),
    commit_sha = character(),
    status = integer(),
    stringsAsFactors = FALSE
  )

  config_dir = tempfile("rayrender-release-configs-")
  dir.create(config_dir, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(releases))) {
    release = releases[i, , drop = FALSE]
    config_path = write_release_config(args$config, config_dir, release)
    command_args = benchmark_args_for_release(
      args,
      args$repo,
      script,
      release,
      config_path
    )
    cat(
      "\nRunning ",
      release$release_version,
      " at ",
      substr(release$commit_sha, 1, 12),
      " with ",
      release$cxx_std,
      "\n",
      sep = ""
    )
    if (isTRUE(args$dry_run)) {
      cat("Rscript ", shell_quote(command_args), "\n", sep = "")
    }

    status = system2("Rscript", command_args)
    if (!is.null(status) && status != 0) {
      failures = rbind(
        failures,
        data.frame(
          release_version = release$release_version,
          commit_sha = release$commit_sha,
          status = status,
          stringsAsFactors = FALSE
        )
      )
      message(
        "Benchmark subprocess failed for ",
        release$release_version,
        " (",
        release$commit_sha,
        ") with status ",
        status
      )
      if (isTRUE(args$stop_on_error)) {
        quit(status = status)
      }
    }
  }

  if (nrow(failures)) {
    message("Release benchmark sweep finished with subprocess failures:")
    for (i in seq_len(nrow(failures))) {
      message(
        "  ",
        failures$release_version[[i]],
        " ",
        failures$commit_sha[[i]],
        " status=",
        failures$status[[i]]
      )
    }
    quit(status = 1)
  }

  cat("\nRelease benchmark sweep complete.\n")
  invisible(TRUE)
}

`%||%` = function(lhs, rhs) {
  if (is.null(lhs)) {
    rhs
  } else {
    lhs
  }
}

run_release_benchmarks()
