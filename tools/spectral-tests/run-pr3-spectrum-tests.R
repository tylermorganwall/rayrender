#!/usr/bin/env Rscript

find_script_path = function() {
  file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    stop("Unable to determine script path", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}

split_command = function(command) {
  parts = strsplit(trimws(command), "[[:space:]]+")[[1]]
  parts[nzchar(parts)]
}

config_value = function(key, fallback = "") {
  value = suppressWarnings(system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "config", key),
    stdout = TRUE
  ))
  value = trimws(paste(value, collapse = " "))
  if (!nzchar(value)) {
    fallback
  } else {
    value
  }
}

run_command = function(command, args, env = character()) {
  status = if (length(env) == 0) {
    system2(command, args)
  } else {
    system2(command, args, env = env)
  }
  if (!identical(status, 0L)) {
    stop(
      "Command failed: ",
      paste(c(command, args), collapse = " "),
      call. = FALSE
    )
  }
  invisible(status)
}

compile_and_run = function(
  label,
  repo_root,
  source,
  extra_flags = character(),
  run_env = character()
) {
  cxx_parts = split_command(config_value("CXX20", Sys.getenv("CXX20", "c++")))
  cxx = cxx_parts[[1]]
  cxx_prefix_args = cxx_parts[-1]

  cxx_flags = split_command(config_value("CXX20FLAGS", ""))
  cpp_flags = split_command(config_value("CPPFLAGS", ""))

  build_dir = tempfile(paste0("rayrender-pr3-", label, "-"))
  dir.create(build_dir)
  executable = file.path(build_dir, paste0("pr3-spectrum-tests-", label))

  args = c(
    cxx_prefix_args,
    cpp_flags,
    cxx_flags,
    "-std=gnu++20",
    "-Wall",
    "-Wextra",
    "-Wpedantic",
    "-I",
    repo_root,
    extra_flags,
    source,
    "-o",
    executable
  )

  message("Compiling PR3 sampled spectrum tests (", label, ")")
  run_command(cxx, args)
  message("Running PR3 sampled spectrum tests (", label, ")")
  run_command(executable, character(), env = run_env)
  invisible(TRUE)
}

script_path = find_script_path()
repo_root = normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)
source = file.path(
  repo_root,
  "tools",
  "spectral-tests",
  "pr3-spectrum-tests.cpp"
)

compile_and_run("normal", repo_root, source)

if (!isTRUE(as.logical(Sys.getenv("RAYRENDER_PR3_SKIP_SANITIZER", "false")))) {
  leak_detection = if (identical(Sys.info()[["sysname"]], "Darwin")) {
    "0"
  } else {
    "1"
  }
  compile_and_run(
    "sanitizer",
    repo_root,
    source,
    extra_flags = c("-fsanitize=address,undefined", "-fno-omit-frame-pointer"),
    run_env = c(
      paste0("ASAN_OPTIONS=detect_leaks=", leak_detection, ":abort_on_error=1"),
      "UBSAN_OPTIONS=halt_on_error=1:abort_on_error=1"
    )
  )
}
