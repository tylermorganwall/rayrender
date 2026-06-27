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

compile_test = function(label, repo_root, extra_flags = character()) {
  cxx_parts = split_command(config_value("CXX20", Sys.getenv("CXX20", "c++")))
  cxx = cxx_parts[[1]]
  cxx_prefix_args = cxx_parts[-1]

  cxx_flags = split_command(config_value("CXX20FLAGS", ""))
  cpp_flags = split_command(config_value("CPPFLAGS", ""))

  build_dir = tempfile(paste0("rayrender-pr19-", label, "-"))
  dir.create(build_dir)
  executable = file.path(build_dir, paste0("pr19-material-tests-", label))

  sources = file.path(
    repo_root,
    c(
      "tools/spectral-tests/pr19-material-tests.cpp",
      "src/materials/spectral_material.cpp",
      "src/materials/spectral_texture.cpp",
      "src/render/spectral_bsdf.cpp",
      "src/render/spectral_dielectric.cpp",
      "src/render/spectral_scene.cpp",
      "src/math/perlin.cpp",
      "src/math/rng.cpp"
    )
  )
  dependency_includes = c(
    file.path(repo_root, "src"),
    system.file("include", package = "spacefillr")
  )
  dependency_includes = dependency_includes[nzchar(dependency_includes)]
  dependency_include_flags = unlist(lapply(dependency_includes, function(path) {
    c("-I", path)
  }))

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
    dependency_include_flags,
    extra_flags,
    sources,
    "-o",
    executable
  )

  message("Compiling PR19 Material tests (", label, ")")
  run_command(cxx, args)
  executable
}

run_test = function(label, executable, extra_env = character()) {
  message("Running PR19 Material tests (", label, ")")
  run_command(executable, character(), env = extra_env)
}

script_path = find_script_path()
repo_root = normalizePath(
  file.path(dirname(script_path), "..", ".."),
  mustWork = TRUE
)

normal_executable = compile_test("normal", repo_root)
run_test("normal", normal_executable)

if (!isTRUE(as.logical(Sys.getenv("RAYRENDER_PR19_SKIP_SANITIZER", "false")))) {
  leak_detection = if (identical(Sys.info()[["sysname"]], "Darwin")) {
    "0"
  } else {
    "1"
  }
  sanitizer_executable = compile_test(
    "sanitizer",
    repo_root,
    extra_flags = c("-fsanitize=address,undefined", "-fno-omit-frame-pointer")
  )
  run_test(
    "sanitizer",
    sanitizer_executable,
    extra_env = c(
      paste0("ASAN_OPTIONS=detect_leaks=", leak_detection, ":abort_on_error=1"),
      "UBSAN_OPTIONS=halt_on_error=1:abort_on_error=1"
    )
  )
}
