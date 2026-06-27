#!/usr/bin/env Rscript

`%||%` = function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    y
  } else {
    x
  }
}

parse_args = function(args) {
  values = list()
  i = 1
  while (i <= length(args)) {
    key = args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    name = gsub("-", "_", sub("^--", "", key))
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      values[[name]] = TRUE
      i = i + 1
    } else {
      values[[name]] = args[[i + 1]]
      i = i + 2
    }
  }
  values
}

script_path = function() {
  file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
  }
  normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE)
}

repo_root = function() {
  normalizePath(file.path(dirname(script_path()), "..", ".."), mustWork = TRUE)
}

sha256_text = function(lines) {
  path = tempfile(fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(lines, path, useBytes = TRUE)
  unname(tools::sha256sum(path))
}

sha256_file = function(path) {
  unname(tools::sha256sum(path))
}

hash_array = function(value) {
  lines = c(
    paste0("type=", typeof(value)),
    paste0("dim=", paste(dim(value), collapse = "x")),
    format(as.numeric(value), digits = 17, scientific = TRUE)
  )
  sha256_text(lines)
}

scene_specs = function() {
  data.frame(
    scene = c(
      "cornell_diffuse_still",
      "material_mix_still",
      "dielectric_still",
      "orbit_animation_frames"
    ),
    type = c("still", "still", "still", "animation"),
    width = c(32L, 32L, 32L, 24L),
    height = c(32L, 32L, 32L, 24L),
    samples = c(4L, 4L, 4L, 2L),
    seed = c(1001L, 1002L, 1003L, 1004L),
    stringsAsFactors = FALSE
  )
}

render_common = function(
  scene,
  spec,
  lookfrom,
  lookat,
  fov = 35,
  clamp_value = 8
) {
  rayrender::render_scene(
    scene = scene,
    width = spec$width,
    height = spec$height,
    samples = spec$samples,
    lookfrom = lookfrom,
    lookat = lookat,
    fov = fov,
    aperture = 0,
    denoise = FALSE,
    preview = FALSE,
    interactive = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    verbose = FALSE,
    parallel = FALSE,
    ambient_light = FALSE,
    clamp_value = clamp_value,
    sample_method = "random",
    integrator_type = "nee",
    tonemap = "raw",
    bloom = FALSE,
    new_page = FALSE
  )
}

run_cornell_diffuse_still = function(spec, artifact_dir) {
  scene = rayrender::generate_cornell()
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = 278,
      y = 130,
      z = 278,
      radius = 90,
      material = rayrender::diffuse(color = "dodgerblue")
    )
  )
  image = render_common(
    scene,
    spec,
    lookfrom = c(278, 278, -800),
    lookat = c(278, 220, 278),
    fov = 40,
    clamp_value = 5
  )
  list(
    output_hash = hash_array(image),
    frame_hashes = NA_character_,
    image_sum = sum(image),
    render_options = "render_scene; cornell diffuse sphere; random; nee; denoise=FALSE; parallel=FALSE"
  )
}

run_material_mix_still = function(spec, artifact_dir) {
  scene = rayrender::generate_ground(
    depth = -0.55,
    material = rayrender::diffuse(color = "white", checkercolor = "grey65")
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = -0.9,
      y = 0,
      z = 0,
      radius = 0.45,
      material = rayrender::diffuse(color = "tomato")
    )
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = 0.15,
      y = 0,
      z = 0.1,
      radius = 0.45,
      material = rayrender::metal(color = "gold", fuzz = 0.08)
    )
  )
  scene = rayrender::add_object(
    scene,
    rayrender::cube(
      x = 1.1,
      y = -0.08,
      z = 0.15,
      width = 0.65,
      angle = c(0, 25, 0),
      material = rayrender::diffuse(color = "seagreen")
    )
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = 0,
      y = 5,
      z = -3,
      radius = 1,
      material = rayrender::light(intensity = 25, importance_sample = TRUE)
    )
  )
  image = render_common(
    scene,
    spec,
    lookfrom = c(3.5, 2.2, 6),
    lookat = c(0, 0, 0),
    fov = 35,
    clamp_value = 8
  )
  list(
    output_hash = hash_array(image),
    frame_hashes = NA_character_,
    image_sum = sum(image),
    render_options = "render_scene; mixed diffuse/metal/cube; random; nee; denoise=FALSE; parallel=FALSE"
  )
}

run_dielectric_still = function(spec, artifact_dir) {
  scene = rayrender::generate_ground(
    depth = -0.5,
    material = rayrender::diffuse(color = "white", checkercolor = "grey70")
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = -0.35,
      y = 0.2,
      z = 0,
      radius = 0.7,
      material = rayrender::dielectric(refraction = 1.45)
    )
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = 1.2,
      y = 0.05,
      z = 0.2,
      radius = 0.55,
      material = rayrender::diffuse(color = "purple")
    )
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = -3,
      y = 6,
      z = -4,
      radius = 1,
      material = rayrender::light(intensity = 35, importance_sample = TRUE)
    )
  )
  image = render_common(
    scene,
    spec,
    lookfrom = c(3.2, 1.5, 5.5),
    lookat = c(0, 0.1, 0),
    fov = 35,
    clamp_value = 10
  )
  list(
    output_hash = hash_array(image),
    frame_hashes = NA_character_,
    image_sum = sum(image),
    render_options = "render_scene; dielectric plus diffuse sphere; random; nee; denoise=FALSE; parallel=FALSE"
  )
}

run_orbit_animation_frames = function(spec, artifact_dir) {
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  }
  scene = rayrender::generate_ground(
    depth = -0.5,
    material = rayrender::diffuse(color = "white", checkercolor = "grey75")
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = -0.4,
      y = 0,
      z = 0,
      radius = 0.5,
      material = rayrender::diffuse(color = "steelblue")
    )
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = 0.8,
      y = 0,
      z = 0.2,
      radius = 0.45,
      material = rayrender::metal(color = "gold", fuzz = 0.1)
    )
  )
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = 0,
      y = 5,
      z = -2,
      radius = 1,
      material = rayrender::light(intensity = 25, importance_sample = TRUE)
    )
  )

  camera_motion = data.frame(
    x = c(3.5, 2.5),
    y = c(1.4, 1.6),
    z = c(5.5, 5.8),
    dx = c(0, 0),
    dy = c(0, 0),
    dz = c(0, 0),
    aperture = c(0, 0),
    fov = c(36, 36),
    focal = c(NA_real_, NA_real_),
    orthox = c(1, 1),
    orthoy = c(1, 1),
    upx = c(0, 0),
    upy = c(1, 1),
    upz = c(0, 0)
  )

  prefix = file.path(artifact_dir, "pr0-orbit-animation-")
  rayrender::render_animation(
    scene = scene,
    camera_motion = camera_motion,
    start_frame = 1,
    end_frame = 2,
    width = spec$width,
    height = spec$height,
    samples = spec$samples,
    filename = prefix,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    verbose = FALSE,
    parallel = FALSE,
    ambient_light = FALSE,
    clamp_value = 8,
    sample_method = "random",
    integrator_type = "nee",
    tonemap = "raw",
    bloom = FALSE
  )
  frames = paste0(prefix, 1:2, ".png")
  if (!all(file.exists(frames))) {
    stop(
      "Animation did not write expected frame files: ",
      paste(frames, collapse = ", ")
    )
  }
  frame_hashes = vapply(frames, sha256_file, character(1))
  list(
    output_hash = sha256_text(frame_hashes),
    frame_hashes = paste(frame_hashes, collapse = ";"),
    image_sum = NA_real_,
    render_options = "render_animation; two orbit frames; random; nee; denoise=FALSE; parallel=FALSE"
  )
}

run_scene = function(scene_name, spec, artifact_dir) {
  fun = switch(
    scene_name,
    "cornell_diffuse_still" = run_cornell_diffuse_still,
    "material_mix_still" = run_material_mix_still,
    "dielectric_still" = run_dielectric_still,
    "orbit_animation_frames" = run_orbit_animation_frames,
    stop("Unknown baseline scene: ", scene_name)
  )
  fun(spec, artifact_dir)
}

write_worker_row = function(args) {
  lib_path = Sys.getenv("R_LIBS_USER", unset = "")
  if (nzchar(lib_path)) {
    .libPaths(c(lib_path, .libPaths()))
  }
  suppressPackageStartupMessages(library(rayrender))

  specs = scene_specs()
  scene_name = args$scene
  spec = specs[specs$scene == scene_name, , drop = FALSE]
  if (nrow(spec) != 1) {
    stop("Unknown scene: ", scene_name)
  }

  artifact_dir = args$artifact_dir %||% tempfile("pr0-baseline-artifacts-")
  if (!dir.exists(artifact_dir)) {
    dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
  }
  set.seed(spec$seed)
  options(cores = 1, Ncpus = 1)

  timing = system.time({
    result = run_scene(scene_name, spec, artifact_dir)
  })

  row = data.frame(
    scene = scene_name,
    type = spec$type,
    width = spec$width,
    height = spec$height,
    samples = spec$samples,
    seed = spec$seed,
    output_hash = result$output_hash,
    frame_hashes = result$frame_hashes,
    image_sum = result$image_sum,
    user_cpu_seconds = unname(timing[["user.self"]]),
    system_cpu_seconds = unname(timing[["sys.self"]]),
    elapsed_seconds = unname(timing[["elapsed"]]),
    peak_rss_kb = NA_real_,
    render_options = result$render_options,
    stringsAsFactors = FALSE
  )
  write.csv(row, args$output, row.names = FALSE, na = "")
}

time_command_args = function(command, args) {
  sysname = Sys.info()[["sysname"]]
  if (file.exists("/usr/bin/time")) {
    if (identical(sysname, "Darwin")) {
      return(list(command = "/usr/bin/time", args = c("-l", command, args)))
    }
    return(list(command = "/usr/bin/time", args = c("-v", command, args)))
  }
  list(command = command, args = args)
}

parse_peak_rss_kb = function(path) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  lines = readLines(path, warn = FALSE)
  darwin = grep("maximum resident set size", lines, value = TRUE)
  if (length(darwin) > 0) {
    value = suppressWarnings(as.numeric(strsplit(
      trimws(darwin[[1]]),
      "[[:space:]]+"
    )[[1]][[1]]))
    return(value / 1024)
  }
  linux = grep("Maximum resident set size", lines, value = TRUE)
  if (length(linux) > 0) {
    value = suppressWarnings(as.numeric(sub(".*:[[:space:]]*", "", linux[[1]])))
    return(value)
  }
  NA_real_
}

run_worker_process = function(scene_name, artifact_dir) {
  row_path = tempfile(fileext = ".csv")
  time_log = tempfile(fileext = ".log")
  worker_args = c(
    script_path(),
    "--worker",
    "--scene",
    scene_name,
    "--output",
    row_path,
    "--artifact-dir",
    artifact_dir
  )
  command = file.path(R.home("bin"), "Rscript")
  timed = time_command_args(command, worker_args)
  status = system2(timed$command, timed$args, stdout = TRUE, stderr = time_log)
  status_code = attr(status, "status")
  if (is.null(status_code)) {
    status_code = 0
  }
  if (status_code != 0 && !file.exists(row_path)) {
    stop(
      "Baseline worker failed for ",
      scene_name,
      ":\n",
      paste(readLines(time_log, warn = FALSE), collapse = "\n")
    )
  }
  if (status_code != 0) {
    warning(
      "Timing wrapper returned status ",
      status_code,
      " for ",
      scene_name,
      "; continuing because the render row was written."
    )
  }
  row = read.csv(row_path, stringsAsFactors = FALSE)
  row$peak_rss_kb = parse_peak_rss_kb(time_log)
  row
}

write_environment = function(path, rows) {
  root = repo_root()
  r_version = R.version.string
  platform = R.version$platform
  sys = Sys.info()
  rayrender_version = as.character(utils::packageVersion("rayrender"))
  rayrender_commit = system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
  pbrt_commit = system2(
    "git",
    c("-C", file.path(root, "mmp", "pbrt-v4"), "rev-parse", "HEAD"),
    stdout = TRUE
  )
  cxx = system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "config", "CXX"),
    stdout = TRUE,
    stderr = TRUE
  )
  cxxflags = system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "config", "CXXFLAGS"),
    stdout = TRUE,
    stderr = TRUE
  )
  cppflags = system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "config", "CPPFLAGS"),
    stdout = TRUE,
    stderr = TRUE
  )

  lines = c(
    "# PR 0 Legacy Baseline Environment",
    "",
    paste0("- Captured at: `", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), "`"),
    paste0("- rayrender package version: `", rayrender_version, "`"),
    paste0("- rayrender commit: `", rayrender_commit, "`"),
    paste0("- pbrt-v4 commit: `", pbrt_commit, "`"),
    paste0("- R: `", r_version, "`"),
    paste0("- R platform: `", platform, "`"),
    paste0(
      "- System: `",
      sys[["sysname"]],
      " ",
      sys[["release"]],
      " ",
      sys[["machine"]],
      "`"
    ),
    paste0("- CXX: `", paste(cxx, collapse = " "), "`"),
    paste0("- CXXFLAGS: `", paste(cxxflags, collapse = " "), "`"),
    paste0("- CPPFLAGS: `", paste(cppflags, collapse = " "), "`"),
    "",
    "## Captured Scenes",
    "",
    paste0("- `", rows$scene, "`: ", rows$render_options)
  )
  writeLines(lines, path)
}

compare_baselines = function(current, baseline_path) {
  baseline = read.csv(baseline_path, stringsAsFactors = FALSE)
  stable_columns = c(
    "scene",
    "type",
    "width",
    "height",
    "samples",
    "seed",
    "output_hash",
    "frame_hashes",
    "render_options"
  )
  missing_columns = setdiff(stable_columns, names(baseline))
  if (length(missing_columns) > 0) {
    stop(
      "Baseline is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  current_key = current[stable_columns]
  baseline_key = baseline[stable_columns]
  for (name in stable_columns) {
    if (is.character(current_key[[name]])) {
      current_key[[name]][is.na(current_key[[name]])] = ""
    }
    if (is.character(baseline_key[[name]])) {
      baseline_key[[name]][is.na(baseline_key[[name]])] = ""
    }
  }
  current_key = current_key[order(current_key$scene), ]
  baseline_key = baseline_key[order(baseline_key$scene), ]
  row.names(current_key) = NULL
  row.names(baseline_key) = NULL
  if (!identical(current_key, baseline_key)) {
    print(current_key)
    print(baseline_key)
    stop("PR 0 legacy baseline hashes do not match recorded values.")
  }
  message("PR 0 legacy baseline hashes match recorded values.")
}

args = parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$worker)) {
  write_worker_row(args)
  quit(status = 0)
}

root = repo_root()
default_output = file.path(
  root,
  "docs",
  "spectral",
  "baselines",
  "pr0-legacy-baselines.csv"
)
default_environment = file.path(
  root,
  "docs",
  "spectral",
  "baselines",
  "pr0-legacy-environment.md"
)
output = args$output %||% default_output
environment_output = args$environment_output %||% default_environment
artifact_dir = args$artifact_dir %||%
  file.path(tempdir(), "rayrender-pr0-baselines")

if (!dir.exists(dirname(output))) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(dirname(environment_output))) {
  dir.create(
    dirname(environment_output),
    recursive = TRUE,
    showWarnings = FALSE
  )
}
if (!dir.exists(artifact_dir)) {
  dir.create(artifact_dir, recursive = TRUE, showWarnings = FALSE)
}

specs = scene_specs()
rows = do.call(
  rbind,
  lapply(specs$scene, run_worker_process, artifact_dir = artifact_dir)
)
write.csv(rows, output, row.names = FALSE, na = "")
write_environment(environment_output, rows)

if (!is.null(args$compare)) {
  compare_baselines(rows, normalizePath(args$compare, mustWork = TRUE))
}

message(
  "Wrote PR 0 legacy baselines to ",
  normalizePath(output, mustWork = TRUE)
)
