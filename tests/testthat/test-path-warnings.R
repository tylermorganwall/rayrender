test_that("transport warnings preserve rendered images and log diagnostics", {
  withr::local_envvar(RAYRENDER_DEBUG_PATHS = "true")
  medium = homogeneous_medium(sigma_a = 0, sigma_s = 0)
  # Each box is closed, but their ordinary medium interiors are not nested.
  # Every camera ray enters both before exiting the first: a path-local failure.
  scene = set_medium(cube(z = .5, xwidth = 4, ywidth = 4), medium) |>
    add_object(set_medium(cube(z = 0, xwidth = 4, ywidth = 4), medium))
  render = function(scene, filename = NA_character_, parallel = FALSE) {
    render_scene(
      scene,
      width = 4,
      height = 4,
      samples = 1,
      lookfrom = c(0, 0, 4),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(1, 1),
      ambient_light = TRUE,
      filename = filename,
      parallel = parallel,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = "raw",
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE
    )
  }
  old = options(cores = 2)
  on.exit(options(old), add = TRUE)
  for (parallel in c(FALSE, TRUE)) {
    file = tempfile(fileext = ".png")
    on.exit(unlink(c(file, paste0(file, ".warnings.log"))), add = TRUE)
    warnings = list()
    result = withCallingHandlers(
      render(scene, file, parallel),
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    )
    expect_length(warnings, 1)
    expect_s3_class(warnings[[1]], "rayrender_path_warning")
    expect_true(file.exists(file))
    expect_true(all(is.finite(png::readPNG(file))))
    expect_true(all(is.finite(result)))
    diagnostics = attr(result, "path_warnings")
    expect_gt(diagnostics$terminated_paths, 0)
    expect_equal(sum(diagnostics$counts), diagnostics$terminated_paths)
    expect_gt(diagnostics$counts[["invalid_exit"]], 0)
    log = attr(result, "path_warning_log")
    expect_equal(log, normalizePath(paste0(file, ".warnings.log")))
    expect_match(paste(readLines(log), collapse = "\n"), "camera origin")
    expect_identical(warnings[[1]]$diagnostics, diagnostics)
  }
  # R's explicit warnings-as-errors policy must not lose the rendered file.
  file = tempfile(fileext = ".png")
  on.exit(unlink(c(file, paste0(file, ".warnings.log"))), add = TRUE)
  old_warn = options(warn = 2)
  expect_error(render(scene, file), "Terminated")
  options(old_warn)
  expect_true(file.exists(file))
  expect_true(file.exists(paste0(file, ".warnings.log")))
  clean = render(sphere())
  expect_null(attr(clean, "path_warnings"))
  expect_null(attr(clean, "path_warning_log"))
})

test_that("warning logs fall back to temp files and work for animation frames", {
  withr::local_envvar(RAYRENDER_DEBUG_PATHS = "1")
  diagnostics = list(
    terminated_paths = 2,
    counts = c(repeated_entry = 2),
    examples = list(repeated_entry = "example failure"),
    examples_per_kind = 8L
  )
  record = save_path_warnings(diagnostics, file.path(tempfile(), "image.png"))
  on.exit(unlink(record$log), add = TRUE)
  expect_true(file.exists(record$log))
  expect_null(save_path_warnings(list(terminated_paths = 0), NA_character_))
  rgb = list(r = matrix(.1, 4, 4), g = matrix(.2, 4, 4), b = matrix(.3, 4, 4))
  attr(rgb, "path_warnings") = diagnostics
  file = tempfile(fileext = ".png")
  on.exit(unlink(c(file, paste0(file, ".warnings.log"))), add = TRUE)
  expect_warning(
    frame <- post_process_frame(
      rgb,
      0,
      file,
      "raw",
      bloom = FALSE,
      plot_scene = FALSE
    ),
    "Terminated 2 path"
  )
  expect_true(file.exists(file))
  expect_identical(attr(frame, "path_warnings"), diagnostics)
  expect_true(file.exists(attr(frame, "path_warning_log")))
})

test_that("failed paths are silent and attached without files by default", {
  withr::local_envvar(RAYRENDER_DEBUG_PATHS = NA_character_)
  withr::local_options(warn = 2, cores = 2)
  medium = homogeneous_medium(sigma_a = 0, sigma_s = 0)
  scene = set_medium(cube(z = .5, xwidth = 4, ywidth = 4), medium) |>
    add_object(set_medium(cube(z = 0, xwidth = 4, ywidth = 4), medium))
  for (parallel in c(FALSE, TRUE)) {
    file = tempfile(fileext = ".png")
    on.exit(unlink(file), add = TRUE)
    expect_no_warning(
      result <- render_scene(
        scene,
        width = 4,
        height = 4,
        samples = 1,
        lookfrom = c(0, 0, 4),
        lookat = c(0, 0, 0),
        fov = 0,
        ortho_dimensions = c(1, 1),
        ambient_light = TRUE,
        filename = file,
        parallel = parallel,
        denoise = FALSE,
        bloom = FALSE,
        tonemap = "raw",
        preview = FALSE,
        plot_scene = FALSE,
        progress = FALSE
      )
    )
    expect_true(file.exists(file))
    expect_true(all(is.finite(result)))
    diagnostics = attr(result, "path_warnings")
    expect_gt(diagnostics$terminated_paths, 0)
    expect_equal(sum(diagnostics$counts), diagnostics$terminated_paths)
    expect_gt(length(diagnostics$examples$invalid_exit), 0)
    expect_null(attr(result, "path_warning_log"))
    expect_false(file.exists(paste0(file, ".warnings.log")))
  }
})

test_that("disabled debugging avoids disk writes for frames and unsaved images", {
  diagnostics = list(terminated_paths = 2, counts = c(invalid_exit = 2))
  rgb = list(r = matrix(.1, 4, 4), g = matrix(.2, 4, 4), b = matrix(.3, 4, 4))
  attr(rgb, "path_warnings") = diagnostics
  for (setting in c(NA_character_, "false", "0", "unrecognized")) {
    withr::local_envvar(RAYRENDER_DEBUG_PATHS = setting)
    before = list.files(tempdir(), all.files = TRUE)
    record = save_path_warnings(diagnostics, NA_character_)
    expect_identical(record$diagnostics, diagnostics)
    expect_null(record$log)
    expect_no_warning(warn_path_failures(record))
    expect_identical(list.files(tempdir(), all.files = TRUE), before)
  }
  file = tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  expect_no_warning(
    frame <- post_process_frame(
      rgb,
      0,
      file,
      "raw",
      bloom = FALSE,
      plot_scene = FALSE
    )
  )
  expect_identical(attr(frame, "path_warnings"), diagnostics)
  expect_null(attr(frame, "path_warning_log"))
  expect_false(file.exists(paste0(file, ".warnings.log")))
})
