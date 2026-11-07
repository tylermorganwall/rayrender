test_that("async output bounds outstanding work and drains the final image", {
  directory = tempfile("async-frames-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  writer = animation_output_writer(TRUE, FALSE)
  on.exit(writer$close(), add = TRUE)
  pixels = list(
    r = matrix(.2, 16, 12),
    g = matrix(.3, 16, 12),
    b = matrix(.4, 16, 12)
  )
  attr(pixels, "volume_statistics") = list(paths = 12L)
  first = file.path(directory, "first.png")
  second = file.path(directory, "second.png")
  writer$process(pixels, 0, first, "raw", bloom = FALSE)
  writer$process(pixels, 0, second, "raw", bloom = FALSE)
  # Accepting the second image requires completion of the first.
  expect_true(file.exists(first))
  result = writer$finish(list(first, second))
  expect_true(file.exists(second))
  expect_identical(attr(result[[1]], "volume_statistics"), list(paths = 12L))
  expect_identical(attr(result[[2]], "volume_statistics"), list(paths = 12L))
  expect_lt(as.numeric(object.size(result)), 4096)
  writer$close()
  expect_silent(writer$close())
})

test_that("async output reports last-frame write failures", {
  writer = animation_output_writer(TRUE, FALSE)
  on.exit(writer$close(), add = TRUE)
  filename = file.path(tempfile("missing-directory-"), "failed.png")
  pixels = list(
    r = matrix(.2, 4, 4),
    g = matrix(.3, 4, 4),
    b = matrix(.4, 4, 4)
  )
  writer$process(pixels, 0, filename, "raw", bloom = FALSE)
  expect_error(writer$finish(list(filename)), "Image output failed.*failed.png")
})

test_that("async output preserves warnings, logs, and warn-as-error behavior", {
  withr::local_envvar(RAYRENDER_DEBUG_PATHS = "true")
  directory = tempfile("async-warnings-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  pixels = list(
    r = matrix(.2, 4, 4),
    g = matrix(.3, 4, 4),
    b = matrix(.4, 4, 4)
  )
  attr(pixels, "path_warnings") = list(terminated_paths = 1L)
  writer = animation_output_writer(TRUE, FALSE)
  on.exit(writer$close(), add = TRUE)
  filename = file.path(directory, "warning.png")
  writer$process(pixels, 0, filename, "raw", bloom = FALSE)
  expect_warning(result <- writer$finish(list(filename)), "Terminated 1 path")
  expect_true(file.exists(filename))
  expect_true(file.exists(attr(result[[1]], "path_warning_log")))
  expect_identical(
    attr(result[[1]], "path_warnings"),
    attr(pixels, "path_warnings")
  )
  withr::local_options(warn = 2)
  writer$process(pixels, 0, filename, "raw", bloom = FALSE)
  expect_error(writer$finish(list(filename, filename)), "Terminated 1 path")
})

test_that("async and synchronous file output match including alpha", {
  directory = tempfile("async-parity-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  motion = generate_camera_motion(
    list(c(0, 0, 5), c(.3, 0, 5)),
    frames = 2,
    type = "linear",
    progress = FALSE
  )
  args = list(
    scene = sphere(material = diffuse("red")),
    camera_motion = motion,
    width = 16,
    height = 12,
    samples = 2,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    parallel = FALSE,
    bloom = TRUE,
    transparent_background = TRUE,
    integrator_type = "nee"
  )
  set.seed(72)
  sync = do.call(
    render_animation,
    c(
      args,
      list(
        async_write = FALSE,
        filename = file.path(directory, "sync_%02d.png")
      )
    )
  )
  set.seed(72)
  async = do.call(
    render_animation,
    c(
      args,
      list(
        async_write = TRUE,
        filename = file.path(directory, "async_%02d.png")
      )
    )
  )
  for (i in 1:2) {
    expect_identical(
      unname(tools::md5sum(sync[[i]])),
      unname(tools::md5sum(async[[i]]))
    )
    expect_null(attr(async[[i]], "rayrender_async_pending"))
  }
  expect_error(render_animation(NULL, async_write = NA), "async_write must be")
})

test_that("closing async output stops pending work and removes transport files", {
  directory = tempfile("async-close-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  writer = animation_output_writer(TRUE, FALSE)
  on.exit(writer$close(), add = TRUE)
  pixels = list(
    r = matrix(.2, 16, 12),
    g = matrix(.3, 16, 12),
    b = matrix(.4, 16, 12)
  )
  writer$process(
    pixels,
    0,
    file.path(directory, "frame.png"),
    "raw",
    bloom = FALSE
  )
  worker = environment(writer$process)$worker
  # Inspect transport ownership only here to verify cleanup without waiting for GC.
  transport = worker$.__enclos_env__$private$options$tmp_files
  expect_true(worker$is_alive())
  expect_true(length(transport) > 0L)
  writer$close()
  expect_false(worker$is_alive())
  expect_false(any(file.exists(transport)))
  expect_silent(writer$close())
})
