test_that("file animations stream pixels unless retention is requested", {
  motion = generate_camera_motion(
    list(c(0, 0, 5), c(0.2, 0, 5)),
    frames = 3,
    type = "linear",
    progress = FALSE
  )
  output_dir = tempfile("animation-stream-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE))
  args = list(
    scene = sphere(material = diffuse("red")),
    camera_motion = motion,
    width = 20,
    height = 16,
    samples = 1,
    parallel = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    bloom = FALSE
  )

  # Cover each native output branch and numerical debug post-processing.
  cases = list(
    list(debug_channel = "none", denoise = FALSE, integrator_type = "nee"),
    list(debug_channel = "none", denoise = FALSE, integrator_type = "auto"),
    list(debug_channel = "normals", denoise = FALSE),
    list(debug_channel = "position", denoise = FALSE)
  )
  if (has_denoiser()) {
    cases = c(cases, list(list(debug_channel = "none", denoise = TRUE)))
  }
  for (settings in cases) {
    args = utils::modifyList(args, settings)
    streamed_pattern = file.path(output_dir, "stream_%03d.png")
    retained_pattern = file.path(output_dir, "retain_%03d.png")
    set.seed(15)
    streamed = do.call(
      render_animation,
      c(args, list(filename = streamed_pattern))
    )
    set.seed(15)
    retained = do.call(
      render_animation,
      c(
        args,
        list(
          filename = retained_pattern,
          return_frames = TRUE
        )
      )
    )
    expect_length(streamed, 3L)
    for (i in seq_along(streamed)) {
      expect_identical(
        as.character(streamed[[i]]),
        sprintf(streamed_pattern, i)
      )
      expect_true(is.array(retained[[i]]))
      expect_identical(
        unname(tools::md5sum(sprintf(streamed_pattern, i))),
        unname(tools::md5sum(sprintf(retained_pattern, i)))
      )
    }
    expect_lt(
      as.numeric(object.size(streamed)),
      as.numeric(object.size(retained))
    )
  }

  # Omitting a filename still returns images; explicit FALSE drops them.
  args$denoise = FALSE
  args$debug_channel = "none"
  frames = do.call(render_animation, args)
  expect_true(all(vapply(frames, is.array, logical(1))))
  discarded = do.call(render_animation, c(args, list(return_frames = FALSE)))
  expect_true(all(vapply(
    discarded,
    function(x) identical(as.character(x), NA_character_),
    logical(1)
  )))

  # A subset reports original frame numbers, not positions within the subset.
  subset = do.call(
    render_animation,
    c(
      args,
      list(
        filename = file.path(output_dir, "subset_%03d.png"),
        start_frame = 2,
        end_frame = 3
      )
    )
  )
  expect_identical(
    unlist(subset),
    file.path(output_dir, c("subset_002.png", "subset_003.png"))
  )
})

test_that("streamed frames preserve diagnostics without image attributes", {
  original = post_process_frame
  diagnostics = list(
    path_warnings = list(terminated_paths = 1L),
    path_warning_log = "test.warnings.log",
    volume_statistics = list(paths = 12L),
    wavefront = list(used = FALSE),
    animation_timing = list(compile_seconds = 1)
  )
  local_mocked_bindings(
    post_process_frame = function(...) {
      frame = original(...)
      for (name in names(diagnostics)) {
        attr(frame, name) = diagnostics[[name]]
      }
      frame
    },
    .package = "rayrender"
  )
  filename = tempfile(fileext = ".png")
  on.exit(unlink(filename))
  frames = render_animation(
    sphere(),
    camera = camera(lookfrom = c(0, 0, 5)),
    filename = filename,
    integrator_type = "nee",
    async_write = FALSE,
    width = 8,
    height = 8,
    samples = 1,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    parallel = FALSE,
    bloom = FALSE
  )
  expect_type(frames[[1]], "character")
  expect_null(dim(frames[[1]]))
  expect_identical(attributes(frames[[1]]), diagnostics)
})

test_that("mixed camera batches resolve retention per output file", {
  output_dir = tempfile("animation-cameras-")
  dir.create(output_dir)
  on.exit(unlink(output_dir, recursive = TRUE))
  motion = generate_camera_motion(
    list(c(0, 0, 5), c(0.2, 0, 5)),
    frames = 2,
    type = "linear",
    progress = FALSE
  )
  for (second_iso in c(100, 200)) {
    # Equal metadata uses one native batch; differing ISO uses the fallback.
    scene = sphere() |>
      add_camera(
        camera(
          motion = motion,
          filename = file.path(output_dir, "saved_%02d.png")
        ),
        name = "saved"
      ) |>
      add_camera(camera(motion = motion, iso = second_iso), name = "memory")
    args = list(
      scene = scene,
      camera = "all",
      width = 8,
      height = 8,
      samples = 1,
      denoise = FALSE,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      parallel = FALSE,
      bloom = FALSE
    )
    if (second_iso == 100) {
      output = do.call(render_animation, args)
    } else {
      expect_warning(output <- do.call(render_animation, args), "one at a time")
    }
    expect_true(all(vapply(output$saved, is.character, logical(1))))
    expect_true(all(vapply(output$memory, is.array, logical(1))))
  }
})

test_that("invalid retention policies fail before scene construction", {
  for (value in list(NA, logical(), c(TRUE, FALSE), 1, "true")) {
    expect_error(
      render_animation(NULL, return_frames = value),
      "return_frames must be"
    )
  }
})
