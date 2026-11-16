test_that("animation overlays follow each frame camera and match still occlusion", {
  motion = generate_camera_motion(
    list(c(0, 0, 6), c(3, 0, 6)),
    lookats = list(c(0, 0, 0), c(0, 0, 0)),
    fovs = 40,
    frames = 2,
    type = "linear",
    progress = FALSE
  )
  scene = sphere(radius = 0.5, material = diffuse("black"))
  labels = screen_text(
    "A",
    x = -0.8,
    y = 0.8,
    color = "green",
    size = 14,
    occlusion = TRUE,
    occlusion_mode = "partial"
  )
  lines = screen_line(
    start = matrix(c(-1, -0.8, 0), nrow = 1),
    end = matrix(c(1, -0.8, 0), nrow = 1),
    color = "red",
    width = 3,
    occlusion = TRUE,
    occlusion_mode = "partial"
  )
  args = list(
    scene = scene,
    width = 80,
    height = 64,
    samples = 1,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    integrator_type = "nee",
    backgroundhigh = "black",
    backgroundlow = "black",
    ambient_light = TRUE,
    screen_text = labels,
    screen_line = lines
  )
  frames = do.call(render_animation, c(args, list(camera_motion = motion)))
  expect_length(frames, 2)
  for (i in seq_along(frames)) {
    still = do.call(
      render_scene,
      c(
        args,
        list(
          lookfrom = unlist(motion[i, c("x", "y", "z")]),
          lookat = unlist(motion[i, c("dx", "dy", "dz")]),
          camera_up = unlist(motion[i, c("upx", "upy", "upz")]),
          fov = motion$fov[i],
          aperture = 0,
          iso = 1
        )
      )
    )
    expect_gt(sum(frames[[i]][,, 1] > 0.5), 10)
    expect_gt(sum(frames[[i]][,, 2] > 0.5), 10)
    expect_equal(as.numeric(frames[[i]]), as.numeric(still), tolerance = 1e-6)
  }
  expect_false(identical(as.numeric(frames[[1]]), as.numeric(frames[[2]])))

  # Scene defaults are how rayshader hands the extracted axis to render_movie_hq.
  saved = set_render_settings(scene, screen_text = labels, screen_line = lines)
  args$scene = saved
  args$screen_text = NULL
  args$screen_line = NULL
  inherited = do.call(
    render_animation,
    c(args, list(camera_motion = motion, start_frame = 2))
  )
  expect_equal(as.numeric(inherited[[1]]), as.numeric(frames[[2]]))
  cleared = do.call(
    render_animation,
    c(
      args,
      list(
        camera_motion = motion,
        start_frame = 2,
        screen_text = NULL,
        screen_line = NULL
      )
    )
  )
  expect_lt(max(cleared[[1]][,, 1:3]), 1e-6)

  # Compatible cameras share a native batch; incompatible metadata uses the
  # single-camera fallback. Both must inherit the same scene overlays.
  for (second_iso in c(100, 200)) {
    args$scene = saved |>
      add_camera(camera(motion = motion, name = "first")) |>
      add_camera(camera(motion = motion, name = "second", iso = second_iso))
    if (second_iso == 100) {
      batch = do.call(render_animation, c(args, list(camera = "all")))
    } else {
      expect_warning(
        batch <- do.call(render_animation, c(args, list(camera = "all"))),
        "camera metadata differs"
      )
    }
    expect_named(batch, c("first", "second"))
    for (camera_frames in batch) {
      for (i in seq_along(frames)) {
        expect_equal(as.numeric(camera_frames[[i]]), as.numeric(frames[[i]]))
      }
    }
  }
  args$scene = saved
  transparent = do.call(
    render_animation,
    c(
      args,
      list(
        camera_motion = motion,
        transparent_background = TRUE
      )
    )
  )
  expect_equal(dim(transparent[[1]])[3], 4)
  expect_gt(sum(transparent[[1]][,, 4] > 0), 10)
  expect_gt(sum(transparent[[1]][,, 1] > 0.5), 10)
  expect_gt(sum(transparent[[1]][,, 2] > 0.5), 10)

  # The background writer receives the same composited output as retained frames.
  directory = tempfile("screen-overlay-frames-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  files = do.call(
    render_animation,
    c(
      args,
      list(
        camera_motion = motion,
        filename = file.path(directory, "frame_%02d.png"),
        async_write = TRUE
      )
    )
  )
  expect_length(files, 2)
  for (i in seq_along(files)) {
    expect_true(file.exists(files[[i]]))
    expect_equal(
      as.numeric(png::readPNG(files[[i]])),
      as.numeric(frames[[i]]),
      tolerance = 1 / 255
    )
  }
})

test_that("animation labels and lines honor anchor and partial occlusion", {
  motion = generate_camera_motion(
    list(c(0, 0, 5), c(0, 0, -5)),
    lookats = list(c(0, 0, 0), c(0, 0, 0)),
    fovs = 40,
    frames = 2,
    type = "linear",
    progress = FALSE
  )
  args = list(
    scene = sphere(material = diffuse("black")),
    camera_motion = motion,
    width = 64,
    height = 64,
    samples = 1,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    integrator_type = "nee",
    ambient_light = TRUE,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  for (mode in c("anchor", "partial")) {
    args$screen_text = screen_text(
      "A",
      z = -1.5,
      color = "green",
      size = 12,
      occlusion = TRUE,
      occlusion_mode = mode
    )
    args$screen_line = screen_line(
      start = matrix(c(-0.3, -0.3, -1.5), nrow = 1),
      end = matrix(c(0.3, -0.3, -1.5), nrow = 1),
      color = "red",
      width = 3,
      occlusion = TRUE,
      occlusion_mode = mode
    )
    frames = do.call(render_animation, args)
    expect_lt(max(frames[[1]][,, 1:3]), 1e-6)
    expect_gt(sum(frames[[2]][,, 1] > 0.5), 10)
    expect_gt(sum(frames[[2]][,, 2] > 0.5), 10)
  }
  args$screen_text$occlusion = FALSE
  args$screen_line$occlusion = FALSE
  visible = do.call(render_animation, args)
  expect_gt(sum(visible[[1]][,, 1] > 0.5), 10)
  expect_gt(sum(visible[[1]][,, 2] > 0.5), 10)
})
