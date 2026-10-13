test_that("Metal animations retain resources but reset independent frame state", {
  options_before = options(cores = 2L)
  on.exit(options(options_before))
  scene = sphere(material = diffuse("steelblue")) |>
    add_object(sphere(
      x = 1.6,
      radius = .6,
      material = light("orange", intensity = 2)
    ))
  motion = generate_camera_motion(
    list(c(0, 1, 5), c(2, 1, 5)),
    frames = 4,
    type = "linear",
    progress = FALSE
  )
  motion$fov = c(40, 0, 360, 40)
  motion$orthox = 5
  motion$orthoy = 4
  args = list(
    scene = scene,
    camera_motion = motion,
    width = 24,
    height = 20,
    samples = 8,
    integrator_type = "metal",
    preview = FALSE,
    denoise = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    bloom = FALSE,
    max_depth = 4,
    transparent_background = TRUE
  )
  set.seed(601)
  frames = suppressWarnings(do.call(render_animation, args))
  if (!isTRUE(attr(frames[[1]], "wavefront")$used)) {
    skip("Metal unavailable")
  }
  expect_length(frames, 4)
  for (i in seq_along(frames)) {
    report = attr(frames[[i]], "wavefront")
    expect_identical(report$scene_reused, i > 1L)
    expect_equal(report$completed_samples, 8)
    expect_equal(report$discarded_paths, 0)
    if (i > 1L) {
      expect_equal(report$upload_seconds, 0)
    }
    set.seed(601)
    separate = do.call(
      render_animation,
      c(args, list(start_frame = i, end_frame = i))
    )
    expect_equal(
      as.numeric(frames[[i]]),
      as.numeric(separate[[1]]),
      tolerance = 0
    )
  }
  # A second animation must start with its own resources, regardless of pointer
  # reuse by the native allocator after the preceding call is destroyed.
  args$scene = sphere(material = light("green"))
  set.seed(601)
  changed = do.call(render_animation, args)
  expect_false(attr(changed[[1]], "wavefront")$scene_reused)
  expect_gt(max(abs(as.numeric(changed[[1]]) - as.numeric(frames[[1]]))), .01)
})

test_that("retained Metal sessions handle lens cameras and motion", {
  scene = sphere(material = diffuse("coral")) |>
    animate_objects(end_position = c(.7, 0, 0))
  motion = generate_camera_motion(
    list(c(0, 1, 6), c(1, 1, 6)),
    frames = 3,
    type = "linear",
    progress = FALSE
  )
  args = list(
    scene = scene,
    camera_motion = motion,
    camera_description_file = "50mm",
    camera_motion_blur = TRUE,
    width = 20,
    height = 16,
    samples = 8,
    integrator_type = "metal",
    preview = FALSE,
    denoise = FALSE,
    parallel = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    bloom = FALSE,
    max_depth = 3
  )
  set.seed(602)
  frames = suppressWarnings(do.call(render_animation, args))
  if (!isTRUE(attr(frames[[1]], "wavefront")$used)) {
    skip("Metal unavailable")
  }
  for (i in seq_along(frames)) {
    expect_identical(attr(frames[[i]], "wavefront")$scene_reused, i > 1L)
    set.seed(602)
    separate = do.call(
      render_animation,
      c(args, list(start_frame = i, end_frame = i))
    )
    expect_equal(
      as.numeric(frames[[i]]),
      as.numeric(separate[[1]]),
      tolerance = 0
    )
  }
})

test_that("unsupported Metal animation sampling retains CPU fallback", {
  motion = generate_camera_motion(
    list(c(0, 1, 4), c(1, 1, 4)),
    frames = 2,
    type = "linear",
    progress = FALSE
  )
  args = list(
    scene = sphere(),
    camera_motion = motion,
    width = 8,
    height = 8,
    samples = c(2, 2),
    sample_method = "stratified",
    preview = FALSE,
    denoise = FALSE,
    parallel = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    bloom = FALSE
  )
  set.seed(603)
  frames = suppressWarnings(do.call(
    render_animation,
    c(args, list(integrator_type = "metal"))
  ))
  set.seed(603)
  reference = do.call(render_animation, c(args, list(integrator_type = "nee")))
  for (i in seq_along(frames)) {
    expect_false(attr(frames[[i]], "wavefront")$used)
    expect_false(attr(frames[[i]], "wavefront")$scene_reused)
    expect_equal(
      as.numeric(frames[[i]]),
      as.numeric(reference[[i]]),
      tolerance = 0
    )
  }
})
