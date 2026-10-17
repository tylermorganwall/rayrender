test_that("camera timestamps sample without seek history and reject ambiguous exposure", {
  motion = generate_camera_motion(
    positions = rbind(c(0, 0, 8), c(2, 1, 8)),
    lookats = rbind(c(0, 0, 0), c(0, 0, 0)),
    frames = 4
  )
  cam = camera(motion = motion)
  sampled = camera_at_scene_times(cam, c(.1, 0, .05, .1), 30, c(0, 0))
  expect_equal(sampled$motion$x[c(1, 4)], rep(tail(cam$motion$x, 1), 2))
  expect_equal(sampled$motion$x[2], cam$motion$x[1])
  expect_equal(
    camera_at_scene_times(cam, 0, 30, c(0, .1))$shutter_offsets,
    c(0, .1)
  )
  expect_error(
    camera_at_scene_times(cam, c(0, 1e12, 1e12 + .25), 30, c(0, 0)),
    "precision"
  )
  expect_error(render_animation(sphere(), times = numeric()), "nonempty")
  expect_error(render_animation(sphere(), times = c(0, NA)), "finite")
  expect_error(render_animation(sphere(), times = 0, shutteropen = 0), "legacy")
  expect_error(
    resolve_scene_timing(0, NULL, 30, c(.1, .1 + 1e-10), FALSE),
    "precision"
  )
  expect_error(
    camera_at_scene_times(camera(), c(0, 100), 30, c(0, 1e-7)),
    "collapses"
  )
})

test_that("explicit offline frames and stills agree at nonmonotonic requested times", {
  withr::local_options(list(cores = 2L))
  path = transform_track(
    c(-1, .5, 2, 4),
    position = rbind(c(-1, 0, 0), c(0, 1, 0), c(1, 0, 0), c(-1, 0, 0)),
    angle = rbind(c(0, 0, 0), c(0, 30, 0), c(0, 80, 0), c(0, 0, 0)),
    position_interpolation = "cubic",
    rotation_interpolation = "smooth"
  )
  scene = animate_keyframes(cube(width = 1.2), path)
  cam = camera(
    lookfrom = c(0, 0, 8),
    lookat = c(0, 0, 0),
    aperture = 0,
    fov = 0,
    ortho_dimensions = c(6, 6)
  )
  scene = add_camera(scene, cam)
  args = list(
    scene = scene,
    width = 20,
    height = 20,
    samples = 1,
    max_depth = 1,
    integrator_type = "nee",
    debug_channel = "position",
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  times = c(2, .25, -1, 4, .25)
  frames = do.call(render_animation, c(args, list(times = times)))
  expect_length(frames, length(times))
  compare_geometry = function(a, b) {
    # Existing offline outputs include alpha/metadata; background debug hits
    # lie on the enclosing environment sphere, whose radius covers the sequence.
    # Compare the actual object's hit mask and numerical world-space positions.
    mask = abs(b[,, 3]) < 2
    expect_gt(sum(mask), 0)
    expect_identical(abs(a[,, 3]) < 2, mask)
    for (k in 1:3) {
      expect_equal(
        as.vector(a[,, k][mask]),
        as.vector(b[,, k][mask]),
        tolerance = 2e-5
      )
    }
  }
  for (i in seq_along(times)) {
    still = do.call(render_scene, c(args, list(time = times[i])))
    compare_geometry(frames[[i]], still)
  }
  expect_identical(frames[[2]], frames[[5]])
  subset = do.call(
    render_animation,
    c(args, list(times = times, start_frame = 2, end_frame = 3))
  )
  expect_equal(subset, frames[2:3])
  # Scene-time batches preserve the existing per-camera return structure.
  args$scene = add_camera(scene, cam, name = "second")
  batch = do.call(render_scene, c(args, list(camera = "all", time = .25)))
  compare_geometry(batch[[1]], frames[[2]])
  expect_equal(batch[[1]], batch[[2]])
  animation_batch = do.call(
    render_animation,
    c(args, list(camera = "all", times = times))
  )
  expect_length(animation_batch, 2L)
  expect_equal(animation_batch[[1]], animation_batch[[2]])
  compare_geometry(animation_batch[[1]][[2]], frames[[2]])

  args$scene = scene
  exposure_times = c(.7, 1.3, .7)
  exposed = do.call(
    render_animation,
    c(args, list(times = exposure_times, shutter_offsets = c(-.3, .3)))
  )
  for (i in seq_along(exposure_times)) {
    still = do.call(
      render_scene,
      c(args, list(time = exposure_times[i], shutter_offsets = c(-.3, .3)))
    )
    compare_geometry(exposed[[i]], still)
  }
  prefix = file.path(withr::local_tempdir(), "frame")
  do.call(
    render_animation,
    c(
      args,
      list(times = times, start_frame = 2, end_frame = 3, filename = prefix)
    )
  )
  expect_identical(
    sort(basename(list.files(dirname(prefix), full.names = TRUE))),
    c("frame2.png", "frame3.png")
  )
})

test_that("closed motion exposures include interior keys and secondary-ray timing", {
  withr::local_options(list(cores = 2L))
  scene = animate_keyframes(
    cube(xwidth = .4, ywidth = 4, zwidth = .1, material = light(intensity = 1)),
    transform_track(c(0, 1, 2, 3), position = cbind(c(-1, 1, -1, -1), 0, 0))
  )
  args = list(
    scene = scene,
    width = 32,
    height = 24,
    samples = 64,
    time = 0,
    shutter_offsets = c(0, 2),
    lookfrom = c(0, 0, 6),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(4, 4),
    aperture = 0,
    integrator_type = "nee",
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE
  )
  set.seed(122)
  image = do.call(render_scene, args)
  expect_gt(mean(image[12:20, 8:16, 1]), .05)
  expect_lt(mean(image[12:20, 8:16, 1]), .5)
  args$integrator_type = "metal"
  set.seed(122)
  device = suppressWarnings(do.call(render_scene, args))
  if (!isTRUE(attr(device, "wavefront")$used)) {
    skip("Metal runtime unavailable; host packing is tested natively")
  }
  expect_equal(mean(device), mean(image), tolerance = .025)
})

test_that("explicit Metal output times reuse geometry and reset each frame", {
  withr::local_options(list(cores = 2L))
  scene = animate_keyframes(
    sphere(radius = .6, material = light(intensity = 1)),
    transform_track(
      c(0, .7, 2.3, 4),
      position = rbind(c(-1, 0, 0), c(0, 1, 0), c(1, 0, 0), c(-1, 0, 0)),
      position_interpolation = "cubic"
    )
  ) |>
    add_camera(camera(
      lookfrom = c(0, 0, 6),
      lookat = c(0, 0, 0),
      aperture = 0,
      fov = 0,
      ortho_dimensions = c(4, 4)
    ))
  args = list(
    scene = scene,
    times = c(2, .25, -1, .25),
    width = 20,
    height = 20,
    samples = 32,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  set.seed(509)
  gpu = suppressWarnings(do.call(
    render_animation,
    c(args, list(integrator_type = "metal"))
  ))
  if (!isTRUE(attr(gpu[[1]], "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  set.seed(509)
  cpu = do.call(render_animation, c(args, list(integrator_type = "nee")))
  for (i in seq_along(gpu)) {
    report = attr(gpu[[i]], "wavefront")
    expect_true(report$used)
    if (i > 1) {
      expect_true(report$scene_reused)
    }
    expect_lt(mean(abs(as.numeric(gpu[[i]]) - as.numeric(cpu[[i]]))), .015)
  }
  expect_identical(as.numeric(gpu[[2]]), as.numeric(gpu[[4]]))
})
