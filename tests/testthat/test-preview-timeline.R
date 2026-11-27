test_that("a headless timeline snapshot agrees with the ordinary explicit still", {
  withr::local_options(list(cores = 2L))
  scene = animate_keyframes(
    sphere(radius = .4, material = light(intensity = 1)),
    transform_track(
      c(0, .7, 2.3, 4),
      position = rbind(c(-1, 0, 0), c(0, 1, 0), c(1, 0, 0), c(-1, 0, 0)),
      position_interpolation = "cubic"
    )
  )
  args = list(
    scene = scene,
    time = 1.5,
    width = 16,
    height = 16,
    samples = 4,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    integrator_type = "nee",
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    lookfrom = c(0, 0, 6),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(4, 4),
    aperture = 0
  )
  set.seed(27)
  still = do.call(render_scene, args)
  set.seed(27)
  bounded = do.call(render_scene, c(args, list(timeline = c(0, 4))))
  expect_identical(as.numeric(still), as.numeric(bounded))
  expect_equal(attr(bounded, "timeline")$published_time, 1.5)
  expect_equal(attr(bounded, "timeline")$displayed_time, 1.5)
  expect_error(render_scene(scene, time = 1, timeline = c(2, 3)), "containing")
  expect_error(
    render_scene(scene, time = 1, timeline = c(0, 3), mode = "animation"),
    "offline"
  )
})

test_that("authored cameras share scene-time sampling in preview and stills", {
  withr::local_options(list(cores = 2L))
  motion = generate_camera_motion(
    positions = rbind(c(-1, 0, 6), c(1, 1, 6)),
    lookats = rbind(c(0, 0, 0), c(0, 0, 0)),
    frames = 3,
    fovs = 40,
    apertures = 0
  )
  motion$time = c(0, .5, 4)
  scene = add_camera(
    sphere(material = light(intensity = 1)),
    camera(motion = motion)
  )
  args = list(
    scene = scene,
    width = 16,
    height = 16,
    samples = 4,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    integrator_type = "nee",
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  for (time in c(3, 0, 1, 4, 1)) {
    set.seed(28)
    still = do.call(render_scene, c(args, list(time = time)))
    set.seed(28)
    bounded = do.call(
      render_scene,
      c(args, list(time = time, timeline = c(0, 4)))
    )
    expect_equal(as.numeric(still), as.numeric(bounded), tolerance = 1e-6)
  }
  set.seed(29)
  exposed = do.call(
    render_scene,
    c(args, list(time = 0, shutter_offsets = c(0, .1)))
  )
  set.seed(29)
  bounded = do.call(
    render_scene,
    c(args, list(time = 0, timeline = c(0, 4), shutter_offsets = c(0, .1)))
  )
  expect_equal(as.numeric(exposed), as.numeric(bounded), tolerance = 1e-6)
})
