test_that("lens selection is retained in authored physical-camera tracks", {
  motion = generate_camera_motion(rbind(c(0, 0, 6), c(1, 0, 6)), frames = 3)
  physical = camera(motion = motion, camera_description_file = "50mm")
  sampled = camera_at_scene_times(physical, c(0, .05), 30, c(0, 0))
  expect_true(all(sampled$motion$fov == -1))
  expect_true(all(sampled$authored_motion$fov == -1))
  expect_error(
    camera_at_scene_times(physical, 0, 30, c(0, .1)),
    "Physical-camera"
  )
  expect_no_error(camera_at_scene_times(physical, 0, 30, c(.01, .01)))
})

test_that("zero-width optics use exposure time without changing the scene clock", {
  motion = generate_camera_motion(
    rbind(c(0, 0, 6), c(1, 0, 6)),
    frames = 2,
    apertures = c(4, 10),
    focal_distances = c(6, 12)
  )
  motion$time = c(2, 4)
  authored = camera(motion = motion, camera_description_file = "50mm")
  sampled = camera_at_scene_times(authored, c(3, 2, 3), 30, c(.5, .5))
  expect_equal(sampled$motion$aperture, c(8.5, 5.5, 8.5))
  expect_equal(sampled$motion$focal, c(10.5, 7.5, 10.5))
  expect_equal(sampled$motion$time, c(3, 2, 3))
  expect_equal(sampled$motion$scene_offset, c(1, 0, 1))
  expect_equal(sampled$scene_time, 2)
  expect_identical(sampled$authored_motion$aperture, motion$aperture)
  expect_error(
    camera_at_scene_times(authored, 3, 30, c(0, .1)),
    "constant aperture"
  )
})

test_that("variable physical snapshots match timeline, still and offline paths", {
  withr::local_options(list(cores = 2L))
  motion = generate_camera_motion(
    rbind(c(0, 0, 6), c(0, 0, 6)),
    lookats = rbind(c(0, 0, 0), c(0, 0, 0)),
    frames = 2,
    apertures = c(4, 10),
    focal_distances = c(6, 12)
  )
  motion$time = c(2, 9)
  settings = list(
    scene = add_camera(
      sphere(radius = .6, material = light(color = "#356F95", intensity = 1)),
      camera(motion = motion, camera_description_file = "50mm", iso = 100)
    ),
    width = 32,
    height = 24,
    samples = 64,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    integrator_type = "nee"
  )
  set.seed(119)
  still = do.call(render_scene, c(settings, list(time = 5.5)))
  set.seed(119)
  timeline = do.call(
    render_scene,
    c(settings, list(time = 5.5, timeline = c(2, 9)))
  )
  expect_gt(mean(still[,, 1:3]), .001)
  expect_equal(as.numeric(timeline), as.numeric(still), tolerance = 1e-7)
  set.seed(119)
  shifted = do.call(
    render_scene,
    c(settings, list(time = 5, timeline = c(2, 9), shutter_offsets = c(.5, .5)))
  )
  expect_equal(as.numeric(shifted), as.numeric(still), tolerance = 1e-7)
  times = c(5.5, 9, 2, 5.5)
  set.seed(119)
  frames = do.call(render_animation, c(settings, list(times = times)))
  expect_identical(as.numeric(frames[[1]]), as.numeric(frames[[4]]))
  expect_lt(mean(abs(frames[[1]][,, 1:3] - still[,, 1:3])), .005)
  expect_gt(mean(abs(frames[[2]][,, 1:3] - frames[[3]][,, 1:3])), .001)
  set.seed(119)
  shifted_frames = do.call(
    render_animation,
    c(settings, list(times = times - .5, shutter_offsets = c(.5, .5)))
  )
  for (i in seq_along(times)) {
    expect_equal(
      as.numeric(shifted_frames[[i]]),
      as.numeric(frames[[i]]),
      tolerance = 1e-7
    )
  }
  expect_error(
    do.call(
      render_scene,
      c(
        settings,
        list(time = 5.5, timeline = c(2, 9), shutter_offsets = c(-.1, .1))
      )
    ),
    "constant aperture and focus"
  )
  settings$integrator_type = "metal"
  set.seed(119)
  device = suppressWarnings(do.call(
    render_scene,
    c(settings, list(time = 5.5, timeline = c(2, 9)))
  ))
  if (!isTRUE(attr(device, "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  expect_lt(mean(abs(device[,, 1:3] - still[,, 1:3])), .005)
  set.seed(119)
  device_frames = do.call(render_animation, c(settings, list(times = times)))
  expect_true(all(vapply(
    device_frames,
    function(x) isTRUE(attr(x, "wavefront")$used),
    logical(1)
  )))
  expect_identical(
    as.numeric(device_frames[[1]]),
    as.numeric(device_frames[[4]])
  )
  for (i in seq_along(times)) {
    expect_lt(mean(abs(device_frames[[i]][,, 1:3] - frames[[i]][,, 1:3])), .005)
  }
})

test_that("fixed physical optics allow authored exposures and preserve random seeks", {
  withr::local_options(list(cores = 2L))
  motion = generate_camera_motion(
    rbind(c(-1, 0, 6), c(1, .5, 6), c(-1, 0, 6)),
    lookats = rbind(c(-1, 0, 0), c(1, .5, 0), c(-1, 0, 0)),
    frames = 3,
    type = "linear",
    apertures = 4,
    focal_distances = 6
  )
  motion$time = c(2, 2.2, 2.4)
  authored = camera(
    motion = motion,
    camera_description_file = "50mm",
    iso = 100
  )
  expect_no_error(camera_at_scene_times(authored, 2, 30, c(0, .4)))
  settings = list(
    scene = add_camera(
      sphere(radius = .4, material = light(intensity = 1)),
      authored
    ),
    width = 32,
    height = 24,
    samples = 64,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    integrator_type = "nee"
  )
  set.seed(117)
  frozen = do.call(render_scene, c(settings, list(time = 2)))
  set.seed(117)
  exposed = do.call(
    render_scene,
    c(settings, list(time = 2, shutter_offsets = c(0, .4)))
  )
  expect_gt(sum(exposed[,, 1] > .002), sum(frozen[,, 1] > .002))
  set.seed(117)
  timeline = do.call(
    render_scene,
    c(settings, list(time = 2, timeline = c(2, 9), shutter_offsets = c(0, .4)))
  )
  expect_equal(as.numeric(timeline), as.numeric(exposed), tolerance = 1e-7)
  set.seed(117)
  frames = do.call(
    render_animation,
    c(settings, list(times = c(2, 2.2, 2.1, 2), shutter_offsets = c(0, .4)))
  )
  expect_identical(as.numeric(frames[[1]]), as.numeric(frames[[4]]))
  expect_lt(mean(abs(frames[[1]][,, 1:3] - exposed[,, 1:3])), .005)
  settings$integrator_type = "metal"
  set.seed(117)
  device = suppressWarnings(do.call(
    render_scene,
    c(settings, list(time = 2, timeline = c(2, 9), shutter_offsets = c(0, .4)))
  ))
  if (!isTRUE(attr(device, "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  expect_lt(mean(abs(device[,, 1:3] - exposed[,, 1:3])), .005)
  set.seed(117)
  device_frames = do.call(
    render_animation,
    c(settings, list(times = c(2, 2.2, 2.1, 2), shutter_offsets = c(0, .4)))
  )
  expect_true(all(vapply(
    device_frames,
    function(x) isTRUE(attr(x, "wavefront")$used),
    logical(1)
  )))
  expect_identical(
    as.numeric(device_frames[[1]]),
    as.numeric(device_frames[[4]])
  )
  expect_lt(mean(abs(device_frames[[1]][,, 1:3] - frames[[1]][,, 1:3])), .005)
})

test_that("explicit physical-camera snapshots retain brightness independently of pose time", {
  withr::local_options(list(cores = 2L))
  args = list(
    scene = add_camera(
      sphere(material = light(intensity = 1)),
      camera(
        lookfrom = c(0, 0, 6),
        lookat = c(0, 0, 0),
        camera_description_file = "50mm",
        aperture = 4,
        focal_distance = 6,
        iso = 100
      )
    ),
    width = 16,
    height = 16,
    samples = 32,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    integrator_type = "nee"
  )
  set.seed(115)
  legacy = do.call(render_scene, args)
  set.seed(115)
  frozen = do.call(render_scene, c(args, list(time = 2)))
  expect_gt(mean(frozen[,, 1:3]), .001)
  expect_equal(as.numeric(frozen), as.numeric(legacy), tolerance = 1e-7)
  set.seed(115)
  exposed = do.call(
    render_scene,
    c(args, list(time = 2, shutter_offsets = c(-.1, .2)))
  )
  expect_equal(as.numeric(exposed), as.numeric(frozen), tolerance = 1e-7)
  set.seed(115)
  frames = do.call(render_animation, c(args, list(times = c(2, 1, 2))))
  # Still and offline entry points initialize different Monte Carlo streams.
  # Their images agree within sampling error; repeated offline times are exact.
  expect_lt(mean(abs(frames[[1]][,, 1:3] - frozen[,, 1:3])), .005)
  expect_identical(as.numeric(frames[[1]]), as.numeric(frames[[3]]))
  args$integrator_type = "metal"
  set.seed(115)
  device = suppressWarnings(do.call(render_scene, c(args, list(time = 2))))
  if (!isTRUE(attr(device, "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  expect_gt(mean(device[,, 1:3]), .001)
  expect_lt(mean(abs(device[,, 1:3] - frozen[,, 1:3])), .005)
})

test_that("authored camera exposure visits interior keys on CPU and Metal", {
  withr::local_options(list(cores = 2L))
  motion = generate_camera_motion(
    positions = rbind(c(-1, 0, 6), c(1, 0, 6), c(-1, 0, 6)),
    lookats = rbind(c(-1, 0, 0), c(1, 0, 0), c(-1, 0, 0)),
    frames = 3,
    type = "linear",
    damp_motion = FALSE,
    fovs = 0,
    ortho_dims = matrix(c(4, 4), nrow = 3, ncol = 2, byrow = TRUE),
    apertures = 0
  )
  motion$time = c(2, 2.2, 2.4)
  cam = camera(motion = motion)
  scene = add_camera(sphere(radius = .25, material = light(intensity = 1)), cam)
  args = list(
    scene = scene,
    width = 32,
    height = 32,
    samples = 64,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    integrator_type = "nee"
  )
  set.seed(511)
  still = do.call(render_scene, c(args, list(time = 2)))
  set.seed(511)
  exposed = do.call(
    render_scene,
    c(args, list(time = 2, shutter_offsets = c(0, .4)))
  )
  expect_gt(sum(exposed[,, 1] > .02), 1.3 * sum(still[,, 1] > .02))
  times = c(2, 2.1, 2)
  set.seed(511)
  frames = do.call(
    render_animation,
    c(args, list(times = times, shutter_offsets = c(0, .4)))
  )
  expect_identical(as.numeric(frames[[1]]), as.numeric(frames[[3]]))
  expect_lt(mean(abs(frames[[1]][,, 1:3] - exposed[,, 1:3])), .01)
  # Batched cameras retain separate authored tracks rather than borrowing the
  # first camera's exposure. Requested times still index output, not playback.
  second = motion
  second$x = second$x + .5
  second$dx = second$dx + .5
  args$scene = add_camera(scene, camera(motion = second), name = "second")
  set.seed(511)
  batch = do.call(
    render_animation,
    c(args, list(camera = "all", times = times, shutter_offsets = c(0, .4)))
  )
  expect_length(batch, 2)
  expect_equal(as.numeric(batch[[1]][[1]]), as.numeric(frames[[1]]))
  expect_gt(mean(abs(batch[[1]][[1]][,, 1:3] - batch[[2]][[1]][,, 1:3])), .001)
  set.seed(511)
  batch_still = do.call(
    render_scene,
    c(args, list(camera = "all", time = 2, shutter_offsets = c(0, .4)))
  )
  expect_lt(mean(abs(batch_still[[2]][,, 1:3] - batch[[2]][[1]][,, 1:3])), .01)
  args$scene = scene
  args$integrator_type = "metal"
  set.seed(511)
  device = suppressWarnings(do.call(
    render_scene,
    c(args, list(time = 2, shutter_offsets = c(0, .4)))
  ))
  if (!isTRUE(attr(device, "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  expect_lt(mean(abs(device[,, 1:3] - exposed[,, 1:3])), .01)
  set.seed(511)
  device_frames = do.call(
    render_animation,
    c(args, list(times = times, shutter_offsets = c(0, .4)))
  )
  expect_true(all(vapply(
    device_frames,
    function(x) isTRUE(attr(x, "wavefront")$used),
    logical(1)
  )))
  expect_true(isTRUE(attr(device_frames[[2]], "wavefront")$scene_reused))
  expect_identical(
    as.numeric(device_frames[[1]]),
    as.numeric(device_frames[[3]])
  )
})

test_that("projection cuts share snapshot and exposure semantics across render paths", {
  withr::local_options(list(cores = 2L))
  motion = generate_camera_motion(
    positions = rbind(c(-1, 0, 6), c(1, 0, 6)),
    lookats = rbind(c(0, 0, 0), c(0, 0, 0)),
    frames = 4,
    fovs = 40,
    apertures = 0
  )
  motion$time = 0:3
  motion$fov = c(40, 0, 360, 40)
  motion$orthox = 4
  motion$orthoy = 3
  scene = add_camera(
    sphere(material = light(intensity = 1)),
    camera(motion = motion)
  )
  args = list(
    scene = scene,
    width = 24,
    height = 18,
    samples = 32,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    integrator_type = "nee"
  )
  times = c(0, 1, 2, 3, 1)
  for (offsets in list(c(0, 0), c(-.5, .5))) {
    set.seed(614)
    frames = do.call(
      render_animation,
      c(args, list(times = times, shutter_offsets = offsets))
    )
    expect_identical(as.numeric(frames[[2]]), as.numeric(frames[[5]]))
    for (index in seq_along(times)) {
      set.seed(614)
      still = do.call(
        render_scene,
        c(args, list(time = times[index], shutter_offsets = offsets))
      )
      set.seed(614)
      preview = do.call(
        render_scene,
        c(
          args,
          list(
            time = times[index],
            timeline = c(0, 3),
            shutter_offsets = offsets
          )
        )
      )
      expect_equal(as.numeric(still), as.numeric(preview), tolerance = 1e-6)
      expect_lt(mean(abs(frames[[index]][,, 1:3] - still[,, 1:3])), .01)
    }
    args$integrator_type = "metal"
    set.seed(614)
    device = suppressWarnings(do.call(
      render_animation,
      c(args, list(times = times, shutter_offsets = offsets))
    ))
    if (
      !all(vapply(
        device,
        function(x) isTRUE(attr(x, "wavefront")$used),
        logical(1)
      ))
    ) {
      skip("Metal runtime unavailable")
    }
    expect_identical(as.numeric(device[[2]]), as.numeric(device[[5]]))
    for (index in seq_along(times)) {
      expect_lt(
        mean(abs(device[[index]][,, 1:3] - frames[[index]][,, 1:3])),
        .02
      )
    }
    args$integrator_type = "nee"
  }
})
