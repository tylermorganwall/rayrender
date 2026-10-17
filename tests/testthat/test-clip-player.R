test_that("players map independent absolute clocks with documented endpoints", {
  path = transform_track(c(2, 3, 4), position = cbind(c(2, 3, 4), 0, 0))
  times = c(10, 9.75, 10.75, 9.5, 10)
  scene = animate_keyframes(
    sphere(),
    path,
    player = clip_player(
      start_time = 10,
      offset = .5,
      speed = 2,
      loop = "repeat"
    )
  )
  poses = animation_evaluate_transform(
    scene$animation_info[[1]]$keyframes,
    times
  )
  expect_equal(
    vapply(poses, function(x) x[1, 4], numeric(1)),
    c(2.5, 2, 2, 3.5, 2.5)
  )
  expect_identical(poses[[1]], poses[[5]])
  expect_null(path$player)

  for (mode in c("hold", "repeat", "pingpong")) {
    for (speed in c(-2, 0, .5, 2)) {
      scene = animate_keyframes(
        sphere(),
        path,
        player = clip_player(
          start_time = .3,
          offset = .7,
          speed = speed,
          loop = mode
        )
      )
      times = c(-4, 8, .3, -1, 2, .3)
      x = speed * (times - .3) + .7
      expected = switch(
        mode,
        hold = pmax(0, pmin(2, x)),
        "repeat" = x %% 2,
        pingpong = 2 - abs((x %% 4) - 2)
      ) +
        2
      poses = animation_evaluate_transform(
        scene$animation_info[[1]]$keyframes,
        times
      )
      expect_equal(
        vapply(poses, function(x) x[1, 4], numeric(1)),
        expected,
        tolerance = 1e-6
      )
    }
  }
  constant = animate_keyframes(
    sphere(),
    transform_track(7, position = c(1, 2, 3)),
    player = clip_player(loop = "repeat")
  )
  poses = animation_evaluate_transform(
    constant$animation_info[[1]]$keyframes,
    c(-1e6, 1e6)
  )
  expect_identical(poses[[1]], poses[[2]])
  expect_error(clip_player(speed = NA_real_), "finite")
  expect_error(clip_player(offset = 1:2), "finite")
  expect_error(clip_player(root_motion = "extract"), "preserve")
  expect_error(
    animate_keyframes(sphere(), path, clip_player(clip = "missing")),
    "named clips"
  )
  invalid = clip_player()
  invalid$speed = c(1, 2)
  expect_error(animate_keyframes(sphere(), path, invalid), "scalars")
  expect_false(animation_loop_seam(path)$pose_continuous)
  closed = transform_track(0:2, position = cbind(c(0, 1, 0), 0, 0))
  expect_true(animation_loop_seam(closed)$pose_continuous)
})

test_that("copying and combining looped attachments does not couple their clocks", {
  path = transform_track(0:2, position = cbind(0:2, 0, 0))
  first = animate_keyframes(sphere(), path, clip_player(loop = "repeat"))
  second = animate_keyframes(
    sphere(),
    path,
    clip_player(offset = .5, speed = -1, loop = "repeat")
  )
  combined = add_object(first, second) |>
    add_camera(camera(name = "loop_camera"))
  copy = unserialize(serialize(combined, NULL))
  expect_identical(combined$animation_info, copy$animation_info)
  expect_identical(ray_scene_cameras(combined), ray_scene_cameras(copy))
  times = c(7, -.1, 1, 0, -.1)
  one = animation_evaluate_transform(copy$animation_info[[1]]$keyframes, times)
  two = animation_evaluate_transform(copy$animation_info[[2]]$keyframes, times)
  expect_identical(one[[2]], one[[5]])
  expect_identical(two[[2]], two[[5]])
  expect_false(identical(one, two))
  expect_identical(first$animation_info[[1]]$keyframes$player$offset, 0)
})

test_that("looped offline frames and seam exposures agree on CPU and Metal", {
  withr::local_options(list(cores = 2L))
  path = transform_track(
    c(0, .5, 1.5, 2),
    position = rbind(c(-1, 0, 0), c(-.5, .7, 0), c(.5, -.3, 0), c(1, 0, 0)),
    position_interpolation = "cubic"
  )
  settings = list(
    width = 18,
    height = 18,
    samples = 64,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  for (mode in c("repeat", "pingpong")) {
    scene = animate_keyframes(
      sphere(radius = .35, material = light(intensity = 1)),
      path,
      clip_player(loop = mode)
    ) |>
      add_camera(camera(
        lookfrom = c(0, 0, 6),
        lookat = c(0, 0, 0),
        aperture = 0,
        fov = 0,
        ortho_dimensions = c(4, 4)
      ))
    times = c(-.001, 0, .001, 1, 1.999, 2, 2.001, 4, -.001)
    set.seed(812)
    cpu = do.call(
      render_animation,
      c(list(scene = scene, times = times, integrator_type = "nee"), settings)
    )
    expect_identical(as.numeric(cpu[[1]]), as.numeric(cpu[[9]]))
    expect_identical(as.numeric(cpu[[2]]), as.numeric(cpu[[8]]))
    expect_gt(mean(abs(as.numeric(cpu[[2]]) - as.numeric(cpu[[4]]))), .01)
    set.seed(812)
    gpu = suppressWarnings(do.call(
      render_animation,
      c(list(scene = scene, times = times, integrator_type = "metal"), settings)
    ))
    if (!isTRUE(attr(gpu[[1]], "wavefront")$used)) {
      skip("Metal runtime unavailable")
    }
    for (i in seq_along(cpu)) {
      expect_true(attr(gpu[[i]], "wavefront")$used)
      expect_lt(mean(abs(as.numeric(cpu[[i]]) - as.numeric(gpu[[i]]))), .02)
    }
    for (backend in c("nee", "metal")) {
      set.seed(813)
      exposure = do.call(
        render_scene,
        c(
          list(
            scene = scene,
            time = 2,
            shutter_offsets = c(-.15, .15),
            integrator_type = backend
          ),
          settings
        )
      )
      expect_true(all(is.finite(exposure)))
      expect_gt(mean(exposure), .005)
      if (backend == "nee") {
        reference = exposure
      } else {
        expect_true(attr(exposure, "wavefront")$used)
        expect_lt(mean(abs(as.numeric(reference) - as.numeric(exposure))), .02)
      }
    }
  }
  # A tiny negative time must remain just before the seam. If float packing
  # would round it to the beginning, the documented warning/fallback is exact.
  scene$animation_info[[1]]$keyframes$player = clip_player(loop = "repeat")
  set.seed(814)
  expected = do.call(
    render_scene,
    c(list(scene = scene, time = -1e-10, integrator_type = "nee"), settings)
  )
  set.seed(814)
  expect_warning(
    actual <- do.call(
      render_scene,
      c(list(scene = scene, time = -1e-10, integrator_type = "metal"), settings)
    ),
    "rounding crosses a loop boundary"
  )
  expect_false(isTRUE(attr(actual, "wavefront")$used))
  expect_identical(as.numeric(actual), as.numeric(expected))
})
