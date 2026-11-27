test_that("track validation catches malformed and singular curves", {
  expect_error(
    keyframe_track(numeric(), matrix(numeric(), 0, 3)),
    "at least one"
  )
  expect_error(keyframe_track(c(0, 0), matrix(0, 2, 3)), "strictly increasing")
  expect_error(keyframe_track(c(0, Inf), matrix(0, 2, 3)), "finite")
  expect_error(keyframe_track(c(0, 1), 1:3), "never recycled")
  expect_error(keyframe_track(0, c(0, NA, 0)), "finite")
  expect_error(keyframe_track(0, c(0, 0, 0, 0), "rotation"), "norm")
  expect_error(
    keyframe_track(c(0, 1), rbind(c(1, 1, 1), c(-1, 1, 1)), "scale"),
    "zero"
  )
  expect_error(
    keyframe_track(
      c(0, 1),
      matrix(1, 2, 3),
      "scale",
      "cubic",
      incoming = matrix(8, 2, 3),
      outgoing = matrix(-8, 2, 3)
    ),
    "zero"
  )
  expect_error(
    keyframe_track(
      c(0, 1),
      rbind(c(0, 0, 0, 1), c(0, 0, 0, -1)),
      "rotation",
      "cubic",
      matrix(0, 2, 4),
      matrix(0, 2, 4)
    ),
    "singularity"
  )
  expect_error(
    transform_track(position = keyframe_track(0, 1, "scalar")),
    "target"
  )
  expect_error(
    transform_track(0, angle = c(0, 0, 0), rotation = c(0, 0, 0, 1)),
    "only one"
  )
})

test_that("keys, STEP sides, nonuniform tangents and SLERP are deterministic", {
  track = keyframe_track(c(-1, .5, 3, 8), c(-2, 1, 6, 16), "scalar", "cubic")
  times = c(8, 0, 1.5, -10, .5, 7, 0)
  expect_equal(
    as.vector(animation_evaluate_track(track, times)),
    pmin(16, pmax(-2, 2 * times))
  )
  step = keyframe_track(c(0, 1, 3), c(2, 4, 8), "scalar", "step")
  expect_equal(
    as.vector(animation_evaluate_track(step, c(-1, .99999, 1, 2, 3))),
    c(2, 2, 4, 4, 8)
  )
  constant = keyframe_track(7, c(-1, 2, 3), "scale")
  expect_equal(
    animation_evaluate_track(constant, c(-20, 7, 100)),
    matrix(rep(c(-1, 2, 3), each = 3), 3)
  )
  q = keyframe_track(c(0, 2), rbind(c(0, 0, 0, 1), c(0, 1, 0, 0)), "rotation")
  expect_equal(
    as.vector(animation_evaluate_track(q, 1)),
    c(0, sqrt(.5), 0, sqrt(.5))
  )
  opposite = keyframe_track(
    c(0, 2),
    rbind(c(0, 0, 0, 1), c(0, 0, 0, -1)),
    "rotation"
  )
  expect_equal(as.vector(animation_evaluate_track(opposite, 1)), c(0, 0, 0, 1))
  # Derivatives are per second: a two-second Hermite segment peaks at one.
  cubic = keyframe_track(
    c(0, 2),
    c(0, 0),
    "scalar",
    "cubic",
    matrix(c(0, -2), 2),
    matrix(c(2, 0), 2)
  )
  expect_equal(as.vector(animation_evaluate_track(cubic, 1)), 1)
})

test_that("smooth rotation has shared derivatives at nonuniform keys", {
  angles = c(0, .6, 1.4, 2.5)
  track = keyframe_track(
    c(0, .7, 2.8, 4),
    cbind(0, sin(angles / 2), 0, cos(angles / 2)),
    "rotation",
    "smooth"
  )
  h = 1e-5
  for (t in c(.7, 2.8)) {
    q = animation_evaluate_track(track, c(t - h, t, t + h))
    expect_equal((q[2, ] - q[1, ]) / h, (q[3, ] - q[2, ]) / h, tolerance = 1e-4)
  }
})

test_that("independent schedules, pivot order and large times survive serialization", {
  path = transform_track(
    position = keyframe_track(c(0, 3), rbind(c(0, 0, 0), c(3, 0, 0))),
    rotation = keyframe_track(
      c(0, 1, 4),
      rbind(c(0, 0, 0, 1), c(0, 0, 1, 0), c(0, 0, 0, 1)),
      "rotation"
    ),
    pivot = c(1, 0, 0)
  )
  at = animation_evaluate_transform(path, 1)[[1]]
  expect_equal(as.vector(at %*% c(2, 0, 0, 1)), c(1, 0, 0, 1), tolerance = 1e-6)
  times = c(3, 1, -1, 2, 1)
  expect_identical(
    animation_evaluate_transform(unserialize(serialize(path, NULL)), times),
    animation_evaluate_transform(path, times)
  )
  big = transform_track(
    1e12 + c(0, .25, .5, 1),
    position = cbind(c(0, 1, 2, 0), 0, 0)
  )
  expect_equal(animation_evaluate_transform(big, 1e12 + .125)[[1]][1, 4], .5)
  angles = rbind(c(20, 30, 40), c(-15, 45, 120))
  euler = transform_track(c(0, 1), angle = angles)
  for (i in 1:2) {
    expect_equal(
      animation_evaluate_transform(euler, i - 1)[[1]],
      generate_rotation_matrix(angles[i, ], 1:3),
      tolerance = 1e-6
    )
  }
})

test_that("scene composition preserves row-owned tracks and camera metadata", {
  path = transform_track(c(0, 1), position = rbind(c(0, 0, 0), c(1, 0, 0)))
  moving = animate_keyframes(sphere(), path)
  moving = add_camera(moving, camera(name = "animation_camera"))
  scene = add_object(moving, cube())
  expect_identical(ray_scene_cameras(scene), ray_scene_cameras(moving))
  expect_identical(scene$animation_info[[1]]$keyframes, path)
  expect_null(scene$animation_info[[2]]$keyframes)
  expect_identical(scene[1, ]$animation_info[[1]]$keyframes, path)
  expect_identical(
    group_objects(scene, translate = c(1, 0, 0))$animation_info,
    scene$animation_info
  )
  expect_identical(
    unserialize(serialize(scene, NULL))$animation_info,
    scene$animation_info
  )
  expect_error(animate_keyframes(moving, path), "already")
  expect_error(animate_objects(moving), "already")
  expect_error(animate_keyframes(animate_objects(sphere()), path), "already")
})

test_that("explicit timing rejects conflicts and unsupported paths", {
  expect_null(resolve_scene_timing(NULL, NULL, 30, NULL, FALSE))
  expect_equal(
    resolve_scene_timing(4, c(4, 4), 30, NULL, FALSE)$offsets,
    c(0, 0)
  )
  expect_error(resolve_scene_timing(4, NULL, 30, NULL, TRUE), "legacy")
  expect_error(resolve_scene_timing(NULL, c(0, 1), 30, NULL, FALSE), "explicit")
  expect_error(resolve_scene_timing(0, NULL, 30, c(1, 0), FALSE), "ordered")
  expect_error(resolve_scene_timing(0, NULL, 0, NULL, FALSE), "fps")
  expect_error(render_scene(sphere(), time = 0, shutteropen = 0), "legacy")
})

test_that("CPU snapshots agree after unrelated seeks and with static geometry", {
  withr::local_options(list(cores = 2L))
  render = function(scene, time = NULL, ...) {
    set.seed(117)
    render_scene(
      scene,
      time = time,
      width = 16,
      height = 16,
      samples = 4,
      lookfrom = c(0, 0, 8),
      lookat = c(0, 0, 0),
      aperture = 0,
      fov = 0,
      ortho_dimensions = c(6, 6),
      integrator_type = "nee",
      ambient_light = TRUE,
      backgroundhigh = "black",
      backgroundlow = "black",
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      min_variance = 0,
      ...
    )
  }
  times = c(0, .7, 2.3, 4)
  path = transform_track(
    times,
    position = cbind(c(-1, 0, 1, -1), c(0, 1, 0, 0), 0),
    position_interpolation = "cubic"
  )
  primitive = sphere(x = .3, radius = .45, material = light())
  moving = animate_keyframes(primitive, path)
  first = render(moving, 1.4)
  expect_gt(max(first), 0)
  expect_gt(sum(abs(first - render(moving, 4))), 0)
  for (t in c(4, -1, 2.3, .1)) {
    render(moving, t)
  }
  expect_identical(render(moving, 1.4), first)
  p = animation_evaluate_transform(path, 1.4)[[1]][1:3, 4]
  expect_equal(
    first,
    render(group_objects(primitive, translate = p), 1.4),
    tolerance = 1e-6
  )
  # Adapted legacy evaluation must match the old two-matrix shutter route.
  legacy = animate_objects(
    primitive,
    end_position = c(1, 0, 0),
    end_angle = c(0, 30, 0)
  )
  expect_equal(
    render(legacy, .4),
    render(legacy, shutteropen = .4, shutterclose = .4),
    tolerance = 1e-6
  )
  expect_error(render(moving), "explicit scene time")
})
