test_that("loop metadata survives camera construction and RDS", {
  for (type in c("linear", "spline", "bezier")) {
    for (damped in c(FALSE, TRUE)) {
      motion = generate_camera_motion(
        list(c(-3, 0, 6), c(3, 0, 6), c(3, 2, 6)),
        frames = 12,
        closed = TRUE,
        type = type,
        damp_motion = damped,
        progress = FALSE
      )
      expect_equal(nrow(motion), 12)
      expect_true(attr(camera(motion = motion)$motion, "closed"))
      endpoints = camera_shutter_endpoints(motion)
      expect_equal(
        unname(endpoints[12, ]),
        as.numeric(motion[1, colnames(endpoints)])
      )
      expect_gt(
        sum((endpoints[12, 1:3] - as.numeric(motion[12, c("x", "y", "z")]))^2),
        0
      )
      path = tempfile(fileext = ".rds")
      saveRDS(motion, path)
      expect_identical(camera_shutter_endpoints(readRDS(path)), endpoints)
      unlink(path)
    }
  }
})

test_that("open and manual paths retain their row counts and endpoint semantics", {
  positions = list(c(-3, 0, 6), c(3, 0, 6), c(3, 2, 6))
  open = generate_camera_motion(
    positions,
    frames = 12,
    type = "linear",
    progress = FALSE
  )
  endpoints = camera_shutter_endpoints(open)
  expect_equal(nrow(open), 12)
  expect_false(isTRUE(attr(open, "closed")))
  expect_equal(
    unname(endpoints[12, ]),
    as.numeric(open[12, colnames(endpoints)])
  )
  manual = generate_camera_motion(
    positions,
    closed = TRUE,
    type = "manual",
    progress = FALSE
  )
  expect_equal(nrow(manual), 3)
  expect_equal(
    unname(camera_shutter_endpoints(manual)[3, ]),
    as.numeric(manual[1, colnames(endpoints)])
  )
  single = generate_camera_motion(
    list(c(0, 0, 6)),
    closed = TRUE,
    type = "manual",
    progress = FALSE
  )
  expect_equal(nrow(camera_shutter_endpoints(single)), 1)
})

test_that("batch selection preserves each camera's actual shutter neighbors", {
  motion = generate_camera_motion(
    list(c(-3, 0, 6), c(3, 0, 6), c(3, 2, 6)),
    frames = 12,
    closed = TRUE,
    type = "linear",
    progress = FALSE
  )
  other = motion
  other$x = other$x + 100
  cameras = list(
    camera(motion = motion, name = "a"),
    camera(motion = other, name = "b")
  )
  plan = camera_batch_plan(
    cameras,
    mode = "animation",
    start_frame = 4,
    end_frame = 6
  )
  expect_equal(plan$shutter_endpoints[c(1, 3), "x"], motion$x[c(5, 7)])
  expect_equal(plan$shutter_endpoints[c(4, 6), "x"], other$x[c(5, 7)])
  plan = camera_batch_plan(
    cameras,
    mode = "animation",
    start_frame = 12,
    end_frame = 12
  )
  expect_equal(plan$shutter_endpoints[, "x"], c(motion$x[1], other$x[1]))
})

test_that("CPU and Metal render the closed final shutter like an explicit wrap", {
  skip_on_cran()
  previous = options(cores = 2L)
  on.exit(options(previous))
  motion = generate_camera_motion(
    list(c(-3, 0, 6), c(3, 0, 6), c(3, 2, 6)),
    frames = 6,
    closed = TRUE,
    type = "linear",
    damp_motion = TRUE,
    progress = FALSE
  )
  scene = sphere(material = light("red")) |>
    add_object(sphere(x = 2, radius = .5, material = light("green")))
  args = list(
    scene = scene,
    camera_motion = motion,
    start_frame = 6,
    end_frame = 6,
    camera_motion_blur = TRUE,
    width = 24,
    height = 24,
    samples = 32,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    bloom = FALSE,
    max_depth = 3,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  for (integrator in c("nee", "metal")) {
    args$integrator_type = integrator
    set.seed(314)
    loop = suppressWarnings(do.call(render_animation, args))[[1]]
    if (integrator == "metal" && !isTRUE(attr(loop, "wavefront")$used)) {
      next
    }
    explicit = motion[c(6, 1), ]
    attr(explicit, "closed") = NULL
    reference_args = modifyList(
      args,
      list(camera_motion = explicit, start_frame = 1, end_frame = 1)
    )
    set.seed(314)
    reference = do.call(render_animation, reference_args)[[1]]
    expect_equal(as.numeric(loop), as.numeric(reference), tolerance = 0)
    set.seed(314)
    still = do.call(
      render_animation,
      modifyList(args, list(camera_motion_blur = FALSE))
    )[[1]]
    expect_gt(max(abs(as.numeric(loop) - as.numeric(still))), .001)
    # A one-frame selection from each camera must use its original loop's
    # neighbor, rather than the next camera or the end of the selected batch.
    second_motion = motion
    second_motion$x = second_motion$x + 3
    batch_args = args
    batch_args$camera_motion = NULL
    batch_args$camera_motion_blur = NULL
    batch_args$cameras = list(
      camera(motion = motion, name = "first", camera_motion_blur = TRUE),
      camera(motion = second_motion, name = "second", camera_motion_blur = TRUE)
    )
    batch_args$mode = "animation"
    set.seed(314)
    batch = do.call(render_camera_batch, batch_args)
    expect_equal(as.numeric(batch$first), as.numeric(loop), tolerance = 0)
    set.seed(314)
    second = do.call(
      render_animation,
      modifyList(args, list(camera_motion = second_motion))
    )[[1]]
    expect_equal(as.numeric(batch$second), as.numeric(second), tolerance = 0)
  }
})
