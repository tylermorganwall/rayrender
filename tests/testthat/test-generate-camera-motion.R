test_that("closed damped camera motion closes the endpoint", {
  positions = list(
    c(0, 0, 0),
    c(10, 0, 0),
    c(10, 10, 0),
    c(0, 10, 0)
  )

  motion = generate_camera_motion(
    positions = positions,
    type = "cubic",
    frames = 12,
    closed = TRUE,
    damp_motion = TRUE,
    progress = FALSE
  )

  expect_equal(
    unname(as.numeric(motion[1, c("x", "y", "z")])),
    unname(as.numeric(motion[nrow(motion), c("x", "y", "z")]))
  )

  bezier_motion = generate_camera_motion(
    positions = positions,
    type = "bezier",
    frames = 12,
    closed = TRUE,
    damp_motion = TRUE,
    progress = FALSE
  )

  expect_equal(
    unname(as.numeric(bezier_motion[1, c("x", "y", "z")])),
    unname(as.numeric(bezier_motion[nrow(bezier_motion), c("x", "y", "z")]))
  )
})

test_that("open damped camera motion keeps the one-way recurrence", {
  positions = list(
    c(0, 0, 0),
    c(10, 0, 0),
    c(10, 10, 0),
    c(0, 10, 0)
  )

  undamped = generate_camera_motion(
    positions = positions,
    type = "linear",
    frames = 8,
    closed = FALSE,
    damp_motion = FALSE,
    progress = FALSE
  )
  damped = generate_camera_motion(
    positions = positions,
    type = "linear",
    frames = 8,
    closed = FALSE,
    damp_motion = TRUE,
    progress = FALSE
  )

  expected = as.matrix(undamped)
  damp_magnitude = 1 - 0.1
  current_pos = expected[1, ]
  for (i in seq_len(nrow(expected))[-1]) {
    current_pos = current_pos *
      damp_magnitude +
      expected[i, ] * (1 - damp_magnitude)
    expected[i, ] = current_pos
  }

  expect_equal(as.matrix(damped), expected)
})

test_that("saved keyframes remove only sequential duplicate camera states", {
  env = rayrender:::ray_environment
  old_keyframes = get("keyframes", envir = env)
  on.exit(assign("keyframes", old_keyframes, envir = env), add = TRUE)

  keyframes = data.frame(
    x = c(0, 1, 1, 2, 1, 1 + 1e-12),
    y = c(0, 0, 0, 0, 0, 0),
    z = c(0, 0, 0, 0, 0, 0),
    dx = c(10, 10, 10, 20, 10, 10),
    dy = c(0, 0, 0, 0, 0, 0),
    dz = c(0, 0, 0, 0, 0, 0),
    aperture = c(0, 0.1, 0.1, 0.2, 0.1, 0.1),
    fov = c(40, 35, 35, 30, 35, 35),
    focal = c(10, 9, 9, 8, 9, 9),
    exposure = c(1, 0.5, 0.5, 0.25, 0.5, 0.5),
    orthox = c(1, 1, 1, 1, 1, 1),
    orthoy = c(1, 1, 1, 1, 1, 1),
    upx = c(0, 0, 0, 0, 0, 0),
    upy = c(1, 1, 1, 1, 1, 1),
    upz = c(0, 0, 0, 0, 0, 0)
  )
  assign("keyframes", keyframes, envir = env)

  expected = keyframes[c(1, 2, 4, 5, 6), , drop = FALSE]
  rownames(expected) = NULL
  expect_equal(get_saved_keyframes(), expected)
})

test_that("bezier camera motion handles repeated scalar keyframe values", {
  keyframes = data.frame(
    x = c(0, 1, 2, 3),
    y = c(0, 0, 0, 0),
    z = c(0, 0, 0, 0),
    dx = c(0, 1, 2, 3),
    dy = c(0, 0, 0, 0),
    dz = c(0, 0, 0, 0),
    aperture = c(0, 1, 1, 2),
    fov = c(40, 35, 35, 30),
    focal = c(1, 1, 1, 1),
    orthox = c(1, 1, 1, 1),
    orthoy = c(1, 1, 1, 1),
    upx = c(0, 0, 0, 0),
    upy = c(1, 1, 1, 1),
    upz = c(0, 0, 0, 0)
  )

  motion = generate_camera_motion(
    keyframes,
    type = "bezier",
    frames = 4,
    progress = FALSE
  )

  expect_equal(nrow(motion), 4)
  expect_true(all(is.finite(as.matrix(motion))))
})
