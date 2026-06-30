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
