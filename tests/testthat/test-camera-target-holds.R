test_that("orbit targets stay exact across smoothing, damping and closed loops", {
  for (target in list(c(0, 0, 0), c(2, 1, -3), c(1e8, -2e8, 3e8))) {
    angles = c(-35, 25, -15) * pi / 180
    positions = sweep(
      cbind(10 * sin(angles), rep(4, 3), 10 * cos(angles)),
      2,
      target,
      "+"
    )
    lookats = matrix(target, 3, 3, byrow = TRUE)
    for (type in c("spline", "linear", "quad", "cubic", "exp")) {
      for (closed in c(FALSE, TRUE)) {
        for (damp in c(FALSE, TRUE)) {
          motion = generate_camera_motion(
            positions,
            lookats = lookats,
            type = type,
            frames = 93,
            closed = closed,
            damp_motion = damp,
            progress = FALSE
          )
          actual = unname(as.matrix(motion[, c("dx", "dy", "dz")]))
          expect_identical(
            actual,
            matrix(target, nrow(motion), 3, byrow = TRUE)
          )
          view = actual - as.matrix(motion[, c("x", "y", "z")])
          view = view / sqrt(rowSums(view^2))
          up = as.matrix(motion[, c("upx", "upy", "upz")])
          expect_equal(
            rowSums(view * up),
            rep(0, nrow(motion)),
            tolerance = 1e-12
          )
          expect_equal(rowSums(up^2), rep(1, nrow(motion)), tolerance = 1e-12)
        }
      }
    }
  }
})

test_that("a stationary target interval stays locked next to moving targets", {
  positions = rbind(c(-5, 3, 10), c(3, 4, 10), c(-2, 3, 9), c(4, 4, 8))
  targets = rbind(c(-2, 0, 0), c(1, 0, 0), c(1, 0, 0), c(3, 1, 0))
  for (damp in c(FALSE, TRUE)) {
    motion = generate_camera_motion(
      positions,
      lookats = targets,
      frames = 121,
      damp_motion = damp,
      progress = FALSE
    )
    expect_identical(
      unname(as.matrix(motion[41:81, c("dx", "dy", "dz")])),
      matrix(c(1, 0, 0), 41, 3, byrow = TRUE)
    )
    expect_gt(max(abs(motion$dx[82:121] - 1)), 0.1)
  }
})

test_that("target holds preserve authored roll and do not snap small pans", {
  positions = rbind(c(-5, 2, 10), c(0, 2, 10), c(5, 2, 10))
  targets = matrix(c(1, 0, 0), 3, 3, byrow = TRUE)
  ups = rbind(c(0, 1, 0), c(1, 1, 0), c(-1, 1, 0))
  motion = generate_camera_motion(
    positions,
    lookats = targets,
    camera_ups = ups,
    frames = 81,
    damp_motion = FALSE,
    progress = FALSE
  )
  forward = targets - positions
  forward = forward / sqrt(rowSums(forward^2))
  expected_up = ups - rowSums(ups * forward) * forward
  expected_up = expected_up / sqrt(rowSums(expected_up^2))
  expect_equal(
    unname(as.matrix(motion[c(1, 41, 81), c("upx", "upy", "upz")])),
    unname(expected_up),
    tolerance = 1e-12
  )
  targets[, 1] = targets[, 1] + c(0, 1e-10, 2e-10)
  orientation = tween_camera_orientation(positions, targets, ups, positions)
  expect_false(any(orientation$fixed_lookat_rows))
})
