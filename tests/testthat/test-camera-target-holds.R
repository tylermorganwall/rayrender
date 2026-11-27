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
          expect_equal(
            actual,
            matrix(target, nrow(motion), 3, byrow = TRUE)
          )
          view = actual - as.matrix(motion[, c("x", "y", "z")])
          view = view / sqrt(rowSums(view^2))
          up = as.matrix(motion[, c("upx", "upy", "upz")])
          if (!damp) {
            expect_equal(
              rowSums(view * up),
              rep(0, nrow(motion)),
              tolerance = 1e-12
            )
            expect_equal(rowSums(up^2), rep(1, nrow(motion)), tolerance = 1e-12)
          } else {
            expect_true(all(is.finite(up)))
            expect_true(all(rowSums(up^2) > 0))
          }
        }
      }
    }
  }
})

test_that("damping filters the complete motion across both target-hold boundaries", {
  positions = rbind(c(-5, 3, 10), c(3, 4, 10), c(-2, 3, 9), c(4, 4, 8))
  targets = rbind(c(-2, 0, 0), c(1, 0, 0), c(1, 0, 0), c(3, 1, 0))
  for (closed in c(FALSE, TRUE)) {
    args = list(
      positions = positions,
      lookats = targets,
      frames = if (closed) 160L else 121L,
      closed = closed,
      progress = FALSE
    )
    undamped = do.call(
      generate_camera_motion,
      c(args, list(damp_motion = FALSE))
    )
    damped = do.call(generate_camera_motion, c(args, list(damp_motion = TRUE)))
    expect_identical(
      unname(as.matrix(undamped[41:81, c("dx", "dy", "dz")])),
      matrix(c(1, 0, 0), 41, 3, byrow = TRUE)
    )
    # Independently apply the documented low-pass recurrence to every column.
    # For a loop, the previous state is its periodic steady-state final row.
    expected = as.matrix(undamped)
    previous = if (closed) {
      weights = .9^(nrow(expected) - seq_len(nrow(expected)))
      colSums(expected * weights) * .1 / (1 - .9^nrow(expected))
    } else {
      expected[1, ]
    }
    for (row in seq_len(nrow(expected))) {
      previous = .9 * previous + .1 * expected[row, ]
      expected[row, ] = previous
    }
    expect_equal(as.matrix(damped), expected, tolerance = 1e-12)
    # Entering a hold must preserve damping lag, not snap the target to the key.
    expect_gt(abs(damped$dx[41] - 1), 1e-3)
    expect_gt(max(abs(damped$dx[82:121] - 1)), .1)
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
