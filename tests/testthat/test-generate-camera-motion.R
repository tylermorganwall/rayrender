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

test_that("spline camera motion changes speed smoothly without curving the path", {
  expect_identical(formals(generate_camera_motion)$type, "spline")

  positions = list(c(0, 0, 0), c(1, 0, 0), c(1, 10, 0))

  spline_motion = generate_camera_motion(
    positions = positions,
    type = "spline",
    frames = 200,
    damp_motion = FALSE,
    progress = FALSE
  )
  linear_motion = generate_camera_motion(
    positions = positions,
    type = "linear",
    frames = 200,
    damp_motion = FALSE,
    progress = FALSE
  )

  keyframe_rows = c(1, 101, 200)
  expect_equal(
    unname(as.matrix(spline_motion[keyframe_rows, c("x", "y", "z")])),
    do.call(rbind, positions)
  )
  expect_equal(spline_motion$y[1:101], rep(0, 101))
  expect_equal(spline_motion$x[101:200], rep(1, 100))

  spline_speed = sqrt(rowSums(
    diff(as.matrix(spline_motion[, c("x", "y", "z")]))^2
  ))
  linear_speed = sqrt(rowSums(
    diff(as.matrix(linear_motion[, c("x", "y", "z")]))^2
  ))
  spline_speed_change = abs(
    spline_speed[101] - spline_speed[100]
  )
  linear_speed_change = abs(
    linear_speed[101] - linear_speed[100]
  )
  expect_lt(spline_speed_change, linear_speed_change / 10)
})

test_that("closed spline camera motion has periodic boundary speed", {
  positions = list(
    c(0, 0, 0),
    c(1, 0, 0),
    c(1, 10, 0),
    c(0, 10, 0)
  )

  motion = generate_camera_motion(
    positions = positions,
    apertures = c(0, 0.1, 0.2, 0.3),
    fovs = c(40, 45, 50, 55),
    ortho_dims = list(c(1, 1), c(2, 1), c(3, 2), c(4, 3)),
    type = "spline",
    frames = 400,
    closed = TRUE,
    damp_motion = FALSE,
    progress = FALSE
  )
  position = as.matrix(motion[, c("x", "y", "z")])
  final_frame = nrow(position)
  start_velocity = (-3 * position[1, ] + 4 * position[2, ] - position[3, ]) / 2
  end_velocity = (3 *
    position[final_frame, ] -
    4 * position[final_frame - 1, ] +
    position[final_frame - 2, ]) /
    2

  expect_equal(
    unname(as.numeric(motion[1, ])),
    unname(as.numeric(motion[final_frame, ]))
  )
  expect_equal(
    sqrt(sum(start_velocity^2)),
    sqrt(sum(end_velocity^2)),
    tolerance = 3e-3
  )
})

test_that("orientation spline preserves keyed states and smooths angular velocity", {
  expect_identical(formals(generate_camera_motion)$smooth_orientation, TRUE)

  positions = rbind(c(0, 0, 0), c(1, 0, 0), c(1, 10, 0))
  angles = c(0, 10, 100) * pi / 180
  lookats = positions + cbind(sin(angles), 0, cos(angles))
  camera_ups = rbind(c(0, 1, 0), c(0, 1, 0), c(0, 1, 0))

  smoothed = generate_camera_motion(
    positions = positions,
    lookats = lookats,
    camera_ups = camera_ups,
    type = "spline",
    frames = 201,
    damp_motion = FALSE,
    progress = FALSE
  )
  direct = generate_camera_motion(
    positions = positions,
    lookats = lookats,
    camera_ups = camera_ups,
    type = "spline",
    smooth_orientation = FALSE,
    frames = 201,
    damp_motion = FALSE,
    progress = FALSE
  )

  position_columns = c("x", "y", "z")
  lookat_columns = c("dx", "dy", "dz")
  up_columns = c("upx", "upy", "upz")
  keyframe_rows = c(1, 101, 201)
  expect_equal(smoothed[, position_columns], direct[, position_columns])
  expect_equal(
    unname(as.matrix(smoothed[keyframe_rows, position_columns])),
    positions
  )
  expect_equal(
    unname(as.matrix(smoothed[keyframe_rows, lookat_columns])),
    lookats
  )
  expect_equal(
    unname(as.matrix(smoothed[keyframe_rows, up_columns])),
    camera_ups
  )

  unwrap_yaw = function(motion) {
    forward = as.matrix(motion[, lookat_columns]) -
      as.matrix(motion[, position_columns])
    yaw = atan2(forward[, 1], forward[, 3])
    cumsum(c(yaw[1], atan2(sin(diff(yaw)), cos(diff(yaw)))))
  }
  angular_velocity_jump = function(motion) {
    yaw = unwrap_yaw(motion)
    keyframe = 101
    left_velocity = (3 *
      yaw[keyframe] -
      4 * yaw[keyframe - 1] +
      yaw[keyframe - 2]) /
      2
    right_velocity = (-3 *
      yaw[keyframe] +
      4 * yaw[keyframe + 1] -
      yaw[keyframe + 2]) /
      2
    abs(left_velocity - right_velocity)
  }

  expect_lt(
    angular_velocity_jump(smoothed),
    angular_velocity_jump(direct) / 50
  )
})

test_that("closed orientation spline has periodic angular velocity", {
  positions = matrix(0, nrow = 4, ncol = 3)
  rotate_x = function(angle) {
    matrix(
      c(
        1,
        0,
        0,
        0,
        cos(angle),
        -sin(angle),
        0,
        sin(angle),
        cos(angle)
      ),
      nrow = 3,
      byrow = TRUE
    )
  }
  rotate_y = function(angle) {
    matrix(
      c(
        cos(angle),
        0,
        sin(angle),
        0,
        1,
        0,
        -sin(angle),
        0,
        cos(angle)
      ),
      nrow = 3,
      byrow = TRUE
    )
  }
  rotate_z = function(angle) {
    matrix(
      c(
        cos(angle),
        -sin(angle),
        0,
        sin(angle),
        cos(angle),
        0,
        0,
        0,
        1
      ),
      nrow = 3,
      byrow = TRUE
    )
  }
  angles = rbind(
    c(0, 0, 0),
    c(20, 30, 10),
    c(-35, 100, 45),
    c(50, 220, -30)
  ) *
    pi /
    180
  orientations = lapply(seq_len(nrow(angles)), function(i) {
    rotate_z(angles[i, 3]) %*%
      rotate_y(angles[i, 2]) %*%
      rotate_x(angles[i, 1])
  })
  lookats = t(vapply(orientations, function(value) value[, 3], numeric(3)))
  camera_ups = t(vapply(orientations, function(value) value[, 2], numeric(3)))

  motion = generate_camera_motion(
    positions = positions,
    lookats = lookats,
    camera_ups = camera_ups,
    type = "spline",
    frames = 401,
    closed = TRUE,
    damp_motion = FALSE,
    progress = FALSE
  )

  forward = as.matrix(motion[, c("dx", "dy", "dz")]) -
    as.matrix(motion[, c("x", "y", "z")])
  forward = forward / sqrt(rowSums(forward^2))
  up = as.matrix(motion[, c("upx", "upy", "upz")])
  up = up / sqrt(rowSums(up^2))
  periodic_derivative_difference = function(value) {
    final_frame = nrow(value)
    start_velocity = (-3 * value[1, ] + 4 * value[2, ] - value[3, ]) / 2
    end_velocity = (3 *
      value[final_frame, ] -
      4 * value[final_frame - 1, ] +
      value[final_frame - 2, ]) /
      2
    max(abs(start_velocity - end_velocity))
  }

  expect_equal(
    unname(as.numeric(motion[1, c("dx", "dy", "dz")])),
    unname(as.numeric(motion[nrow(motion), c("dx", "dy", "dz")]))
  )
  expect_equal(
    unname(as.numeric(motion[1, c("upx", "upy", "upz")])),
    unname(as.numeric(motion[nrow(motion), c("upx", "upy", "upz")]))
  )
  expect_lt(periodic_derivative_difference(forward), 1e-5)
  expect_lt(periodic_derivative_difference(up), 1e-5)
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

test_that("linear camera motion interpolates lookats directly", {
  motion = generate_camera_motion(
    positions = list(c(0, 0, 0), c(0, 0, 0)),
    lookats = list(c(1, 0, 0), c(0, 0, 1)),
    camera_ups = list(c(0, 1, 0), c(0, 1, 0)),
    type = "linear",
    smooth_orientation = FALSE,
    frames = 3,
    progress = FALSE
  )

  middle_lookat = unname(as.numeric(motion[2, c("dx", "dy", "dz")]))
  expect_equal(middle_lookat, c(0.5, 0, 0.5))
})

test_that("linear camera motion interpolates up vectors directly", {
  motion = generate_camera_motion(
    positions = list(c(0, 0, 0), c(0, 0, 0)),
    lookats = list(c(0, 0, 1), c(0, 0, 1)),
    camera_ups = list(c(0, 1, 0), c(0, 0, 1)),
    type = "linear",
    smooth_orientation = FALSE,
    frames = 3,
    progress = FALSE
  )

  middle_up = unname(as.numeric(motion[2, c("upx", "upy", "upz")]))
  expect_equal(middle_up, c(0, 0.5, 0.5))
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
