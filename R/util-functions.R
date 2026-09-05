#' Calculate Control Points
#'
#' @keywords internal
calculate_control_points = function(s_mat) {
  nr = nrow(s_mat)
  nr2 = nr - 2
  if (nr == 1) {
    stop("Only one point passed, unable to draw curve.")
  }
  if (nr == 2) {
    vec = s_mat[2, ] - s_mat[1, ]
    return(list(matrix(
      c(
        s_mat[1, ],
        s_mat[1, ] + 1 / 3 * vec,
        s_mat[1, ] + 2 / 3 * vec,
        s_mat[2, ]
      ),
      ncol = 3,
      byrow = TRUE
    )))
  }
  spline_matrix = diag(nr2) *
    4 +
    diag(nr2 + 1)[-1, -(nr2 + 1)] +
    t(diag(nr2 + 1)[-1, -(nr2 + 1)])
  inv_spline_matrix = solve(spline_matrix)
  if (nr == 3) {
    new_b = 1 / 4 * (6 * s_mat[2, ] - s_mat[1, ] - s_mat[3, ])
    b_vec = rbind(s_mat[1, ], new_b, s_mat[3, ])
    return_points = list()
    for (i in seq_len(nrow(s_mat) - 1)) {
      vec = b_vec[i + 1, ] - b_vec[i, ]
      return_points[[i]] = matrix(
        c(
          s_mat[i, ],
          b_vec[i, ] + 1 / 3 * vec,
          b_vec[i, ] + 2 / 3 * vec,
          s_mat[i + 1, ]
        ),
        ncol = 3,
        byrow = TRUE
      )
    }
    return(return_points)
  }
  s_vec = matrix(0, nrow = nr - 2, ncol = 3)
  s_vec[1, ] = 6 * s_mat[2, ] - s_mat[1, ]
  for (i in seq_len(nr - 4)) {
    s_vec[i + 1, ] = 6 * s_mat[i + 2, ]
  }
  s_vec[nr - 2, ] = 6 * s_mat[nr - 1, ] - s_mat[nr, ]
  b_vec = inv_spline_matrix %*% s_vec
  b_vec = rbind(s_mat[1, ], b_vec, s_mat[nr, ])
  return_points = list()
  for (i in seq_len(nrow(s_mat) - 1)) {
    vec = b_vec[i + 1, ] - b_vec[i, ]
    return_points[[i]] = matrix(
      c(
        s_mat[i, ],
        b_vec[i, ] + 1 / 3 * vec,
        b_vec[i, ] + 2 / 3 * vec,
        s_mat[i + 1, ]
      ),
      ncol = 3,
      byrow = TRUE
    )
  }
  return(return_points)
}

#' Calculate Control Points (straight)
#'
#' @keywords internal
calculate_control_points_straight = function(s_mat) {
  nr = nrow(s_mat)
  if (nr == 1) {
    stop("Only one point passed, unable to draw curve.")
  }
  return_points = list()
  for (i in seq_len(nrow(s_mat) - 1)) {
    vec = s_mat[i + 1, ] - s_mat[i, ]
    return_points[[i]] = matrix(
      c(
        s_mat[i, ],
        s_mat[i, ] + 1 / 3 * vec,
        s_mat[i, ] + 2 / 3 * vec,
        s_mat[i + 1, ]
      ),
      ncol = 3,
      byrow = TRUE
    )
  }
  return(return_points)
}

#' Clamp Values
#'
#' @keywords internal
clamp = function(v, min = 0, max = Inf) {
  v[v < min] = min
  v[v > max] = max
  v
}

#' Lerp
#'
#' @param t Interpolation distance
#' @param v1 Value 1
#' @param v2 Value 2
#' @return Linearly interpolated value
#'
#' @keywords internal
lerp = function(t, v1, v2) {
  return((1 - t) * v1 + t * v2)
}

#' Quad-in-out
#'
#' @param t Value
#' @return number
#'
#' @keywords internal
quadInOut = function(t) {
  ifelse(t * 2 <= 1, (2 * t)^2 / 2, (2 - (2 * t - 2)^2) / 2)
}

#' Cubic-in-out
#'
#' @param t Value
#' @return number
#'
#' @keywords internal
cubicInOut = function(t) {
  ifelse(t * 2 <= 1, (2 * t)^3 / 2, ((2 * t - 2)^3 + 2) / 2)
}

#' Cubic-in-out
#'
#' @param t Value
#' @return number
#'
#' @keywords internal
expInOut = function(t) {
  ifelse(t * 2 <= 1, 2^(-10 * (1 - 2 * t)) / 2, (2 - 2^(-10 * (2 * t - 1))) / 2)
}

#' Slerp
#'
#' @param vec1 Value
#' @param vec2 Value
#' @param n Value
#' @return number
#'
#' @keywords internal
slerp = function(vec1, vec2, n) {
  t = seq(0, 1, length.out = n + 2)
  subtended_angle = acos(sum(vec1 * vec2))
  return_vecs = list()
  for (i in 1:(n + 2)) {
    return_vecs[[i]] = sin((1 - t[i]) * subtended_angle) /
      sin(subtended_angle) *
      vec1 +
      sin(t[i] * subtended_angle) / sin(subtended_angle) * vec2
  }
  return(return_vecs)
}

#' Tween
#'
#' @param vals Numeric values.
#' @param n Number of frames.
#' @param ease Default `"cubic"`. Interpolation type.
#' @param closed Default `FALSE`. Whether spline interpolation should match speed
#' at the path boundaries.
#' @return number
#'
#' @keywords internal
tween = function(vals, n, ease = "cubic", closed = FALSE) {
  if (length(vals) == 1) {
    return(rep(vals, n))
  }
  if (ease == "spline") {
    return(tween_spline(vals, n, closed))
  }
  len_vals = rep(0, length(vals) - 1)
  free_vals = n - length(vals)
  counter = 1
  for (i in seq_len(free_vals)) {
    len_vals[counter] = len_vals[counter] + 1
    counter = counter + 1
    if (counter > (length(vals) - 1)) {
      counter = 1
    }
  }
  tlist = list()
  for (i in seq_len((length(vals) - 1))) {
    if (ease == "cubic") {
      tlist[[i]] = cubicInOut(seq(0, 1, length.out = len_vals[i] + 2)[
        -(len_vals[i] + 2)
      ])
    } else if (ease == "exp") {
      tlist[[i]] = expInOut(seq(0, 1, length.out = len_vals[i] + 2)[
        -(len_vals[i] + 2)
      ])
    } else if (ease == "quad") {
      tlist[[i]] = quadInOut(seq(0, 1, length.out = len_vals[i] + 2)[
        -(len_vals[i] + 2)
      ])
    } else {
      tlist[[i]] = seq(0, 1, length.out = len_vals[i] + 2)[-(len_vals[i] + 2)]
    }
  }
  final_vals = list()
  counter = 1
  for (i in seq_len(length(tlist))) {
    seg = tlist[[i]]
    for (j in seq_len(length(seg))) {
      final_vals[[counter]] = lerp(seg[j], vals[i], vals[i + 1])
      counter = counter + 1
    }
  }
  final_vals[[counter]] = vals[length(vals)]
  if (length(final_vals) != n) {
    stop(
      "Length of interpolated sequence (",
      length(final_vals),
      ") not equal to frames (",
      n,
      ")"
    )
  }
  return(unlist(final_vals))
}

#' Interpolate Values with Smooth Speed
#'
#' @param vals Numeric values.
#' @param n Number of frames.
#' @param closed Default `FALSE`. Whether to match speed at the path boundaries.
#' @return Numeric vector of interpolated values.
#'
#' @keywords internal
tween_spline = function(vals, n, closed = FALSE) {
  interpolated = tween_spline_path(
    matrix(vals, ncol = 1),
    n = n,
    closed = closed
  )
  return(as.numeric(interpolated[, 1]))
}

#' Interpolate Speed Along a Piecewise Linear Path
#'
#' @param points Matrix of path points.
#' @param n Number of frames.
#' @param closed Default `FALSE`. Whether to match speed at the path boundaries.
#' @return Matrix of interpolated points.
#'
#' @keywords internal
tween_spline_path = function(points, n, closed = FALSE) {
  points = as.matrix(points)
  if (nrow(points) == 1) {
    return(matrix(
      points[1, ],
      nrow = n,
      ncol = ncol(points),
      byrow = TRUE,
      dimnames = list(NULL, colnames(points))
    ))
  }

  len_vals = rep(0, nrow(points) - 1)
  free_vals = n - nrow(points)
  counter = 1
  for (i in seq_len(free_vals)) {
    len_vals[counter] = len_vals[counter] + 1
    counter = counter + 1
    if (counter > (nrow(points) - 1)) {
      counter = 1
    }
  }

  segment_frames = len_vals + 1
  keyframe_frames = c(1, 1 + cumsum(segment_frames))
  segment_delta = points[-1, , drop = FALSE] -
    points[-nrow(points), , drop = FALSE]
  segment_distance = sqrt(rowSums(segment_delta^2))
  segment_speed = segment_distance / segment_frames

  weighted_speed = function(
    previous_speed,
    next_speed,
    previous_frames,
    next_frames
  ) {
    if (previous_speed <= 0 || next_speed <= 0) {
      return(0)
    }
    previous_weight = 2 * next_frames + previous_frames
    next_weight = next_frames + 2 * previous_frames
    (previous_weight + next_weight) /
      (previous_weight / previous_speed + next_weight / next_speed)
  }
  endpoint_speed = function(
    primary_speed,
    adjacent_speed,
    primary_frames,
    adjacent_frames
  ) {
    speed = ((2 * primary_frames + adjacent_frames) *
      primary_speed -
      primary_frames * adjacent_speed) /
      (primary_frames + adjacent_frames)
    if (speed <= 0) {
      return(0)
    }
    if (adjacent_speed == 0 && speed > 3 * primary_speed) {
      return(3 * primary_speed)
    }
    speed
  }

  knot_speed = numeric(nrow(points))
  if (length(segment_speed) == 1) {
    knot_speed[] = segment_speed
  } else {
    for (i in 2:(nrow(points) - 1)) {
      knot_speed[i] = weighted_speed(
        segment_speed[i - 1],
        segment_speed[i],
        segment_frames[i - 1],
        segment_frames[i]
      )
    }
    if (closed) {
      boundary_speed = weighted_speed(
        segment_speed[length(segment_speed)],
        segment_speed[1],
        segment_frames[length(segment_frames)],
        segment_frames[1]
      )
      knot_speed[1] = boundary_speed
      knot_speed[length(knot_speed)] = boundary_speed
    } else {
      knot_speed[1] = endpoint_speed(
        segment_speed[1],
        segment_speed[2],
        segment_frames[1],
        segment_frames[2]
      )
      final_segment = length(segment_speed)
      knot_speed[length(knot_speed)] = endpoint_speed(
        segment_speed[final_segment],
        segment_speed[final_segment - 1],
        segment_frames[final_segment],
        segment_frames[final_segment - 1]
      )
    }
  }

  interpolated = matrix(
    0,
    nrow = n,
    ncol = ncol(points),
    dimnames = list(NULL, colnames(points))
  )
  for (i in seq_along(segment_distance)) {
    frame_rows = seq.int(keyframe_frames[i], keyframe_frames[i + 1])
    time = (frame_rows - keyframe_frames[i]) / segment_frames[i]
    if (segment_distance[i] < sqrt(.Machine$double.eps)) {
      fraction = cubicInOut(time)
    } else {
      time2 = time^2
      time3 = time^3
      start_slope = segment_frames[i] * knot_speed[i] / segment_distance[i]
      end_slope = segment_frames[i] * knot_speed[i + 1] / segment_distance[i]
      fraction = (time3 - 2 * time2 + time) *
        start_slope +
        (-2 * time3 + 3 * time2) +
        (time3 - time2) * end_slope
      fraction = pmax(0, pmin(1, fraction))
    }

    interpolated[frame_rows, ] = matrix(
      points[i, ],
      nrow = length(frame_rows),
      ncol = ncol(points),
      byrow = TRUE
    ) +
      fraction %o% segment_delta[i, ]
  }

  return(interpolated)
}

#' Interpolate Camera Orientation with a Quaternion Spline
#'
#' @param positions Keyframe camera positions.
#' @param lookats Keyframe lookat positions.
#' @param camera_ups Keyframe camera up vectors.
#' @param output_positions Interpolated camera positions.
#' @param closed Default `FALSE`. Whether to use periodic orientation tangents.
#' @return List containing interpolated lookat positions and camera up vectors.
#'
#' @keywords internal
tween_camera_orientation = function(
  positions,
  lookats,
  camera_ups,
  output_positions,
  closed = FALSE
) {
  positions = as.matrix(positions)
  lookats = as.matrix(lookats)
  camera_ups = as.matrix(camera_ups)
  output_positions = as.matrix(output_positions)

  normalize_vector = function(value, fallback) {
    magnitude = sqrt(sum(value^2))
    if (!is.finite(magnitude) || magnitude < sqrt(.Machine$double.eps)) {
      return(fallback)
    }
    return(value / magnitude)
  }
  cross_vector = function(a, b) {
    c(
      a[2] * b[3] - a[3] * b[2],
      a[3] * b[1] - a[1] * b[3],
      a[1] * b[2] - a[2] * b[1]
    )
  }
  quaternion_normalize = function(value) {
    value / sqrt(sum(value^2))
  }
  quaternion_multiply = function(a, b) {
    c(
      a[1] * b[1] - sum(a[2:4] * b[2:4]),
      a[1] * b[2:4] + b[1] * a[2:4] + cross_vector(a[2:4], b[2:4])
    )
  }
  quaternion_inverse = function(value) {
    c(value[1], -value[2:4]) / sum(value^2)
  }
  quaternion_log = function(value) {
    value = quaternion_normalize(value)
    vector_magnitude = sqrt(sum(value[2:4]^2))
    if (vector_magnitude < sqrt(.Machine$double.eps)) {
      return(c(0, 0, 0, 0))
    }
    half_angle = atan2(vector_magnitude, value[1])
    c(0, value[2:4] * half_angle / vector_magnitude)
  }
  quaternion_exp = function(value) {
    vector_magnitude = sqrt(sum(value[2:4]^2))
    if (vector_magnitude < sqrt(.Machine$double.eps)) {
      return(c(1, 0, 0, 0))
    }
    c(
      cos(vector_magnitude),
      value[2:4] * sin(vector_magnitude) / vector_magnitude
    )
  }
  quaternion_slerp = function(a, b, amount) {
    cosine = sum(a * b)
    if (cosine < 0) {
      b = -b
      cosine = -cosine
    }
    cosine = pmax(-1, pmin(1, cosine))
    if (cosine > 0.9995) {
      return(quaternion_normalize((1 - amount) * a + amount * b))
    }
    angle = acos(cosine)
    quaternion_normalize(
      sin((1 - amount) * angle) /
        sin(angle) *
        a +
        sin(amount * angle) / sin(angle) * b
    )
  }
  matrix_to_quaternion = function(value) {
    trace = sum(diag(value))
    if (trace > 0) {
      scale = 2 * sqrt(trace + 1)
      quaternion = c(
        0.25 * scale,
        (value[3, 2] - value[2, 3]) / scale,
        (value[1, 3] - value[3, 1]) / scale,
        (value[2, 1] - value[1, 2]) / scale
      )
    } else if (value[1, 1] > value[2, 2] && value[1, 1] > value[3, 3]) {
      scale = 2 * sqrt(1 + value[1, 1] - value[2, 2] - value[3, 3])
      quaternion = c(
        (value[3, 2] - value[2, 3]) / scale,
        0.25 * scale,
        (value[1, 2] + value[2, 1]) / scale,
        (value[1, 3] + value[3, 1]) / scale
      )
    } else if (value[2, 2] > value[3, 3]) {
      scale = 2 * sqrt(1 + value[2, 2] - value[1, 1] - value[3, 3])
      quaternion = c(
        (value[1, 3] - value[3, 1]) / scale,
        (value[1, 2] + value[2, 1]) / scale,
        0.25 * scale,
        (value[2, 3] + value[3, 2]) / scale
      )
    } else {
      scale = 2 * sqrt(1 + value[3, 3] - value[1, 1] - value[2, 2])
      quaternion = c(
        (value[2, 1] - value[1, 2]) / scale,
        (value[1, 3] + value[3, 1]) / scale,
        (value[2, 3] + value[3, 2]) / scale,
        0.25 * scale
      )
    }
    quaternion_normalize(quaternion)
  }
  quaternion_to_matrix = function(value) {
    value = quaternion_normalize(value)
    w = value[1]
    x = value[2]
    y = value[3]
    z = value[4]
    matrix(
      c(
        1 - 2 * (y^2 + z^2),
        2 * (x * y - z * w),
        2 * (x * z + y * w),
        2 * (x * y + z * w),
        1 - 2 * (x^2 + z^2),
        2 * (y * z - x * w),
        2 * (x * z - y * w),
        2 * (y * z + x * w),
        1 - 2 * (x^2 + y^2)
      ),
      nrow = 3,
      byrow = TRUE
    )
  }
  camera_quaternion = function(position, lookat, camera_up) {
    forward = normalize_vector(lookat - position, c(0, 0, 1))
    up = normalize_vector(camera_up, c(0, 1, 0))
    right = cross_vector(up, forward)
    if (sum(right^2) < 1e-12) {
      up = if (abs(forward[2]) < 0.999) c(0, 1, 0) else c(1, 0, 0)
      right = cross_vector(up, forward)
      if (sum(right^2) < 1e-12) {
        right = cross_vector(c(0, 0, 1), forward)
      }
    }
    right = normalize_vector(right, c(1, 0, 0))
    up = normalize_vector(cross_vector(forward, right), c(0, 1, 0))
    matrix_to_quaternion(cbind(right, up, forward))
  }
  orientation_tangent = function(
    current,
    previous,
    next_value,
    previous_frames,
    next_frames
  ) {
    if (sum(current * previous) < 0) {
      previous = -previous
    }
    if (sum(current * next_value) < 0) {
      next_value = -next_value
    }
    previous_log = quaternion_log(
      quaternion_multiply(quaternion_inverse(current), previous)
    )
    next_log = quaternion_log(
      quaternion_multiply(quaternion_inverse(current), next_value)
    )
    tangent_log = -0.5 *
      (next_frames * previous_log + previous_frames * next_log) /
      (previous_frames + next_frames)
    quaternion_normalize(
      quaternion_multiply(current, quaternion_exp(tangent_log))
    )
  }
  quaternion_squad = function(start, start_tangent, end_tangent, end, amount) {
    quaternion_slerp(
      quaternion_slerp(start, end, amount),
      quaternion_slerp(start_tangent, end_tangent, amount),
      2 * amount * (1 - amount)
    )
  }

  number_keyframes = nrow(positions)
  number_segments = number_keyframes - 1
  number_frames = nrow(output_positions)
  if (number_segments < 1) {
    stop("At least two camera orientation keyframes are required.")
  }
  if (number_frames < number_keyframes) {
    stop("`frames` must be at least the number of camera keyframes.")
  }

  extra_frames = number_frames - number_keyframes
  segment_frames = rep(1, number_segments)
  if (extra_frames > 0) {
    segment_frames = segment_frames +
      tabulate(
        rep(seq_len(number_segments), length.out = extra_frames),
        nbins = number_segments
      )
  }
  keyframe_frames = c(1, 1 + cumsum(segment_frames))

  quaternions = matrix(0, nrow = number_keyframes, ncol = 4)
  for (i in seq_len(number_keyframes)) {
    quaternions[i, ] = camera_quaternion(
      positions[i, ],
      lookats[i, ],
      camera_ups[i, ]
    )
    if (i > 1 && sum(quaternions[i - 1, ] * quaternions[i, ]) < 0) {
      quaternions[i, ] = -quaternions[i, ]
    }
  }

  tangents = quaternions
  if (number_keyframes > 2) {
    for (i in 2:(number_keyframes - 1)) {
      tangents[i, ] = orientation_tangent(
        quaternions[i, ],
        quaternions[i - 1, ],
        quaternions[i + 1, ],
        segment_frames[i - 1],
        segment_frames[i]
      )
    }
  }
  if (closed && number_keyframes > 2) {
    tangents[1, ] = orientation_tangent(
      quaternions[1, ],
      quaternions[number_keyframes - 1, ],
      quaternions[2, ],
      segment_frames[number_segments],
      segment_frames[1]
    )
    tangent_sign = if (
      sum(quaternions[number_keyframes, ] * quaternions[1, ]) < 0
    ) {
      -1
    } else {
      1
    }
    tangents[number_keyframes, ] = tangent_sign * tangents[1, ]
  }

  output_quaternions = matrix(0, nrow = number_frames, ncol = 4)
  for (i in seq_len(number_segments)) {
    frame_rows = seq.int(keyframe_frames[i], keyframe_frames[i + 1])
    amount = (frame_rows - keyframe_frames[i]) / segment_frames[i]
    for (j in seq_along(frame_rows)) {
      output_quaternions[frame_rows[j], ] = quaternion_squad(
        quaternions[i, ],
        tangents[i, ],
        tangents[i + 1, ],
        quaternions[i + 1, ],
        amount[j]
      )
    }
  }

  lookat_distance = sqrt(rowSums((lookats - positions)^2))
  output_lookat_distance = tween_spline(
    lookat_distance,
    n = number_frames,
    closed = closed
  )
  output_lookats = matrix(0, nrow = number_frames, ncol = 3)
  output_camera_ups = matrix(0, nrow = number_frames, ncol = 3)
  for (i in seq_len(number_frames)) {
    camera_basis = quaternion_to_matrix(output_quaternions[i, ])
    output_lookats[i, ] = output_positions[i, ] +
      camera_basis[, 3] * output_lookat_distance[i]
    output_camera_ups[i, ] = camera_basis[, 2]
  }

  output_lookats[keyframe_frames, ] = lookats
  output_camera_ups[keyframe_frames, ] = camera_ups
  list(lookats = output_lookats, camera_ups = output_camera_ups)
}

#' Generate Translation Matrix
#'
#' @param delta Distance
#' @return number
#'
#' @keywords internal
#'
generate_translation_matrix = function(delta) {
  m = matrix(
    c(1, 0, 0, delta[1], 0, 1, 0, delta[2], 0, 0, 1, delta[3], 0, 0, 0, 1),
    4,
    4,
    byrow = T
  )
  return(m)
}

#' Generate Rotation Matrix (order)
#'
#' @param angles Angles
#' @param order_rotation Order of Rotation
#' @return number
#'
#' @keywords internal
#'
generate_rotation_matrix = function(angles, order_rotation) {
  M = diag(4)
  for (i in 1:3) {
    if (order_rotation[i] == 1) {
      if (angles[1] != 0) {
        M = RotateX(angles[1]) %*% M
      }
    }
    if (order_rotation[i] == 2) {
      if (angles[2] != 0) {
        M = RotateY(angles[2]) %*% M
      }
    }
    if (order_rotation[i] == 3) {
      if (angles[3] != 0) {
        M = RotateZ(angles[3]) %*% M
      }
    }
  }
  return(M)
}

#' Generate Rotation Matrix X
#'
#' @param theta Angle
#' @return number
#'
#' @keywords internal
#'
RotateX = function(theta) {
  sinTheta = sinpi(theta / 180)
  cosTheta = cospi(theta / 180)
  M = matrix(
    c(
      1,
      0,
      0,
      0,
      0,
      cosTheta,
      -sinTheta,
      0,
      0,
      sinTheta,
      cosTheta,
      0,
      0,
      0,
      0,
      1
    ),
    4,
    4,
    byrow = T
  )
  return(M)
}

#' Generate Rotation Matrix Y
#'
#' @param theta Angle
#' @return number
#'
#' @keywords internal
#'
RotateY = function(theta) {
  sinTheta = sinpi(theta / 180)
  cosTheta = cospi(theta / 180)
  M = matrix(
    c(
      cosTheta,
      0,
      sinTheta,
      0,
      0,
      1,
      0,
      0,
      -sinTheta,
      0,
      cosTheta,
      0,
      0,
      0,
      0,
      1
    ),
    4,
    4,
    byrow = T
  )
  return(M)
}

#' Generate Rotation Matrix Z
#'
#' @param theta Angle
#' @return number
#'
#' @keywords internal
#'
RotateZ = function(theta) {
  sinTheta = sinpi(theta / 180)
  cosTheta = cospi(theta / 180)
  M = matrix(
    c(
      cosTheta,
      -sinTheta,
      0,
      0,
      sinTheta,
      cosTheta,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      1
    ),
    4,
    4
  )
  return(M)
}

#' Generate Rotation Matrix Axis
#'
#' @param theta Angle
#' @param axis The rotation axis
#' @return matrix
#'
#' @keywords internal
#'
RotateAxis = function(theta, axis) {
  a = axis / sqrt(sum(axis * axis))
  sinTheta = sinpi(theta / 180)
  cosTheta = cospi(theta / 180)
  m = diag(4)
  # Compute rotation of first basis vector
  m[1, 1] = a[1] * a[1] + (1 - a[1] * a[1]) * cosTheta
  m[1, 2] = a[1] * a[2] * (1 - cosTheta) - a[3] * sinTheta
  m[1, 3] = a[1] * a[3] * (1 - cosTheta) + a[2] * sinTheta
  m[1, 4] = 0

  # Compute rotations of second and third basis vectors
  m[2, 1] = a[1] * a[2] * (1 - cosTheta) + a[3] * sinTheta
  m[2, 2] = a[2] * a[2] + (1 - a[2] * a[2]) * cosTheta
  m[2, 3] = a[2] * a[3] * (1 - cosTheta) - a[1] * sinTheta
  m[2, 4] = 0

  m[3, 1] = a[1] * a[3] * (1 - cosTheta) - a[2] * sinTheta
  m[3, 2] = a[2] * a[3] * (1 - cosTheta) + a[1] * sinTheta
  m[3, 3] = a[3] * a[3] + (1 - a[3] * a[3]) * cosTheta
  m[3, 4] = 0
  return(m)
}

#' Calculate Final Angle
#'
#' @keywords internal
calculate_final_twist = function(
  full_control_points,
  breaks,
  t_vals,
  t_vec,
  s_vec,
  r_vec
) {
  r_vec0 = r_vec
  s_vec0 = s_vec
  t_vec0 = t_vec
  for (i in seq_len(breaks - 1)) {
    t_val0 = t_vals[i]
    if (t_val0 < 0) {
      t_val0 = 0
    }
    t_val1 = t_vals[i + 1]
    if (t_val1 < 0) {
      t_val1 = 0
    }

    rot_mat = matrix(c(s_vec, r_vec, t_vec), 3, 3)
    t_temp0 = t_val0 - floor(t_val0)
    if (i != breaks - 1) {
      t_temp1 = t_val1 - floor(t_val1)
    } else {
      t_temp1 = 1
    }

    i0 = floor(t_val0) + 1
    if (i != breaks - 1) {
      i1 = floor(t_val1) + 1
    } else {
      i1 = max(c(1, floor(t_val1 + 1e-8)))
    }

    cp0 = full_control_points[[i0]]
    if (i1 <= length(full_control_points)) {
      cp1 = full_control_points[[i1]]
    } else {
      cp1 = cp0
    }

    x0 = eval_bezier(cp0, t_temp0)
    x1 = eval_bezier(cp1, t_temp1)

    #Evaluate next set of vectors
    v1 = x1 - x0
    c1 = sum(v1 * v1)
    rl = r_vec - (2 / c1) * sum(v1 * r_vec) * v1
    tl = t_vec - (2 / c1) * sum(v1 * t_vec) * v1

    next_deriv = eval_bezier_deriv(cp1, t_temp1)
    t_vec_prev = next_deriv / sqrt(sum(next_deriv * next_deriv))

    v2 = t_vec_prev - tl
    c2 = sum(v2 * v2)
    if (c2 != 0) {
      t_vec = t_vec_prev
      r_vec = rl - (2 / c2) * sum(v2 * rl) * v2
      s_vec = cross_prod(t_vec, r_vec)
    }
  }
  angle_r = acos(sum(r_vec0 * r_vec))
  angle_s = acos(sum(s_vec0 * s_vec))
  angle_t = acos(sum(t_vec0 * t_vec))

  return(c(angle_r, angle_s, angle_t))
}

#' Add Points to Polygon
#'
#' @param polygon Polygon
#' @param added_points Default `0`
#' @return matrix
#'
#' @keywords internal
add_points_polygon = function(polygon, added_points = 0L) {
  existing_verts = nrow(polygon)
  total_verts = existing_verts + (existing_verts - 1L) * added_points
  return_polygon = matrix(0, ncol = 3, nrow = total_verts)
  return_polygon[, 1] = tween(polygon[, 1], n = total_verts, ease = "linear")
  return_polygon[, 2] = tween(polygon[, 2], n = total_verts, ease = "linear")
  return(return_polygon)
}

#' Print time
#'
#' @return Nothing
#' @keywords internal
init_time = function() {
  assign("init_time", proc.time()[3], envir = ray_environment)
  assign("prev_time", proc.time()[3], envir = ray_environment)
}

#' Get time
#'
#' @return Nothing
#' @keywords internal
get_time = function(init = TRUE) {
  if (init) {
    get("init_time", envir = ray_environment)
  } else {
    get("prev_time", envir = ray_environment)
  }
}

#' Print time
#'
#' @return Nothing
#' @keywords internal
print_time = function(verbose = FALSE, message_text = "") {
  if (verbose) {
    time_now = proc.time()[3]
    message(sprintf(
      "%-27s: %0.1f secs (Total: %0.1f secs)",
      message_text,
      time_now - get_time(FALSE),
      time_now - get_time(TRUE)
    ))
    assign("prev_time", time_now, envir = ray_environment)
  }
}
