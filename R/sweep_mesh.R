#' Validate a scalar sweep argument
#' @keywords internal
#' @noRd
sweep_scalar = function(x, name, lower = -Inf, upper = Inf, integer = FALSE) {
  if (
    !is.numeric(x) ||
      length(x) != 1L ||
      !is.finite(x) ||
      x < lower ||
      x > upper ||
      (integer && x != floor(x))
  ) {
    stop(
      sprintf(
        "`%s` must be a finite %s in [%s, %s].",
        name,
        if (integer) "integer" else "number",
        lower,
        upper
      ),
      call. = FALSE
    )
  }
  x
}

#' Normalize a sweep profile
#' @keywords internal
#' @noRd
sweep_profile = function(polygon, name, added = 0) {
  p = tryCatch(grDevices::xy.coords(polygon), error = function(e) NULL)
  if (is.null(p)) {
    stop(
      sprintf("`%s` must specify numeric x/y coordinates.", name),
      call. = FALSE
    )
  }
  p = cbind(p$x, p$y)
  if (nrow(p) < 3L || any(!is.finite(p))) {
    stop(
      sprintf("`%s` must have at least three finite vertices.", name),
      call. = FALSE
    )
  }
  if (all(p[1, ] == p[nrow(p), ])) {
    p = p[-nrow(p), , drop = FALSE]
  }
  next_i = c(seq_len(nrow(p))[-1L], 1L)
  lengths = sqrt(rowSums((p[next_i, , drop = FALSE] - p)^2))
  centered = sweep(p, 2, colMeans(p))
  area = sum(
    centered[, 1] * centered[next_i, 2] - centered[next_i, 1] * centered[, 2]
  )
  if (
    nrow(p) < 3L ||
      any(lengths == 0) ||
      !is.finite(area) ||
      abs(area) <= .Machine$double.eps * sum(lengths)^2
  ) {
    stop(
      sprintf(
        "`%s` must have nonzero area and distinct consecutive vertices.",
        name
      ),
      call. = FALSE
    )
  }
  # Preserve the first vertex so reversing winding preserves correspondence.
  if (area > 0) {
    p = p[c(1L, nrow(p):2L), , drop = FALSE]
  }
  if (added > 0) {
    next_i = c(seq_len(nrow(p))[-1L], 1L)
    p = do.call(
      rbind,
      lapply(seq_len(nrow(p)), function(i) {
        t = (0:added) / (added + 1)
        (1 - t) %o% p[i, ] + t %o% p[next_i[i], ]
      })
    )
  }
  unname(p)
}

#' Prepare sweep control points
#' @keywords internal
#' @noRd
sweep_controls = function(
  points,
  closed,
  closed_smooth,
  straight,
  precomputed
) {
  if (precomputed) {
    if (
      !is.list(points) ||
        !length(points) ||
        any(vapply(
          points,
          function(p) {
            !is.matrix(p) ||
              !is.numeric(p) ||
              !identical(dim(p), c(4L, 3L)) ||
              any(!is.finite(p))
          },
          TRUE
        ))
    ) {
      stop(
        "Precomputed `points` must be a list of finite numeric 4-by-3 matrices.",
        call. = FALSE
      )
    }
    cps = points
    if (length(cps) > 1L) {
      for (i in seq_len(length(cps) - 1L)) {
        if (
          !isTRUE(all.equal(
            cps[[i]][4, ],
            cps[[i + 1L]][1, ],
            check.attributes = FALSE
          ))
        ) {
          stop(
            "Precomputed Bezier segments must meet at their endpoints.",
            call. = FALSE
          )
        }
      }
    }
    if (
      closed &&
        !isTRUE(all.equal(
          cps[[1]][1, ],
          cps[[length(cps)]][4, ],
          check.attributes = FALSE
        ))
    ) {
      stop(
        "Precomputed closed paths must return to their starting point.",
        call. = FALSE
      )
    }
  } else {
    if (is.data.frame(points)) {
      points = as.matrix(points)
    }
    if (is.list(points)) {
      if (
        !length(points) ||
          any(vapply(
            points,
            function(p) !is.numeric(p) || length(p) != 3L,
            TRUE
          ))
      ) {
        stop(
          "`points` must contain numeric length-three vectors.",
          call. = FALSE
        )
      }
      points = do.call(rbind, points)
    }
    if (
      !is.matrix(points) ||
        !is.numeric(points) ||
        ncol(points) != 3L ||
        nrow(points) < 2L ||
        any(!is.finite(points))
    ) {
      stop(
        "`points` must specify at least two finite three-dimensional points.",
        call. = FALSE
      )
    }
    if (
      closed && nrow(points) > 2L && all(points[1, ] == points[nrow(points), ])
    ) {
      points = points[-nrow(points), , drop = FALSE]
    }
    closed = closed && nrow(points) > 2L
    if (
      any(
        rowSums(
          (points[-1, , drop = FALSE] - points[-nrow(points), , drop = FALSE])^2
        ) ==
          0
      )
    ) {
      stop("Consecutive path points must be distinct.", call. = FALSE)
    }
    n = nrow(points)
    if (closed && closed_smooth && !straight) {
      # Periodic cubic interpolation: matching positions, first and second derivatives.
      previous = c(n, seq_len(n - 1L))
      following = c(2:n, 1L)
      system = diag(4, n)
      system[cbind(seq_len(n), previous)] = system[cbind(
        seq_len(n),
        previous
      )] +
        1
      system[cbind(seq_len(n), following)] = system[cbind(
        seq_len(n),
        following
      )] +
        1
      derivatives = solve(
        system,
        3 * (points[following, ] - points[previous, ])
      )
      cps = lapply(seq_len(n), function(i) {
        rbind(
          points[i, ],
          points[i, ] + derivatives[i, ] / 3,
          points[following[i], ] - derivatives[following[i], ] / 3,
          points[following[i], ]
        )
      })
    } else if (straight) {
      if (closed) {
        points = rbind(points, points[1, ])
      }
      cps = calculate_control_points_straight(points)
    } else {
      cps = calculate_control_points(points)
      if (closed) {
        first = cps[[1]]
        last = cps[[length(cps)]]
        cps[[length(cps) + 1L]] = rbind(
          last[4, ],
          2 * last[4, ] - last[3, ],
          2 * first[1, ] - first[2, ],
          first[1, ]
        )
      }
    }
  }
  list(cps = cps, closed = closed)
}

#' Evaluate a sweep curve and tangent
#' @keywords internal
#' @noRd
sweep_curve = function(cps, t, closed) {
  s = min(1, max(0, t)) * length(cps)
  index = min(length(cps), floor(s) + 1L)
  local = s - index + 1
  point = eval_bezier(cps[[index]], local)
  tangent = eval_bezier_deriv(cps[[index]], local)
  unit = function(x) {
    len = sqrt(sum(x^2))
    if (!is.finite(len) || len == 0) {
      stop(
        "The path has a zero or nonfinite tangent (a cusp or stationary point).",
        call. = FALSE
      )
    }
    x / len
  }
  tangent = unit(tangent)
  # At polyline joints, use the tangent bisector. For smooth curves this is unchanged.
  if (abs(s - round(s)) < 1e-12 && (s > 0 && s < length(cps) || closed)) {
    before = if (s == 0) length(cps) else as.integer(round(s))
    after = if (s == length(cps)) 1L else as.integer(round(s)) + 1L
    tangent = unit(
      unit(eval_bezier_deriv(cps[[before]], 1)) +
        unit(eval_bezier_deriv(cps[[after]], 0))
    )
  }
  list(point = point, tangent = tangent)
}

#' Adaptive arc-length table for a sweep
#' @keywords internal
#' @noRd
sweep_arc_table = function(cps, tolerance) {
  size = sum(vapply(
    cps,
    function(p) sum(sqrt(rowSums((p[-1, ] - p[-4, ])^2))),
    0.0
  ))
  if (!is.finite(size) || size == 0) {
    stop("The path must have finite nonzero length.", call. = FALSE)
  }
  samples = list(c(0, cps[[1]][1, ]))
  append_segment = function(p, a, b, depth = 0L) {
    chord = p[4, ] - p[1, ]
    error = max(
      sqrt(sum((p[2, ] - p[1, ] - chord / 3)^2)),
      sqrt(sum((p[3, ] - p[1, ] - 2 * chord / 3)^2))
    )
    if (depth >= 4L && error <= tolerance * size) {
      samples[[length(samples) + 1L]] <<- c(b, p[4, ])
      return(invisible(NULL))
    }
    if (depth >= 20L) {
      stop("Arc-length tolerance is too small for this path.", call. = FALSE)
    }
    a01 = (p[1, ] + p[2, ]) / 2
    a12 = (p[2, ] + p[3, ]) / 2
    a23 = (p[3, ] + p[4, ]) / 2
    b01 = (a01 + a12) / 2
    b12 = (a12 + a23) / 2
    middle = (b01 + b12) / 2
    mid_t = (a + b) / 2
    append_segment(rbind(p[1, ], a01, b01, middle), a, mid_t, depth + 1L)
    append_segment(rbind(middle, b12, a23, p[4, ]), mid_t, b, depth + 1L)
  }
  for (i in seq_along(cps)) {
    append_segment(cps[[i]], (i - 1) / length(cps), i / length(cps))
  }
  samples = do.call(rbind, samples)
  distance = sqrt(rowSums(
    (samples[-1, 2:4, drop = FALSE] -
      samples[-nrow(samples), 2:4, drop = FALSE])^2
  ))
  if (any(!is.finite(distance)) || any(distance == 0)) {
    stop(
      "The path contains a zero-length sampled segment; remove repeated points or cusps.",
      call. = FALSE
    )
  }
  list(t = samples[, 1], distance = c(0, cumsum(distance)) / sum(distance))
}

#' Transport a sweep frame by double reflection
#' @keywords internal
#' @noRd
sweep_transport = function(normal, from, to) {
  delta = to$point - from$point
  length2 = sum(delta^2)
  if (length2 == 0) {
    result = normal
  } else {
    reflected = normal - 2 * sum(delta * normal) / length2 * delta
    reflected_tangent = from$tangent -
      2 * sum(delta * from$tangent) / length2 * delta
    correction = to$tangent - reflected_tangent
    correction2 = sum(correction^2)
    result = if (correction2 < 1e-24) {
      reflected
    } else {
      reflected - 2 * sum(correction * reflected) / correction2 * correction
    }
  }
  result = result - sum(result * to$tangent) * to$tangent
  len = sqrt(sum(result^2))
  if (!is.finite(len) || len < 1e-12) {
    stop("Cannot transport a frame through a path cusp.", call. = FALSE)
  }
  result / len
}

#' Construct sweep width interpolation
#' @keywords internal
#' @noRd
sweep_width = function(width, width_end, ease) {
  if (is.numeric(width) && is.null(dim(width))) {
    if (length(width) == 1L) {
      if (length(width_end) == 1L && is.na(width_end)) {
        width_end = width
      }
      width = c(width, width_end)
    }
    width = list(x = seq(0, 1, length.out = length(width)), y = width)
  }
  values = tryCatch(grDevices::xy.coords(width), error = function(e) NULL)
  if (
    is.null(values) ||
      length(values$x) < 2L ||
      any(!is.finite(c(values$x, values$y))) ||
      any(values$y < 0) ||
      any(values$x < 0 | values$x > 1) ||
      any(diff(values$x) <= 0)
  ) {
    stop(
      "`width` must give nonnegative finite widths at strictly increasing positions in [0, 1].",
      call. = FALSE
    )
  }
  if (values$x[1] > 0) {
    values$x = c(0, values$x)
    values$y = c(values$y[1], values$y)
  }
  if (tail(values$x, 1) < 1) {
    values$x = c(values$x, 1)
    values$y = c(values$y, tail(values$y, 1))
  }
  if (ease == "spline") {
    return(stats::splinefun(values$x, values$y, method = "monoH.FC"))
  }
  function(u) {
    i = pmin(
      length(values$x) - 1L,
      findInterval(u, values$x, all.inside = TRUE)
    )
    t = (u - values$x[i]) / (values$x[i + 1L] - values$x[i])
    t = switch(
      ease,
      linear = t,
      quad = quadInOut(t),
      cubic = cubicInOut(t),
      exp = (expInOut(t) - expInOut(0)) / (expInOut(1) - expInOut(0))
    )
    (1 - t) * values$y[i] + t * values$y[i + 1L]
  }
}

#' Build indexed sweep geometry without scene or material state
#' @keywords internal
#' @noRd
sweep_mesh_data = function(
  points,
  polygon = NA,
  polygon_end = NA,
  breaks = NA,
  closed = FALSE,
  closed_smooth = TRUE,
  polygon_add_points = 0,
  twists = 0,
  texture_repeats = 1,
  straight = FALSE,
  precomputed_control_points = FALSE,
  width = 1,
  width_end = NA,
  width_ease = "spline",
  smooth_normals = FALSE,
  u_min = 0,
  u_max = 1,
  linear_step = FALSE,
  end_caps = c(TRUE, TRUE),
  initial_normal = NULL,
  smooth_angle = 180,
  arc_tolerance = 1e-5
) {
  for (name in c(
    "closed",
    "closed_smooth",
    "straight",
    "precomputed_control_points",
    "smooth_normals",
    "linear_step"
  )) {
    value = get(name)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
    }
  }
  if (!is.logical(end_caps) || length(end_caps) != 2L || anyNA(end_caps)) {
    stop("`end_caps` must contain two logical values.", call. = FALSE)
  }
  sweep_scalar(polygon_add_points, "polygon_add_points", 0, 10000, TRUE)
  sweep_scalar(twists, "twists")
  sweep_scalar(texture_repeats, "texture_repeats", 0)
  sweep_scalar(u_min, "u_min", 0)
  sweep_scalar(u_max, "u_max", 0)
  sweep_scalar(smooth_angle, "smooth_angle", 0, 180)
  sweep_scalar(arc_tolerance, "arc_tolerance", .Machine$double.eps, 1)
  width_ease = match.arg(
    width_ease,
    c("spline", "linear", "quad", "cubic", "exp")
  )
  missing_profile = function(x) {
    is.null(x) || (is.atomic(x) && length(x) == 1L && is.na(x))
  }
  if (missing_profile(polygon)) {
    theta = 2 * pi * (0:29) / 30
    polygon = cbind(sin(theta), cos(theta)) / 2
  }
  p = sweep_profile(polygon, "polygon", polygon_add_points)
  q = if (missing_profile(polygon_end)) {
    p
  } else {
    sweep_profile(polygon_end, "polygon_end", polygon_add_points)
  }
  if (nrow(p) != nrow(q)) {
    stop(
      "`polygon` and `polygon_end` must have the same number of vertices after closing duplicates are removed.",
      call. = FALSE
    )
  }
  controls = sweep_controls(
    points,
    closed,
    closed_smooth,
    straight,
    precomputed_control_points
  )
  cps = controls$cps
  closed = controls$closed
  if (is.atomic(breaks) && length(breaks) == 1L && is.na(breaks)) {
    breaks = 20L * length(cps)
  }
  sweep_scalar(breaks, "breaks", 2, .Machine$integer.max, TRUE)
  width_at = sweep_width(width, width_end, width_ease)
  if (!closed && (u_min > 1 || u_max > 1 || u_max < u_min)) {
    stop("Open paths require 0 <= u_min <= u_max <= 1.", call. = FALSE)
  }
  if (u_min == u_max) {
    return(list())
  }
  span = u_max - u_min
  full = closed && span >= 1
  if (closed) {
    if (span < 0) {
      span = span %% 1
    }
    u_min = u_min %% 1
    if (full) {
      u_min = 0
      span = 1
      end_caps = c(FALSE, FALSE)
    }
    u_max = u_min + span
  }
  if (span == 0) {
    return(list())
  }
  if (full && breaks < 4) {
    stop(
      "Full closed sweeps require at least four `breaks` (three distinct rings).",
      call. = FALSE
    )
  }
  intervals = if (u_max > 1) {
    list(c(u_min, 1), c(0, u_max - 1))
  } else {
    list(c(u_min, u_max))
  }
  arc = sweep_arc_table(cps, arc_tolerance)
  curve_samples = lapply(arc$t, function(t) sweep_curve(cps, t, closed))
  first = curve_samples[[1]]
  if (is.null(initial_normal)) {
    second = eval_bezier_2nd_deriv(cps[[1]], 0)
    normal = cross_prod(first$tangent, second)
    if (sqrt(sum(normal^2)) < 1e-10 * max(1, sqrt(sum(second^2)))) {
      normal = diag(3)[, which.min(abs(first$tangent))]
    }
  } else {
    if (
      !is.numeric(initial_normal) ||
        length(initial_normal) != 3L ||
        any(!is.finite(initial_normal))
    ) {
      stop("`initial_normal` must be a finite three-vector.", call. = FALSE)
    }
    normal = initial_normal
  }
  normal = normal - sum(normal * first$tangent) * first$tangent
  if (sqrt(sum(normal^2)) < 1e-12) {
    stop(
      "`initial_normal` must not be zero or parallel to the initial tangent.",
      call. = FALSE
    )
  }
  frame_normals = vector("list", length(curve_samples))
  frame_normals[[1]] = normal / sqrt(sum(normal^2))
  for (i in seq_len(length(curve_samples) - 1L)) {
    frame_normals[[i + 1L]] = sweep_transport(
      frame_normals[[i]],
      curve_samples[[i]],
      curve_samples[[i + 1L]]
    )
  }
  correction = 0
  if (closed) {
    last = tail(frame_normals, 1)[[1]]
    correction = atan2(
      sum(first$tangent * cross_prod(last, frame_normals[[1]])),
      sum(last * frame_normals[[1]])
    )
  }
  parameter_at = if (linear_step) {
    stats::approxfun(arc$distance, arc$t, rule = 2)
  } else {
    function(u) u
  }
  distance_at = stats::approxfun(arc$t, arc$distance, rule = 2)
  ring_at = function(u) {
    t = parameter_at(u)
    i = findInterval(t, arc$t, all.inside = TRUE)
    curve = sweep_curve(cps, t, closed)
    s = sweep_transport(frame_normals[[i]], curve_samples[[i]], curve)
    # Rotate around the tangent using a signed closure correction.
    angle = correction * distance_at(t)
    s = cos(angle) * s + sin(angle) * cross_prod(curve$tangent, s)
    r = cross_prod(s, curve$tangent)
    angle = u * twists * 2 * pi
    rotation = matrix(
      c(cos(angle), -sin(angle), sin(angle), cos(angle)),
      2,
      2,
      byrow = TRUE
    )
    profile = ((1 - u) * p + u * q) %*% t(rotation)
    basis = cbind(s, r)
    offsets = profile %*% t(basis)
    list(
      vertices = sweep(offsets * width_at(u), 2, curve$point, "+"),
      offsets = offsets,
      profile = profile,
      tangent = curve$tangent,
      width = width_at(u)
    )
  }
  seam_map = NULL
  if (closed && (full || length(intervals) == 2L)) {
    a = ring_at(0)$vertices
    b = ring_at(1)$vertices
    tolerance = 1e-8 * max(1, max(abs(a)), max(abs(b)))
    for (shift in 0:(nrow(p) - 1L)) {
      candidate = (seq_len(nrow(p)) + shift - 1L) %% nrow(p) + 1L
      if (max(abs(a[candidate, ] - b)) < tolerance) {
        seam_map = candidate
        break
      }
    }
    if (is.null(seam_map)) {
      stop(
        "Closed sweep profiles do not meet: widths and profiles must match at the seam, and twists must respect the profile's rotational symmetry.",
        call. = FALSE
      )
    }
  }
  grid = seq(0, 1, length.out = breaks)
  # Include polyline corners, even when the requested sample grid misses them.
  if (straight) {
    knots = seq(0, 1, length.out = length(cps) + 1L)
    grid = sort(unique(c(grid, if (linear_step) distance_at(knots) else knots)))
  }
  result = list()
  for (part in seq_along(intervals)) {
    interval = intervals[[part]]
    samples = sort(unique(c(
      interval,
      grid[grid > interval[1] & grid < interval[2]]
    )))
    rings = lapply(samples, ring_at)
    if (!is.null(seam_map) && tail(samples, 1) == 1) {
      rings[[length(rings)]]$vertices = ring_at(0)$vertices[
        seam_map,
        ,
        drop = FALSE
      ]
    }
    caps = c(
      part == 1L && end_caps[1],
      part == length(intervals) && end_caps[2]
    )
    result[[part]] = sweep_assemble(
      rings,
      samples,
      ring_at,
      smooth_normals,
      smooth_angle,
      texture_repeats,
      caps,
      if (full) seam_map else NULL,
      seam_map
    )
  }
  result
}

#' Assemble indexed sweep rings and caps
#' @keywords internal
#' @noRd
sweep_assemble = function(
  rings,
  samples,
  ring_at,
  smooth,
  smooth_angle,
  texture_repeats,
  caps,
  join_map,
  periodic_map
) {
  n = nrow(rings[[1]]$vertices)
  nr = length(rings)
  next_i = c(2:n, 1L)
  previous = c(n, seq_len(n - 1L))
  vertices = texcoords = normals = vector("list", nr)
  indices = vector("list", nr)
  offset = 0L
  for (i in seq_len(nr)) {
    ring = rings[[i]]
    if (i == nr && !is.null(join_map)) {
      indices[[i]] = indices[[1]][
        if (length(indices[[1]]) == 1L) 1L else join_map
      ]
      vertices[[i]] = matrix(numeric(), ncol = 3)
    } else {
      vertices[[i]] = if (ring$width == 0) {
        ring$vertices[1, , drop = FALSE]
      } else {
        ring$vertices
      }
      indices[[i]] = offset + seq_len(nrow(vertices[[i]]))
      offset = offset + nrow(vertices[[i]])
    }
    lengths = sqrt(rowSums((ring$profile[next_i, ] - ring$profile)^2))
    if (any(lengths == 0)) {
      stop(
        "Interpolated profiles must have distinct consecutive vertices.",
        call. = FALSE
      )
    }
    texcoords[[i]] = cbind(
      c(0, cumsum(lengths)) / sum(lengths),
      samples[i] * texture_repeats
    )
    if (smooth) {
      h = 1e-5
      lo = max(0, samples[i] - h)
      hi = min(1, samples[i] + h)
      a = ring_at(lo)$vertices
      b = ring_at(hi)$vertices
      if (!is.null(periodic_map) && samples[i] == 0) {
        a = ring_at(1 - h)$vertices[order(periodic_map), ]
        lo = -h
      }
      if (!is.null(periodic_map) && samples[i] == 1) {
        b = ring_at(h)$vertices[periodic_map, ]
        hi = 1 + h
      }
      longitudinal = (b - a) / (hi - lo)
      incoming = ring$offsets - ring$offsets[previous, ]
      outgoing = ring$offsets[next_i, ] - ring$offsets
      unit_rows = function(x) {
        sizes = sqrt(rowSums(x^2))
        if (any(!is.finite(sizes)) || any(sizes == 0)) {
          stop(
            "Cannot compute normals at a degenerate sweep surface.",
            call. = FALSE
          )
        }
        x / sizes
      }
      incoming = unit_rows(incoming)
      outgoing = unit_rows(outgoing)
      cross_rows = function(a, b) {
        cbind(
          a[, 2] * b[, 3] - a[, 3] * b[, 2],
          a[, 3] * b[, 1] - a[, 1] * b[, 3],
          a[, 1] * b[, 2] - a[, 2] * b[, 1]
        )
      }
      ni = unit_rows(cross_rows(incoming, longitudinal))
      no = unit_rows(cross_rows(outgoing, longitudinal))
      blend = rowSums(incoming * outgoing) >= cospi(smooth_angle / 180) - 1e-12
      if (any(blend)) {
        averaged = unit_rows(
          ni[blend, , drop = FALSE] + no[blend, , drop = FALSE]
        )
        ni[blend, ] = no[blend, ] = averaged
      }
      normals[[i]] = rbind(ni, no)
    }
  }
  faces = texture_faces = normal_faces = list()
  for (i in seq_len(nr - 1L)) {
    lower = indices[[i]]
    upper = indices[[i + 1L]]
    lu = (i - 1L) * (n + 1L) + seq_len(n)
    uu = i * (n + 1L) + seq_len(n)
    ln = (i - 1L) * (2L * n)
    un = i * (2L * n)
    if (length(lower) == 1L && length(upper) == 1L) {
      next
    }
    if (length(lower) == 1L) {
      f = cbind(lower, upper[next_i], upper)
      tf = cbind(lu, uu + 1L, uu)
      nf = cbind(ln + n + seq_len(n), un + next_i, un + n + seq_len(n))
    } else if (length(upper) == 1L) {
      f = cbind(lower, lower[next_i], upper)
      tf = cbind(lu, lu + 1L, uu)
      nf = cbind(ln + n + seq_len(n), ln + next_i, un + n + seq_len(n))
    } else {
      f = rbind(
        cbind(lower, lower[next_i], upper),
        cbind(upper, lower[next_i], upper[next_i])
      )
      tf = rbind(cbind(lu, lu + 1L, uu), cbind(uu, lu + 1L, uu + 1L))
      nf = rbind(
        cbind(ln + n + seq_len(n), ln + next_i, un + n + seq_len(n)),
        cbind(un + n + seq_len(n), ln + next_i, un + next_i)
      )
    }
    faces[[length(faces) + 1L]] = f
    texture_faces[[length(texture_faces) + 1L]] = tf
    normal_faces[[length(normal_faces) + 1L]] = nf
  }
  surface = NULL
  if (length(faces)) {
    surface = list(
      vertices = unname(do.call(rbind, vertices)),
      indices = do.call(rbind, faces) - 1L,
      texcoords = do.call(rbind, texcoords),
      tex_indices = do.call(rbind, texture_faces) - 1L
    )
    if (smooth) {
      surface$normals = do.call(rbind, normals)
      surface$norm_indices = do.call(rbind, normal_faces) - 1L
    }
  }
  cap_meshes = list()
  for (end in which(caps)) {
    ring = rings[[if (end == 1L) 1L else nr]]
    if (ring$width == 0) {
      next
    }
    triangles = matrix(decido::earcut(ring$profile), ncol = 3L, byrow = TRUE)
    if (nrow(triangles) != n - 2L) {
      stop(
        "Cannot triangulate sweep cap: the profile must be a simple polygon without holes.",
        call. = FALSE
      )
    }
    desired = ring$tangent * if (end == 1L) -1 else 1
    v = ring$vertices[triangles[1, ], ]
    if (sum(cross_prod(v[2, ] - v[1, ], v[3, ] - v[1, ]) * desired) < 0) {
      triangles = triangles[, c(1, 3, 2)]
    }
    uv = sweep(ring$profile, 2, apply(ring$profile, 2, min))
    uv = sweep(uv, 2, apply(uv, 2, max), "/")
    cap_meshes[[length(cap_meshes) + 1L]] = list(
      vertices = ring$vertices,
      indices = triangles - 1L,
      normals = matrix(desired, nrow = 1L),
      norm_indices = matrix(0L, nrow(triangles), 3L),
      texcoords = uv,
      tex_indices = triangles - 1L
    )
  }
  list(surface = surface, caps = cap_meshes)
}
