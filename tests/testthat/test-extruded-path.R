sweep_test_line = function() rbind(c(0, 0, 0), c(0, 0, 10))
sweep_test_square = function() {
  rbind(c(-.5, -.5), c(-.5, .5), c(.5, .5), c(.5, -.5))
}
sweep_test_loop = function() {
  rbind(c(0, 0, 0), c(2, 0, 0), c(2, 2, 1), c(0, 2, 0), c(-1, 1, 2))
}
sweep_test_edges = function(mesh) {
  f = mesh$indices
  edges = rbind(f[, c(1, 2)], f[, c(2, 3)], f[, c(3, 1)])
  keys = paste(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  list(
    counts = table(keys),
    orientation = tapply(sign(edges[, 2] - edges[, 1]), keys, sum)
  )
}
sweep_test_crosses = function(mesh) {
  f = mesh$indices + 1L
  a = mesh$vertices[f[, 1], , drop = FALSE]
  b = mesh$vertices[f[, 2], , drop = FALSE] - a
  c = mesh$vertices[f[, 3], , drop = FALSE] - a
  cbind(
    b[, 2] * c[, 3] - b[, 3] * c[, 2],
    b[, 3] * c[, 1] - b[, 1] * c[, 3],
    b[, 1] * c[, 2] - b[, 2] * c[, 1]
  )
}

test_that("closed sweeps weld their seam for both curve modes and handednesses", {
  loop = sweep_test_loop()
  for (points in list(
    loop,
    loop[nrow(loop):1, ],
    sweep(loop, 2, c(-1, 1, 1), "*")
  )) {
    for (smooth in c(TRUE, FALSE)) {
      m = sweep_mesh_data(
        points,
        polygon = sweep_test_square(),
        closed = TRUE,
        closed_smooth = smooth,
        breaks = 40
      )[[1]]$surface
      edges = sweep_test_edges(m)
      expect_true(all(edges$counts == 2L))
      expect_true(all(edges$orientation == 0))
      expect_true(all(rowSums(sweep_test_crosses(m)^2) > 0))
      expect_equal(nrow(m$vertices), 39 * 4)
    }
  }
  for (twists in c(.25, .5, 1, -1)) {
    m = sweep_mesh_data(
      loop,
      polygon = sweep_test_square(),
      closed = TRUE,
      twists = twists,
      breaks = 40
    )[[1]]$surface
    expect_true(all(sweep_test_edges(m)$counts == 2L))
  }
  expect_error(
    sweep_mesh_data(
      loop,
      polygon = sweep_test_square(),
      closed = TRUE,
      twists = .1
    ),
    "rotational symmetry"
  )
  expect_error(
    sweep_mesh_data(loop, closed = TRUE, width_end = 2),
    "must match"
  )
})

test_that("periodic control points have matching first and second derivatives", {
  cps = sweep_controls(sweep_test_loop(), TRUE, TRUE, FALSE, FALSE)$cps
  for (i in seq_along(cps)) {
    next_i = i %% length(cps) + 1L
    expect_equal(cps[[i]][4, ], cps[[next_i]][1, ])
    expect_equal(
      eval_bezier_deriv(cps[[i]], 1),
      eval_bezier_deriv(cps[[next_i]], 0)
    )
    expect_equal(
      eval_bezier_2nd_deriv(cps[[i]], 1),
      eval_bezier_2nd_deriv(cps[[next_i]], 0)
    )
  }
})

test_that("trimming preserves full-path frames and evaluates exact endpoints", {
  args = list(
    points = sweep_test_loop(),
    polygon = sweep_test_square(),
    breaks = 101,
    smooth_normals = TRUE
  )
  full = do.call(sweep_mesh_data, args)[[1]]$surface
  cropped = do.call(sweep_mesh_data, c(args, list(u_min = .4)))[[1]]$surface
  expect_equal(
    cropped$vertices,
    full$vertices[161:nrow(full$vertices), ],
    tolerance = 1e-10
  )
  expect_equal(
    cropped$normals,
    full$normals[321:nrow(full$normals), ],
    tolerance = 1e-8
  )
  m = sweep_mesh_data(
    sweep_test_line(),
    breaks = 11,
    u_min = .23,
    u_max = .61
  )[[1]]$surface
  expect_equal(range(m$vertices[, 3]), c(2.3, 6.1))
  expect_equal(range(m$texcoords[, 2]), c(.23, .61))
})

test_that("the final ring uses its exact twist and morphed profile", {
  p = sweep_test_square()
  m = sweep_mesh_data(
    sweep_test_line(),
    polygon = p,
    twists = .25,
    breaks = 10
  )[[1]]$surface
  a = m$vertices[1, 1:2]
  b = tail(m$vertices, 4)[1, 1:2]
  expect_equal(acos(sum(a * b) / sqrt(sum(a^2) * sum(b^2))) * 180 / pi, 90)
  m = sweep_mesh_data(
    sweep_test_line(),
    polygon = p,
    polygon_end = 2 * p,
    u_max = .6,
    breaks = 11
  )[[1]]$surface
  expect_equal(
    unname(sqrt(rowSums(tail(m$vertices, 4)[, 1:2]^2))),
    rep(sqrt(.5) * 1.6, 4)
  )
})

test_that("profile closure and representation do not change winding", {
  p = sweep(sweep_test_square(), 2, c(0, -10), "+")
  args = list(points = sweep_test_line(), breaks = 3)
  a = do.call(sweep_mesh_data, c(args, list(polygon = p)))[[1]]$surface
  for (profile in list(
    rbind(p, p[1, ]),
    data.frame(p),
    list(x = p[, 1], y = p[, 2]),
    p[c(1, 4, 3, 2), ]
  )) {
    b = do.call(sweep_mesh_data, c(args, list(polygon = profile)))[[1]]$surface
    expect_equal(a$vertices, b$vertices)
    expect_equal(a$indices, b$indices)
  }
  center = colMeans(a$vertices[1:4, ])
  f = a$indices + 1L
  radial = (a$vertices[f[, 1], ] +
    a$vertices[f[, 2], ] +
    a$vertices[f[, 3], ]) /
    3
  radial = sweep(radial, 2, center)
  radial[, 3] = 0
  expect_true(all(rowSums(sweep_test_crosses(a) * radial) > 0))
  expect_no_error(sweep_mesh_data(
    sweep_test_line(),
    polygon = p,
    polygon_end = rbind(p, p[1, ])
  ))
})

test_that("taper normals include radius and remain correct at the final ring", {
  m = sweep_mesh_data(
    sweep_test_line(),
    width = 1,
    width_end = 2,
    smooth_normals = TRUE,
    breaks = 11
  )[[1]]$surface
  expect_equal(
    m$normals[, 3] / sqrt(rowSums(m$normals[, 1:2]^2)),
    rep(-.05, nrow(m$normals)),
    tolerance = 1e-7
  )
  expect_equal(rowSums(m$normals^2), rep(1, nrow(m$normals)), tolerance = 1e-10)
})

test_that("normals include profile morphing and twist derivatives and retain creases", {
  p = sweep_test_square()
  m = sweep_mesh_data(
    sweep_test_line(),
    polygon = p,
    polygon_end = 2 * p,
    twists = .25,
    initial_normal = c(1, 0, 0),
    smooth_normals = TRUE,
    smooth_angle = 30,
    breaks = 11
  )[[1]]$surface
  u = .5
  angle = u * pi / 2
  theta = pi / 2
  x = p[1, 1]
  y = p[1, 2]
  tangent = c(
    cos(angle) *
      x -
      sin(angle) * y +
      (1 + u) * theta * (-sin(angle) * x - cos(angle) * y),
    -sin(angle) *
      x -
      cos(angle) * y -
      (1 + u) * theta * (cos(angle) * x - sin(angle) * y),
    10
  )
  expect_equal(sum(m$normals[45, ] * tangent), 0, tolerance = 1e-7)
  expect_gt(sqrt(sum((m$normals[41, ] - m$normals[45, ])^2)), .5)
})

test_that("wrapped intervals honor cap choices and share all material IDs", {
  for (caps in list(
    c(FALSE, FALSE),
    c(TRUE, FALSE),
    c(FALSE, TRUE),
    c(TRUE, TRUE)
  )) {
    object = extruded_path(
      sweep_test_loop(),
      closed = TRUE,
      breaks = 20,
      u_min = .8,
      u_max = 1.2,
      end_caps = caps
    )
    expect_equal(nrow(object), 1L)
    expect_length(object$shape_info[[1]]$mesh_info[[1]]$shapes, 2L + sum(caps))
    ids = vapply(object$shape_info, function(s) s$material_id, 0.0)
    expect_length(unique(ids), 1)
  }
  object = extruded_path(
    sweep_test_loop(),
    closed = TRUE,
    breaks = 20,
    u_min = .8,
    u_max = 1.2,
    material_caps = diffuse("red")
  )
  ids = vapply(object$shape_info, function(s) s$material_id, 0.0)
  expect_length(ids, 2)
  expect_false(ids[1] == ids[2])
  expect_length(object$shape_info[[1]]$mesh_info[[1]]$shapes, 2)
  expect_length(object$shape_info[[2]]$mesh_info[[1]]$shapes, 2)
})

test_that("caps use trimmed widths, profiles and outward orientation", {
  part = sweep_mesh_data(
    sweep_test_line(),
    width = 0,
    width_end = 1,
    u_min = .2,
    u_max = .6,
    breaks = 11
  )[[1]]
  expect_length(part$caps, 2)
  for (i in 1:2) {
    expect_equal(unique(part$caps[[i]]$vertices[, 3]), c(2, 6)[i])
    expect_true(all(sweep_test_crosses(part$caps[[i]])[, 3] * c(-1, 1)[i] > 0))
  }
  expect_length(sweep_mesh_data(sweep_test_line(), closed = TRUE)[[1]]$caps, 2)
  p = rbind(c(0, 0), c(0, 2), c(1, 1), c(2, 2), c(2, 0))
  expect_no_error(sweep_mesh_data(
    sweep_test_line(),
    polygon = p,
    polygon_add_points = 2
  ))
})

test_that("zero-width tips use triangle fans with no degenerate triangles", {
  part = sweep_mesh_data(
    sweep_test_line(),
    width = 0,
    width_end = 1,
    smooth_normals = TRUE,
    breaks = 11
  )[[1]]
  m = part$surface
  expect_equal(nrow(m$vertices), 301)
  expect_equal(nrow(m$indices), 570)
  expect_true(all(rowSums(sweep_test_crosses(m)^2) > 0))
  expect_length(part$caps, 1)
  expect_true(all(is.finite(m$normals)))
  expect_null(extruded_path(sweep_test_line(), width = 0, breaks = 5))
})

test_that("adaptive arc sampling handles nonuniform speed and polyline corners", {
  cp = list(rbind(c(0, 0, 0), c(0, 0, .2), c(0, 0, .6), c(0, 0, 1)))
  m = sweep_mesh_data(
    cp,
    precomputed_control_points = TRUE,
    linear_step = TRUE,
    breaks = 11
  )[[1]]$surface
  expect_equal(
    unique(m$vertices[, 3]),
    seq(0, 1, length.out = 11),
    tolerance = 2e-5
  )
  expect_no_error(sweep_mesh_data(
    cp,
    precomputed_control_points = TRUE,
    u_min = .2
  ))
  points = rbind(c(0, 0, 0), c(0, 0, 1), c(3, 0, 1))
  m = sweep_mesh_data(points, straight = TRUE, breaks = 4)[[1]]$surface
  centers = t(vapply(
    split(
      seq_len(nrow(m$vertices)),
      rep(seq_len(nrow(m$vertices) / 30), each = 30)
    ),
    function(i) colMeans(m$vertices[i, ]),
    numeric(3)
  ))
  expect_true(any(rowSums((sweep(centers, 2, points[2, ]))^2) < 1e-20))
})

test_that("invalid sweep inputs fail early without unbounded searches", {
  line = sweep_test_line()
  expect_error(
    sweep_mesh_data(
      line,
      polygon = rbind(c(0, 0), c(1, 0), c(2, 0)),
      smooth_normals = TRUE
    ),
    "nonzero area"
  )
  expect_error(
    sweep_mesh_data(rbind(line[1, ], line[1, ], line[2, ])),
    "distinct"
  )
  expect_error(sweep_mesh_data(rbind(c(0, 0, 0), c(0, 0, Inf))), "finite")
  expect_error(
    sweep_mesh_data(list(matrix(0, 3, 3)), precomputed_control_points = TRUE),
    "4-by-3"
  )
  expect_error(sweep_mesh_data(line, width = -1), "nonnegative")
  expect_error(
    sweep_mesh_data(line, width = list(x = c(0, 0, 1), y = c(1, 2, 1))),
    "strictly increasing"
  )
  expect_error(sweep_mesh_data(line, initial_normal = c(0, 0, 1)), "parallel")
  for (breaks in list(1, 3.2, Inf, c(3, 4), "3")) {
    expect_error(sweep_mesh_data(line, breaks = breaks), "breaks")
  }
  expect_error(sweep_mesh_data(line, u_max = Inf), "u_max")
  expect_error(sweep_mesh_data(line, end_caps = TRUE), "two logical")
  expect_error(sweep_mesh_data(line, u_min = .8, u_max = .2), "Open paths")
  expect_length(sweep_mesh_data(line, u_min = .2, u_max = .2), 0)
})

test_that("profile orientation, texture indexing and materials reach the renderer", {
  mesh = sweep_mesh_data(
    sweep_test_line(),
    polygon = sweep_test_square(),
    initial_normal = c(0, 1, 0),
    breaks = 3
  )[[1]]$surface
  expect_equal(as.numeric(mesh$vertices[1, ]), c(-.5, -.5, 0))
  uv = mesh$texcoords[mesh$tex_indices + 1L, 1]
  expect_true(all(
    apply(matrix(uv, ncol = 3), 1, function(x) diff(range(x))) <= .25 + 1e-12
  ))
  mat = diffuse("steelblue")
  scene = extruded_path(
    sweep_test_line(),
    breaks = 5,
    material = mat,
    x = 1,
    angle = c(0, 30, 0),
    scale = 2
  )
  expect_true(all(scene$x == 1))
  expect_equal(scene$material[[1]], mat[[1]])
  set.seed(1)
  image = render_scene(
    scene,
    width = 16,
    height = 16,
    samples = 2,
    lookfrom = c(4, 4, 25),
    lookat = c(1, 0, 10),
    fov = 60,
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_gt(diff(range(image)), 0)
})

test_that("sweep surfaces and caps form one absorbing SSS boundary", {
  scene = extruded_path(
    sweep_test_line(),
    breaks = 5,
    material = subsurface(sigma_a = .1, sigma_s = 0, refraction = 1)
  )
  set.seed(42)
  image = render_scene(
    scene,
    width = 4,
    height = 4,
    samples = 2048,
    sample_method = "random",
    lookfrom = c(0, 0, -2),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.25, .25),
    aperture = 0,
    ambient_light = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white",
    min_variance = 0,
    tonemap = "raw",
    denoise = FALSE,
    bloom = FALSE,
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_equal(mean(image[,, 1:3]), exp(-1), tolerance = .015)
})

test_that("smooth zero-width tips render with sparse normal indices", {
  scene = extruded_path(
    sweep_test_line(),
    width = 0,
    width_end = 2,
    breaks = 8,
    smooth_normals = TRUE
  )
  mesh = scene$shape_info[[1]]$mesh_info[[1]]
  expect_lt(
    length(unique(as.vector(mesh$shapes[[1]]$norm_indices))),
    nrow(mesh$normals[[1]])
  )
  set.seed(2)
  image = render_scene(
    scene,
    width = 8,
    height = 8,
    samples = 2,
    lookfrom = c(4, 2, 12),
    lookat = c(0, 0, 6),
    fov = 60,
    debug_channel = "normals",
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot = FALSE,
    denoise = FALSE,
    bloom = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_gt(diff(range(image)), 0)
})

test_that("wrapped capped sweeps weld exactly across their shared seam", {
  parts = sweep_mesh_data(
    sweep_test_loop(),
    polygon = sweep_test_square(),
    closed = TRUE,
    twists = .25,
    breaks = 40,
    u_min = .8,
    u_max = 1.2
  )
  a = tail(parts[[1]]$surface$vertices, 4)
  b = head(parts[[2]]$surface$vertices, 4)
  key = function(x) {
    apply(x, 1, function(row) paste(sprintf("%.17g", row), collapse = ","))
  }
  expect_setequal(key(a), key(b))
  expect_length(parts[[1]]$caps, 1)
  expect_length(parts[[2]]$caps, 1)
  expect_error(
    sweep_mesh_data(sweep_test_loop(), closed = TRUE, breaks = 2),
    "at least four"
  )
})
