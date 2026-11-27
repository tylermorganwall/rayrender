test_that("curve subdivision is bounded and preserves one R scene row", {
  curve = bezier_curve()
  expect_equal(nrow(curve), 1)
  expect_identical(curve$shape_info[[1]]$shape_properties$split_depth, 3L)
  expect_identical(
    bezier_curve(split_depth = 0)$shape_info[[1]]$shape_properties$split_depth,
    0L
  )
  for (value in list(-1, 11, 1.5, NA, Inf, c(1, 2), "3")) {
    expect_error(bezier_curve(split_depth = value), "split_depth")
  }
  expect_error(bezier_curve(width = -1))
  expect_error(bezier_curve(width = 0, width_end = 0))
  expect_error(bezier_curve(u_min = -1))
  expect_error(bezier_curve(type = "ribbon", normal = c(0, 0, 0)))
  expect_error(bezier_curve(p4 = c(0, Inf, 0)))
  expect_equal(nrow(bezier_curve(width_end = 0)), 1)
})

test_that("hair rejects singular material parameters", {
  for (value in list(0, -1, 2, Inf, NA, c(.2, .3))) {
    expect_error(hair(beta_m = value), "roughness")
    expect_error(hair(beta_n = value), "roughness")
  }
  expect_error(hair(eta = 1), "eta")
  expect_error(hair(alpha = Inf), "alpha")
  expect_silent(hair(sigma_a = 0, beta_m = .1, beta_n = 1))
})

test_that("PBRT splitdepth and ribbon-chain normals are preserved", {
  path = tempfile(fileext = ".pbrt")
  on.exit(unlink(path))
  writeLines(
    c(
      'WorldBegin',
      'Shape "curve" "string type" "ribbon" "integer splitdepth" 2',
      '"point3 P" [0 0 0  0 1 0  0 2 0  0 3 0  0 4 0  0 5 0  0 6 0]',
      '"normal N" [0 0 1  1 0 1  0 0 1]'
    ),
    path
  )
  scene = read_pbrt(path)$scene
  expect_identical(read_pbrt(path)$render_args$roulette_active_depth, 5L)
  expect_equal(nrow(scene), 2)
  props = lapply(scene$shape_info, function(x) x$shape_properties)
  expect_identical(props[[1]]$split_depth, 2L)
  expect_equal(props[[1]]$normal_end, c(1, 0, 1))
  expect_equal(props[[2]]$normal, c(1, 0, 1))
  expect_equal(props[[2]]$normal_end, c(0, 0, 1))
  writeLines(
    'WorldBegin Shape "curve" "point3 P" [0 0 0 0 1 0 0 2 0 0 3 0]',
    path
  )
  expect_identical(
    read_pbrt(path)$scene$shape_info[[1]]$shape_properties$split_depth,
    3L
  )
})
