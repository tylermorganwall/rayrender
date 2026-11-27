test_that("uniform cubic B-splines preserve positions, joins and widths", {
  filename = tempfile(fileext = '.pbrt')
  on.exit(unlink(filename))
  points = rbind(c(0, 0, 0), c(1, 2, 0), c(2, -1, 1), c(4, 3, 0), c(5, 0, 2))
  writeLines(
    c(
      'WorldBegin Shape "curve" "string basis" "bspline"',
      paste('"point3 P" [', paste(t(points), collapse = ' '), ']'),
      '"float width0" .2 "float width1" .8 "integer splitdepth" 1'
    ),
    filename
  )
  scene = read_pbrt(filename)$scene
  expect_equal(nrow(scene), 2)
  spans = lapply(scene$shape_info, function(x) x$shape_properties)
  expect_equal(spans[[1]]$width, .2)
  expect_equal(spans[[1]]$width_end, .5)
  expect_equal(spans[[2]]$width, .5)
  expect_equal(spans[[2]]$width_end, .8)
  expect_equal(spans[[1]]$p4, spans[[2]]$p1)
  expect_equal(spans[[1]]$p4 - spans[[1]]$p3, spans[[2]]$p2 - spans[[2]]$p1)
  for (i in 1:2) {
    control = do.call(rbind, spans[[i]][c('p1', 'p2', 'p3', 'p4')])
    for (u in c(0, .17, .5, .83, 1)) {
      bezier = c((1 - u)^3, 3 * u * (1 - u)^2, 3 * u^2 * (1 - u), u^3)
      bspline = c(
        (1 - u)^3,
        3 * u^3 - 6 * u^2 + 4,
        -3 * u^3 + 3 * u^2 + 3 * u + 1,
        u^3
      ) /
        6
      expect_equal(
        as.vector(bezier %*% control),
        as.vector(bspline %*% points[i + 0:3, ]),
        tolerance = 1e-12
      )
    }
  }
})
