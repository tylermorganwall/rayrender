polygon_test_square = function() rbind(c(0, 0), c(2, 0), c(2, 1), c(0, 1))
polygon_test_close = function(p) rbind(p, p[1, ])

test_that("polygon wrappers preserve rayvertex geometry and rayrender placement", {
  p = polygon_test_square()
  geometry = list(polygon = p, plane = "zx", scale = c(-2, 3, 1), top = -2)
  material = glossy("steelblue")
  expected = do.call(rayvertex::extruded_polygon_mesh, geometry)
  scene = do.call(
    extruded_polygon,
    c(
      geometry,
      list(
        x = 5,
        y = 2,
        z = -3,
        angle = c(10, 20, 30),
        order_rotation = c(3, 1, 2),
        material = material
      )
    )
  )
  expect_identical(scene$shape_info[[1]]$mesh_info[[1]], expected)
  expect_equal(c(scene$x, scene$y, scene$z), c(5, 2, -3))
  expect_equal(scene$transforms[[1]]$angle[[1]], c(10, 20, 30))
  expect_equal(scene$transforms[[1]]$order_rotation[[1]], c(3, 1, 2))
  expect_equal(scene$material[[1]], material[[1]])
  expect_true(scene$shape_info[[1]]$shape_properties$override_material)
  for (value in list(NA, Inf, "one", numeric(), c(1, 2))) {
    for (name in c("x", "y", "z")) {
      expect_error(
        do.call(
          extruded_polygon,
          c(list(polygon = p), setNames(list(value), name))
        ),
        name
      )
    }
  }
})

test_that("fixed extrusions render as closed SSS boundaries", {
  p = polygon_test_square()
  cases = lapply(c("xz", "zx", "xy", "yx", "yz", "zy"), function(plane) {
    list(polygon = p, plane = plane)
  })
  cases = c(
    cases,
    list(
      list(polygon = p[c(1, 2, 2, 3, 4), ]),
      list(polygon = p, scale = c(-1, 1, 1), top = -1),
      list(
        polygon = rbind(p * 5, p + 1, p + cbind(rep(6, 4), rep(2, 4))),
        holes = c(5, 9)
      )
    )
  )
  if (requireNamespace("sf", quietly = TRUE)) {
    close = polygon_test_close
    s = sf::st_sf(
      roof = 2,
      geometry = sf::st_sfc(sf::st_multipolygon(list(
        list(close(p * 5), close(p + 1)),
        list(close(p * 5 + 20), close(p + 21))
      )))
    )
    cases = c(
      cases,
      list(list(polygon = s, top = 0, bottom = 0, data_column_top = "roof"))
    )
  }
  for (args in cases) {
    obj = do.call(
      extruded_polygon,
      c(args, list(material = subsurface(sigma_a = .1, sigma_s = 1)))
    )
    v = obj$shape_info[[1]]$mesh_info[[1]]$vertices[[1]]
    center = colMeans(v)
    extent = max(apply(v, 2, function(x) diff(range(x))))
    set.seed(314159)
    image = render_scene(
      obj,
      width = 6,
      height = 6,
      samples = 2,
      lookfrom = center + extent * c(2, 1.5, 2.5),
      lookat = center,
      fov = 35,
      ambient_light = TRUE,
      backgroundhigh = "white",
      backgroundlow = "white",
      preview = FALSE,
      plot_scene = FALSE,
      parallel = FALSE,
      progress = FALSE,
      tonemap = "raw",
      denoise = FALSE,
      bloom = FALSE
    )
    expect_true(all(is.finite(image)))
  }
})
