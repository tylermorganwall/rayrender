polygon_test_square = function() rbind(c(0, 0), c(2, 0), c(2, 1), c(0, 1))
polygon_test_mesh = function(...) {
  extruded_polygon(...)$shape_info[[1]]$mesh_info[[1]]
}
polygon_test_close = function(p) rbind(p, p[1, ])
polygon_test_metrics = function(m) {
  f = m$indices + 1L
  v = m$vertices
  a = v[f[, 1], , drop = FALSE]
  b = v[f[, 2], , drop = FALSE] - a
  c = v[f[, 3], , drop = FALSE] - a
  normals = cbind(
    b[, 2] * c[, 3] - b[, 3] * c[, 2],
    b[, 3] * c[, 1] - b[, 1] * c[, 3],
    b[, 1] * c[, 2] - b[, 2] * c[, 1]
  )
  edges = rbind(f[, c(1, 2)], f[, c(2, 3)], f[, c(3, 1)])
  keys = paste(pmin(edges[, 1], edges[, 2]), pmax(edges[, 1], edges[, 2]))
  list(
    counts = table(keys),
    balance = tapply(sign(edges[, 2] - edges[, 1]), keys, sum),
    areas = sqrt(rowSums(normals^2)) / 2,
    volume = sum(rowSums(sweep(a, 2, v[1, ]) * normals)) / 6,
    normals = normals
  )
}
polygon_test_closed = function(mesh, volume) {
  stats = polygon_test_metrics(mesh)
  expect_true(all(stats$counts == 2L))
  expect_true(all(stats$balance == 0))
  expect_true(all(stats$areas > 0))
  expect_equal(stats$volume, volume, tolerance = 1e-10)
}

test_that("every extrusion plane and reflection has outward caps and walls", {
  p = polygon_test_square()
  for (plane in c("xz", "zx", "xy", "yx", "yz", "zy")) {
    for (height in c(-3, 3)) {
      for (scale in list(1, c(-2, 3, 1), c(-2, -3, 1), c(-2, -3, -1))) {
        m = polygon_test_mesh(p, plane = plane, top = height, scale = scale)
        polygon_test_closed(m, 6 * abs(prod(rep(scale, length.out = 3))))
        expect_equal(nrow(m$vertices), 8L)
        center = colMeans(m$vertices)
        f = m$indices + 1L
        centroids = (m$vertices[f[, 1], ] +
          m$vertices[f[, 2], ] +
          m$vertices[f[, 3], ]) /
          3
        expect_true(all(
          rowSums(
            polygon_test_metrics(m)$normals * sweep(centroids, 2, center)
          ) >
            0
        ))
      }
    }
    expect_equal(
      polygon_test_mesh(p, plane = toupper(plane)),
      polygon_test_mesh(p, plane = plane)
    )
  }
})

test_that("ring representations, winding, duplicates and collinearity preserve the solid", {
  p = polygon_test_square()
  for (ring in list(
    p,
    p[4:1, ],
    polygon_test_close(p),
    p[c(1, 1, 2, 2, 3, 4, 1, 1), ],
    rbind(p[1, ], c(1, 0), p[2:4, ]),
    data.frame(p),
    list(x = p[, 1], y = p[, 2])
  )) {
    for (horizontal in c(FALSE, TRUE)) {
      for (vertical in c(FALSE, TRUE)) {
        m = polygon_test_mesh(
          ring,
          flip_horizontal = horizontal,
          flip_vertical = vertical
        )
        polygon_test_closed(m, 2)
        expect_equal(nrow(m$vertices), 8L)
      }
    }
  }
  concave = rbind(c(0, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 2), c(0, 2))
  polygon_test_closed(polygon_test_mesh(concave), 3)
})

test_that("multiple holes retain caps and inward cavity walls for either ring winding", {
  outer = polygon_test_square() * 5
  hole = polygon_test_square() + 1
  hole2 = hole + cbind(rep(5, 4), rep(1, 4))
  for (closed in c(TRUE, FALSE)) {
    for (reverse in c(TRUE, FALSE)) {
      rings = list(outer, hole, hole2)
      if (reverse) {
        rings = lapply(rings, function(p) p[4:1, ])
      }
      if (closed) {
        rings = lapply(rings, polygon_test_close)
      }
      sizes = vapply(rings, nrow, 0L)
      m = polygon_test_mesh(
        do.call(rbind, rings),
        holes = cumsum(sizes)[1:2] + 1L
      )
      polygon_test_closed(m, 46)
      expect_equal(nrow(m$vertices), 24L)
      # Both caps must cover the outer area minus the two holes.
      f = m$indices + 1L
      on_cap = apply(matrix(m$vertices[f, 2], ncol = 3), 1, function(y) {
        length(unique(y)) == 1L
      })
      expect_equal(sum(polygon_test_metrics(m)$areas[on_cap]), 92)
    }
  }
})

test_that("large translations and centering preserve geometry and orientation", {
  p = polygon_test_square()
  base = polygon_test_mesh(p)
  for (shift in c(1e6, 1e9, -1e9)) {
    m = polygon_test_mesh(p + shift)
    polygon_test_closed(m, 2)
    expect_equal(sweep(m$vertices, 2, c(-shift, 0, shift)), base$vertices)
    expect_equal(m$indices, base$indices)
    expect_equal(
      polygon_test_mesh(p + shift, center = TRUE),
      polygon_test_mesh(p, center = TRUE)
    )
  }
  object = extruded_polygon(p, x = 5, y = 2, z = -3)
  expect_equal(c(object$x, object$y, object$z), c(5, 2, -3))
  expect_equal(object$shape_info[[1]]$mesh_info[[1]]$vertices, base$vertices)
})

test_that("hole indices refer to the input before redundant vertices are removed", {
  p = polygon_test_square()
  rings = list(
    (p * 5)[c(1, 2, 2, 3, 4, 1), ],
    (p + 1)[c(1, 1, 2, 3, 4, 1), ],
    (p + cbind(rep(6, 4), rep(2, 4)))[c(1, 2, 3, 3, 4), ]
  )
  starts = cumsum(vapply(rings, nrow, 0L))[1:2] + 1L
  polygon_test_closed(
    polygon_test_mesh(do.call(rbind, rings), holes = starts),
    46
  )
})

test_that("equal heights emit exactly one nondegenerate cap", {
  for (plane in c("xz", "zx", "xy", "yx", "yz", "zy")) {
    m = polygon_test_mesh(
      polygon_test_square(),
      top = 2,
      bottom = 2,
      plane = plane
    )
    expect_equal(nrow(m$vertices), 4L)
    expect_equal(nrow(m$indices), 2L)
    expect_equal(sum(polygon_test_metrics(m)$areas), 2)
  }
})

test_that("invalid polygon inputs fail before triangulation", {
  p = polygon_test_square()
  for (holes in list(
    c(5, 9, 8),
    c(5, 5),
    5.9,
    c(0, 5),
    -1,
    NA,
    Inf,
    "5",
    integer(),
    99,
    3
  )) {
    expect_error(
      extruded_polygon(rbind(p * 5, p + 1, p + 2), holes = holes),
      "holes"
    )
  }
  for (value in list(NA, Inf, "one", numeric(), c(1, 2))) {
    for (name in c("top", "bottom", "scale_data", "x", "y", "z")) {
      expect_error(
        do.call(
          extruded_polygon,
          c(list(polygon = p), setNames(list(value), name))
        ),
        name
      )
    }
  }
  for (scale in list(0, c(1, 0, 1), c(1, 2), Inf, NA, numeric(), "one")) {
    expect_error(extruded_polygon(p, scale = scale), "scale")
  }
  for (plane in list("x", NA, 1, c("xy", "xz"))) {
    expect_error(extruded_polygon(p, plane = plane), "plane")
  }
  for (flag in c("center", "flip_horizontal", "flip_vertical")) {
    expect_error(
      do.call(extruded_polygon, c(list(polygon = p), setNames(list(NA), flag))),
      flag
    )
  }
  for (ring in list(
    NULL,
    matrix(numeric(), ncol = 2),
    p[1:2, ],
    p[c(1, 1, 1), ],
    rbind(p, c(NA, 0)),
    rbind(p, c(Inf, 0))
  )) {
    expect_error(extruded_polygon(ring), "vertices|coordinates")
  }
  expect_error(extruded_polygon(cbind(0:3, 0)), "area|retrace")
  expect_error(extruded_polygon(p[c(1, 3, 2, 4), ]), "self-intersect")
  # Nonzero signed area does not make a self-intersecting ring valid.
  expect_error(
    extruded_polygon(rbind(c(0, 0), c(4, 0), c(0, 3), c(3, 3), c(1, -1))),
    "self-intersect"
  )
  expect_error(extruded_polygon(rbind(p * 5, p + 20), holes = 5), "inside")
  expect_error(extruded_polygon(rbind(p * 5, p), holes = 5), "touch")
  expect_error(
    extruded_polygon(rbind(p * 5, p + 1, p + 1.5), holes = c(5, 9)),
    "cross"
  )
  expect_error(
    extruded_polygon(rbind(p * 5, p * 2 + 1, p * .5 + 1.5), holes = c(5, 9)),
    "contain"
  )
  expect_error(extruded_polygon(p * 1e-200, scale = 1e-200), "collapsed")
})

test_that("sf preserves multipart holes, heights and XYZ geometry", {
  skip_if_not_installed("sf")
  p = polygon_test_close(polygon_test_square())
  multi = sf::st_multipolygon(list(
    list(p * 5, p + 1),
    list(p * 5 + 20, p + 21)
  ))
  features = sf::st_sf(
    roof = c(2, 3),
    base = c(-1, 0),
    geometry = sf::st_sfc(multi, sf::st_polygon(list(p + 40)))
  )
  m = polygon_test_mesh(
    features,
    data_column_top = "roof",
    data_column_bottom = "base"
  )
  polygon_test_closed(m, 2 * 48 * 3 + 2 * 3)
  for (part in list(c(0, 10, -1, 2), c(20, 30, -1, 2), c(40, 42, 0, 3))) {
    vertices = m$vertices[
      -m$vertices[, 1] >= part[1] & -m$vertices[, 1] <= part[2],
      ,
      drop = FALSE
    ]
    expect_equal(range(vertices[, 2]), part[3:4])
  }
  xyz = sf::st_sf(
    geometry = sf::st_sfc(sf::st_multipolygon(list(list(cbind(p, 100)))))
  )
  expect_equal(polygon_test_mesh(xyz), polygon_test_mesh(p))
  for (g in list(
    sf::st_polygon(),
    sf::st_multipolygon(),
    sf::st_point(c(0, 0))
  )) {
    expect_error(
      extruded_polygon(sf::st_sf(geometry = sf::st_sfc(g))),
      "Empty|vertices|POLYGON"
    )
  }
  expect_error(extruded_polygon(sf::st_sf(geometry = sf::st_sfc())), "Empty")
})

test_that("sf extrusion decisions use each resolved height pair", {
  skip_if_not_installed("sf")
  p = polygon_test_close(polygon_test_square())
  s = sf::st_sf(
    high = c(0, 2, -2),
    low = c(0, 0, 0),
    geometry = sf::st_sfc(lapply(0:2, function(i) {
      sf::st_polygon(list(p + i * 5))
    }))
  )
  a = polygon_test_mesh(s, data_column_top = "high", data_column_bottom = "low")
  b = polygon_test_mesh(
    s,
    top = 0,
    bottom = 0,
    data_column_top = "high",
    data_column_bottom = "low"
  )
  expect_equal(a, b)
  expect_equal(nrow(a$indices), 26L)
  expect_equal(nrow(a$vertices), 20L)
  expect_true(all(polygon_test_metrics(a)$areas > 0))
  flat = polygon_test_mesh(s[1, ], data_column_top = "high")
  expect_equal(nrow(flat$indices), 2L)
  polygon_test_closed(
    polygon_test_mesh(s[2, ], data_column_top = "high", scale_data = -2),
    8
  )
  for (values in list(c(0, NA, 1), c(0, Inf, 1), c("a", "b", "c"))) {
    s$bad = values
    expect_error(
      extruded_polygon(s, data_column_top = "bad"),
      "finite numeric height"
    )
  }
  expect_warning(
    extruded_polygon(s, data_column_bottom = "missing"),
    "data_column_bottom"
  )
  expect_error(
    extruded_polygon(s, data_column_top = c("high", "low")),
    "single column name"
  )
})

test_that("SpatialPolygons work without a data frame or ring-order assumptions", {
  skip_if_not_installed("sf")
  p = polygon_test_close(polygon_test_square())
  s = sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(sf::st_multipolygon(list(
      list(p * 5, p + 1),
      list(p * 5 + 20, p + 21)
    )))
  )
  spatial = sf::as_Spatial(s)
  plain = methods::as(spatial, "SpatialPolygons")
  for (object in list(spatial, plain)) {
    polygon_test_closed(polygon_test_mesh(object), 96)
  }
  # Swap the ring slots without changing their explicit outer/hole status.
  plain@polygons[[1]]@Polygons = plain@polygons[[1]]@Polygons[c(4, 1, 2, 3)]
  polygon_test_closed(polygon_test_mesh(plain), 96)
  # An island inside a larger hole owns its own smaller hole.
  nested = sf::st_sf(
    geometry = sf::st_sfc(sf::st_multipolygon(list(
      list(p * 10, p * 8 + 1),
      list(p * 4 + 2, p + 3)
    )))
  )
  nested_sp = methods::as(sf::as_Spatial(nested), "SpatialPolygons")
  nested_sp@polygons[[1]]@Polygons = nested_sp@polygons[[1]]@Polygons[c(
    4,
    2,
    3,
    1
  )]
  polygon_test_closed(polygon_test_mesh(nested_sp), 200 - 128 + 32 - 2)
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
    v = obj$shape_info[[1]]$mesh_info[[1]]$vertices
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
