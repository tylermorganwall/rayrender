prague_test_light = function(...) {
  sky_light(
    40.7,
    -74,
    as.POSIXct('2026-06-21 20:35:00', tz = 'America/New_York'),
    atmosphere = TRUE,
    ...
  )
}

test_that('atmospheric descriptions validate physical units and supported model settings', {
  light = prague_test_light(
    meters_per_unit = 10,
    atmosphere_origin = c(20, 30, 40)
  )
  expect_true(light$atmosphere)
  expect_equal(light$meters_per_unit, 10)
  expect_identical(unserialize(serialize(light, NULL)), light)
  for (x in list(0, -1, NA_real_, Inf, c(1, 2), 'meters')) {
    expect_error(prague_test_light(meters_per_unit = x), 'meters_per_unit')
  }
  for (x in list(c(1, 2), c(0, NA, 0), Inf, 'origin')) {
    expect_error(prague_test_light(atmosphere_origin = x), 'atmosphere_origin')
  }
  for (x in list(
    list(hosek = TRUE),
    list(moon = TRUE),
    list(stars = TRUE),
    list(visibility = 10),
    list(altitude = 20000),
    list(albedo = 2),
    list(resolution = 16.5),
    list(unknown = 1),
    list(render_mode = 'infrared')
  )) {
    expect_error(
      prague_test_light(sky_args = x),
      'Atmospheric|atmospheric|sky_args'
    )
  }
  time = light$datetime
  expect_error(sky_light(0, 0, time, atmosphere = NA), 'atmosphere')
  expect_error(sky_light(0, 0, time, atmosphere = NULL), 'atmosphere')
  expect_false(sky_light(0, 0, time)$atmosphere)
})

test_that('atmosphere requires nee and is unique in a scene', {
  scene = sphere() |> add_infinite_light(prague_test_light())
  expect_error(prepare_scene_list(scene, integrator_type = 'rtiow'), 'nee')
  expect_error(prepare_scene_list(scene, integrator_type = 'basic'), 'nee')
  duplicate = scene |> add_infinite_light(prague_test_light(name = 'second'))
  expect_error(
    prepare_scene_list(duplicate, integrator_type = 'nee'),
    'only one sky_light'
  )
})

test_that('native preparation uses public metadata without generating an image', {
  skip_if_not_installed('skymodelr')
  skip_if_not('get_prague_sky_metadata' %in% getNamespaceExports('skymodelr'))
  seen = NULL
  testthat::local_mocked_bindings(
    get_prague_sky_metadata = function(...) {
      seen <<- list(...)
      list(
        filename = 'dataset.dat',
        elevation_deg = -1.7,
        azimuth_deg = 303,
        angular_diameter_deg = .53,
        rgb_gain = c(1, .9, 1.1)
      )
    },
    generate_sky_latlong = function(...) stop('unexpected image generation'),
    generate_sun_disk = function(...) stop('unexpected disk generation'),
    .package = 'skymodelr'
  )
  light = prague_test_light(
    sky_args = list(altitude = 250, visibility = 40),
    meters_per_unit = 5,
    atmosphere_origin = c(0, -50, 0)
  )
  result = prepare_infinite_light(light)
  expect_equal(result$type, 'prague')
  expect_equal(result$altitude, 250)
  expect_equal(result$meters_per_unit, 5)
  expect_equal(result$origin, c(0, -50, 0))
  expect_equal(result$rgb_gain, c(1, .9, 1.1))
  expect_equal(seen$altitude, 250)
  expect_equal(seen$lon, -74)
  scene = sphere() |> add_infinite_light(light)
  info = prepare_scene_list(scene, integrator_type = 'nee')$render_info
  expect_true(info$has_atmosphere)
  expect_true(info$hasbackground)
})

# These checks need the large external coefficient dataset. Opt in explicitly;
# the image-review script remains the primary validation of the resulting light.
prague_test_description = function() {
  skip_if_not(identical(Sys.getenv('RAYRENDER_PRAGUE_TESTS'), 'true'))
  skip_if_not_installed('skymodelr')
  skip_if_not('get_prague_sky_metadata' %in% getNamespaceExports('skymodelr'))
  prepare_infinite_light(prague_test_light(sky_args = list(resolution = 16)))
}

test_that('Prague native Sun sees the depressed horizon from high clouds', {
  d = prague_test_description()
  d$render_mode = 'sun'
  a = d$azimuth * pi / 180
  e = d$elevation * pi / 180
  w = c(-sin(a) * cos(e), sin(e), cos(a) * cos(e))
  q = query_prague_atmosphere(
    d,
    rbind(c(0, 0, 0), c(0, 5000, 0)),
    rbind(w, w),
    c(0, 0)
  )
  expect_equal(q$radiance[1, ], c(0, 0, 0))
  expect_gt(sum(q$radiance[2, ]), 1)
  expect_equal(q$transmission, matrix(1, 2, 3))
  expect_equal(q$inscatter, matrix(0, 2, 3))
  # The solar proposal retains its support at both heights, including u=1.
  q = query_prague_atmosphere(
    d,
    rbind(c(0, 5000, 0), c(0, 0, 0)),
    rbind(c(.5, .5), c(1, 1)),
    c(0, 0),
    TRUE
  )
  expect_true(all(is.finite(q$direction)))
  expect_true(all(q$pdf > 0))
})

test_that('short atmospheric segments converge to identity and physical scaling agrees', {
  d = prague_test_description()
  p = matrix(rep(c(0, 200, 0), 4), ncol = 3, byrow = TRUE)
  w = matrix(rep(c(0, .1, 1), 4), ncol = 3, byrow = TRUE)
  distance = c(0, 1e-6, 100, 3000)
  q = query_prague_atmosphere(d, p, w, distance)
  expect_equal(q$transmission[1, ], rep(1, 3))
  expect_equal(q$inscatter[1, ], rep(0, 3))
  expect_lt(max(abs(q$transmission[2, ] - 1)), 1e-6)
  expect_lt(max(abs(q$inscatter[2, ])), 1e-6)
  expect_true(all(q$transmission >= 0 & q$transmission <= 1))
  scaled = d
  scaled$meters_per_unit = 10
  scaled$origin = c(3, 4, 5)
  r = query_prague_atmosphere(
    scaled,
    sweep(p / 10, 2, scaled$origin, '+'),
    w,
    distance / 10
  )
  expect_equal(q$radiance, r$radiance, tolerance = 1e-5)
  expect_equal(q$transmission, r$transmission, tolerance = 1e-5)
  expect_equal(q$inscatter, r$inscatter, tolerance = 1e-5)
})

test_that('atmospheric sky PDFs match sampling at interpolated altitudes', {
  d = prague_test_description()
  d$render_mode = 'atmosphere'
  grid = expand.grid(u = (seq_len(64) - 0.5) / 64, v = (seq_len(64) - 0.5) / 64)
  u = as.matrix(grid)
  n = nrow(u)
  p = cbind(rep(0, n), rep(c(0, 700, 5000, 12000), length.out = n), rep(0, n))
  sample = query_prague_atmosphere(d, p, u, rep(0, n), sample = TRUE)
  evaluated = query_prague_atmosphere(
    d,
    p,
    sample$direction,
    rep(0, n),
    build_sampler = TRUE
  )
  expect_equal(sample$pdf, evaluated$pdf, tolerance = 1e-5)
  expect_true(all(sample$pdf > 0 & is.finite(sample$pdf)))
  # Equal-solid-angle quadrature independently checks PDF normalization.
  z = 1 - 2 * u[, 1]
  radius = sqrt(1 - z^2)
  w = cbind(radius * cos(2 * pi * u[, 2]), z, radius * sin(2 * pi * u[, 2]))
  p[, 2] = 700
  uniform = query_prague_atmosphere(d, p, w, rep(0, n), build_sampler = TRUE)
  expect_equal(mean(uniform$pdf) * 4 * pi, 1, tolerance = .025)
})

test_that('atmospheric preparation requests unattenuated disks and replaces the built-in Sun', {
  seen = list()
  local_mocked_bindings(
    prepare_prague_sky_light = function(light) {
      list(type = 'prague', include_sun = TRUE)
    },
    prepare_celestial_light = function(light, atmospheric_attenuation = TRUE) {
      seen[[light$type]] <<- list(
        attenuation = atmospheric_attenuation,
        altitude = light$sky_args$altitude
      )
      list(type = 'disk', body = light$type)
    }
  )
  sky = prague_test_light(sky_args = list(altitude = 1200))
  time = sky$datetime
  sun = sun_light(40.7, -74, time)
  moon = moon_light(40.7, -74, time)
  result = prepare_scene_infinite_lights(list(sky, sun, moon))
  expect_false(result[[1]]$include_sun)
  expect_false(seen$sun$attenuation)
  expect_false(seen$moon$attenuation)
  expect_equal(seen$sun$altitude, 1200)
  expect_null(sun$sky_args$altitude) # The reusable description is unchanged.
  result = prepare_scene_infinite_lights(list(sky, moon))
  expect_true(result[[1]]$include_sun)
  moon$sky_args$altitude = 500
  prepare_scene_infinite_lights(list(sky, moon))
  expect_equal(seen$moon$altitude, 500)
  prepare_celestial_light(sun)
  expect_true(seen$sun$attenuation)
})
