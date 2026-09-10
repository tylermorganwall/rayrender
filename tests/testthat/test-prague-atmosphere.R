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
    list(turbidity = 3),
    list(moon_atmosphere = TRUE),
    list(render_mode = 'infrared')
  )) {
    expect_error(
      do.call(prague_test_light, x),
      names(x)
    )
  }
  time = light$datetime
  expect_error(sky_light(0, 0, time, atmosphere = NA), 'atmosphere')
  expect_error(sky_light(0, 0, time, atmosphere = NULL), 'atmosphere')
  expect_false(sky_light(0, 0, time)$atmosphere)
})

test_that('native lighting switches validate and retain serialized defaults', {
  for (field in c(
    'attenuation',
    'query_altitude',
    'haze_in_volumes',
    'deferred_haze',
    'cache_spectra',
    'transmission_table'
  )) {
    for (value in list(NULL, NA, 1, 'yes', c(TRUE, FALSE))) {
      expect_error(
        do.call(prague_test_light, setNames(list(value), field)),
        field
      )
    }
  }
  expect_error(
    prague_test_light(query_altitude = FALSE),
    'requires query_altitude'
  )
  light = prague_test_light(attenuation = FALSE, query_altitude = FALSE)
  expect_false(light$attenuation)
  expect_false(light$query_altitude)
  old = prague_test_light()
  old$attenuation = old$query_altitude = NULL
  old$haze_in_volumes = NULL
  old$deferred_haze = NULL
  old$cache_spectra = old$transmission_table = old$haze_correction_probability = NULL
  old$transmission_table_max_mb = NULL
  expect_silent(validate_infinite_light(old))
})

test_that('transmission table memory limits accept zero, fractions, and infinity', {
  for (value in list(NULL, NA_real_, NaN, -Inf, -1, TRUE, '512', c(1, 2))) {
    expect_error(
      prague_test_light(transmission_table_max_mb = value),
      'transmission_table_max_mb'
    )
  }
  for (value in c(0, .5, 512, 1024, Inf)) {
    light = prague_test_light(transmission_table_max_mb = value)
    expect_identical(light$transmission_table_max_mb, value)
    expect_identical(unserialize(serialize(light, NULL)), light)
  }
})

test_that('sampled haze correction validates probability and requires deferred transport', {
  default = prague_test_light()
  expect_true(default$deferred_haze)
  expect_equal(default$haze_correction_probability, .5)
  eager = prague_test_light(deferred_haze = FALSE)
  expect_false(eager$deferred_haze)
  expect_equal(eager$haze_correction_probability, 1)
  for (value in list(NULL, NA, NaN, Inf, 0, -1, 1.1, TRUE, '0.5', c(.5, 1))) {
    expect_error(
      prague_test_light(
        deferred_haze = TRUE,
        haze_correction_probability = value
      ),
      'haze_correction_probability'
    )
  }
  expect_error(
    prague_test_light(deferred_haze = FALSE, haze_correction_probability = .5),
    'requires deferred_haze'
  )
  light = prague_test_light(
    deferred_haze = TRUE,
    haze_correction_probability = .25
  )
  expect_equal(light$haze_correction_probability, .25)
  expect_identical(unserialize(serialize(light, NULL)), light)
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
    altitude = 250,
    visibility = 40,
    meters_per_unit = 5,
    atmosphere_origin = c(0, -50, 0)
  )
  result = prepare_infinite_light(light)
  expect_equal(result$type, 'prague')
  expect_equal(result$altitude, 250)
  expect_equal(result$meters_per_unit, 5)
  expect_equal(result$origin, c(0, -50, 0))
  expect_true(result$haze_in_volumes)
  expect_true(result$deferred_haze)
  expect_equal(result$haze_correction_probability, .5)
  expect_true(result$cache_spectra)
  expect_true(result$transmission_table)
  expect_equal(result$transmission_table_max_mb, 512)
  expect_equal(result$rgb_gain, c(1, .9, 1.1))
  expect_equal(seen$altitude, 250)
  expect_equal(seen$lon, -74)
  scene = sphere() |> add_infinite_light(light)
  info = prepare_scene_list(scene, integrator_type = 'nee')$render_info
  expect_true(info$has_atmosphere)
  expect_true(info$hasbackground)
  light$attenuation = light$query_altitude = FALSE
  light$haze_in_volumes = FALSE
  light$deferred_haze = TRUE
  light$haze_correction_probability = .5
  light$cache_spectra = light$transmission_table = FALSE
  light$transmission_table_max_mb = Inf
  result = prepare_infinite_light(light)
  expect_false(result$attenuation)
  expect_false(result$query_altitude)
  expect_false(result$haze_in_volumes)
  expect_true(result$deferred_haze)
  expect_equal(result$haze_correction_probability, .5)
  expect_false(result$cache_spectra)
  expect_false(result$transmission_table)
  expect_equal(result$transmission_table_max_mb, Inf)
  light$attenuation = light$query_altitude = NULL
  light$haze_in_volumes = NULL
  light$deferred_haze = NULL
  light$haze_correction_probability = light$cache_spectra = light$transmission_table = NULL
  light$transmission_table_max_mb = NULL
  result = prepare_infinite_light(light)
  expect_true(result$attenuation)
  expect_true(result$query_altitude)
  expect_true(result$haze_in_volumes)
  expect_false(result$deferred_haze)
  expect_equal(result$haze_correction_probability, 1)
  expect_true(result$cache_spectra)
  expect_true(result$transmission_table)
  expect_equal(result$transmission_table_max_mb, 512)
})

# These checks need the large external coefficient dataset. Opt in explicitly;
# the image-review script remains the primary validation of the resulting light.
prague_test_description = function() {
  skip_if_not(identical(Sys.getenv('RAYRENDER_PRAGUE_TESTS'), 'true'))
  skip_if_not_installed('skymodelr')
  skip_if_not('get_prague_sky_metadata' %in% getNamespaceExports('skymodelr'))
  prepare_infinite_light(prague_test_light(resolution = 16))
}

test_that('native exact cache controls preserve spectra when switched between queries', {
  d = prague_test_description()
  positions = rbind(c(0, 0, 0), c(0, 100, 0), c(0, 5000, 0))
  directions = rbind(c(1, 0, 0), c(1, .01, 0), c(1, -.03, 0))
  distance = c(.01, 100, 10000)
  d$cache_spectra = d$transmission_table = FALSE
  original = query_prague_atmosphere(d, positions, directions, distance)
  for (cache in c(TRUE, FALSE)) {
    for (table in c(TRUE, FALSE)) {
      d$cache_spectra = cache
      d$transmission_table = table
      actual = query_prague_atmosphere(d, positions, directions, distance)
      expect_identical(actual, original)
    }
  }
  d$cache_spectra = d$transmission_table = TRUE
  for (limit in c(0, 1, 256, Inf, 0)) {
    d$transmission_table_max_mb = limit
    expect_identical(
      query_prague_atmosphere(d, positions, directions, distance),
      original
    )
  }
  for (limit in list(-1, NA_real_, NaN, -Inf, TRUE, c(1, 2))) {
    d$transmission_table_max_mb = limit
    expect_error(
      query_prague_atmosphere(d, positions, directions, distance),
      'transmission_table_max_mb'
    )
  }
  d$transmission_table_max_mb = NULL
  d$haze_correction_probability = .5
  d$deferred_haze = FALSE
  expect_error(
    query_prague_atmosphere(d, positions, directions, distance),
    'requires deferred_haze'
  )
  d$deferred_haze = TRUE
  # Diagnostic segment queries remain deterministic even when rendering samples
  # the horizon correction. The native C++ test checks both sampled outcomes.
  expect_identical(
    query_prague_atmosphere(d, positions, directions, distance),
    original
  )
})

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

test_that('lighting without finite haze retains altitude-dependent sunset', {
  d = prague_test_description()
  d$render_mode = 'sun'
  a = d$azimuth * pi / 180
  e = d$elevation * pi / 180
  w = c(-sin(a) * cos(e), sin(e), cos(a) * cos(e))
  p = rbind(c(0, 0, 0), c(0, 5000, 0))
  directions = rbind(w, w)
  reference = query_prague_atmosphere(d, p, directions, c(1000, 1000))
  d$attenuation = FALSE
  q = query_prague_atmosphere(d, p, directions, c(1000, 1000))
  expect_equal(q$radiance, reference$radiance)
  expect_equal(q$radiance[1, ], rep(0, 3))
  expect_gt(sum(q$radiance[2, ]), 1)
  expect_equal(q$transmission, matrix(1, 2, 3))
  expect_equal(q$inscatter, matrix(0, 2, 3))
  # Freezing the observer also freezes the horizon. A high reference observer
  # sees the Sun everywhere; a ground reference observer sees it nowhere.
  d$query_altitude = FALSE
  fixed = query_prague_atmosphere(d, p, directions, c(1000, 1000))
  expect_equal(fixed$radiance, matrix(0, 2, 3))
  d$altitude = 5000
  high = query_prague_atmosphere(d, p, directions, c(1000, 1000))
  expect_equal(high$radiance[1, ], q$radiance[2, ])
  expect_equal(high$radiance[2, ], q$radiance[2, ])
  d$attenuation = TRUE
  expect_error(
    query_prague_atmosphere(d, p, directions, c(0, 0)),
    'requires query_altitude'
  )
})

test_that('fixed-altitude native sampling uses one consistent proposal', {
  d = prague_test_description()
  d$attenuation = d$query_altitude = FALSE
  d$altitude = 700
  d$render_mode = 'atmosphere'
  u = as.matrix(expand.grid(
    u = (seq_len(32) - .5) / 32,
    v = (seq_len(32) - .5) / 32
  ))
  n = nrow(u)
  p = cbind(rep(0, n), rep(c(0, 5000), length.out = n), rep(0, n))
  sampled = query_prague_atmosphere(d, p, u, rep(0, n), sample = TRUE)
  evaluated = query_prague_atmosphere(
    d,
    p,
    sampled$direction,
    rep(0, n),
    build_sampler = TRUE
  )
  expect_equal(sampled$pdf, evaluated$pdf, tolerance = 1e-5)
  expect_true(all(is.finite(sampled$pdf) & sampled$pdf > 0))
  at_reference = p
  at_reference[, 2] = d$altitude
  fixed = query_prague_atmosphere(d, at_reference, sampled$direction, rep(0, n))
  expect_equal(sampled$radiance, fixed$radiance, tolerance = 1e-5)
  z = 1 - 2 * u[, 1]
  radius = sqrt(1 - z^2)
  w = cbind(radius * cos(2 * pi * u[, 2]), z, radius * sin(2 * pi * u[, 2]))
  uniform = query_prague_atmosphere(d, p, w, rep(0, n), build_sampler = TRUE)
  expect_equal(mean(uniform$pdf) * 4 * pi, 1, tolerance = .04)
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
  sky = prague_test_light(altitude = 1200)
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

test_that('lighting-only queries retain their altitude through empty cloud boundaries', {
  prague_test_description()
  light = prague_test_light(
    attenuation = FALSE,
    resolution = 16
  )
  empty = sphere(x = 50000, radius = 1) |> add_infinite_light(light)
  boundary = cube(y = 3000, xwidth = 10000, ywidth = 6000, zwidth = 6000) |>
    set_medium(homogeneous_medium(sigma_s = 0))
  settings = list(
    lookfrom = c(0, 300, 0),
    lookat = c(0, 6000, 3000),
    width = 16,
    height = 12,
    samples = 4,
    sample_method = 'sobol',
    integrator_type = 'nee',
    max_depth = 2,
    min_variance = 0,
    aperture = 0,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = 'raw',
    auto_exposure = FALSE,
    ambient_light = FALSE,
    preview = FALSE,
    interactive = FALSE,
    plot_scene = FALSE,
    progress = FALSE
  )
  set.seed(42)
  reference = do.call(render_scene, c(list(scene = empty), settings))
  set.seed(42)
  enclosed = do.call(
    render_scene,
    c(list(scene = add_object(empty, boundary)), settings)
  )
  expect_equal(enclosed, reference, tolerance = 1e-5)
  set.seed(42)
  transparent = do.call(
    render_scene,
    c(list(scene = empty, transparent_background = TRUE), settings)
  )
  expect_equal(as.numeric(transparent[,, 4]), rep(0, 16 * 12))
})

# Small orthographic images make every ray follow effectively the same physical
# line, allowing the renderer to be checked against explicit air-interval queries.
prague_region_settings = function() {
  list(
    lookfrom = c(0, 200, 0),
    lookat = c(0, 200, 6000),
    fov = 0,
    ortho_dimensions = c(.002, .002),
    width = 4,
    height = 4,
    samples = 4,
    sample_method = 'sobol',
    integrator_type = 'nee',
    max_depth = 1,
    min_variance = 0,
    aperture = 0,
    denoise = FALSE,
    bloom = FALSE,
    clamp_value = Inf,
    tonemap = 'raw',
    auto_exposure = FALSE,
    ambient_light = FALSE,
    preview = FALSE,
    interactive = FALSE,
    plot_scene = FALSE,
    progress = FALSE
  )
}

test_that('outside-only haze composes air intervals around nested empty volumes', {
  prague_test_description()
  # Compare the analytic integral with eager transport, then deferred transport
  # with an unsampled horizon correction. The default is stochastic.
  sky = prague_test_light(
    haze_in_volumes = FALSE,
    deferred_haze = FALSE,
    resolution = 16
  )
  panel = xy_rect(
    y = 200,
    z = 6000,
    xwidth = 20000,
    ywidth = 20000,
    material = light('white', intensity = 1),
    flipped = TRUE
  )
  volume = cube(
    y = 200,
    z = 3000,
    xwidth = 20000,
    ywidth = 20000,
    zwidth = 4000
  ) |>
    set_medium(homogeneous_medium(sigma_s = 0))
  settings = prague_region_settings()
  render = function(scene, render_settings = settings) {
    set.seed(42)
    do.call(
      render_scene,
      c(list(scene = add_infinite_light(scene, sky)), render_settings)
    )
  }
  description = prepare_infinite_light(sky)
  intervals = query_prague_atmosphere(
    description,
    rbind(c(0, 200, 0), c(0, 200, 5000)),
    rbind(c(0, 0, 1), c(0, 0, 1)),
    c(1000, 1000)
  )
  expected = intervals$inscatter[1, ] +
    intervals$transmission[1, ] *
      (intervals$inscatter[2, ] + intervals$transmission[2, ])
  image = render(add_object(panel, volume))
  observed = vapply(1:3, function(i) mean(image[,, i]), numeric(1))
  expect_equal(observed, expected, tolerance = 3e-4)
  sky$deferred_haze = TRUE
  deferred = render(add_object(panel, volume))
  observed = vapply(1:3, function(i) mean(deferred[,, i]), numeric(1))
  expect_equal(observed, expected, tolerance = 3e-4)
  inner = cube(
    y = 200,
    z = 3000,
    xwidth = 10000,
    ywidth = 10000,
    zwidth = 2000
  ) |>
    set_medium(homogeneous_medium(sigma_s = 0))
  nested = render(add_object(add_object(panel, volume), inner))
  expect_equal(nested, image, tolerance = 3e-4)
  settings$lookfrom = c(0, 200, 3000)
  inside = render(add_object(panel, volume), settings)
  observed = vapply(1:3, function(i) mean(inside[,, i]), numeric(1))
  expect_equal(
    observed,
    intervals$inscatter[2, ] + intervals$transmission[2, ],
    tolerance = 3e-4
  )
})

test_that('outside-only background opacity excludes exactly the volume interval', {
  prague_test_description()
  sky = prague_test_light(
    haze_in_volumes = FALSE,
    resolution = 16
  )
  empty = sphere(x = 50000, radius = 1)
  volume = cube(
    y = 200,
    z = 3000,
    xwidth = 20000,
    ywidth = 20000,
    zwidth = 4000
  ) |>
    set_medium(homogeneous_medium(sigma_s = 0))
  settings = prague_region_settings()
  settings$lookat = c(0, 3800, 6000)
  settings$transparent_background = TRUE
  w = c(0, .6, 1) / sqrt(1 + .6^2)
  queries = query_prague_atmosphere(
    prepare_infinite_light(sky),
    rbind(c(0, 200, 0), c(0, 3200, 5000)),
    rbind(w, w),
    c(1000 / w[3], Inf)
  )
  expected_alpha = 1 -
    mean(queries$transmission[1, ] * queries$transmission[2, ])
  set.seed(42)
  image = do.call(
    render_scene,
    c(
      list(scene = add_object(empty, volume) |> add_infinite_light(sky)),
      settings
    )
  )
  expect_equal(mean(image[,, 4]), expected_alpha, tolerance = 3e-4)
  # Without a volume boundary both modes have the same sky, radiance, and alpha.
  set.seed(42)
  outside = do.call(
    render_scene,
    c(list(scene = add_infinite_light(empty, sky)), settings)
  )
  sky$haze_in_volumes = TRUE
  set.seed(42)
  all_air = do.call(
    render_scene,
    c(list(scene = add_infinite_light(empty, sky)), settings)
  )
  expect_equal(outside, all_air, tolerance = 1e-5)
})
