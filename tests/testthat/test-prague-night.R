test_that("nighttime Prague radiance is black while transmission remains active", {
  skip_if_not(identical(Sys.getenv("RAYRENDER_PRAGUE_TESTS"), "true"))
  light = sky_light(
    40.7,
    -74,
    as.POSIXct("2026-06-21 00:00:00", tz = "America/New_York"),
    moon = FALSE,
    sun = FALSE,
    sampling_resolution = 16,
    transmission_table_max_mb = 0
  )
  d = prepare_infinite_light(light)
  expect_lt(d$elevation, -4.2)
  positions = rbind(c(0, 0, 0), c(0, 1000, 0), c(0, 15000, 0), c(0, 10, 0))
  directions = rbind(c(0, 1, 0), c(0, 1, 1), c(0, 1, 0), c(0, 1, 0))
  distances = c(1, 100, 1000, Inf)
  twilight = d
  twilight$elevation = -4.2
  reference = query_prague_atmosphere(
    twilight,
    positions,
    directions,
    distances
  )
  expect_gt(max(reference$radiance), 0)
  d$include_sun = TRUE
  for (elevation in c(-4.200001, -30, -90)) {
    d$elevation = elevation
    for (mode in c("all", "atmosphere", "sun")) {
      d$render_mode = mode
      result = query_prague_atmosphere(d, positions, directions, distances)
      expect_true(all(result$radiance == 0))
      expect_true(all(result$inscatter == 0))
      expect_equal(result$transmission, reference$transmission)
      expect_true(all(result$transmission > 0 & result$transmission <= 1))
    }
  }
  # A black light still needs a finite fallback PDF if sampled on its own.
  sampled = query_prague_atmosphere(
    d,
    positions,
    matrix(c(.2, .7), nrow(positions), 2, byrow = TRUE),
    distances,
    sample = TRUE
  )
  expect_true(all(is.finite(sampled$pdf) & sampled$pdf > 0))
  expect_true(all(sampled$radiance == 0))
  below = query_prague_atmosphere(
    d,
    matrix(c(0, 10, 0), 1),
    matrix(c(0, -1, 0), 1),
    Inf
  )
  expect_true(all(below$transmission == 0))
})

test_that("nighttime skies retain independently enabled Moon, stars, and planets", {
  skip_if_not(identical(Sys.getenv("RAYRENDER_PRAGUE_TESTS"), "true"))
  # Known source images isolate each switch from ephemeris visibility, lunar
  # phase, and whether a small render happens to hit a star or planet.
  local_mocked_bindings(
    generate_moon_disk = function(...) {
      list(
        image = array(1, c(16, 16, 3)),
        azimuth_deg = 0,
        elevation_deg = 45,
        angular_diameter_deg = .53,
        projection = "rectilinear",
        atmospheric_attenuation = FALSE
      )
    },
    generate_stars = function(...) array(2, c(16, 32, 4)),
    generate_planets = function(...) array(3, c(16, 32, 4)),
    .package = "skymodelr"
  )
  light = sky_light(
    40.7,
    -74,
    as.POSIXct("2026-06-21 00:00:00", tz = "America/New_York") +
      as.numeric(Sys.time()) / 1e6,
    sun = FALSE,
    moon = FALSE,
    sampling_resolution = 16,
    moon_resolution = 16,
    celestial_resolution = 16,
    transmission_table_max_mb = 0
  )
  render = function(light, lookat = c(0, 11, 1)) {
    render_scene(
      sphere(x = 100) |> add_infinite_light(light),
      width = 4,
      height = 4,
      samples = 1,
      lookfrom = c(0, 10, 0),
      lookat = lookat,
      fov = 0,
      ortho_dimensions = c(.001, .001),
      aperture = 0,
      integrator_type = "nee",
      parallel = FALSE,
      preview = FALSE,
      interactive = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      auto_exposure = FALSE,
      ambient_light = FALSE,
      tonemap = "raw"
    )[,, 1:3]
  }
  expect_true(all(render(light) == 0))
  contributions = list()
  for (body in c("moon", "stars", "planets")) {
    enabled = light
    enabled[[body]] = TRUE
    contributions[[body]] = render(enabled)
    expect_true(all(is.finite(contributions[[body]])))
    expect_gt(min(contributions[[body]]), 0)
  }
  light$moon = light$stars = light$planets = TRUE
  combined = render(light)
  expect_equal(combined, Reduce(`+`, contributions), tolerance = 1e-5)
  light$haze = FALSE
  expect_equal(render(light), combined, tolerance = 1e-5)
  # Earth still blocks the otherwise full-sphere star and planet maps.
  expect_true(all(render(light, lookat = c(0, 9, 1)) == 0))
})
