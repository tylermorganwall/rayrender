sky_celestial_test_light = function(...) {
  sky_light(
    40.7,
    -74,
    as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York"),
    ...
  )
}

test_that("sky celestial controls validate before preparation", {
  light = sky_celestial_test_light()
  expect_true(light$sun)
  expect_true(light$moon)
  expect_false(light$stars)
  expect_false(light$planets)
  expect_identical(unserialize(serialize(light, NULL)), light)
  for (field in c("sun", "moon", "stars", "planets", "earthshine")) {
    for (value in list(NULL, NA, 1, "yes", c(TRUE, FALSE))) {
      expect_error(
        do.call(sky_celestial_test_light, setNames(list(value), field)),
        field
      )
    }
  }
  for (field in c(
    "sun_resolution",
    "moon_resolution",
    "celestial_resolution",
    "number_cores"
  )) {
    for (value in list(NULL, NA, Inf, 0, -1, 16.5, "32", c(32, 32))) {
      expect_error(
        do.call(sky_celestial_test_light, setNames(list(value), field)),
        field
      )
    }
  }
  expect_error(sky_celestial_test_light(star_width = 0), "star_width")
  expect_error(
    sky_celestial_test_light(stars_exposure = 1024),
    "stars_exposure"
  )
  expect_error(
    sky_celestial_test_light(earthshine_albedo = -1),
    "earthshine_albedo"
  )
  expect_error(
    sky_celestial_test_light(solar_irradiance_w_m2 = 0),
    "solar_irradiance"
  )
  expect_error(
    sky_celestial_test_light(sun = FALSE, render_mode = "sun"),
    "requires sun"
  )
})

test_that("automatic disks inherit the sky and explicit disks replace them", {
  seen = list()
  local_mocked_bindings(
    prepare_prague_sky_light = function(light) {
      list(type = "prague", include_sun = light$sun)
    },
    prepare_celestial_light = function(light, atmospheric_attenuation = TRUE) {
      seen[[length(seen) + 1L]] <<- light
      list(
        type = "disk",
        body = light$type,
        name = light$name,
        atmospheric_attenuation = atmospheric_attenuation
      )
    }
  )
  sky = sky_celestial_test_light(
    altitude = 1200,
    visibility = 40,
    albedo = .3,
    rotation = 25,
    intensity = 2,
    name = "evening",
    sun_resolution = 32,
    moon_resolution = 64,
    earthshine = FALSE,
    earthshine_albedo = .2,
    solar_irradiance_w_m2 = 1400,
    number_cores = 2,
    prague_rgb_correction_gain = c(1, .95, .9)
  )
  scene = sphere() |> add_infinite_light(sky)
  result = prepare_scene_infinite_lights(ray_scene_infinite_lights(scene))
  expect_length(result, 3)
  expect_false(result[[1]]$include_sun)
  expect_equal(vapply(seen, function(x) x$type, character(1)), c("sun", "moon"))
  for (disk in seen) {
    expect_identical(disk$datetime, sky$datetime)
    expect_equal(c(disk$lat, disk$long), c(sky$lat, sky$long))
    expect_equal(disk$sky_args$altitude, 1200)
    expect_equal(disk$sky_args$visibility, 40)
    expect_equal(disk$sky_args$albedo, .3)
    expect_equal(disk$sky_args$number_cores, 2)
    expect_equal(disk$sky_args$prague_rgb_correction_gain, c(1, .95, .9))
    expect_equal(disk$rotation, 25)
    expect_equal(disk$intensity, 2)
  }
  expect_equal(seen[[1]]$resolution, 32)
  expect_equal(seen[[2]]$resolution, 64)
  expect_false(seen[[2]]$moon_args$earthshine)
  expect_equal(seen[[2]]$moon_args$earthshine_albedo, .2)
  expect_equal(seen[[2]]$moon_args$solar_irradiance_w_m2, 1400)
  expect_false(result[[2]]$atmospheric_attenuation)
  expect_false(result[[3]]$atmospheric_attenuation)
  expect_named(list_infinite_lights(scene), "evening")
  expect_identical(get_infinite_light(scene, "evening"), sky)
  expect_length(
    ray_scene_infinite_lights(remove_infinite_light(scene, "evening")),
    0
  )

  explicit = moon_light(
    sky$lat,
    sky$long,
    sky$datetime,
    sky_args = list(altitude = 500),
    intensity = .5,
    rotation = 80,
    name = "custom"
  )
  seen = list()
  result = prepare_scene_infinite_lights(list(sky, explicit))
  expect_length(result, 3)
  expect_equal(vapply(seen, function(x) x$type, character(1)), c("moon", "sun"))
  expect_identical(seen[[1]], explicit)

  seen = list()
  sky$sun = sky$moon = FALSE
  result = prepare_scene_infinite_lights(list(sky))
  expect_length(result, 1)
  expect_length(seen, 0)
  expect_false(result[[1]]$include_sun)
  sky$sun = TRUE
  sky$sky_args$render_mode = "atmosphere"
  expect_length(prepare_scene_infinite_lights(list(sky)), 1)
})

test_that("native star and planet maps are unfiltered and cached independently of sky transport", {
  skip_if_not_installed("skymodelr")
  seen = list()
  generator = function(body, ...) {
    seen[[body]] <<- list(...)
    array(if (body == "stars") 2 else 3, c(16, 32, 4))
  }
  local_mocked_bindings(
    generate_stars = function(...) generator("stars", ...),
    generate_planets = function(...) generator("planets", ...),
    generate_sky_latlong = function(...) stop("unexpected atmosphere image"),
    .package = "skymodelr"
  )
  sky = sky_celestial_test_light(
    stars = TRUE,
    planets = TRUE,
    star_width = 2,
    stars_exposure = 1,
    celestial_resolution = 16
  )
  sky$datetime = sky$datetime + as.numeric(Sys.time()) / 1e6
  result = prepare_sky_celestial_background(sky)
  withr::defer(unlink(result$filename))
  for (args in seen) {
    expect_false(args$atmosphere_effects)
    expect_false(args$upper_hemisphere_only)
    expect_identical(args$datetime, sky$datetime)
    expect_equal(args$resolution, 16)
  }
  expect_equal(seen$stars$star_width, 2)
  expect_equal(seen$planets$planet_width, 2)
  pixels = rayimage::ray_read_image(result$filename)
  expect_equal(as.numeric(pixels[,, 1:3]), rep(7, 16 * 32 * 3))
  seen = list()
  sky$intensity = .5
  sky$rotation = 40
  sky$name = "renamed"
  sky$sky_args$altitude = NULL
  sky$sky_args$visibility = 120
  sky$attenuation = FALSE
  cached = prepare_sky_celestial_background(sky)
  expect_length(seen, 0)
  expect_identical(result$filename, cached$filename)
  expect_equal(cached$intensity, .5)
  expect_equal(cached$rotation, 40)
  expect_equal(cached$name, "renamed::background")
})

test_that("image skies expose the same Sun and Moon switches", {
  skip_if_not_installed("skymodelr")
  seen = NULL
  local_mocked_bindings(
    generate_sky_latlong = function(filename, ...) {
      seen <<- list(...)
      writeBin(as.raw(1), filename)
    },
    .package = "skymodelr"
  )
  sky = sky_celestial_test_light()
  time = sky$datetime + as.numeric(Sys.time()) / 1e6
  light = sky_light_image(sky$lat, sky$long, time, sun = FALSE, moon = TRUE)
  result = prepare_infinite_light(light)
  withr::defer(unlink(result$filename))
  expect_equal(seen$render_mode, "atmosphere")
  expect_true(seen$moon)
  expect_null(seen$sun)
  for (value in list(NULL, NA, 1, "yes", c(TRUE, FALSE))) {
    expect_error(sky_light_image(0, 0, time, sun = value), "sun")
    expect_error(sky_light_image(0, 0, time, moon = value), "moon")
  }
})
