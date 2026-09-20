test_that("location and UTC edits rebuild descriptions without mutating the scene", {
  original = sky_light(40, -74, as.POSIXct("2026-06-21 16:00:00", tz = "UTC"))
  captured = NULL
  local_mocked_bindings(prepare_scene_infinite_lights = function(lights) {
    captured <<- lights
    list("prepared")
  })
  controls = native_sky_controls(list(original))
  expect_identical(controls$datetime, "2026-06-21 16:00:00")
  result = controls$update(-20, 120, "2026-12-21 03:04:05")
  expect_identical(result$error, "")
  expect_identical(result$lights, list("prepared"))
  expect_identical(captured[[1]]$lat, -20)
  expect_identical(captured[[1]]$long, 120)
  expect_identical(
    format(captured[[1]]$datetime, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    "2026-12-21 03:04:05"
  )
  expect_identical(original$lat, 40)
  expect_identical(controls$latitude, 40)
})

test_that("invalid dates and coordinates do not rebuild or normalize silently", {
  controls = native_sky_controls(list(sky_light(
    0,
    0,
    as.POSIXct("2026-01-01", tz = "UTC")
  )))
  local_mocked_bindings(prepare_scene_infinite_lights = function(...) {
    stop("Unexpected rebuild")
  })
  for (time in c(
    "bad",
    "2026-02-30 00:00:00",
    "2026-01-01 25:00:00",
    "2026-01-01 00:00:00suffix"
  )) {
    expect_match(controls$update(0, 0, time)$error, "valid UTC")
  }
  expect_match(controls$update(91, 0, "2026-01-01 00:00:00")$error, "Latitude")
  expect_null(native_sky_controls(list()))
})

test_that("image sky controls default to Hosek and switch models without changing other lights", {
  sky = sky_light_image(
    40,
    -74,
    as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    resolution = 16,
    moon = FALSE
  )
  other = list(type = "image", name = "unchanged")
  captured = NULL
  local_mocked_bindings(prepare_scene_infinite_lights = function(lights) {
    captured <<- lights
    list("prepared")
  })
  controls = native_sky_controls(list(other, sky))
  expect_identical(controls$model, 0L)
  expect_identical(controls$index, 1L)
  expect_true(is.finite(controls$elevation))
  result = controls$update(20, 10, "2026-06-21 12:00:00", 1L)
  expect_identical(result$error, "")
  expect_identical(captured[[2]]$type, "sky")
  expect_true(captured[[2]]$atmosphere)
  expect_true(captured[[2]]$haze)
  expect_true(captured[[2]]$query_altitude)
  expect_identical(captured[[1]], other)
  expect_true(sky$sky_args$hosek)
  result = controls$update(20, 10, "2026-06-21 12:00:00", 0L)
  expect_identical(result$error, "")
  expect_true(captured[[2]]$sky_args$hosek)
})

test_that("native Prague settings survive temporary Hosek selection", {
  sky = sky_light(
    40,
    -74,
    as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    haze = FALSE,
    altitude = 100,
    moon = FALSE
  )
  captured = NULL
  local_mocked_bindings(prepare_scene_infinite_lights = function(lights) {
    captured <<- lights[[1]]
    list("prepared")
  })
  controls = native_sky_controls(list(sky))
  expect_identical(controls$model, 1L)
  expect_identical(controls$update(40, -74, controls$datetime, 0L)$error, "")
  expect_identical(captured$type, "sky_image")
  expect_true(captured$sky_args$hosek)
  expect_identical(controls$update(40, -74, controls$datetime, 1L)$error, "")
  expect_true(captured$atmosphere)
  expect_false(captured$haze)
  expect_equal(captured$sky_args$altitude, 100)
})

test_that("manual image Sun direction and load failures return reviewable results", {
  sky = sky_light_image(
    0,
    0,
    as.POSIXct("2026-06-21 12:00:00", tz = "UTC"),
    moon = FALSE
  )
  controls = native_sky_controls(list(sky))
  direction = NULL
  local_mocked_bindings(
    native_sky_manual_image = function(sky, elevation, azimuth) {
      direction <<- c(elevation, azimuth)
      list(type = "image", name = sky$name)
    },
    prepare_scene_infinite_lights = identity
  )
  result = controls$update(0, 0, controls$datetime, 0L, 30, 210)
  expect_identical(result$error, "")
  expect_equal(direction, c(30, 210))
  expect_equal(result$elevation, 30)
  expect_equal(result$azimuth, 210)
  local_mocked_bindings(prepare_scene_infinite_lights = function(...) {
    stop("Missing sky data")
  })
  expect_match(
    controls$update(0, 0, controls$datetime, 1L)$error,
    "Missing sky data"
  )
  expect_true(sky$sky_args$hosek)
})

test_that("Prague promotion preserves shared sky settings and Hosek restores its image options", {
  sky = sky_light_image(
    40,
    -74,
    as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    altitude = 120,
    visibility = 50,
    albedo = 0.3,
    intensity = 0.8,
    rotation = 25,
    name = "editor-sky",
    resolution = 32,
    turbidity = 4,
    moon = FALSE,
    sun = FALSE,
    prague_rgb_correction_strength = 0.6,
    environment_light_bake_white = TRUE
  )
  local_mocked_bindings(
    prepare_scene_infinite_lights = identity,
    native_sky_manual_image = function(...) {
      stop("Native Prague must not bake an image")
    }
  )
  controls = native_sky_controls(list(sky))
  result = controls$update(35, -100, controls$datetime, 1L, 20, 210)
  expect_identical(result$error, "")
  native = result$lights[[1]]
  expect_true(native$atmosphere)
  expect_true(native$haze)
  expect_true(native$query_altitude)
  expect_false(native$sun)
  expect_false(native$moon)
  expect_equal(native$lat, 35)
  expect_equal(native$long, -100)
  expect_equal(native$intensity, sky$intensity)
  expect_equal(native$rotation, sky$rotation)
  expect_identical(native$name, sky$name)
  for (field in c(
    "altitude",
    "visibility",
    "albedo",
    "prague_rgb_correction_strength"
  )) {
    expect_equal(native$sky_args[[field]], sky$sky_args[[field]])
  }
  expect_equal(native$sky_args$resolution, 64)
  expect_equal(c(result$elevation, result$azimuth), c(20, 210))
  restored = controls$update(40, -74, controls$datetime, 0L)
  expect_identical(restored$error, "")
  expect_identical(restored$lights[[1]], sky)

  sky$sky_args$hosek = FALSE
  controls = native_sky_controls(list(sky))
  expect_identical(controls$model, 1L)
  result = controls$update(40, -74, controls$datetime)
  expect_identical(result$error, "")
  expect_true(result$lights[[1]]$atmosphere)
})

test_that("Sun updates return the clamped direction used to generate either sky model", {
  sky = sky_light_image(
    0,
    0,
    as.POSIXct("2026-06-21 12:00:00", tz = "UTC"),
    moon = FALSE
  )
  controls = native_sky_controls(list(sky))
  generated = NULL
  local_mocked_bindings(
    native_sky_manual_image = function(sky, elevation, azimuth) {
      generated <<- c(elevation, azimuth)
      list(type = "image")
    },
    prepare_scene_infinite_lights = identity
  )
  for (model in 0:1) {
    for (elevation in c(90, 89.99, 89.9, 45, -10)) {
      result = controls$update(0, 0, controls$datetime, model, elevation, 135)
      expect_identical(result$error, "")
      expect_equal(result$elevation, min(elevation, 89.9))
      expect_equal(result$azimuth, 135)
      if (model == 0L) expect_equal(generated, c(result$elevation, 135))
    }
  }
  expect_match(
    controls$update(0, 0, controls$datetime, 0L, 91, 135)$error,
    "Invalid Sun"
  )
})

test_that("a zenith image request produces the same finite, lit sky as the capped angle", {
  sky = sky_light_image(
    0,
    0,
    as.POSIXct("2026-06-21 12:00:00", tz = "UTC"),
    resolution = 64,
    moon = FALSE
  )
  capped = native_sky_manual_image(sky, 89.9, 135)
  zenith = native_sky_manual_image(sky, 90, 135)
  expect_identical(zenith$filename, capped$filename)
  pixels = rayimage::ray_read_image(zenith$filename)
  expect_true(all(is.finite(pixels)))
  expect_gt(max(pixels[,, 1:3]), 0)
})
