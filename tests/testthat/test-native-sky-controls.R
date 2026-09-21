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

test_that("Prague numeric inputs update metadata before preparing lights and reject invalid drafts", {
  sky = sky_light(
    40,
    -74,
    as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    altitude = 100,
    meters_per_unit = 5,
    moon = FALSE
  )
  controls = native_sky_controls(list(sky))
  expect_equal(controls$base_altitude, 100)
  expect_equal(controls$meters_per_unit, 5)
  local_mocked_bindings(prepare_scene_infinite_lights = identity)
  result = controls$update(
    40,
    -74,
    controls$datetime,
    1L,
    base_altitude = 800,
    meters_per_unit = 200
  )
  expect_identical(result$error, "")
  updated = result$lights[[1]]
  expect_equal(updated$sky_args$altitude, 800)
  expect_equal(updated$meters_per_unit, 200)
  expect_equal(result$elevation, native_sky_position(updated)$elevation)
  expect_equal(sky$sky_args$altitude, 100)
  expect_equal(sky$meters_per_unit, 5)
  local_mocked_bindings(prepare_scene_infinite_lights = function(...) {
    stop("Unexpected rebuild")
  })
  for (bad in list(-1, 15001, Inf, NaN, NA, c(1, 2), "10")) {
    expect_match(
      controls$update(40, -74, controls$datetime, base_altitude = bad)$error,
      "Base altitude"
    )
  }
  for (bad in list(0, -1, 1e13, Inf, NaN, NA, c(1, 2), "10")) {
    expect_match(
      controls$update(40, -74, controls$datetime, meters_per_unit = bad)$error,
      "Atmosphere scale"
    )
  }
})

test_that("export replay retains numeric Prague settings and accepts old sky state", {
  sky = sky_light(
    40,
    -74,
    as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    altitude = 100,
    meters_per_unit = 5,
    moon = FALSE
  )
  local_mocked_bindings(prepare_scene_infinite_lights = function(lights) {
    list(list(
      type = "prague",
      altitude = lights[[1]]$sky_args$altitude,
      meters_per_unit = lights[[1]]$meters_per_unit
    ))
  })
  state = list(
    model = 1L,
    latitude = 40,
    longitude = -74,
    datetime = "2026-06-21 16:00:00",
    manual = FALSE,
    elevation = 30,
    azimuth = 135,
    haze = FALSE,
    altitude = TRUE,
    base_altitude = 750,
    meters_per_unit = 300
  )
  restored = native_sky_restore(list(sky), state)
  expect_equal(restored$lights[[1]]$altitude, 750)
  expect_equal(restored$lights[[1]]$meters_per_unit, 300)
  expect_equal(restored$controls$base_altitude, 750)
  expect_equal(restored$controls$meters_per_unit, 300)
  expect_false(restored$lights[[1]]$haze)
  state$base_altitude = state$meters_per_unit = NULL
  restored = native_sky_restore(list(sky), state)
  expect_equal(restored$lights[[1]]$altitude, 100)
  expect_equal(restored$lights[[1]]$meters_per_unit, 5)
  expect_equal(restored$controls$base_altitude, 100)
  expect_equal(restored$controls$meters_per_unit, 5)
})

test_that("sun dragging uses temporary sky resolution and release restores scene quality", {
  sky = sky_light_image(
    40,
    -74,
    as.POSIXct("2026-06-21 16:00:00", tz = "UTC"),
    resolution = 512,
    moon = FALSE
  )
  controls = native_sky_controls(list(sky))
  local_mocked_bindings(
    prepare_scene_infinite_lights = identity,
    native_sky_manual_image = function(sky, elevation, azimuth) {
      list(type = "image", resolution = sky$sky_args$resolution)
    }
  )
  for (model in 0:1) {
    fast = controls$update(
      40,
      -74,
      controls$datetime,
      model,
      30,
      135,
      fast = TRUE
    )
    full = controls$update(40, -74, controls$datetime, model, 30, 135)
    expect_identical(fast$error, "")
    expect_identical(full$error, "")
    if (model == 0) {
      expect_equal(fast$lights[[1]]$resolution, 128)
      expect_equal(full$lights[[1]]$resolution, 512)
    } else {
      expect_equal(fast$lights[[1]]$sky_args$resolution, 16)
      expect_equal(full$lights[[1]]$sky_args$resolution, 64)
    }
    expect_equal(fast$elevation, full$elevation)
    expect_equal(fast$azimuth, full$azimuth)
  }
  expect_equal(sky$sky_args$resolution, 512)
  sky$sky_args$resolution = 16
  controls = native_sky_controls(list(sky))
  result = controls$update(40, -74, controls$datetime, 0L, 30, 135, fast = TRUE)
  expect_equal(result$lights[[1]]$resolution, 16)
})


test_that("Hosek Fast updates keep a separate full-resolution Sun and no baked duplicate", {
  sky = sky_light_image(
    0,
    0,
    as.POSIXct("2026-06-21 12:00:00", tz = "UTC"),
    resolution = 512,
    render_mode = "sun",
    moon = FALSE,
    rotation = 20,
    intensity = 0.7
  )
  controls = native_sky_controls(list(sky))
  for (elevation in c(30, 45, 90)) {
    fast = controls$update(
      0,
      0,
      controls$datetime,
      0L,
      elevation,
      135,
      fast = TRUE
    )
    full = controls$update(0, 0, controls$datetime, 0L, elevation, 135)
    expect_identical(fast$error, "")
    expect_identical(full$error, "")
    expect_equal(fast$elevation, min(89.9, elevation))
    expect_equal(vapply(fast$lights, `[[`, "", "type"), c("image", "disk"))
    sun = fast$lights[[2]]
    expect_identical(sun$filename, full$lights[[2]]$filename)
    expect_equal(sun$rotation, 20)
    expect_equal(sun$intensity, 0.7)
    angle = min(89.9, elevation) * pi / 180
    expect_equal(
      sun$direction,
      c(
        -sinpi(135 / 180) * cos(angle),
        sin(angle),
        cospi(135 / 180) * cos(angle)
      )
    )
    pixels = rayimage::ray_read_image(sun$filename, normalize = FALSE)
    expect_true(all(is.finite(pixels)))
    expect_gt(sum(pixels[,, 1:3]), 0)
    expect_equal(
      max(abs(rayimage::ray_read_image(fast$lights[[1]]$filename)[,, 1:3])),
      0
    )
  }
  expect_true(sky$sun)
  expect_equal(sky$sky_args$resolution, 512)
})

test_that("initial Hosek skies, explicit suns and exported edits use consistent light splitting", {
  sky = sky_light_image(
    0,
    0,
    as.POSIXct("2026-06-21 12:00:00", tz = "UTC"),
    resolution = 16,
    moon = FALSE
  )
  initial = prepare_scene_infinite_lights(list(sky))
  expect_equal(vapply(initial, `[[`, "", "type"), c("image", "disk"))
  controls = native_sky_controls(list(sky))
  replay = native_sky_restore(
    list(sky),
    list(
      model = 0L,
      latitude = 0,
      longitude = 0,
      datetime = controls$datetime,
      manual = TRUE,
      elevation = 90,
      azimuth = 135,
      haze = TRUE,
      altitude = TRUE
    )
  )
  expect_equal(vapply(replay$lights, `[[`, "", "type"), c("image", "disk"))
  expect_equal(replay$controls$elevation, 89.9)
  explicit = sun_light(
    0,
    0,
    sky$datetime,
    sky_args = list(hosek = TRUE),
    resolution = 16
  )
  replaced = prepare_scene_infinite_lights(list(sky, explicit))
  expect_length(replaced, 2)
  expect_equal(replaced[[2]]$name, "sun")
  sky$sun = FALSE
  expect_length(prepare_scene_infinite_lights(list(sky)), 1)
  sky$sun = TRUE
  sky$sky_args$render_mode = "atmosphere"
  expect_length(prepare_scene_infinite_lights(list(sky)), 1)
})

test_that("Hosek disks retain white adaptation and invalid drafts are rejected", {
  sky = sky_light_image(
    0,
    0,
    as.POSIXct("2026-06-21 12:00:00", tz = "UTC"),
    resolution = 16,
    moon = FALSE,
    environment_light_bake_white = TRUE,
    environment_light_bake_white_target = "D50"
  )
  lights = prepare_scene_infinite_lights(list(sky))
  expect_length(lights, 2)
  for (light in lights) {
    pixels = rayimage::ray_read_image(light$filename)
    expect_true(all(is.finite(pixels)))
    expect_gt(sum(pixels[,, 1:3]), 0)
  }
  for (args in list(
    list(hosek = TRUE, turbidity = 0),
    list(hosek = TRUE, elevation = 45),
    list(hosek = TRUE, azimuth = 0)
  )) {
    expect_error(sun_light(0, 0, sky$datetime, sky_args = args), "sky_args")
  }
})
