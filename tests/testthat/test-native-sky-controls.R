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
