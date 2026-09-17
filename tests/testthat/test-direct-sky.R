test_that("direct skies validate angles and keep location/time defaults", {
  for (constructor in list(sky_light, sky_light_image)) {
    light = constructor(elevation = 25, azimuth = 135)
    expect_equal(light$elevation, 25)
    expect_equal(light$azimuth, 135)
    expect_null(light$datetime)
    expect_false(if (light$atmosphere) light$moon else light$sky_args$moon)
    expect_output(print(light), "elevation: 25 degrees")
    expect_identical(unserialize(serialize(light, NULL)), light)
    for (bad in list(NA, Inf, "25", numeric(), c(10, 20), -91, 91)) {
      expect_error(constructor(elevation = bad, azimuth = 135), "elevation")
    }
    for (bad in list(NA, Inf, "135", numeric(), c(10, 20), -1, 360)) {
      expect_error(constructor(elevation = 25, azimuth = bad), "azimuth")
    }
    expect_error(constructor(elevation = 25), "azimuth")
    expect_error(constructor(azimuth = 135), "elevation")
    expect_error(
      constructor(lat = 0, elevation = 25, azimuth = 135),
      "not both"
    )
    for (field in c("moon", "stars", "planets")) {
      expect_error(
        do.call(
          constructor,
          c(
            list(elevation = 25, azimuth = 135),
            setNames(list(TRUE), field)
          )
        ),
        "requires lat, long, and datetime"
      )
    }
  }
})

test_that("direct image skies use the direct generator and cache by angles", {
  skip_if_not_installed("skymodelr")
  seen = list()
  local_mocked_bindings(
    generate_sky = function(filename, ...) {
      seen[[length(seen) + 1L]] <<- list(...)
      writeBin(as.raw(1), filename)
    },
    generate_sky_latlong = function(...) stop("unexpected ephemeris"),
    .package = "skymodelr"
  )
  light = sky_light_image(
    elevation = 23.45678,
    azimuth = 123.45678,
    resolution = 32,
    turbidity = 4,
    sun = FALSE
  )
  image = prepare_infinite_light(light)
  withr::defer(unlink(image$filename))
  expect_length(seen, 1)
  expect_equal(seen[[1]]$elevation, light$elevation)
  expect_equal(seen[[1]]$azimuth, light$azimuth)
  expect_equal(seen[[1]]$turbidity, 4)
  expect_equal(seen[[1]]$render_mode, "atmosphere")
  expect_null(seen[[1]]$moon)
  expect_null(seen[[1]]$datetime)
  light$intensity = 2
  expect_equal(prepare_infinite_light(light)$filename, image$filename)
  expect_length(seen, 1)
  light$azimuth = 234.56789
  other = prepare_infinite_light(light)
  withr::defer(unlink(other$filename))
  expect_length(seen, 2)
  expect_false(identical(image$filename, other$filename))
})

test_that("direct native skies use supplied geometry and the native solar sampler", {
  skip_if_not_installed("skymodelr")
  local_mocked_bindings(
    get_prague_sky_metadata = function(...) {
      list(
        filename = "model.dat",
        elevation_deg = 80,
        azimuth_deg = 180,
        angular_diameter_deg = 0.54,
        rgb_gain = c(1, 0.9, 0.8)
      )
    },
    .package = "skymodelr"
  )
  light = sky_light(elevation = 25, azimuth = 135, visibility = 40)
  expect_length(sky_light_celestial_lights(light), 0)
  prepared = prepare_scene_infinite_lights(list(light))
  expect_length(prepared, 1)
  expect_equal(prepared[[1]]$elevation, 25)
  expect_equal(prepared[[1]]$azimuth, 135)
  expect_equal(prepared[[1]]$angular_diameter, 0.533)
  expect_equal(prepared[[1]]$visibility, 40)
  expect_equal(prepared[[1]]$rgb_gain, c(1, 0.9, 0.8))
  expect_true(prepared[[1]]$include_sun)
  light$sun = FALSE
  expect_false(prepare_scene_infinite_lights(list(light))[[1]]$include_sun)
})

test_that("only pkgdown preparation downloads missing sky datasets", {
  skip_if_not_installed("skymodelr")
  downloaded = list()
  installed = character()
  local_mocked_bindings(
    list_sky_data = function() data.frame(file = installed),
    download_sky_data = function(sea_level, wide_spectrum) {
      downloaded[[length(downloaded) + 1L]] <<- c(sea_level, wide_spectrum)
    },
    .package = "skymodelr"
  )
  withr::local_envvar(IN_PKGDOWN = "false")
  prepare_pkgdown_sky_data(native = TRUE)
  expect_length(downloaded, 0)
  withr::local_envvar(IN_PKGDOWN = "true")
  prepare_pkgdown_sky_data(native = TRUE)
  prepare_pkgdown_sky_data()
  prepare_pkgdown_sky_data(wide_spectrum = TRUE)
  prepare_pkgdown_sky_data(altitude = 1000)
  expect_equal(
    downloaded,
    list(c(FALSE, FALSE), c(TRUE, FALSE), c(TRUE, TRUE), c(FALSE, FALSE))
  )
  installed = c(
    "SkyModelDataset.dat",
    "SkyModelDatasetGround.dat",
    "PragueSkyModelDatasetGroundInfra.dat"
  )
  prepare_pkgdown_sky_data(native = TRUE)
  prepare_pkgdown_sky_data()
  prepare_pkgdown_sky_data(wide_spectrum = TRUE)
  expect_length(downloaded, 4)
})
