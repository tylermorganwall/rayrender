test_that("environment light white balance bakes EXR white_current", {
  skip_if_not_installed("rayimage")
  skip_if_not_installed("libopenexr")

  environment_array = array(1, dim = c(2, 4, 3))
  environment_image = rayimage::ray_read_image(
    environment_array,
    source_linear = TRUE,
    assume_colorspace = rayimage::CS_SRGB,
    assume_white = "D60"
  )
  environment_path = tempfile(fileext = ".exr")
  rayimage::ray_write_image(environment_image, environment_path)

  environment_info = prepare_environment_light_white_balance(
    environment_path,
    environment_light_bake_white = TRUE,
    environment_light_bake_white_target = "D65"
  )
  on.exit(unlink(environment_info$cleanup), add = TRUE)

  expect_true(file.exists(environment_info$environment_light))

  baked_environment = rayimage::ray_read_image(
    environment_info$environment_light,
    normalize = FALSE
  )
  expect_equal(
    unname(attr(baked_environment, "white_current")),
    environment_light_white_xyz("D65"),
    tolerance = 1e-6
  )
  expect_false(isTRUE(all.equal(
    as.numeric(baked_environment[1, 1, 1:3]),
    as.numeric(environment_image[1, 1, 1:3])
  )))
})

test_that("environment light white balance no-ops when disabled", {
  environment_path = tempfile(fileext = ".exr")
  environment_info = prepare_environment_light_white_balance(
    environment_path,
    environment_light_bake_white = FALSE
  )

  expect_identical(environment_info$environment_light, environment_path)
  expect_length(environment_info$cleanup, 0L)
})

test_that("image skies cache adapted pixels without regenerating the sky", {
  skip_if_not_installed("skymodelr")
  skip_if_not_installed("rayimage")
  skip_if_not_installed("libopenexr")

  source_image = rayimage::ray_read_image(
    array(1, dim = c(2, 4, 3)),
    source_linear = TRUE,
    assume_colorspace = rayimage::CS_SRGB,
    assume_white = "D60"
  )
  calls = 0L
  forwarded = NULL
  local_mocked_bindings(
    generate_sky_latlong = function(filename, ...) {
      calls <<- calls + 1L
      forwarded <<- list(...)
      rayimage::ray_write_image(source_image, filename)
    },
    .package = "skymodelr"
  )
  time = as.POSIXct("2026-06-21 18:00:00", tz = "UTC") +
    as.numeric(Sys.time()) / 1e6
  files = character()
  withr::defer(unlink(files))
  raw_light = sky_light_image(0, 0, time, resolution = 2, moon = FALSE)
  raw = prepare_infinite_light(raw_light)
  files = c(files, raw$filename)

  light = sky_light_image(
    0,
    0,
    time,
    resolution = 2,
    moon = FALSE,
    environment_light_bake_white = TRUE,
    environment_light_bake_white_target = "D65"
  )
  baked = prepare_infinite_light(light)
  files = c(files, baked$filename)
  expect_equal(calls, 1L)
  expect_false(any(grepl("bake_white", names(forwarded))))
  expect_false(identical(baked$filename, raw$filename))

  image = rayimage::ray_read_image(baked$filename, normalize = FALSE)
  expect_equal(
    unname(attr(image, "white_current")),
    environment_light_white_xyz("D65"),
    tolerance = 1e-6
  )
  expect_false(isTRUE(all.equal(
    as.numeric(image[1, 1, 1:3]),
    as.numeric(source_image[1, 1, 1:3])
  )))

  light$intensity = 2
  light$rotation = 45
  light$name = "adapted sky"
  reused = prepare_infinite_light(light)
  expect_identical(reused$filename, baked$filename)
  expect_equal(reused$intensity, 2)
  expect_equal(reused$rotation, 45)
  expect_identical(reused$name, "adapted sky")

  light$environment_light_bake_white_target = "D50"
  warmer = prepare_infinite_light(light)
  files = c(files, warmer$filename)
  expect_false(identical(warmer$filename, baked$filename))
  expect_equal(
    unname(attr(
      rayimage::ray_read_image(warmer$filename, normalize = FALSE),
      "white_current"
    )),
    environment_light_white_xyz("D50"),
    tolerance = 1e-6
  )
  light$environment_light_bake_white_target = 2 *
    environment_light_white_xyz("D50")
  expect_identical(prepare_infinite_light(light)$filename, warmer$filename)
  expect_equal(calls, 1L)

  # Baking leaves the original image available for unadapted lights.
  expect_identical(prepare_infinite_light(raw_light)$filename, raw$filename)
  expect_equal(
    unname(attr(
      rayimage::ray_read_image(raw$filename, normalize = FALSE),
      "white_current"
    )),
    environment_light_white_xyz("D60"),
    tolerance = 1e-6
  )
})

test_that("white balance belongs to the image sky constructor", {
  time = as.POSIXct("2026-06-21 18:00:00", tz = "UTC")
  for (value in list(NA, NULL, 1, c(TRUE, FALSE))) {
    expect_error(
      sky_light_image(0, 0, time, environment_light_bake_white = value),
      "environment_light_bake_white must be TRUE or FALSE"
    )
  }
  expect_error(
    sky_light_image(
      0,
      0,
      time,
      environment_light_bake_white = TRUE,
      environment_light_bake_white_target = "invalid"
    ),
    "Unknown.*environment_light_bake_white_target"
  )
  expect_error(
    sky_light_image(
      0,
      0,
      time,
      environment_light_bake_white = TRUE,
      exr_metadata = FALSE
    ),
    "requires exr_metadata = TRUE"
  )
  expect_error(
    render_scene(sphere(), environment_light_bake_white = TRUE),
    "unused argument"
  )
  expect_error(
    render_animation(sphere(), environment_light_bake_white_target = "D65"),
    "unused argument"
  )
})
