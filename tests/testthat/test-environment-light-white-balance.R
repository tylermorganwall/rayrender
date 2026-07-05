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
