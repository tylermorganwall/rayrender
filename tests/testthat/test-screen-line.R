test_that("screen_line recycles endpoints and style", {
  lines = screen_line(
    start = matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE),
    end = matrix(c(1, 0, 0, 2, 1, 1), ncol = 3, byrow = TRUE),
    offset = c(1, 2),
    width = 4,
    color = "red",
    alpha = 0.5,
    occlusion = TRUE,
    occlusion_mode = "partial"
  )

  expect_s3_class(lines, "ray_screen_line")
  expect_equal(lines$x, c(0, 1))
  expect_equal(lines$xend, c(1, 2))
  expect_equal(lines$x_offset, c(1, 1))
  expect_equal(lines$yend_offset, c(2, 2))
  expect_equal(lines$width, c(4, 4))
  expect_equal(lines$occlusion_mode, c("line", "line"))
})

test_that("add_screen_line composites antialiased line into an image", {
  image_array = array(1, c(40, 40, 3))
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 0,
    ortho_dimensions = c(4, 4),
    nx = 40,
    ny = 40
  )

  output = add_screen_line(
    image_array,
    screen_line(
      x = -1,
      y = 0,
      z = 0,
      xend = 1,
      yend = 0,
      zend = 0,
      width = 3,
      color = "black"
    ),
    camera_info
  )

  expect_equal(dim(output), dim(image_array))
  expect_lt(min(output[,, 1]), 1)
  expect_true(any(output[,, 1] > 0 & output[,, 1] < 1))
})

test_that("prepare_screen_line_preview marks partial line occlusion overlays", {
  overlays = prepare_screen_line_preview(
    screen_line(occlusion = TRUE, occlusion_mode = "line")
  )

  expect_true(overlays$active)
  expect_true(overlays$lines[[1]]$partial_occlusion)
  expect_false(overlays$lines[[1]]$occlusion)
  expect_true(screen_line_needs_native_overlay(
    screen_line(occlusion = TRUE, occlusion_mode = "partial")
  ))
})

test_that("prepare_screen_line_occlusion uses midpoint for anchor occlusion", {
  occlusion = prepare_screen_line_occlusion(
    screen_line(
      x = 0,
      y = 0,
      z = 0,
      xend = 2,
      yend = 2,
      zend = 2,
      occlusion = TRUE
    )
  )

  expect_true(occlusion$active)
  expect_equal(occlusion$x, 1)
  expect_equal(occlusion$y, 1)
  expect_equal(occlusion$z, 1)
  expect_true(occlusion$occlusion)
})
