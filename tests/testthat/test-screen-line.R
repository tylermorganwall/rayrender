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

test_that("screen_line inputs can be supplied as a list", {
  lines = normalize_screen_line(list(
    screen_line(x = 1, y = 2, z = 3, xend = 4, yend = 5, zend = 6),
    screen_line(x = 7, y = 8, z = 9, xend = 10, yend = 11, zend = 12)
  ))

  expect_equal(lines$x, c(1, 7))
  expect_equal(lines$y, c(2, 8))
  expect_equal(lines$z, c(3, 9))
  expect_equal(lines$xend, c(4, 10))
  expect_equal(lines$yend, c(5, 11))
  expect_equal(lines$zend, c(6, 12))
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

test_that("post_process_scene treats native screen line overlay as complete stack", {
  width = 40
  height = 30
  rgb_mat = list(
    r = matrix(1, width, height),
    g = matrix(1, width, height),
    b = matrix(1, width, height),
    a = matrix(1, width, height)
  )
  screen_line_overlay = array(0, c(height, width, 4))
  screen_line_overlay[,, 3] = 1
  screen_line_overlay[,, 4] = 1
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 90,
    nx = width,
    ny = height
  )

  output = post_process_scene(
    rgb_mat,
    iso = 1,
    use_iso = FALSE,
    tonemap = "raw",
    debug_channel = 0,
    filename = NA,
    plot_scene = FALSE,
    bloom = FALSE,
    screen_line = screen_line(
      x = -1,
      y = 0,
      z = 0,
      xend = 1,
      yend = 0,
      zend = 0,
      width = 12,
      color = "red"
    ),
    camera_info = camera_info,
    screen_line_overlay = screen_line_overlay
  )

  expect_equal(max(output[,, 1]), 0)
  expect_equal(max(output[,, 2]), 0)
  expect_equal(min(output[,, 3]), 1)
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
