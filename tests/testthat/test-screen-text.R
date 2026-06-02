test_that("screen_text recycles labels and anchor data", {
  labels = screen_text(
    label = c("a", "b"),
    point = matrix(c(0, 0, 0, 1, 1, 1), ncol = 3, byrow = TRUE),
    offset = c(4, -2),
    size = 20,
    halo_color = "white",
    halo_offset = c(1, 3),
    halo_expand = 2,
    halo_gap_fill = 3,
    halo_gap_fill_alpha_threshold = 0.4,
    occlusion = TRUE,
    occlusion_mode = "label",
    occlusion_tolerance = 0.01
  )

  expect_s3_class(labels, "ray_screen_text")
  expect_equal(labels$x, c(0, 1))
  expect_equal(labels$y_offset, c(-2, -2))
  expect_equal(labels$size, c(20, 20))
  expect_equal(labels$halo_color, c("white", "white"))
  expect_equal(labels$halo_y_offset, c(3, 3))
  expect_equal(labels$halo_expand, c(2, 2))
  expect_equal(labels$halo_gap_fill, c(3, 3))
  expect_equal(labels$halo_gap_fill_alpha_threshold, c(0.4, 0.4))
  expect_equal(labels$occlusion, c(TRUE, TRUE))
  expect_equal(labels$occlusion_mode, c("label", "label"))
  expect_equal(labels$occlusion_tolerance, c(0.01, 0.01))
})

test_that("screen_text inputs can be supplied as a list", {
  labels = normalize_screen_text(list(
    screen_text("a", x = 1, y = 2, z = 3),
    screen_text("b", x = 4, y = 5, z = 6)
  ))

  expect_equal(labels$label, c("a", "b"))
  expect_equal(labels$x, c(1, 4))
  expect_equal(labels$y, c(2, 5))
  expect_equal(labels$z, c(3, 6))
})

test_that("project_points_to_screen matches perspective camera center", {
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 90,
    nx = 101,
    ny = 101
  )

  projected = project_points_to_screen(
    matrix(c(0, 0, 0, 1, 1, 0), ncol = 3, byrow = TRUE),
    camera_info
  )

  expect_equal(projected$x[1], 51)
  expect_equal(projected$y[1], 51)
  expect_true(projected$in_frame[1])
  expect_lt(projected$x[2], projected$x[1])
  expect_lt(projected$y[2], projected$y[1])
})

test_that("project_points_to_screen matches orthographic camera center", {
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 0,
    ortho_dimensions = c(4, 2),
    nx = 101,
    ny = 51
  )

  projected = project_points_to_screen(
    matrix(c(0, 0, 0, 2, 1, 0), ncol = 3, byrow = TRUE),
    camera_info
  )

  expect_equal(projected$x[1], 51)
  expect_equal(projected$y[1], 26)
  expect_equal(projected$x[2], 1)
  expect_equal(projected$y[2], 1)
})

test_that("project_points_to_screen uses final renderer camera orientation", {
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 0,
    ortho_dimensions = c(4, 4),
    nx = 101,
    ny = 101,
    screen_camera_origin = c(10, 0, 0),
    screen_camera_u = c(0, 0, 1),
    screen_camera_v = c(0, 1, 0),
    screen_camera_w = c(-1, 0, 0),
    screen_camera_fov = 0,
    screen_camera_ortho_dimensions = c(4, 4)
  )

  projected = project_points_to_screen(
    matrix(c(0, 0, 1), ncol = 3, byrow = TRUE),
    camera_info
  )

  expect_equal(projected$x, 26)
  expect_equal(projected$y, 51)
  expect_true(projected$in_frame)
})

test_that("orthographic screen projection matches post-processed image orientation", {
  camera_info = list(
    lookfrom = c(1, 10, 0),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 0,
    ortho_dimensions = c(10, 10),
    nx = 400,
    ny = 400
  )

  projected = project_points_to_screen(
    matrix(c(0, 0, 2), ncol = 3, byrow = TRUE),
    camera_info
  )

  expect_lt(projected$x, 200.5)
  expect_equal(projected$y, 200.5)
})

test_that("add_screen_text composites text into an image", {
  image_array = array(1, c(80, 100, 3))
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 90,
    nx = 100,
    ny = 80
  )

  output = add_screen_text(
    image_array,
    screen_text("A", size = 20, color = "black", hjust = 0.5, vjust = 0.5),
    camera_info
  )

  expect_equal(dim(output), dim(image_array))
  expect_true(any(output[,, 1] < 1))
})

test_that("add_screen_text composites text halo into an image", {
  image_array = array(0, c(80, 100, 3))
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 90,
    nx = 100,
    ny = 80
  )

  output = add_screen_text(
    image_array,
    screen_text(
      "A",
      size = 20,
      color = "black",
      hjust = 0.5,
      vjust = 0.5,
      halo_color = "white",
      halo_expand = 4,
      halo_blur = 1
    ),
    camera_info
  )

  expect_equal(dim(output), dim(image_array))
  expect_true(any(output[,, 1] > 0.5))
})

test_that("add_screen_text skips labels hidden by visibility vector", {
  image_array = array(1, c(80, 100, 3))
  camera_info = list(
    lookfrom = c(0, 0, -10),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 90,
    nx = 100,
    ny = 80
  )

  output = add_screen_text(
    image_array,
    screen_text("A", size = 20, color = "black", hjust = 0.5, vjust = 0.5),
    camera_info,
    screen_text_visible = FALSE
  )

  expect_equal(output, image_array)
})

test_that("post_process_scene treats native screen text overlay as complete stack", {
  width = 40
  height = 30
  rgb_mat = list(
    r = matrix(1, width, height),
    g = matrix(1, width, height),
    b = matrix(1, width, height),
    a = matrix(1, width, height)
  )
  screen_text_overlay = array(0, c(height, width, 4))
  screen_text_overlay[,, 3] = 1
  screen_text_overlay[,, 4] = 1
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
    screen_text = screen_text(
      "covered",
      size = 18,
      color = "red",
      background_color = "red",
      background_alpha = 1,
      hjust = 0.5,
      vjust = 0.5
    ),
    camera_info = camera_info,
    screen_text_overlay = screen_text_overlay
  )

  expect_equal(max(output[,, 1]), 0)
  expect_equal(max(output[,, 2]), 0)
  expect_equal(min(output[,, 3]), 1)
})

test_that("prepare_screen_text_preview caches foreground and halo overlays", {
  overlays = prepare_screen_text_preview(
    screen_text(
      "A",
      halo_color = "white",
      halo_expand = 2,
      halo_blur = 1
    )
  )

  expect_true(overlays$active)
  expect_length(overlays$overlays, 2)
  expect_equal(dim(overlays$overlays[[1]]$image)[3], 4)
  expect_equal(dim(overlays$overlays[[2]]$image)[3], 4)
  expect_lt(overlays$overlays[[1]]$x_offset, overlays$overlays[[2]]$x_offset)
  expect_false(overlays$overlays[[1]]$occlusion)
  expect_false(overlays$overlays[[1]]$partial_occlusion)
})

test_that("prepare_screen_text_preview marks partial label occlusion overlays", {
  overlays = prepare_screen_text_preview(
    screen_text("A", occlusion = TRUE, occlusion_mode = "partial")
  )

  expect_true(overlays$active)
  expect_true(overlays$overlays[[1]]$partial_occlusion)
  expect_false(overlays$overlays[[1]]$occlusion)
  expect_true(screen_text_needs_native_overlay(
    screen_text("A", occlusion = TRUE, occlusion_mode = "label")
  ))
})

test_that("prepare_screen_text_occlusion caches per-label occlusion inputs", {
  occlusion = prepare_screen_text_occlusion(
    screen_text(
      c("A", "B"),
      x = c(0, 1),
      occlusion = c(TRUE, FALSE),
      occlusion_tolerance = c(0.01, 0.02)
    )
  )

  expect_true(occlusion$active)
  expect_equal(occlusion$x, c(0, 1))
  expect_equal(occlusion$occlusion, c(TRUE, FALSE))
  expect_equal(occlusion$occlusion_tolerance, c(0.01, 0.02))
})

test_that("prepare_screen_text_occlusion ignores partial label occlusion", {
  occlusion = prepare_screen_text_occlusion(
    screen_text(
      c("A", "B"),
      occlusion = TRUE,
      occlusion_mode = c("label", "anchor")
    )
  )

  expect_true(occlusion$active)
  expect_equal(occlusion$occlusion, c(FALSE, TRUE))
})
