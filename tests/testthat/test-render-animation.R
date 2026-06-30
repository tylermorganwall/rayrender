test_that("animation frame post-processing matches scene orientation", {
  rgb_mat = list(
    r = matrix(1:6, nrow = 2, ncol = 3),
    g = matrix(11:16, nrow = 2, ncol = 3),
    b = matrix(21:26, nrow = 2, ncol = 3)
  )

  frame = post_process_frame(
    rgb_mat = rgb_mat,
    debug_channel = 2,
    filename = "",
    tonemap = "raw",
    bloom = FALSE,
    write_file = FALSE,
    plot_scene = FALSE
  )

  expect_equal(frame[,, 1], fliplr(flipud(t(rgb_mat$r))))
  expect_equal(frame[,, 2], fliplr(flipud(t(rgb_mat$g))))
  expect_equal(frame[,, 3], fliplr(flipud(t(rgb_mat$b))))
})
