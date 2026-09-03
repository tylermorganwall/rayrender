test_that("post_process_scene does not plot cancelled renders", {
  rgb_mat = list(
    r = matrix(0.1, 2, 2),
    g = matrix(0.2, 2, 2),
    b = matrix(0.3, 2, 2)
  )
  attr(rgb_mat, "render_cancelled") = TRUE
  plot_calls = 0
  local_mocked_bindings(
    plot_image = function(...) {
      plot_calls <<- plot_calls + 1
    },
    .package = "rayimage"
  )

  output = post_process_scene(
    rgb_mat,
    iso = 1,
    use_iso = FALSE,
    tonemap = "raw",
    debug_channel = "none",
    filename = NA,
    plot_scene = TRUE,
    bloom = FALSE
  )

  expect_equal(plot_calls, 0)
  expect_equal(dim(output)[1:2], c(2, 2))
})
