test_that("the selected tone map reaches the native renderer", {
  for (method in c("raw", "hbd", "reinhard", "uncharted")) {
    captured = NULL
    local_mocked_bindings(
      render_scene_rcpp = function(
        scene,
        camera_info,
        scene_info,
        render_info
      ) {
        captured <<- camera_info$tonemap
        list(
          r = matrix(0.18, 3, 3),
          g = matrix(0.18, 3, 3),
          b = matrix(0.18, 3, 3)
        )
      }
    )
    image = render_scene(
      sphere(),
      width = 3,
      height = 3,
      samples = 1,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = method
    )
    expect_identical(captured, method)
    expect_true(all(is.finite(image)))
  }
})
