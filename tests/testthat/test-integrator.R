test_that("legacy integrators write a PNG without printing to stdout", {
  scene = generate_cornell()
  integrators = c("rtiow", "nee", "basic")
  for (i in seq_along(integrators)) {
    for (denoise in c(FALSE, TRUE)) {
      filename = withr::local_tempfile(fileext = ".png")
      set.seed(i)
      expect_output(
        render_scene(
          scene,
          samples = 4,
          denoise = denoise,
          preview = FALSE,
          plot_scene = FALSE,
          width = 6,
          height = 6,
          filename = filename,
          parallel = FALSE,
          integrator_type = integrators[i]
        ),
        NA
      )
      expect_equal(dim(png::readPNG(filename))[1:2], c(6, 6))
    }
  }
})
