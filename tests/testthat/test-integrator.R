test_that("legacy integrators write a PNG and respect the color debug setting", {
  scene = generate_cornell()
  debug_build = tolower(Sys.getenv("RAY_COLOR_DEBUG", "false")) == "true"
  expected_output = if (debug_build) "Channel: 0" else NA
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
        expected_output
      )
      expect_equal(dim(png::readPNG(filename))[1:2], c(6, 6))
    }
  }
})
