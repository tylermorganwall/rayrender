test_that("Metal moving closed boundaries preserve analytic priority attenuation", {
  withr::local_options(list(cores = 6L))
  # A unit thickness of each active medium gives optical depth .4 + .3.
  # Linear shutter scaling from 1 to .5 has a known integrated transmittance.
  body = cube(
    xwidth = 10,
    ywidth = 10,
    zwidth = 2,
    material = subsurface(
      sigma_a = .4,
      sigma_s = 0,
      refraction = 1,
      priority = 2
    )
  ) |>
    add_object(cube(
      xwidth = 11,
      ywidth = 11,
      zwidth = 1,
      material = dielectric(
        refraction = 1,
        attenuation = rep(.3, 3),
        priority = 1
      )
    )) |>
    animate_objects(end_position = c(.3, .2, 0), end_scale = c(1, 1, .5))
  expected = (exp(-.35) - exp(-.7)) / .35
  for (refit in c(FALSE, TRUE)) {
    scene = if (refit) add_object(body, sphere(x = 50)) else body
    for (backend in c("metal", "nee")) {
      set.seed(846)
      image = render_scene(
        scene,
        integrator_type = backend,
        width = 12,
        height = 12,
        samples = 256,
        sample_method = "sobol",
        shutter_speed = 1,
        lookfrom = c(0, 0, 4),
        lookat = c(0, 0, 0),
        fov = 0,
        ortho_dimensions = c(.5, .5),
        max_depth = 12,
        ambient_light = TRUE,
        backgroundlow = "white",
        backgroundhigh = "white",
        preview = FALSE,
        plot_scene = FALSE,
        progress = FALSE,
        denoise = FALSE,
        min_variance = 0,
        bloom = FALSE,
        tonemap = "raw"
      )
      if (backend == "metal") {
        report = attr(image, "wavefront")
        if (!isTRUE(report$used)) {
          expect_match(report$fallback, "not enabled|device|macOS 11")
          skip("Metal unavailable")
        }
        expect_equal(report$discarded_paths, 0)
      }
      expect_equal(mean(image[,, 1:3]), expected, tolerance = .02)
    }
  }
})
