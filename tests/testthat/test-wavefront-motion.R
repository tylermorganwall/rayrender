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

test_that("Metal geometry shutter times vary across pixels at one sample", {
  withr::local_options(list(cores = 6L))
  # The slab travels two units. In the interior of its swept silhouette, its
  # .4-unit width gives exactly .2 coverage. Shared whole-image times instead
  # leave fully white/black columns, even after averaging vertically.
  scene = cube(
    x = -1,
    xwidth = .4,
    ywidth = 5,
    zwidth = .1,
    material = light(intensity = 1)
  ) |>
    animate_objects(end_position = c(2, 0, 0))
  render = function(speed) {
    set.seed(61)
    render_scene(
      scene,
      integrator_type = "metal",
      width = 64,
      height = 96,
      samples = 1,
      lookfrom = c(0, 0, 6),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(4, 4),
      shutter_speed = speed,
      ambient_light = FALSE,
      backgroundlow = "black",
      backgroundhigh = "black",
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      bloom = FALSE,
      tonemap = "raw"
    )
  }
  image = render(1)
  report = attr(image, "wavefront")
  if (!isTRUE(report$used)) {
    expect_match(report$fallback, "not enabled|device|macOS 11")
    skip("Metal unavailable")
  }
  expect_equal(report$discarded_paths, 0)
  coverage = colMeans(image[, 23:42, 1])
  expect_lt(max(abs(coverage - .2)), .15)
  expect_lt(abs(mean(coverage) - .2), .04)
  frozen = render(Inf)
  expect_lt(max(frozen[, 25:42, 1]), .001)
  expect_gt(mean(frozen[, 15:17, 1]), .9)
})

test_that("Metal frozen animated shading matches the same static placement", {
  withr::local_options(list(cores = 6L))
  alpha = matrix(1, 16, 16)
  alpha[5:12, 5:12] = 0
  bump = outer(
    seq(0, 2 * pi, length.out = 16),
    seq(0, 2 * pi, length.out = 16),
    function(x, y) .5 + .25 * sin(x) * cos(y)
  )
  paint = diffuse(
    "gold",
    checkercolor = "steelblue",
    alpha_texture = alpha,
    bump_texture = bump,
    bump_intensity = .01
  )
  objects = sphere(x = -.8, y = .8, radius = .7, material = paint) |>
    add_object(cube(x = .8, y = .8, material = paint)) |>
    add_object(sphere(y = 4, radius = .5, material = light(intensity = 6)))
  moving = animate_objects(
    objects,
    end_position = c(.3, .2, 0),
    end_angle = c(0, 37, 0),
    end_scale = c(1.2, .8, 1.1)
  )
  static = group_objects(
    objects,
    translate = c(.3, .2, 0),
    angle = c(0, 37, 0),
    scale = c(1.2, .8, 1.1)
  )
  render = function(objects) {
    set.seed(416)
    render_scene(
      add_object(generate_ground(), objects),
      integrator_type = "metal",
      width = 40,
      height = 40,
      samples = 16,
      max_depth = 5,
      lookfrom = c(0, 2, 7),
      lookat = c(0, 1, 0),
      fov = 35,
      shutteropen = 1,
      shutterclose = 1,
      shutter_speed = Inf,
      ambient_light = FALSE,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      bloom = FALSE,
      tonemap = "raw"
    )
  }
  actual = render(moving)
  report = attr(actual, "wavefront")
  if (!isTRUE(report$used)) {
    expect_match(report$fallback, "not enabled|device|macOS 11")
    skip("Metal unavailable")
  }
  expected = render(static)
  expect_equal(report$discarded_paths, 0)
  expect_lt(sqrt(mean((actual[,, 1:3] - expected[,, 1:3])^2)), .001)
})
