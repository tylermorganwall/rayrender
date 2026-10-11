test_that('invisible volumes preserve the sky observer when finite haze is disabled', {
  if (!any(skymodelr::list_sky_data()$file == 'SkyModelDataset.dat')) {
    skip('Full-altitude Prague data not installed')
  }
  sky = sky_light(
    elevation = 20,
    azimuth = 160,
    sun = FALSE,
    moon = FALSE,
    haze = FALSE,
    query_altitude = TRUE,
    meters_per_unit = 1000,
    sampling_resolution = 16,
    prague_rgb_correction = FALSE
  )
  base = sphere(x = 20000) |> add_infinite_light(sky)
  # This box moves the tracing origin several kilometres above the camera,
  # but crossing vacuum must not change the sky's physical observer position.
  split = base |>
    add_object(set_medium(
      cube(y = 4, z = 6, xwidth = 100, ywidth = 4, zwidth = 8),
      homogeneous_medium(sigma_s = 0)
    ))
  render = function(scene, backend = 'metal') {
    set.seed(2026)
    render_scene(
      scene,
      integrator_type = backend,
      width = 20,
      height = 20,
      samples = 4,
      lookfrom = c(0, .35, -8),
      lookat = c(0, 2.4, 3),
      fov = 47,
      iso = 4,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      min_variance = 0,
      parallel = FALSE
    )
  }
  empty = suppressWarnings(render(base))
  if (!isTRUE(attr(empty, 'wavefront')$used)) {
    skip('Metal unavailable')
  }
  vacuum = render(split)
  reference = render(split, 'nee')
  expect_true(isTRUE(attr(vacuum, 'wavefront')$used))
  expect_equal(
    as.numeric(vacuum[,, 1:3]),
    as.numeric(empty[,, 1:3]),
    tolerance = 1e-5
  )
  expect_equal(mean(vacuum[,, 1:3]), mean(reference[,, 1:3]), tolerance = .005)
  # Direct-light connections through the same box must also retain the surface
  # observer, rather than evaluating incident sky at the box's far boundary.
  ground = xz_rect(xwidth = 160, zwidth = 160, material = diffuse('#469D60'))
  lit = render(add_object(base, ground))
  lit_vacuum = render(add_object(split, ground))
  expect_lt(mean(abs(lit_vacuum[,, 1:3] - lit[,, 1:3])), 1e-5)
})

test_that('Metal evaluates Prague lighting and preserves haze across vacuum boundaries', {
  probe = suppressWarnings(render_scene(
    sphere(),
    integrator_type = 'metal',
    width = 4,
    height = 4,
    samples = 1,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    parallel = FALSE
  ))
  if (!isTRUE(attr(probe, 'wavefront')$used)) {
    skip('Metal unavailable')
  }
  datasets = skymodelr::list_sky_data()
  if (!any(datasets$file == 'SkyModelDataset.dat')) {
    skip('Full-altitude Prague data not installed')
  }
  sky = sky_light(
    elevation = 20,
    azimuth = 160,
    visibility = 50,
    sun = FALSE,
    moon = FALSE,
    haze_in_volumes = TRUE,
    sampling_resolution = 16,
    prague_rgb_correction = FALSE
  )
  base = sphere(x = 20000, radius = 10) |> add_infinite_light(sky)
  vacuum = base |>
    add_object(set_medium(cube(width = 1000), homogeneous_medium(sigma_s = 0)))
  render = function(scene, backend = 'metal', rotation = 0, samples = 16) {
    set.seed(804)
    render_scene(
      scene,
      integrator_type = backend,
      width = 20,
      height = 20,
      samples = samples,
      lookfrom = c(0, 20, 400),
      lookat = c(0, 40, 0),
      fov = 40,
      iso = 3,
      rotate_env = rotation,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      min_variance = 0,
      parallel = FALSE,
      max_depth = 3
    )
  }
  empty = render(base)
  split = render(vacuum)
  expect_true(isTRUE(attr(split, 'wavefront')$used))
  expect_equal(attr(split, 'wavefront')$atmosphere_seconds, 0)
  expect_true(all(is.finite(split)))
  expect_equal(as.numeric(split), as.numeric(empty), tolerance = 1e-5)
  cpu = render(base, 'nee')
  expect_equal(mean(empty[,, 1:3]), mean(cpu[,, 1:3]), tolerance = .005)
  rotated = render(base, rotation = 80)
  expect_gt(mean(abs(rotated[,, 1:3] - empty[,, 1:3])), .01)
  cpu_rotated = render(base, 'nee', 80)
  expect_equal(
    mean(rotated[,, 1:3]),
    mean(cpu_rotated[,, 1:3]),
    tolerance = .005
  )
  wall = xy_rect(
    z = -1000,
    xwidth = 10000,
    ywidth = 10000,
    material = diffuse('grey50')
  ) |>
    add_infinite_light(sky)
  reflected = render(wall, samples = 128)
  reference = render(wall, 'nee', samples = 128)
  expect_equal(
    apply(reflected[,, 1:3], 3, mean),
    apply(reference[,, 1:3], 3, mean),
    tolerance = .004
  )
  fog = wall |>
    add_object(set_medium(
      cube(width = 600),
      homogeneous_medium(sigma_s = .002, sigma_a = .0002, g = .3, haze = TRUE)
    ))
  scattered = render(fog, samples = 128)
  reference = render(fog, 'nee', samples = 128)
  expect_equal(attr(scattered, 'wavefront')$discarded_paths, 0)
  expect_equal(
    apply(scattered[,, 1:3], 3, mean),
    apply(reference[,, 1:3], 3, mean),
    tolerance = .01
  )
  # Density thresholds now split native trilinear cells directly on the GPU.
  threshold = base |>
    add_object(cloud(
      y = 40,
      width = 400,
      height = 400,
      depth = 400,
      resolution = 24,
      optical_depth = 1,
      haze = TRUE,
      haze_density_threshold = .05,
      seed = 27
    ))
  threshold_gpu = render(threshold, samples = 4)
  expect_true(isTRUE(attr(threshold_gpu, 'wavefront')$used))
  expect_equal(attr(threshold_gpu, 'wavefront')$atmosphere_seconds, 0)
  expect_equal(attr(threshold_gpu, 'wavefront')$discarded_paths, 0)
})
