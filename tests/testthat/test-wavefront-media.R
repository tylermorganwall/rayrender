# Compare integrated observables, not individual noise samples from different
# backends. Analytic slabs test the transport independently of the CPU integrator.
wavefront_medium_render = function(
  scene,
  backend = 'metal',
  inside = FALSE,
  samples = 64,
  max_depth = 5
) {
  set.seed(513)
  render_scene(
    scene,
    integrator_type = backend,
    width = 24,
    height = 24,
    samples = samples,
    sample_method = 'sobol',
    parallel = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    min_variance = 0,
    lookfrom = if (inside) c(0, 0, 0) else c(0, 0, 4),
    lookat = c(0, 0, -1),
    fov = 0,
    ortho_dimensions = c(1, 1),
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white',
    max_depth = max_depth
  )
}

test_that('Metal homogeneous media reproduce absorption and absorption-weighted emission', {
  sigma = c(.2, .6, 1)
  slab = set_medium(
    cube(width = 2),
    homogeneous_medium(sigma_a = sigma, sigma_s = 0)
  )
  image = suppressWarnings(wavefront_medium_render(slab))
  if (!isTRUE(attr(image, 'wavefront')$used)) {
    skip('Metal ray tracing unavailable')
  }
  expect_equal(apply(image[,, 1:3], 3, mean), exp(-2 * sigma), tolerance = .015)
  expect_equal(attr(image, 'wavefront')$discarded_paths, 0)
  inside = wavefront_medium_render(slab, inside = TRUE)
  expect_equal(apply(inside[,, 1:3], 3, mean), exp(-sigma), tolerance = .015)
  emission = c(.5, 1, 2)
  glowing = set_medium(
    cube(width = 2),
    homogeneous_medium(sigma_a = sigma, sigma_s = 0, emission = emission)
  )
  image = wavefront_medium_render(glowing)
  expect_equal(
    apply(image[,, 1:3], 3, mean),
    emission + (1 - emission) * exp(-2 * sigma),
    tolerance = .02
  )
  # Reflection spends the one allowed bounce before reaching a volume behind
  # the camera. Absorption still applies, but terminal medium emission does not.
  mirror = xy_rect(xwidth = 20, ywidth = 20, material = metal('white'))
  terminal = wavefront_medium_render(
    mirror |>
      add_object(set_medium(
        cube(z = 6, width = 2),
        homogeneous_medium(sigma_a = sigma, sigma_s = 0, emission = emission)
      )),
    max_depth = 1
  )
  unlit = wavefront_medium_render(
    mirror |>
      add_object(set_medium(
        cube(z = 6, width = 2),
        homogeneous_medium(sigma_a = sigma, sigma_s = 0)
      )),
    max_depth = 1
  )
  expect_gt(mean(unlit[,, 1:3]), .1)
  expect_equal(as.numeric(terminal), as.numeric(unlit), tolerance = 1e-6)
})

test_that('Metal preserves nested media, transforms and native grid interpolation', {
  probe = suppressWarnings(wavefront_medium_render(set_medium(
    cube(),
    homogeneous_medium()
  )))
  if (!isTRUE(attr(probe, 'wavefront')$used)) {
    skip('Metal ray tracing unavailable')
  }
  scene = set_medium(
    cube(width = 4),
    homogeneous_medium(sigma_a = .1, sigma_s = 0)
  ) |>
    add_object(set_medium(
      cube(width = 2),
      homogeneous_medium(sigma_a = .7, sigma_s = 0)
    ))
  image = wavefront_medium_render(scene)
  expect_equal(mean(image[,, 1:3]), exp(-1.6), tolerance = .015)
  transform = diag(c(2, 1.5, .8, 1))
  grid = grid_medium(
    array(seq(.1, 1, length.out = 5 * 7 * 9), c(5, 7, 9)),
    sigma_a = c(.1, .2, .4),
    sigma_s = .7,
    g = .5,
    medium_transform = transform,
    temperature = array(seq(1200, 8500, length.out = 5 * 7 * 9), c(5, 7, 9))
  )
  scene = set_medium(cube(width = 2), grid)
  gpu = wavefront_medium_render(scene, samples = 128)
  cpu = wavefront_medium_render(scene, 'nee', samples = 128)
  expect_true(all(is.finite(gpu)))
  expect_equal(
    apply(gpu[,, 1:3], 3, mean),
    apply(cpu[,, 1:3], 3, mean),
    tolerance = .015
  )
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
})

test_that('Metal sparse OpenVDB and NanoVDB retain transforms, tiles and temperature', {
  probe = suppressWarnings(wavefront_medium_render(set_medium(
    cube(),
    homogeneous_medium()
  )))
  if (!isTRUE(attr(probe, 'wavefront')$used)) {
    skip('Metal ray tracing unavailable')
  }
  for (file in c(
    'transformed.vdb',
    'transformed-zip.nvdb',
    'transformed-blosc.nvdb',
    'tiles.nvdb'
  )) {
    transform = diag(4)
    if (file == 'tiles.nvdb') {
      transform[1:3, 4] = -3.5
    }
    medium = vdb_medium(
      test_path('fixtures', 'volumes', file),
      sigma_a = .3,
      sigma_s = .4,
      temperature_grid = if (file != 'tiles.nvdb') 'temperature' else NULL,
      medium_transform = transform,
      emission_scale = .1
    )
    scene = set_medium(cube(width = 2), medium)
    gpu = wavefront_medium_render(scene)
    cpu = wavefront_medium_render(scene, 'nee')
    expect_true(isTRUE(attr(gpu, 'wavefront')$used))
    expect_true(all(is.finite(gpu)))
    expect_equal(
      apply(gpu[,, 1:3], 3, mean),
      apply(cpu[,, 1:3], 3, mean),
      tolerance = .02
    )
  }
})

test_that('Metal keeps unsupported extreme thermal fields explicit', {
  scene = set_medium(
    cube(),
    grid_medium(array(1, c(2, 2, 2)), temperature = array(25000, c(2, 2, 2)))
  )
  expect_warning(
    image <- wavefront_medium_render(scene),
    '20000 kelvin.*CPU NEE'
  )
  expect_false(isTRUE(attr(image, 'wavefront')$used))
})

test_that('Metal cloud shadow flights tolerate representable rounding', {
  skip_if_not_installed('ambient')
  scene = generate_ground(material = diffuse('grey30')) |>
    add_object(cloud(
      y = 1.5,
      width = 3,
      height = 2,
      depth = 2,
      resolution = 40,
      optical_depth = 4,
      seed = 27,
      g = .4
    )) |>
    add_infinite_light(disk_light(
      direction = c(-1, 2, 1),
      angular_diameter = 25,
      intensity = 5
    ))
  set.seed(157)
  image = suppressWarnings(render_scene(
    scene,
    integrator_type = 'metal',
    width = 192,
    height = 192,
    samples = 32,
    sample_method = 'sobol',
    lookfrom = c(3, 2.8, 6),
    lookat = c(0, 1, 0),
    fov = 38,
    min_variance = 0,
    max_depth = 6,
    denoise = FALSE,
    bloom = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    ambient_light = FALSE,
    parallel = FALSE
  ))
  if (!isTRUE(attr(image, 'wavefront')$used)) {
    skip('Metal unavailable')
  }
  expect_equal(attr(image, 'wavefront')$discarded_paths, 0)
  expect_true(all(is.finite(image)))
})

test_that('Metal shadow grid flights preserve colored absorption and clipped endpoints', {
  probe = suppressWarnings(wavefront_medium_render(set_medium(
    cube(),
    homogeneous_medium()
  )))
  if (!isTRUE(attr(probe, 'wavefront')$used)) {
    skip('Metal unavailable')
  }
  # The camera misses this off-axis slab. Only direct-light shadow rays cross it,
  # isolating ratio tracking from camera free-flight sampling and multiple scattering.
  render = function(medium, seed) {
    scene = xy_rect(xwidth = 8, ywidth = 8, material = diffuse('white')) |>
      add_object(set_medium(cube(x = 2, z = 2, width = 1), medium)) |>
      add_infinite_light(disk_light(
        direction = c(1, 0, 1),
        angular_diameter = 5,
        intensity = 100
      ))
    set.seed(seed)
    image = render_scene(
      scene,
      integrator_type = 'metal',
      width = 24,
      height = 24,
      samples = 256,
      sample_method = 'sobol',
      max_depth = 1,
      lookfrom = c(0, 0, 3),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(.2, .2),
      aperture = 0,
      tonemap = 'raw',
      denoise = FALSE,
      bloom = FALSE,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      min_variance = 0,
      ambient_light = FALSE
    )
    if (!isTRUE(attr(image, 'wavefront')$used)) {
      skip('Metal unavailable')
    }
    expect_equal(attr(image, 'wavefront')$discarded_paths, 0)
    expect_true(all(is.finite(image)))
    apply(image[,, 1:3], 3, mean)
  }
  sigma = c(.2, .6, 1)
  # The transformed grid extends past the actual boundary: flights must stop at
  # the geometry endpoint rather than integrating the rest of the grid.
  grid = grid_medium(
    array(1, c(17, 19, 23)),
    sigma_a = sigma,
    sigma_s = 0,
    medium_transform = diag(c(2, .8, 2, 1))
  )
  homogeneous = homogeneous_medium(sigma_a = sigma, sigma_s = 0)
  results = lapply(c(157, 513, 2026), function(seed) {
    a = render(homogeneous, seed)
    b = render(grid, seed)
    data.frame(
      seed,
      channel = 1:3,
      homogeneous = a,
      grid = b,
      relative_error = (b - a) / a
    )
  })
  results = do.call(rbind, results)
  expect_lt(max(abs(results$relative_error)), .02)
  empty = render(
    grid_medium(array(0, c(17, 19, 23)), sigma_a = sigma, sigma_s = 0),
    513
  )
  vacuum = render(homogeneous_medium(sigma_a = 0, sigma_s = 0), 513)
  expect_equal(empty, vacuum, tolerance = 1e-6)
})
