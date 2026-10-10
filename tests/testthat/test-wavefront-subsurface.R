# Physical invariants and CPU comparisons for the public GPU transport path.
# Tiny images keep tests bounded; aggregate expectations avoid per-pixel noise tests.
render_wavefront_body = function(
  scene,
  backend = 'metal',
  samples = 64,
  lookfrom = c(0, 0, 4),
  lookat = c(0, 0, 0),
  width = 12,
  depth = 4
) {
  set.seed(738)
  image = render_scene(
    scene,
    integrator_type = backend,
    width = width,
    height = width,
    samples = samples,
    lookfrom = lookfrom,
    lookat = lookat,
    fov = 0,
    ortho_dimensions = c(.5, .5),
    max_depth = depth,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white',
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    min_variance = 0,
    bloom = FALSE
  )
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    reason = attr(image, 'wavefront')$fallback
    if (grepl('not enabled|device|macOS 11', reason)) {
      skip('Metal is unavailable')
    }
    stop('Unexpected GPU fallback: ', reason)
  }
  image
}

test_that('Metal SSS conserves white light and supports camera-inside paths', {
  scene = sphere(
    material = subsurface(sigma_a = 0, sigma_s = 3, refraction = 1.3)
  )
  gpu = suppressWarnings(render_wavefront_body(scene, samples = 16, depth = 1))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal is unavailable')
  }
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_gt(attr(gpu, 'wavefront')$scattering_events, 0)
  expect_equal(mean(gpu[,, 1:3]), 1, tolerance = .015)
  # A matched boundary has a constant-radiance interior as well as exterior.
  scene = sphere(
    material = subsurface(sigma_a = 0, sigma_s = 3, refraction = 1)
  )
  inside = render_wavefront_body(
    scene,
    samples = 64,
    lookfrom = c(0, 0, 0),
    lookat = c(0, 0, -1)
  )
  expect_equal(attr(inside, 'wavefront')$discarded_paths, 0)
  expect_equal(mean(inside[,, 1:3]), 1, tolerance = .03)
})

test_that('Metal RGB free flights retain vacuum and absorption channels', {
  sigma = c(0, .4, .8)
  scene = cube(
    xwidth = 10,
    ywidth = 10,
    zwidth = 1,
    material = subsurface(sigma_a = sigma, sigma_s = 0, refraction = 1)
  )
  gpu = suppressWarnings(render_wavefront_body(scene, samples = 256))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal is unavailable')
  }
  expect_equal(apply(gpu[,, 1:3], 3, mean), exp(-sigma), tolerance = .025)
  expect_equal(attr(gpu, 'wavefront')$scattering_events, 0)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
})

test_that('Metal priority displaces SSS with higher priority glass', {
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
  )
  scene = body |>
    add_object(cube(
      xwidth = 11,
      ywidth = 11,
      zwidth = 1,
      material = dielectric(
        refraction = 1,
        attenuation = rep(.3, 3),
        priority = 1
      )
    ))
  gpu = suppressWarnings(render_wavefront_body(scene, samples = 256))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal is unavailable')
  }
  expect_equal(mean(gpu[,, 1:3]), exp(-.7), tolerance = .025)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  # With glass completely enclosing the body, even strong SSS is inactive.
  hidden = sphere(
    radius = 1.1,
    material = dielectric(refraction = 1, priority = 0)
  ) |>
    add_object(sphere(
      material = subsurface(sigma_a = 10, sigma_s = 20, priority = 1)
    ))
  image = render_wavefront_body(hidden, samples = 4)
  expect_equal(attr(image, 'wavefront')$scattering_events, 0)
  expect_equal(mean(image[,, 1:3]), 1, tolerance = 1e-5)
})

test_that('Metal random walks agree with CPU transport for rough anisotropic bodies', {
  scene = sphere(
    material = subsurface(
      sigma_a = c(.1, .2, .4),
      sigma_s = c(2, 3, 4),
      g = .4,
      roughness = .2
    )
  ) |>
    add_infinite_light(disk_light(
      direction = c(-1, 2, 1),
      angular_diameter = 50,
      intensity = 5
    ))
  gpu = suppressWarnings(render_wavefront_body(scene, samples = 256))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal is unavailable')
  }
  cpu = render_wavefront_body(scene, 'nee', samples = 256)
  expect_lt(
    max(abs(apply(gpu[,, 1:3], 3, mean) - apply(cpu[,, 1:3], 3, mean))),
    .06
  )
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
})

test_that('Metal triangulates analytic shapes and cylindrical/ribbon curves', {
  shapes = list(
    sphere(),
    ellipsoid(a = 1, b = .8, c = .6),
    cylinder(),
    disk(),
    bezier_curve(),
    bezier_curve(type = 'ribbon')
  )
  for (scene in shapes) {
    gpu = suppressWarnings(render_wavefront_body(scene, samples = 2, width = 8))
    if (!isTRUE(attr(gpu, 'wavefront')$used)) {
      skip('Metal is unavailable')
    }
    expect_gt(attr(gpu, 'wavefront')$tessellated_primitives, 0)
    expect_true(all(is.finite(gpu)))
  }
})

test_that('closed cylinders support SSS and open cylinders are rejected', {
  material = subsurface(sigma_a = 0, sigma_s = 2, refraction = 1.3)
  scene = cylinder(length = 2, material = material)
  gpu = suppressWarnings(render_wavefront_body(scene, samples = 16, depth = 1))
  if (isTRUE(attr(gpu, 'wavefront')$used)) {
    expect_equal(mean(gpu[,, 1:3]), 1, tolerance = .02)
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  }
  cpu = render_wavefront_body(scene, 'nee', samples = 16, depth = 1)
  expect_equal(mean(cpu[,, 1:3]), 1, tolerance = .02)
  expect_error(
    render_wavefront_body(cylinder(capped = FALSE, material = material), 'nee'),
    'closed shape'
  )
  expect_error(
    render_wavefront_body(cylinder(phi_max = 180, material = material), 'nee'),
    'closed shape'
  )
})
