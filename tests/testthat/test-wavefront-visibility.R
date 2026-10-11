# Coverage is independent of the RGB walk. Test full-frame analytic values as
# well as camera, reflected, and directly sampled views of invisible emitters.
render_wavefront_visibility = function(
  scene,
  backend = 'metal',
  transparent = TRUE,
  samples = 64,
  ...
) {
  set.seed(618)
  image = withCallingHandlers(
    render_scene(
      scene,
      integrator_type = backend,
      transparent_background = transparent,
      width = 24,
      height = 24,
      samples = samples,
      max_depth = 4,
      sample_method = 'sobol',
      min_variance = 0,
      parallel = FALSE,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = 'raw',
      lookfrom = c(0, 0, 4),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(1, 1),
      ...
    ),
    warning = function(w) {
      # Unsupported platforms still exercise the CPU package; omit only the
      # expected availability warning before skipping these GPU comparisons.
      if (
        backend == 'metal' &&
          grepl('not enabled|Metal device|macOS 11', conditionMessage(w))
      ) {
        invokeRestart('muffleWarning')
      }
    }
  )
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    expect_match(
      attr(image, 'wavefront')$fallback,
      'not enabled|device|macOS 11'
    )
    skip('Metal ray tracing unavailable')
  }
  expect_true(all(is.finite(image)))
  if (backend == 'metal') {
    expect_equal(attr(image, 'wavefront')$discarded_paths, 0)
  }
  image
}

test_that('Metal hides invisible area lights but retains their illumination', {
  scene = xy_rect(xwidth = 20, ywidth = 20, material = diffuse('white'))
  emitter = xy_rect(
    z = 2,
    xwidth = 4,
    ywidth = 4,
    angle = c(0, 180, 0),
    material = light(intensity = 2, invisible = TRUE)
  )
  scene = add_object(scene, emitter)
  gpu = render_wavefront_visibility(scene, ambient_light = FALSE)
  cpu = render_wavefront_visibility(scene, 'nee', ambient_light = FALSE)
  expect_true(all(gpu[,, 4] == 1))
  expect_gt(mean(gpu[,, 1:3]), .1)
  expect_equal(mean(gpu[,, 1:3]), mean(cpu[,, 1:3]), tolerance = .06)

  # An invisible emitter by itself leaves the sky visible (or alpha zero).
  clear = render_wavefront_visibility(
    emitter,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  expect_equal(max(abs(clear)), 0)
  sky = render_wavefront_visibility(
    emitter,
    transparent = FALSE,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  expect_equal(mean(sky), 1, tolerance = 1e-6)
  # Analytic emitters have front and back intersections; skipping both must
  # terminate normally without leaving an opaque silhouette in alpha.
  hidden_sphere = render_wavefront_visibility(
    sphere(material = light(invisible = TRUE)),
    samples = 4,
    ambient_light = TRUE
  )
  expect_equal(max(abs(hidden_sphere)), 0)

  # The smooth mirror must not reveal the source behind the camera. A diffuse
  # surface above did see its emission, so hiding it is not unconditional.
  mirror = xy_rect(xwidth = 20, ywidth = 20, material = metal('white')) |>
    add_object(xy_rect(
      z = 6,
      xwidth = 20,
      ywidth = 20,
      angle = c(0, 180, 0),
      material = light(intensity = 4, invisible = TRUE)
    ))
  reflected = render_wavefront_visibility(mirror, ambient_light = FALSE)
  expect_equal(max(abs(reflected[,, 1:3])), 0)
  expect_true(all(reflected[,, 4] == 1))

  # A connection to a different source must cross the invisible area light,
  # rather than stopping at the first geometry returned by an any-hit query.
  through = xy_rect(xwidth = 20, ywidth = 20, material = diffuse('white')) |>
    add_object(xy_rect(
      z = 2,
      xwidth = 20,
      ywidth = 20,
      material = light(intensity = 0, invisible = TRUE)
    )) |>
    add_light(point_light(position = c(0, 0, 6), intensity = 40))
  gpu = render_wavefront_visibility(through, ambient_light = FALSE)
  cpu = render_wavefront_visibility(through, 'nee', ambient_light = FALSE)
  expect_gt(mean(gpu[,, 1:3]), .1)
  expect_equal(mean(gpu[,, 1:3]), mean(cpu[,, 1:3]), tolerance = .02)
})

test_that('Metal camera coverage preserves surfaces, cutouts and environment illumination', {
  wall = xy_rect(xwidth = 20, ywidth = 20, material = diffuse('grey50'))
  opaque = render_wavefront_visibility(
    wall,
    transparent = FALSE,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  rgba = render_wavefront_visibility(
    wall,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  expect_equal(
    as.numeric(rgba[,, 1:3]),
    as.numeric(opaque[,, 1:3]),
    tolerance = 1e-6
  )
  expect_true(all(rgba[,, 4] == 1))
  glass = render_wavefront_visibility(
    cube(width = 2, material = dielectric(refraction = 1)),
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  expect_true(all(glass[,, 4] == 1))

  alpha = tempfile(fileext = '.png')
  on.exit(unlink(alpha))
  texture = array(1, c(2, 2, 4))
  texture[,, 4] = .35
  png::writePNG(texture, alpha)
  cutout = xy_rect(
    xwidth = 20,
    ywidth = 20,
    material = diffuse('red', alpha_texture = alpha)
  )
  gpu = render_wavefront_visibility(cutout, ambient_light = TRUE)
  cpu = render_wavefront_visibility(cutout, 'nee', ambient_light = TRUE)
  expect_equal(mean(gpu[,, 4]), .35, tolerance = .005)
  expect_equal(mean(gpu[,, 4]), mean(cpu[,, 4]), tolerance = .015)
  blocked = render_wavefront_visibility(
    add_object(cutout, cube(z = -2, width = 2)),
    ambient_light = TRUE
  )
  expect_true(all(blocked[,, 4] == 1))
})

test_that('Metal transparent volumes match analytic extinction and priority', {
  sigma = c(.2, .6, 1)
  slab = set_medium(
    cube(width = 2),
    homogeneous_medium(sigma_a = sigma, sigma_s = .3)
  )
  gpu = render_wavefront_visibility(
    slab,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  cpu = render_wavefront_visibility(
    slab,
    'nee',
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  expect_equal(
    mean(gpu[,, 4]),
    1 - mean(exp(-2 * (sigma + .3))),
    tolerance = .001
  )
  expect_equal(mean(gpu[,, 4]), mean(cpu[,, 4]), tolerance = .015)
  nested = set_medium(
    cube(width = 4),
    homogeneous_medium(sigma_a = .1, sigma_s = 0)
  ) |>
    add_object(set_medium(
      cube(width = 2),
      homogeneous_medium(sigma_a = .7, sigma_s = 0)
    ))
  image = render_wavefront_visibility(nested, ambient_light = TRUE)
  expect_equal(mean(image[,, 4]), 1 - exp(-1.6), tolerance = .001)
  grid = set_medium(
    cube(width = 2),
    grid_medium(
      array(.4, c(3, 3, 3)),
      sigma_a = .8,
      sigma_s = .2,
      bounds = rbind(c(-1, -1, -1), c(1, 1, 1))
    )
  )
  image = render_wavefront_visibility(grid, ambient_light = TRUE)
  expect_equal(mean(image[,, 4]), 1 - exp(-.8), tolerance = .015)
})

test_that('Metal transparent Prague skies retain finite haze coverage', {
  if (!any(skymodelr::list_sky_data()$file == 'SkyModelDataset.dat')) {
    skip('Full-altitude Prague data not installed')
  }
  for (haze in c(FALSE, TRUE)) {
    # Geometry is outside the view. The camera sees only the sky, so without
    # finite haze both RGB and alpha must be zero; with haze compare the CPU.
    scene = sphere(x = 10) |>
      add_infinite_light(sky_light(
        elevation = 30,
        azimuth = 160,
        sun = FALSE,
        moon = FALSE,
        haze = haze,
        sampling_resolution = 16,
        prague_rgb_correction = FALSE
      ))
    gpu = render_wavefront_visibility(scene, samples = 4)
    cpu = render_wavefront_visibility(scene, 'nee', samples = 4)
    if (!haze) {
      expect_equal(max(abs(gpu)), 0)
    }
    expect_equal(mean(gpu[,, 4]), mean(cpu[,, 4]), tolerance = .015)
    expect_equal(mean(gpu[,, 1:3]), mean(cpu[,, 1:3]), tolerance = .025)
  }
})
