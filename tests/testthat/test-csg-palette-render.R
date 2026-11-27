# Surface palettes must reach the correct shader queue and light proposal, not
# just replace the fallback's base color.
render_palette_scene = function(scene, backend, ...) {
  set.seed(62)
  image = render_scene(
    scene,
    integrator_type = backend,
    width = 20,
    height = 20,
    samples = 256,
    sample_method = 'sobol',
    lookfrom = c(0, 0, 4),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.15, .15),
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white',
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    min_variance = 0,
    ...
  )
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    reason = attr(image, 'wavefront')$fallback
    if (grepl('not enabled|device|macOS 11', reason)) {
      skip('Metal unavailable')
    }
    stop('Unexpected CSG material fallback: ', reason)
  }
  if (backend == 'metal') {
    expect_equal(attr(image, 'wavefront')$discarded_paths, 0)
    expect_equal(attr(image, 'wavefront')$csg_failed_queries, 0)
  }
  expect_true(all(is.finite(image)))
  apply(image[,, 1:3], 3, mean)
}

test_that('CSG component palettes dispatch OpenPBR and preserve complete shaders', {
  shape = csg_combine(
    csg_sphere(),
    csg_sphere(
      material = openpbr(base_color = '#cc5544', specular_roughness = .3)
    ),
    'mix',
    radius = .6
  )
  scene = csg_object(shape, material = metal('white', fuzz = .1))
  cpu = render_palette_scene(scene, 'nee')
  gpu = render_palette_scene(scene, 'metal')
  expect_equal(gpu, cpu, tolerance = .045)
  expect_gt(cpu[1], cpu[3])
})

test_that('CSG component emission is sampled and inherited textures still load', {
  shape = csg_combine(
    csg_sphere(),
    csg_sphere(material = light('blue', intensity = 2)),
    'mix',
    radius = .5
  )
  scene = csg_object(shape, material = diffuse('black'))
  for (backend in c('nee', 'metal')) {
    expect_equal(
      render_palette_scene(scene, backend),
      c(0, 0, 1),
      tolerance = .05
    )
  }
  image = array(0, c(3, 3, 3))
  image[,, 2] = 1
  shape = csg_sphere(material = diffuse(image_texture = image))
  scene = csg_object(shape, material = diffuse('red')) |>
    add_object(sphere(x = 20, material = diffuse(image_texture = image)))
  # The following ordinary primitive exercises texture indices after a palette.
  for (backend in c('nee', 'metal')) {
    expect_equal(
      render_palette_scene(scene, backend),
      c(0, 1, 0),
      tolerance = .04
    )
  }
})

test_that('a component emitter illuminates other objects through light sampling', {
  emitter = csg_combine(
    csg_sphere(radius = .5),
    csg_sphere(radius = .5, material = light(intensity = 10)),
    'mix'
  )
  scene = sphere(material = diffuse('white')) |>
    add_object(csg_object(
      emitter,
      x = 1.2,
      y = 1.5,
      z = 2,
      material = diffuse('black')
    ))
  # At depth one there is no BSDF bounce to accidentally discover the emitter:
  # this isolates registration and material evaluation in the light proposal.
  cpu = render_palette_scene(scene, 'nee', max_depth = 1)
  gpu = render_palette_scene(scene, 'metal', max_depth = 1)
  expect_equal(gpu, cpu, tolerance = .08)
  expect_gt(mean(cpu), 1.05)
})
