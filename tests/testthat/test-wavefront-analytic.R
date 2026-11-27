render_wavefront_analytic = function(
  scene,
  backend = 'metal',
  lookfrom = c(0, 0, 4),
  ambient = TRUE,
  samples = 128
) {
  set.seed(508)
  image = suppressWarnings(render_scene(
    scene,
    integrator_type = backend,
    width = 32,
    height = 32,
    samples = samples,
    lookfrom = lookfrom,
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(2.4, 2.4),
    max_depth = 12,
    min_variance = 0,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    bloom = FALSE,
    ambient_light = ambient,
    backgroundlow = if (ambient) 'white' else 'black',
    backgroundhigh = if (ambient) 'white' else 'black'
  ))
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    expect_match(
      attr(image, 'wavefront')$fallback,
      'not enabled|device|macOS 11'
    )
    skip('Metal ray tracing is unavailable')
  }
  image
}

test_that('Metal spheres and transformed ellipsoids remain analytic', {
  scenes = list(
    sphere(
      material = diffuse(texture_checker(
        'coral',
        'navy',
        texture_coordinates('uv', scale = 8)
      ))
    ),
    ellipsoid(
      a = .9,
      b = .65,
      c = .45,
      angle = c(20, 35, 12),
      material = diffuse(texture_direction_mix('gold', 'steelblue'))
    ),
    create_instances(
      sphere(material = diffuse('coral')),
      x = c(-.55, .55),
      scale_x = c(.4, -.4),
      scale_y = c(.8, .8),
      scale_z = c(.3, .3)
    )
  )
  for (scene in scenes) {
    gpu = render_wavefront_analytic(scene)
    cpu = render_wavefront_analytic(scene, 'nee')
    expect_equal(attr(gpu, 'wavefront')$triangles, 0)
    expect_gt(attr(gpu, 'wavefront')$analytic_primitives, 0)
    expect_equal(attr(gpu, 'wavefront')$tessellated_primitives, 0)
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
    expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .025)
  }
})

test_that('transmitted rays retain the opposite intersection of the same sphere', {
  scene = sphere(material = dielectric(color = c(.25, .5, .8), refraction = 1))
  gpu = render_wavefront_analytic(scene)
  cpu = render_wavefront_analytic(scene, 'nee')
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .015)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  # An interior camera must initialize membership and cross the outward root.
  scene = sphere(radius = 3, material = dielectric(refraction = 1.4))
  gpu = render_wavefront_analytic(scene, lookfrom = c(0, 0, .5))
  cpu = render_wavefront_analytic(scene, 'nee', lookfrom = c(0, 0, .5))
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .02)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
})

test_that('analytic spherical lights use a consistent direct and BSDF-hit density', {
  scene = sphere(radius = .8, material = diffuse('grey70')) |>
    add_object(sphere(
      x = 2,
      y = 3,
      z = 2,
      radius = .65,
      material = light(intensity = 12)
    ))
  gpu = render_wavefront_analytic(scene, ambient = FALSE)
  cpu = render_wavefront_analytic(scene, 'nee', ambient = FALSE)
  expect_equal(attr(gpu, 'wavefront')$analytic_primitives, 2)
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .025)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
})

test_that('ellipsoidal emitter sampling agrees with a triangulated reference', {
  skip_if_not_installed('rayvertex')
  material = light(intensity = 12)
  object = sphere(radius = .8, material = diffuse('grey70'))
  analytic = add_object(
    object,
    ellipsoid(
      x = 2,
      y = 3,
      z = 2,
      a = .9,
      b = .55,
      c = .4,
      angle = c(10, 20, 30),
      material = material
    )
  )
  mesh = rayvertex::sphere_mesh(scale = c(.9, .55, .4))
  reference = add_object(
    object,
    raymesh_model(
      mesh,
      x = 2,
      y = 3,
      z = 2,
      angle = c(10, 20, 30),
      material = material,
      override_material = TRUE,
      calculate_consistent_normals = FALSE
    )
  )
  gpu = render_wavefront_analytic(analytic, ambient = FALSE, samples = 256)
  triangulated = render_wavefront_analytic(
    reference,
    ambient = FALSE,
    samples = 256
  )
  expect_lt(mean(abs(gpu[,, 1:3] - triangulated[,, 1:3])), .025)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
})

test_that('large analytic ground spheres do not shadow their sampled lights', {
  scene = generate_ground(material = diffuse('grey70')) |>
    add_object(sphere(
      x = -2,
      y = 4,
      z = 2,
      radius = 1.5,
      material = light(intensity = 10)
    ))
  images = lapply(c('metal', 'nee'), function(backend) {
    set.seed(361)
    suppressWarnings(render_scene(
      scene,
      integrator_type = backend,
      width = 48,
      height = 48,
      samples = 128,
      lookfrom = c(3, 2.5, 5),
      lookat = c(0, -1, 0),
      fov = 0,
      ortho_dimensions = c(5, 5),
      ambient_light = FALSE,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      bloom = FALSE
    ))
  })
  if (!isTRUE(attr(images[[1]], 'wavefront')$used)) {
    expect_match(
      attr(images[[1]], 'wavefront')$fallback,
      'not enabled|device|macOS 11'
    )
    skip('Metal ray tracing is unavailable')
  }
  expect_lt(mean(abs(images[[1]][,, 1:3] - images[[2]][,, 1:3])), .025)
  expect_equal(attr(images[[1]], 'wavefront')$discarded_paths, 0)
})
