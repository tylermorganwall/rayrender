render_csg_test = function(scene, backend = 'metal', settings = list()) {
  set.seed(903)
  image = do.call(
    render_scene,
    modifyList(
      list(
        scene = scene,
        integrator_type = backend,
        width = 32,
        height = 32,
        samples = 128,
        sample_method = 'sobol',
        lookfrom = c(0, 0, 4),
        lookat = c(0, 0, 0),
        fov = 0,
        ortho_dimensions = c(2.6, 2.6),
        max_depth = 16,
        ambient_light = TRUE,
        backgroundlow = 'white',
        backgroundhigh = 'white',
        preview = FALSE,
        plot_scene = FALSE,
        progress = FALSE,
        min_variance = 0,
        denoise = FALSE,
        bloom = FALSE
      ),
      settings
    )
  )
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    reason = attr(image, 'wavefront')$fallback
    if (grepl('not enabled|device|macOS 11', reason)) {
      skip('Metal unavailable')
    }
    stop('Unexpected CSG fallback: ', reason)
  }
  if (backend == 'metal') {
    expect_equal(attr(image, 'wavefront')$csg_failed_queries, 0)
    expect_true(all(is.finite(image)))
  }
  image
}

test_that('native Metal CSG preserves procedural textures and Boolean surfaces', {
  shape = csg_combine(csg_sphere(), csg_box(width = c(1.5, 1.5, 3)), 'subtract')
  scene = csg_object(
    shape,
    material = diffuse(texture_checker(
      'coral',
      'steelblue',
      texture_coordinates('object', scale = 4)
    ))
  )
  gpu = render_csg_test(scene)
  cpu = render_csg_test(scene, 'nee')
  expect_equal(attr(gpu, 'wavefront')$triangles, 0)
  expect_equal(attr(gpu, 'wavefront')$analytic_primitives, 1)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .035)
})

test_that('CSG media respect priority, interior cameras, and thin shells', {
  # Exact path lengths test priority independently of either implementation.
  scene = csg_object(
    csg_box(width = c(10, 10, 2)),
    material = subsurface(
      sigma_a = .4,
      sigma_s = 0,
      refraction = 1,
      priority = 2
    )
  ) |>
    add_object(csg_object(
      csg_box(width = c(11, 11, 1)),
      material = dielectric(
        refraction = 1,
        attenuation = rep(.3, 3),
        priority = 1
      )
    ))
  for (backend in c('metal', 'nee')) {
    image = render_csg_test(
      scene,
      backend,
      list(ortho_dimensions = c(.5, .5), samples = 256)
    )
    expect_equal(mean(image[,, 1:3]), exp(-.7), tolerance = .025)
  }
  scene = csg_object(
    csg_onion(csg_sphere(), thickness = .002),
    material = dielectric(refraction = 1, attenuation = rep(40, 3))
  )
  image = render_csg_test(scene, settings = list(ortho_dimensions = c(.1, .1)))
  expect_equal(mean(image[,, 1:3]), exp(-.32), tolerance = .02)
  scene = csg_object(
    csg_sphere(radius = 3),
    material = subsurface(
      sigma_a = 0,
      sigma_s = 2,
      refraction = 1
    )
  )
  image = render_csg_test(
    scene,
    settings = list(
      lookfrom = c(0, 0, 0),
      lookat = c(0, 0, -1),
      samples = 32
    )
  )
  expect_equal(mean(image[,, 1:3]), 1, tolerance = .025)
  expect_equal(attr(image, 'wavefront')$discarded_paths, 0)
})

test_that('CSG diffusion and transformed placements stay on Metal', {
  scene = csg_object(
    csg_combine(
      csg_sphere(x = -.3, radius = .65),
      csg_sphere(x = .3, radius = .65),
      'blend',
      radius = .3
    ),
    material = subsurface_diffusion(color = c(.8, .5, .3), radius = .2)
  )
  gpu = render_csg_test(scene, settings = list(samples = 1024))
  cpu = render_csg_test(scene, 'nee', list(samples = 1024))
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .04)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  ball_material = subsurface_diffusion(color = c(.8, .5, .3), radius = .2)
  csg_ball = render_csg_test(
    csg_object(csg_sphere(), material = ball_material),
    settings = list(samples = 512)
  )
  analytic_ball = render_csg_test(
    sphere(material = ball_material),
    settings = list(samples = 512)
  )
  expect_equal(
    apply(csg_ball[,, 1:3], 3, mean),
    apply(analytic_ball[,, 1:3], 3, mean),
    tolerance = .01
  )
  moving = scene |>
    animate_objects(end_position = c(.3, 0, 0), end_angle = c(0, 35, 0))
  image = render_csg_test(moving, settings = list(samples = 16))
  expect_true(all(is.finite(image)))
  expect_equal(attr(image, 'wavefront')$discarded_paths, 0)
})

test_that('CSG emitters sample holes without biased rejection', {
  scene = sphere(radius = .8, material = diffuse('grey70')) |>
    add_object(csg_object(
      csg_translate(
        csg_torus(radius = .65, minor_radius = .22),
        x = 1.5,
        y = 2,
        z = 1
      ),
      material = light(intensity = 15)
    ))
  gpu = render_csg_test(
    scene,
    settings = list(ambient_light = FALSE, samples = 512)
  )
  cpu = render_csg_test(
    scene,
    'nee',
    list(ambient_light = FALSE, samples = 512)
  )
  expect_gt(mean(gpu[,, 1:3]), .02)
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .035)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
})

test_that('CSG explicit media and unsupported nesting have explicit behavior', {
  scene = csg_object(csg_box(width = c(4, 4, 2))) |>
    set_medium(homogeneous_medium(sigma_a = .4, sigma_s = 0))
  image = render_csg_test(scene, settings = list(ortho_dimensions = c(.5, .5)))
  expect_equal(mean(image[,, 1:3]), exp(-.8), tolerance = .02)
  object = csg_sphere()
  for (i in seq_len(33)) {
    object = csg_round(object, radius = .001)
  }
  # Test fallback itself without the helper's requirement for a used GPU backend.
  expect_warning(
    image <- render_scene(
      csg_object(object),
      integrator_type = 'metal',
      width = 4,
      height = 4,
      samples = 1,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      bloom = FALSE,
      denoise = FALSE
    ),
    'CSG.*nesting.*CPU NEE'
  )
  expect_false(attr(image, 'wavefront')$used)
})
