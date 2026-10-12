test_that('Metal directly intersects all curve types with hair and textured materials', {
  render = function(scene, backend = 'metal', settings = list()) {
    set.seed(123)
    image = suppressWarnings(do.call(
      render_scene,
      modifyList(
        list(
          scene = scene,
          integrator_type = backend,
          width = 40,
          height = 40,
          samples = 512,
          sample_method = 'sobol',
          lookfrom = c(0, .1, 4),
          lookat = c(0, 0, 0),
          fov = 35,
          max_depth = 4,
          min_variance = 0,
          preview = FALSE,
          plot_scene = FALSE,
          progress = FALSE,
          denoise = FALSE,
          bloom = FALSE,
          ambient_light = FALSE,
          tonemap = 'raw'
        ),
        settings
      )
    ))
    if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
      expect_match(
        attr(image, 'wavefront')$fallback,
        'not enabled|device|macOS 11'
      )
      skip('Metal is unavailable')
    }
    image
  }
  for (type in c('flat', 'cylinder', 'ribbon')) {
    scene = bezier_curve(
      p1 = c(-.5, -1, 0),
      p2 = c(.5, -.3, .1),
      p3 = c(-.5, .3, -.1),
      p4 = c(.5, 1, 0),
      width = .3,
      width_end = .12,
      type = type,
      normal = c(0, 0, 1),
      normal_end = c(.2, 0, 1),
      material = hair(color = c(.5, .2, .08), beta_m = .4, beta_n = .5)
    ) |>
      add_infinite_light(disk_light(
        direction = c(-1, 2, 2),
        intensity = 3,
        angular_diameter = 35
      ))
    gpu = render(scene)
    cpu = render(scene, 'nee')
    report = attr(gpu, 'wavefront')
    expect_equal(report$triangles, 0)
    expect_equal(report$tessellated_primitives, 0)
    expect_equal(report$analytic_primitives, 8)
    expect_equal(report$discarded_paths, 0)
    expect_gt(mean(gpu[,, 1:3]), .0001)
    expect_lt(abs(mean(gpu[,, 1:3]) / mean(cpu[,, 1:3]) - 1), .08)
  }
  # Pure emission isolates transformed silhouettes/UV lookup from BSDF sampling.
  scene = bezier_curve(
    p1 = c(-.4, -1, 0),
    p2 = c(.7, -.3, 0),
    p3 = c(-.7, .3, 0),
    p4 = c(.4, 1, 0),
    width = .3,
    material = light(
      color = 'red',
      gradient_color = 'blue',
      gradient_transpose = TRUE,
      gradient_type = 'rgb',
      intensity = 1
    ),
    type = 'flat',
    scale = c(-1.2, .8, 1),
    angle = c(0, 15, 0)
  )
  gpu = render(scene)
  cpu = render(scene, 'nee')
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .004)
  expect_gt(sum(gpu[,, 1]), 1)
  expect_gt(sum(gpu[,, 3]), 1)

  # Exercise both whole-assembly reference traversal and partial-motion refits.
  moving = animate_objects(
    scene,
    end_position = c(.5, 0, 0),
    end_angle = c(0, 25, 0)
  )
  for (mixed in c(FALSE, TRUE)) {
    world = if (mixed) {
      add_object(
        moving,
        sphere(x = .8, y = -.7, radius = .15, material = diffuse())
      )
    } else {
      moving
    }
    gpu = render(world, settings = list(shutter_speed = 1))
    cpu = render(world, 'nee', list(shutter_speed = 1))
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
    expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .006)
    expect_lt(abs(mean(gpu[,, 1:3]) / mean(cpu[,, 1:3]) - 1), .04)
  }

  # A curve emitter must illuminate another surface, not just appear bright
  # when the camera hits it. No environment light supplies the floor here.
  lamp = bezier_curve(
    p1 = c(-.6, 1.3, 0),
    p2 = c(-.2, 1.3, 0),
    p3 = c(.2, 1.3, 0),
    p4 = c(.6, 1.3, 0),
    width = .2,
    type = 'flat',
    material = light(intensity = 20)
  ) |>
    add_object(xz_rect(
      xwidth = 6,
      zwidth = 6,
      y = -.8,
      material = diffuse('white')
    ))
  lit = render(lamp, settings = list(max_depth = 1))
  expect_equal(attr(lit, 'wavefront')$discarded_paths, 0)
  expect_gt(mean(lit[26:40, , 1:3]), .001)
})
