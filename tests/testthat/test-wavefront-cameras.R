test_that('Metal camera projections and moving geometry preserve CPU image energy', {
  render = function(scene, backend = 'metal', settings = list()) {
    set.seed(104)
    do.call(
      render_scene,
      modifyList(
        list(
          scene = scene,
          integrator_type = backend,
          width = 40,
          height = 32,
          samples = 128,
          lookfrom = c(0, 1, 6),
          lookat = c(0, 1, 0),
          iso = 100,
          preview = FALSE,
          plot_scene = FALSE,
          progress = FALSE,
          denoise = FALSE,
          bloom = FALSE,
          ambient_light = FALSE,
          backgroundhigh = 'black',
          backgroundlow = 'black',
          max_depth = 3
        ),
        settings
      )
    )
  }
  base = sphere(x = -1, y = .7, radius = .7, material = light('red')) |>
    add_object(sphere(x = 1, y = 1.2, radius = .8, material = light('blue')))
  probe = suppressWarnings(render(base, settings = list(samples = 1)))
  if (!isTRUE(attr(probe, 'wavefront')$used)) {
    skip('Metal unavailable')
  }
  configurations = list(
    list(fov = 360),
    list(camera_description_file = '50mm', aperture = 4, focal_distance = 6)
  )
  for (settings in configurations) {
    gpu = render(base, settings = settings)
    cpu = render(base, 'nee', settings)
    expect_true(isTRUE(attr(gpu, 'wavefront')$used))
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
    expect_lt(abs(mean(gpu[,, 1:3]) / mean(cpu[,, 1:3]) - 1), .04)
    expect_lt(sqrt(mean((gpu[,, 1:3] - cpu[,, 1:3])^2)), .01)
  }
  # A moving emitter checks both hit geometry and frame accumulation. Infinite
  # shutter speed must reproduce its reference pose, rather than its midpoint.
  moving = animate_objects(
    base,
    end_position = c(2, 0, 0),
    end_angle = c(0, 30, 0)
  )
  frozen = render(moving, settings = list(shutter_speed = Inf))
  still = render(base)
  expect_equal(as.numeric(frozen), as.numeric(still), tolerance = 1e-6)
  gpu = render(moving, settings = list(shutter_speed = 1))
  cpu = render(moving, 'nee', list(shutter_speed = 1))
  expect_true(isTRUE(attr(gpu, 'wavefront')$used))
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_lt(abs(mean(gpu[,, 1:3]) / mean(cpu[,, 1:3]) - 1), .03)
  expect_gt(mean(abs(gpu[,, 1:3] - frozen[,, 1:3])), .005)

  # All primitives share an affine motion: Metal keeps a tight reference BVH
  # and transforms rays. Include a mesh, a quadric, and nonuniform scaling so
  # normalization of the transformed ray direction cannot go unnoticed.
  skip_if_not_installed('rayvertex')
  assembly = sphere(x = -1, y = .7, radius = .7, material = light('red')) |>
    add_object(raymesh_model(
      rayvertex::cube_mesh(),
      x = 1,
      y = 1,
      material = light('blue'),
      override_material = TRUE
    )) |>
    animate_objects(
      end_position = c(.4, .2, 0),
      end_angle = c(0, 155, 0),
      end_scale = c(.8, 1.2, 1.1)
    )
  # Bright moving silhouettes need enough temporal samples to compare independent
  # CPU/GPU sampling streams; keep the same absolute image tolerance below.
  settings = list(shutter_speed = 1, samples = 1024, sample_method = 'sobol')
  gpu = render(assembly, settings = settings)
  cpu = render(assembly, 'nee', settings)
  expect_true(isTRUE(attr(gpu, 'wavefront')$used))
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_lt(abs(mean(gpu[,, 1:3]) / mean(cpu[,, 1:3]) - 1), .04)
  expect_lt(sqrt(mean((gpu[,, 1:3] - cpu[,, 1:3])^2)), .03)
})

test_that('Metal motion refits shared meshes and transformed medium boundaries', {
  skip_if_not_installed('rayvertex')
  mesh = raymesh_model(
    rayvertex::cube_mesh(),
    material = diffuse('steelblue'),
    override_material = TRUE
  )
  copies = create_instances(
    mesh,
    x = c(-1, 1),
    y = c(.5, .7),
    scale_x = c(.6, .8)
  )
  moving = animate_objects(
    copies,
    end_position = c(.5, .4, 0),
    end_angle = c(0, 45, 0)
  )
  fog = set_medium(
    cube(width = 2),
    homogeneous_medium(sigma_s = .2, sigma_a = .02)
  ) |>
    animate_objects(start_position = c(-.5, 1, 0), end_position = c(.5, 1, 0))
  scene = moving |>
    add_object(fog) |>
    add_object(xz_rect(
      y = -1,
      xwidth = 8,
      zwidth = 8,
      material = diffuse('grey40')
    )) |>
    add_infinite_light(disk_light(
      direction = c(-1, 2, 1),
      intensity = 6,
      angular_diameter = 30
    ))
  render = function(backend) {
    set.seed(704)
    render_scene(
      scene,
      integrator_type = backend,
      width = 32,
      height = 32,
      samples = 64,
      lookfrom = c(0, 2, 6),
      lookat = c(0, .5, 0),
      shutter_speed = 1,
      max_depth = 4,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      ambient_light = FALSE
    )
  }
  gpu = suppressWarnings(render('metal'))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal unavailable')
  }
  cpu = render('nee')
  expect_gt(attr(gpu, 'wavefront')$mesh_instances, 0)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_lt(abs(mean(gpu[,, 1:3]) / mean(cpu[,, 1:3]) - 1), .05)
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .03)
})
