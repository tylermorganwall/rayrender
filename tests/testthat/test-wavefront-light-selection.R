test_that('Metal light BVHs preserve lighting from moving mesh instances', {
  skip_if_not_installed('rayvertex')
  source = raymesh_model(
    rayvertex::cube_mesh(),
    material = light(intensity = 3),
    override_material = TRUE
  )
  lamps = create_instances(
    source,
    x = c(-2, 0, 2),
    y = c(2, 3, 2),
    scale_x = c(.3, .5, .2),
    scale_y = c(.2, .3, .4),
    scale_z = c(.4, .3, .2)
  )
  render = function(moving, backend) {
    scene = if (moving) {
      animate_objects(
        lamps,
        end_position = c(.7, 0, 0),
        end_angle = c(0, 20, 0)
      )
    } else {
      lamps
    }
    scene = add_object(
      scene,
      xz_rect(xwidth = 10, zwidth = 10, material = diffuse('grey70'))
    ) |>
      add_infinite_light(disk_light(
        direction = c(-1, 2, 1),
        intensity = .2,
        angular_diameter = 30
      ))
    set.seed(191)
    image = suppressWarnings(render_scene(
      scene,
      integrator_type = backend,
      width = 24,
      height = 24,
      samples = 512,
      sample_method = 'sobol',
      lookfrom = c(0, 4, 7),
      lookat = c(0, 0, 0),
      fov = 40,
      aperture = 0,
      shutter_speed = 1,
      ambient_light = FALSE,
      max_depth = 3,
      min_variance = 0,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = 'raw'
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
  for (moving in c(FALSE, TRUE)) {
    gpu = render(moving, 'metal')
    cpu = render(moving, 'nee')
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
    expect_gt(attr(gpu, 'wavefront')$mesh_instances, 0)
    expect_lt(abs(mean(gpu[,, 1:3]) / mean(cpu[,, 1:3]) - 1), .04)
    # Exclude directly visible lamps from the floor-only transport comparison.
    floor = apply(cpu[,, 1:3], c(1, 2), max) < .8
    expect_lt(abs(mean(gpu[,, 1][floor]) / mean(cpu[,, 1][floor]) - 1), .04)
  }
})
