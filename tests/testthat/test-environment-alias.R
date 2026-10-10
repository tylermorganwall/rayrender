test_that('image alias sampling preserves energy with blue-noise and mixed rotated lights', {
  path = tempfile(fileext = '.png')
  on.exit(unlink(path))
  # More cells than the 256 blue-noise strata: coarse midpoint samples alone
  # repeatedly miss small alias intervals on this map.
  map = array(.002, c(256, 512, 3))
  map[40:55, 100:125, ] = 1
  png::writePNG(map, path)
  scene = sphere(material = diffuse('white')) |>
    add_object(xz_rect(
      y = -1,
      xwidth = 8,
      zwidth = 8,
      material = diffuse('grey70')
    )) |>
    add_infinite_light(infinite_light(path, intensity = 25, rotation = 57)) |>
    add_infinite_light(disk_light(
      direction = c(-1, 2, 3),
      intensity = .2,
      angular_diameter = 30
    ))
  args = list(
    scene = scene,
    integrator_type = 'nee',
    width = 40,
    height = 40,
    lookfrom = c(0, 1, 6),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(4, 4),
    min_variance = 0,
    max_depth = 4,
    ambient_light = FALSE,
    bloom = FALSE,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE
  )
  set.seed(943)
  reference = do.call(
    render_scene,
    c(args, list(samples = 1024, sample_method = 'sobol'))
  )
  set.seed(943)
  blue = do.call(
    render_scene,
    c(args, list(samples = 128, sample_method = 'sobol_blue'))
  )
  expect_true(all(is.finite(blue)))
  expect_lt(abs(mean(blue[,, 1:3]) - mean(reference[,, 1:3])), .01)
  expect_lt(sqrt(mean((blue[,, 1:3] - reference[,, 1:3])^2)), .04)
})
