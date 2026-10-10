# Render comparisons exercise the public scene compiler and complete GPU queue
# pipeline, including reconstruction and film accumulation. No exact-noise
# comparison is expected between two independent sampling implementations.
test_that("Metal fallback preserves the existing CPU result", {
  scene = cube(material = dielectric())
  render = function(backend) {
    set.seed(729)
    render_scene(
      scene,
      integrator_type = backend,
      width = 8,
      height = 8,
      samples = 2,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      lookfrom = c(0, 0, 4),
      max_depth = 3
    )
  }
  cpu = render('nee')
  expect_warning(gpu <- render('metal'), 'MediumBoundary.*CPU NEE')
  expect_false(attr(gpu, 'wavefront')$used)
  expect_identical(as.numeric(gpu), as.numeric(cpu))
})

test_that("Metal wavefront traces diffuse geometry and respects depth", {
  scene = xy_rect(xwidth = 100, ywidth = 100, material = diffuse(c(.3, .5, .7)))
  render = function(backend, depth = 2L, material = scene) {
    set.seed(136)
    render_scene(
      material,
      integrator_type = backend,
      width = 24,
      height = 24,
      samples = 32,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      lookfrom = c(0, 0, 4),
      fov = 0,
      ortho_dimensions = c(1, 1),
      ambient_light = TRUE,
      backgroundlow = 'white',
      backgroundhigh = 'white',
      max_depth = depth,
      bloom = FALSE
    )
  }
  message = NULL
  gpu = withCallingHandlers(render('metal'), warning = function(w) {
    message <<- conditionMessage(w)
    invokeRestart('muffleWarning')
  })
  report = attr(gpu, 'wavefront')
  if (!isTRUE(report$used)) {
    expect_match(message, 'Metal.*(not enabled|device|macOS 11)')
    skip('Metal ray tracing is unavailable on this build/device')
  }
  expect_true(all(is.finite(gpu)))
  expect_identical(report$completed_samples, 32)
  expect_identical(report$triangles, 2)
  expect_equal(
    as.numeric(render('metal', depth = 1L)),
    as.numeric(gpu),
    tolerance = 1e-6
  )
  cpu = render('nee')
  means = function(x) apply(x[,, 1:3, drop = FALSE], 3, mean)
  expect_equal(means(gpu), means(cpu), tolerance = .015)
  expect_identical(as.numeric(render('metal')), as.numeric(gpu))
  rough = xy_rect(
    xwidth = 100,
    ywidth = 100,
    material = diffuse(c(.3, .5, .7), sigma = 55)
  )
  expect_equal(
    means(render('metal', material = rough)),
    means(render('nee', material = rough)),
    tolerance = .02
  )
})

test_that("Metal direct lights and environment rotation agree with CPU", {
  scene = cube(material = diffuse('coral')) |>
    add_object(xz_rect(
      y = -.51,
      xwidth = 5,
      zwidth = 5,
      material = diffuse('grey70')
    )) |>
    add_infinite_light(disk_light(
      intensity = 8,
      angular_diameter = 35,
      direction = c(-1, 2, 1)
    ))
  render = function(backend, rotation = 0) {
    set.seed(791)
    render_scene(
      scene,
      integrator_type = backend,
      width = 32,
      height = 32,
      samples = 64,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      lookfrom = c(3, 2, 5),
      fov = 40,
      ambient_light = FALSE,
      rotate_env = rotation,
      max_depth = 4,
      bloom = FALSE
    )
  }
  gpu = suppressWarnings(render('metal'))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal unavailable; checked in the transport test')
  }
  cpu = render('nee')
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .025)
  rotated = render('metal', 100)
  expect_gt(mean(abs(gpu[,, 1:3] - rotated[,, 1:3])), .025)
  expect_lt(mean(abs(rotated[,, 1:3] - render('nee', 100)[,, 1:3])), .025)
})

test_that("Metal exports owned mesh triangles and transformed instances", {
  skip_if_not_installed('rayvertex')
  mesh = rayvertex::sphere_mesh(low_poly = TRUE)
  count = nrow(mesh$shapes[[1]]$indices)
  # OBJ-like partial normal indices must make those complete faces flat.
  mesh$shapes[[1]]$norm_indices[seq(1L, count, by = 2L), 1] = -1L
  prototype = raymesh_model(
    mesh,
    material = diffuse('steelblue'),
    calculate_consistent_normals = FALSE,
    override_material = TRUE
  )
  scene = create_instances(
    prototype,
    x = c(-1.2, 1.2),
    scale_x = c(-1, 1),
    scale_y = c(1, .7)
  ) |>
    add_object(xz_rect(
      y = -1,
      xwidth = 9,
      zwidth = 9,
      material = diffuse('grey70')
    )) |>
    add_infinite_light(disk_light(
      intensity = 8,
      angular_diameter = 35,
      direction = c(-1, 2, 1)
    ))
  render = function(backend) {
    set.seed(25)
    render_scene(
      scene,
      integrator_type = backend,
      width = 48,
      height = 32,
      samples = 128,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      lookfrom = c(0, 2, 7),
      lookat = c(0, 0, 0),
      fov = 42,
      ambient_light = FALSE,
      max_depth = 4,
      bloom = FALSE
    )
  }
  gpu = suppressWarnings(render('metal'))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal unavailable; checked in the transport test')
  }
  expect_equal(attr(gpu, 'wavefront')$triangles, 2 * count + 2)
  cpu = render('nee')
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .025)
})

test_that("Metal samples area and spot lights with textured surfaces", {
  path = tempfile(fileext = '.png')
  on.exit(unlink(path))
  image = array(0, c(2, 2, 3))
  image[,, 1] = c(.2, .8, .4, .9)
  image[,, 2] = c(.7, .1, .3, .5)
  image[,, 3] = .2
  png::writePNG(image, path)
  scene = cube(
    material = diffuse(image_texture = path, image_offset = c(.2, -.3))
  ) |>
    add_object(xz_rect(
      y = -.51,
      xwidth = 6,
      zwidth = 6,
      material = diffuse('grey70')
    )) |>
    add_object(xz_rect(
      y = 3,
      xwidth = 2,
      zwidth = 2,
      flipped = TRUE,
      material = light(intensity = 4)
    )) |>
    add_light(spot_light(
      position = c(-2, 3, 2),
      direction = c(2, -3, -2),
      intensity = 15
    ))
  render = function(backend) {
    set.seed(309)
    render_scene(
      scene,
      integrator_type = backend,
      width = 32,
      height = 32,
      samples = 128,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      min_variance = 0,
      lookfrom = c(3, 2, 5),
      fov = 40,
      ambient_light = FALSE,
      max_depth = 4,
      bloom = FALSE
    )
  }
  gpu = suppressWarnings(render('metal'))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    skip('Metal unavailable; checked in the transport test')
  }
  cpu = render('nee')
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .025)
})
