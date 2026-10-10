# Compare complete transport, including texture coordinates, mapped scattering,
# transparent continuation, and shadows. Keep assertions per rendered image.
render_wavefront_texture = function(
  scene,
  backend,
  samples = 128L,
  ambient_light = TRUE,
  ...
) {
  set.seed(429)
  render_scene(
    scene,
    integrator_type = backend,
    width = 32,
    height = 32,
    samples = samples,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    min_variance = 0,
    max_depth = 4,
    bloom = FALSE,
    lookfrom = c(0, 0, 5),
    fov = 0,
    ortho_dimensions = c(2, 2),
    ambient_light = ambient_light,
    backgroundlow = 'white',
    backgroundhigh = 'white',
    ...
  )
}

expect_wavefront_texture = function(
  scene,
  tolerance = .03,
  ambient_light = TRUE
) {
  gpu = suppressWarnings(render_wavefront_texture(
    scene,
    'metal',
    ambient_light = ambient_light
  ))
  if (!isTRUE(attr(gpu, 'wavefront')$used)) {
    expect_match(attr(gpu, 'wavefront')$fallback, 'not enabled|device|macOS 11')
    skip('Metal ray tracing is unavailable')
  }
  cpu = render_wavefront_texture(scene, 'nee', ambient_light = ambient_light)
  expect_identical(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_true(all(is.finite(gpu)))
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), tolerance)
  list(gpu = gpu, cpu = cpu)
}

test_that('Metal evaluates composable operations in UV and world coordinates', {
  uv = texture_coordinates(
    'uv',
    scale = c(7, 5, 1),
    offset = c(-.3, .2, 0),
    rotation = 21
  )
  noise = texture_noise(coordinates = uv, octaves = 4, seed = 41)
  images = list(
    texture_checker('coral', 'navy', uv),
    texture_mix('red', 'blue', noise),
    texture_scale(texture_mix('gold', 'green', texture_gradient(0, 1)), .65),
    texture_mix(
      'purple',
      'cyan',
      texture_channel(texture_checker('gold', 'navy', uv), 'luminance')
    ),
    texture_checker('white', 'grey20', texture_coordinates('world', scale = 2)),
    texture_direction_mix('red', 'blue', direction = c(0, 0, 1))
  )
  for (color in images) {
    expect_wavefront_texture(xy_rect(
      xwidth = 4,
      ywidth = 4,
      material = diffuse(color)
    ))
  }
})

test_that('Metal distinguishes object and world mappings through instances', {
  for (space in c('object', 'world')) {
    color = texture_checker(
      'coral',
      'navy',
      texture_coordinates(space, scale = 3, rotation = 17)
    )
    shape = cube(material = diffuse(color), angle = c(13, 27, 8))
    scene = create_instances(
      shape,
      x = c(-.45, .45),
      scale_x = c(.6, -.6),
      scale_y = c(.8, .8),
      scale_z = c(.5, .5)
    )
    expect_wavefront_texture(scene, .035)
  }
})

test_that('Metal uploads graph images with the CPU encoding and addressing', {
  pixels = array(0, c(9, 13, 3))
  pixels[,, 1] = rep(seq(0, 1, length.out = 13), each = 9)
  pixels[,, 2] = seq(0, 1, length.out = 9)
  pixels[,, 3] = .2
  file = tempfile(fileext = '.png')
  on.exit(unlink(file))
  png::writePNG(pixels, file)
  for (wrap in c('repeat', 'clamp')) {
    color = texture_image(
      file,
      coordinates = texture_coordinates(
        'uv',
        scale = c(-2, 3, 1),
        offset = c(.23, -.4, 0)
      ),
      wrap = wrap,
      encoding = 'srgb'
    )
    expect_wavefront_texture(xy_rect(
      xwidth = 4,
      ywidth = 4,
      material = diffuse(color)
    ))
  }
})

test_that('Metal graph roughness reaches both microfacet and OpenPBR', {
  rough = texture_mix(.15, .8, texture_gradient(0, 1))
  for (material in list(
    microfacet('gold', roughness = rough),
    openpbr(base_color = 'coral', specular_roughness = rough)
  )) {
    expect_wavefront_texture(sphere(radius = .8, material = material), .04)
  }
})

test_that('Metal alpha continues through layered cutouts', {
  alpha = array(1, c(8, 8, 4))
  alpha[,, 4] = rep(c(0, .4, 1, .4), length.out = 64)
  file = tempfile(fileext = '.png')
  on.exit(unlink(file))
  png::writePNG(alpha, file)
  foreground = xy_rect(
    z = .3,
    xwidth = 3,
    ywidth = 3,
    material = diffuse(
      'coral',
      alpha_texture = file,
      image_offset = c(.13, -.24)
    )
  )
  scene = foreground |>
    add_object(xy_rect(
      z = -.1,
      xwidth = 3,
      ywidth = 3,
      material = diffuse('steelblue')
    ))
  expect_wavefront_texture(scene, .04)
  expect_wavefront_texture(
    sphere(radius = .8, material = diffuse('gold', alpha_texture = file)),
    .04
  )
})

test_that('Metal alpha coverage modulates shadows under a directional light', {
  file = tempfile(fileext = '.png')
  on.exit(unlink(file))
  pixels = array(1, c(8, 8, 4))
  pixels[,, 4] = rep(c(0, .4, 1, .4), length.out = 64)
  png::writePNG(pixels, file)
  scene = xy_rect(xwidth = 4, ywidth = 4, material = diffuse('grey80')) |>
    add_object(xy_rect(
      x = -.45,
      y = .4,
      z = 1,
      xwidth = 1.2,
      ywidth = 1.2,
      material = diffuse('coral', alpha_texture = file)
    )) |>
    add_light(point_light(c(-3, 2, 4), intensity = 80))
  expect_wavefront_texture(scene, .035, ambient_light = FALSE)
})

test_that('Metal alpha holes do not consume scattering depth', {
  file = tempfile(fileext = '.png')
  on.exit(unlink(file))
  pixels = array(1, c(2, 2, 4))
  pixels[,, 4] = 0
  png::writePNG(pixels, file)
  scene = xy_rect(
    z = -1,
    xwidth = 4,
    ywidth = 4,
    material = diffuse('coral')
  )
  for (z in seq(-.8, .8, length.out = 8)) {
    scene = add_object(
      scene,
      xy_rect(
        z = z,
        xwidth = 4,
        ywidth = 4,
        material = diffuse('blue', alpha_texture = file)
      )
    )
  }
  expect_wavefront_texture(scene)
})

test_that('Metal bump heights preserve filtered diffuse and glossy appearance', {
  file = tempfile(fileext = '.png')
  on.exit(unlink(file))
  xy = seq(0, 2 * pi, length.out = 65)[-65]
  height = .5 + .25 * outer(sin(3 * xy), cos(2 * xy))
  png::writePNG(height, file)
  for (sigma in c(0, 55)) {
    material = diffuse(
      'grey70',
      sigma = sigma,
      bump_texture = file,
      bump_intensity = .12,
      image_repeat = c(-2, 1.5),
      image_offset = c(.1, -.2)
    )
    expect_wavefront_texture(
      xy_rect(xwidth = 3, ywidth = 3, material = material),
      .035
    )
  }
  expect_wavefront_texture(
    sphere(
      radius = .8,
      material = glossy('steelblue', bump_texture = file, bump_intensity = .015)
    ),
    .04
  )
})

test_that('Metal mesh bump and alpha maps survive transformed placements', {
  skip_if_not_installed('rayvertex')
  height = tempfile(fileext = '.png')
  alpha = tempfile(fileext = '.png')
  on.exit(unlink(c(height, alpha)))
  png::writePNG(
    outer(
      seq(0, 1, length.out = 32),
      seq(0, 1, length.out = 32),
      function(x, y) .5 + .2 * sin(16 * x) * cos(20 * y)
    ),
    height
  )
  pixels = array(1, c(16, 16, 4))
  pixels[,, 4] = rep(c(0, 1), each = 4, length.out = 256)
  png::writePNG(pixels, alpha)
  mesh = rayvertex::sphere_mesh(radius = .75)
  for (material in list(
    diffuse('coral', bump_texture = height, bump_intensity = .02),
    diffuse('gold', alpha_texture = alpha)
  )) {
    shape = raymesh_model(
      mesh,
      angle = c(23, 17, 9),
      material = material,
      override_material = TRUE,
      calculate_consistent_normals = FALSE
    )
    expect_wavefront_texture(shape, .04)
  }
})

test_that('Metal bump heights scale with nonuniform and mirrored instances', {
  file = tempfile(fileext = '.png')
  on.exit(unlink(file))
  height = outer(
    seq(0, 1, length.out = 64),
    seq(0, 1, length.out = 64),
    function(x, y) .5 + .25 * sin(24 * x) * cos(21 * y)
  )
  png::writePNG(height, file)
  shape = xy_rect(
    xwidth = 2,
    ywidth = 2,
    material = diffuse('grey70', bump_texture = file, bump_intensity = .15)
  )
  scene = create_instances(
    shape,
    x = c(-.55, .55),
    scale_x = c(.45, -.45),
    scale_y = c(.8, .8),
    scale_z = c(.125, .125)
  ) |>
    add_infinite_light(disk_light(
      direction = c(-1, 1, 2),
      angular_diameter = 10,
      intensity = 100
    ))
  expect_wavefront_texture(scene, .02, ambient_light = FALSE)
})
