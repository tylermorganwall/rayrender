# Render both traversal backends, without denoising masking transport errors.
render_wavefront_sharing = function(
  scene,
  backend = 'metal',
  samples = 128,
  ...
) {
  set.seed(754)
  image = suppressWarnings(render_scene(
    scene,
    integrator_type = backend,
    width = 40,
    height = 40,
    samples = samples,
    lookfrom = c(0, 1, 6),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(4, 4),
    min_variance = 0,
    max_depth = 16,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ...
  ))
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    expect_match(
      attr(image, 'wavefront')$fallback,
      '^(Metal was not enabled|this Metal device does not support|Metal wavefront requires macOS 11)'
    )
    skip('Metal ray tracing is unavailable')
  }
  image
}

test_that('nested mesh instances share geometry while preserving texture and bump frames', {
  skip_if_not_installed('rayvertex')
  mesh = rayvertex::sphere_mesh()
  material = diffuse(
    texture_checker('gold', 'steelblue', texture_coordinates('uv', scale = 5)),
    bump_texture = outer(
      seq(0, 8 * pi, length.out = 64),
      seq(0, 8 * pi, length.out = 64),
      function(x, y) .5 + .25 * sin(x) * cos(y)
    ),
    bump_intensity = .08
  )
  prototype = raymesh_model(mesh, material = material, override_material = TRUE)
  scene = create_instances(
    prototype,
    x = c(-.9, .9),
    scale_x = c(-.65, .65),
    scale_y = c(.6, .8),
    scale_z = .5
  ) |>
    create_instances(y = c(-.8, .8), angle_y = c(-15, 20)) |>
    add_infinite_light(disk_light(
      intensity = 5,
      angular_diameter = 40,
      direction = c(-1, 2, 3)
    ))
  gpu = render_wavefront_sharing(scene, ambient_light = FALSE)
  cpu = render_wavefront_sharing(scene, 'nee', ambient_light = FALSE)
  report = attr(gpu, 'wavefront')
  expect_equal(report$unique_meshes, 1)
  expect_equal(report$mesh_instances, 4)
  expect_equal(report$triangles, nrow(mesh$shapes[[1]]$indices))
  expect_equal(report$discarded_paths, 0)
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .035)
})

test_that('emissive mesh placements retain distinct positions and world-area PDFs', {
  skip_if_not_installed('rayvertex')
  prototype = raymesh_model(
    rayvertex::sphere_mesh(low_poly = TRUE),
    material = light(intensity = 6),
    override_material = TRUE
  )
  scene = create_instances(
    prototype,
    x = c(-1.4, 1.4),
    y = c(1, 1.5),
    scale_x = c(.3, -.4),
    scale_y = c(.5, .35),
    scale_z = .3
  ) |>
    add_object(sphere(radius = .7, material = diffuse('grey70'))) |>
    add_object(xz_rect(
      y = -.7,
      xwidth = 8,
      zwidth = 8,
      material = diffuse('white')
    ))
  gpu = render_wavefront_sharing(scene, samples = 512, ambient_light = FALSE)
  cpu = render_wavefront_sharing(
    scene,
    'nee',
    samples = 512,
    ambient_light = FALSE
  )
  expect_equal(attr(gpu, 'wavefront')$mesh_instances, 2)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .04)
})

test_that('instanced SSS and glass keep separate priority boundary identities', {
  skip_if_not_installed('rayvertex')
  mesh = rayvertex::sphere_mesh()
  for (milk in list(
    subsurface(sigma_a = .05, sigma_s = 3, priority = 2),
    subsurface_diffusion(color = 'white', radius = .2, priority = 2)
  )) {
    prototype = raymesh_model(
      mesh,
      material = milk,
      override_material = TRUE
    ) |>
      add_object(sphere(
        radius = 1.04,
        material = dielectric(refraction = 1.4, priority = 1)
      ))
    scene = create_instances(
      prototype,
      x = c(-.9, .9),
      scale_x = c(.65, -.65),
      scale_y = .8,
      scale_z = .65
    )
    gpu = render_wavefront_sharing(
      scene,
      samples = 256,
      ambient_light = TRUE,
      backgroundlow = 'white',
      backgroundhigh = 'white'
    )
    cpu = render_wavefront_sharing(
      scene,
      'nee',
      samples = 256,
      ambient_light = TRUE,
      backgroundlow = 'white',
      backgroundhigh = 'white'
    )
    expect_equal(attr(gpu, 'wavefront')$mesh_instances, 2)
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
    expect_true(all(is.finite(gpu)))
    expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .045)
  }
})

test_that('image-environment alias lighting agrees with CPU under rotation', {
  # Non-power-of-two dimensions also exercise zero-mass alias padding.
  map = array(.002, c(127, 251, 3))
  map[25:35, 50:65, ] = 1
  map[48:64, 180:200, 3] = .25
  file = tempfile(fileext = '.png')
  on.exit(unlink(file))
  png::writePNG(map, file)
  scene = sphere(material = diffuse('white')) |>
    add_object(xz_rect(
      y = -1,
      xwidth = 8,
      zwidth = 8,
      material = diffuse('grey70')
    ))
  # At 128 spp the default is the eight-bit blue-noise sampler; alias lookup
  # must restore within-stratum precision instead of visiting only 256 columns.
  for (angle in c(0, 73)) {
    gpu = render_wavefront_sharing(
      scene,
      samples = 128,
      ambient_light = FALSE,
      environment_light = file,
      rotate_env = angle,
      intensity_env = 30
    )
    cpu = render_wavefront_sharing(
      scene,
      'nee',
      samples = 512,
      ambient_light = FALSE,
      environment_light = file,
      rotate_env = angle,
      intensity_env = 30
    )
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
    expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .025)
  }
})

test_that('shared mesh alpha masks survive mirrored placements and mixed traversal', {
  skip_if_not_installed('rayvertex')
  opacity = outer(seq_len(32), seq_len(32), function(x, y) {
    as.numeric((x %/% 4 + y %/% 4) %% 2)
  })
  prototype = raymesh_model(
    rayvertex::sphere_mesh(),
    material = diffuse('coral', alpha_texture = opacity),
    override_material = TRUE,
    calculate_consistent_normals = FALSE
  )
  scene = create_instances(
    prototype,
    x = c(-.85, .85),
    scale_x = c(-.7, .7),
    scale_y = .8,
    scale_z = .6
  ) |>
    add_object(sphere(y = -.7, z = -1, radius = .5, material = diffuse('navy')))
  gpu = render_wavefront_sharing(
    scene,
    samples = 256,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  # A separately transformed copy of each mesh is an independent reference
  # for sharing: UVs, alpha continuation and shadow tests must survive the TLAS.
  expanded = raymesh_model(
    rayvertex::sphere_mesh(),
    x = -.85,
    scale = c(-.7, .8, .6),
    material = diffuse('coral', alpha_texture = opacity),
    override_material = TRUE,
    calculate_consistent_normals = FALSE
  ) |>
    add_object(raymesh_model(
      rayvertex::sphere_mesh(),
      x = .85,
      scale = c(.7, .8, .6),
      material = diffuse('coral', alpha_texture = opacity),
      override_material = TRUE,
      calculate_consistent_normals = FALSE
    )) |>
    add_object(sphere(y = -.7, z = -1, radius = .5, material = diffuse('navy')))
  reference = render_wavefront_sharing(
    expanded,
    samples = 256,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  expect_equal(attr(gpu, 'wavefront')$mesh_instances, 2)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_lt(mean(abs(gpu[,, 1:3] - reference[,, 1:3])), .035)
  cpu = render_wavefront_sharing(
    scene,
    'nee',
    samples = 256,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white'
  )
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .015)
})
