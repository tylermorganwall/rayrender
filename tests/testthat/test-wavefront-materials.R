# Compare complete transport, not just acceptance by the scene exporter. A broad
# white environment and fixed framing give stable RGB estimates at modest cost.
render_wavefront_material = function(
  scene,
  backend = 'metal',
  samples = 128,
  lookfrom = c(0, 0, 4),
  depth = 8
) {
  set.seed(137)
  image = render_scene(
    scene,
    integrator_type = backend,
    width = 16,
    height = 16,
    samples = samples,
    lookfrom = lookfrom,
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.5, .5),
    max_depth = depth,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white',
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    min_variance = 0,
    bloom = FALSE
  )
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    reason = attr(image, 'wavefront')$fallback
    if (grepl('not enabled|device|macOS 11', reason)) {
      skip('Metal is unavailable')
    }
    stop('Unexpected GPU fallback: ', reason)
  }
  image
}

wavefront_rgb_mean = function(image) apply(image[,, 1:3, drop = FALSE], 3, mean)

test_that('Metal transports the reflective, transmissive, sheet and fiber materials', {
  materials = list(
    metal = metal(color = c(.8, .6, .3)),
    fuzzy_metal = metal(color = c(.7, .8, .9), fuzz = .15),
    ggx = microfacet(color = c(.8, .6, .3), roughness = .2),
    beckmann = microfacet(
      color = c(.4, .6, .8),
      roughness = c(.15, .4),
      microfacet = 'beckmann'
    ),
    transmission = microfacet(
      color = c(.8, .9, .95),
      roughness = .2,
      transmission = TRUE
    ),
    matched = microfacet(color = c(.3, .5, .8), transmission = TRUE, eta = 1),
    glossy = glossy(color = c(.4, .6, .8), gloss = .7),
    sheet = translucent(
      reflectance = c(.2, .3, .1),
      transmittance = c(.3, .4, .6)
    ),
    fiber = hair(color = c(.6, .3, .1), beta_m = .4, beta_n = .5)
  )
  for (name in names(materials)) {
    scene = xy_rect(xwidth = 100, ywidth = 100, material = materials[[name]])
    gpu = suppressWarnings(render_wavefront_material(scene))
    cpu = render_wavefront_material(scene, 'nee')
    expect_true(all(is.finite(gpu)), info = name)
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0, info = name)
    expect_equal(
      wavefront_rgb_mean(gpu),
      wavefront_rgb_mean(cpu),
      tolerance = .04,
      info = name
    )
  }
})

test_that('Metal OpenPBR keeps its layered surface rather than becoming plain glass', {
  materials = list(
    rough_metal = openpbr(
      base_color = c(.8, .5, .2),
      base_metalness = 1,
      specular_roughness = .3
    ),
    velvet = openpbr(
      base_color = c(.15, .3, .5),
      fuzz_weight = .9,
      fuzz_color = c(.6, .8, 1)
    ),
    coated = openpbr(
      base_color = c(.6, .3, .1),
      coat_weight = 1,
      coat_roughness = .15
    ),
    thin_film = openpbr(
      base_metalness = 1,
      thin_film_weight = 1,
      thin_film_thickness = .45
    ),
    thin_sheet = openpbr(
      base_color = c(.3, .6, .9),
      geometry_thin_walled = TRUE,
      transmission_weight = .5
    )
  )
  for (name in names(materials)) {
    scene = xy_rect(xwidth = 100, ywidth = 100, material = materials[[name]])
    gpu = suppressWarnings(render_wavefront_material(scene))
    cpu = render_wavefront_material(scene, 'nee')
    expect_true(all(is.finite(gpu)), info = name)
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0, info = name)
    expect_equal(
      wavefront_rgb_mean(gpu),
      wavefront_rgb_mean(cpu),
      tolerance = .04,
      info = name
    )
  }
})

test_that('Metal normalized diffusion handles priority and nonlocal exits', {
  material = subsurface_diffusion(
    color = c(.8, .7, .5),
    radius = .12,
    priority = 2
  )
  scenes = list(
    bare = sphere(material = material),
    priority = sphere(material = material) |>
      add_object(sphere(
        x = .35,
        radius = .9,
        material = dielectric(priority = 1)
      )),
    hidden = sphere(material = material) |>
      add_object(sphere(
        radius = 1.1,
        material = dielectric(refraction = 1, priority = 0)
      ))
  )
  for (name in names(scenes)) {
    gpu = suppressWarnings(render_wavefront_material(
      scenes[[name]],
      samples = 256
    ))
    cpu = render_wavefront_material(scenes[[name]], 'nee', samples = 256)
    expect_true(all(is.finite(gpu)), info = name)
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0, info = name)
    expect_equal(
      wavefront_rgb_mean(gpu),
      wavefront_rgb_mean(cpu),
      tolerance = .06,
      info = name
    )
    if (name == 'hidden') {
      expect_equal(attr(gpu, 'wavefront')$scattering_events, 0)
    }
  }
})

test_that('Metal OpenPBR supports refracting and scattering solid interiors', {
  materials = list(
    glass = openpbr(
      transmission_weight = 1,
      transmission_color = c(.8, .9, .95),
      transmission_depth = 2,
      specular_roughness = .12
    ),
    sss = openpbr(
      subsurface_weight = 1,
      subsurface_radius = .4,
      subsurface_color = c(.85, .7, .5),
      specular_roughness = .1
    )
  )
  for (name in names(materials)) {
    scene = sphere(material = materials[[name]])
    gpu = suppressWarnings(render_wavefront_material(scene, samples = 256))
    cpu = render_wavefront_material(scene, 'nee', samples = 256)
    expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0, info = name)
    expect_equal(
      wavefront_rgb_mean(gpu),
      wavefront_rgb_mean(cpu),
      tolerance = .045,
      info = name
    )
  }
})

test_that('Metal diffusion resolves coincident contacts and a camera inside', {
  scene = cube(
    width = 2,
    material = subsurface_diffusion(
      color = c(.8, .7, .6),
      radius = .1,
      priority = 2
    )
  ) |>
    add_object(cube(
      z = .6,
      xwidth = 2,
      ywidth = 2,
      zwidth = .8,
      material = dielectric(refraction = 1.5, priority = 1)
    ))
  gpu = suppressWarnings(render_wavefront_material(scene, samples = 256))
  cpu = render_wavefront_material(scene, 'nee', samples = 256)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_equal(
    wavefront_rgb_mean(gpu),
    wavefront_rgb_mean(cpu),
    tolerance = .06
  )
  inside = cube(width = 4, material = subsurface_diffusion())
  expect_warning(
    image <- render_wavefront_material(
      inside,
      samples = 1,
      lookfrom = c(0, 0, 1)
    ),
    'terminated|discarded'
  )
  expect_equal(attr(image, 'wavefront')$discarded_paths, 256)
  expect_equal(wavefront_rgb_mean(image), c(0, 0, 0))
})

test_that('Metal hair uses the cylindrical and ribbon curve impact parameters', {
  scene = bezier_curve(
    p1 = c(-.12, -1, 0),
    p2 = c(-.12, -.33, 0),
    p3 = c(-.12, .33, 0),
    p4 = c(-.12, 1, 0),
    width = .12,
    material = hair(color = c(.7, .35, .15), beta_m = .4, beta_n = .5)
  ) |>
    add_object(bezier_curve(
      p1 = c(.12, -1, 0),
      p2 = c(.12, -.33, 0),
      p3 = c(.12, .33, 0),
      p4 = c(.12, 1, 0),
      width = .12,
      material = hair(color = c(.1, .3, .6), beta_m = .35, beta_n = .5)
    ))
  gpu = suppressWarnings(render_wavefront_material(scene, samples = 256))
  cpu = render_wavefront_material(scene, 'nee', samples = 256)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_equal(
    wavefront_rgb_mean(gpu),
    wavefront_rgb_mean(cpu),
    tolerance = .04
  )
  ribbon = bezier_curve(
    p1 = c(0, -1, 0),
    p2 = c(0, -.33, 0),
    p3 = c(0, .33, 0),
    p4 = c(0, 1, 0),
    width = .2,
    type = 'ribbon',
    material = hair(color = c(.7, .35, .15), beta_m = .4, beta_n = .5)
  )
  gpu = suppressWarnings(render_wavefront_material(ribbon, samples = 256))
  cpu = render_wavefront_material(ribbon, 'nee', samples = 256)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_equal(
    wavefront_rgb_mean(gpu),
    wavefront_rgb_mean(cpu),
    tolerance = .05
  )
})

test_that('Matching OpenPBR IOR does not erase its layered shadow response', {
  scene = xz_rect(material = diffuse(c(.6, .6, .6))) |>
    add_object(cube(
      y = .6,
      width = .7,
      material = openpbr(
        specular_ior = 1,
        transmission_weight = .7,
        base_color = c(.3, .5, .7)
      )
    )) |>
    add_light(point_light(position = c(0, 3, 0), intensity = 5))
  images = lapply(c('metal', 'nee'), function(backend) {
    set.seed(173)
    image = suppressWarnings(render_scene(
      scene,
      integrator_type = backend,
      width = 24,
      height = 24,
      lookfrom = c(0, 1, 4),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(1, 1),
      samples = 256,
      ambient_light = FALSE,
      preview = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      min_variance = 0,
      max_depth = 8
    ))
    if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
      if (
        grepl('not enabled|device|macOS 11', attr(image, 'wavefront')$fallback)
      ) {
        skip('Metal is unavailable')
      }
      stop(attr(image, 'wavefront')$fallback)
    }
    image
  })
  expect_equal(
    wavefront_rgb_mean(images[[1]]),
    wavefront_rgb_mean(images[[2]]),
    tolerance = .06
  )
})

test_that('Metal preserves image roughness lookup and UV offsets', {
  path = tempfile(fileext = '.png')
  withr::defer(unlink(path))
  roughness = array(rep(c(.05, .3, .65, .9), 3), c(2, 2, 3))
  png::writePNG(roughness, path)
  materials = list(
    microfacet(
      roughness_texture = path,
      roughness_range = c(.05, .9),
      image_offset = c(.3, -.2)
    ),
    glossy(
      color = 'steelblue',
      gloss = .6,
      roughness_texture = path,
      image_offset = c(.3, -.2)
    ),
    openpbr(
      base_metalness = 1,
      base_color = c(.8, .5, .2),
      roughness_texture = path,
      image_repeat = c(2, 3),
      image_offset = c(.3, -.2)
    )
  )
  for (material in materials) {
    scene = xy_rect(xwidth = .5, ywidth = .5, material = material)
    gpu = suppressWarnings(render_wavefront_material(scene))
    cpu = render_wavefront_material(scene, 'nee')
    expect_equal(
      wavefront_rgb_mean(gpu),
      wavefront_rgb_mean(cpu),
      tolerance = .04
    )
  }
})

test_that('Metal supports reflecting OpenPBR emitters', {
  scene = xy_rect(xwidth = 10, ywidth = 10, material = diffuse('grey60')) |>
    add_object(xz_rect(
      y = 2,
      z = 2,
      angle = c(180, 0, 0),
      xwidth = 2,
      zwidth = 2,
      material = openpbr(
        base_color = c(.4, .6, .8),
        emission_luminance = 5,
        emission_color = c(1, .6, .2),
        coat_weight = .4
      )
    ))
  gpu = suppressWarnings(render_wavefront_material(scene, samples = 256))
  cpu = render_wavefront_material(scene, 'nee', samples = 256)
  expect_equal(attr(gpu, 'wavefront')$discarded_paths, 0)
  expect_equal(
    wavefront_rgb_mean(gpu),
    wavefront_rgb_mean(cpu),
    tolerance = .05
  )
})
