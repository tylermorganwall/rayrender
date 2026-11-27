# Open interfaces must remain real refracting surfaces, including on shared meshes.
render_open_dielectric = function(scene, backend = "metal", samples = 64, ...) {
  set.seed(847)
  image = render_scene(
    scene,
    integrator_type = backend,
    width = 20,
    height = 20,
    samples = samples,
    lookfrom = c(0, 0, 4),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.6, .6),
    max_depth = 12,
    ambient_light = FALSE,
    backgroundlow = "black",
    backgroundhigh = "black",
    min_variance = 0,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ...
  )
  if (backend == "metal" && !isTRUE(attr(image, "wavefront")$used)) {
    expect_match(
      attr(image, "wavefront")$fallback,
      "not enabled|device|macOS 11"
    )
    skip("Metal ray tracing is unavailable")
  }
  if (backend == "metal") {
    expect_equal(attr(image, "wavefront")$discarded_paths, 0)
  }
  image
}

test_that("open sheets attenuate the entire segment to opaque surfaces", {
  sigma = c(.1, .3, .7)
  scene = xy_rect(
    z = 1,
    xwidth = 100,
    ywidth = 100,
    material = dielectric(refraction = 1, attenuation = sigma)
  ) |>
    add_object(xy_rect(
      z = -1,
      xwidth = 100,
      ywidth = 100,
      material = light(intensity = 1)
    ))
  gpu = render_open_dielectric(scene, samples = 4)
  cpu = render_open_dielectric(scene, "nee", samples = 4)
  expect_equal(apply(gpu[,, 1:3], 3, mean), exp(-2 * sigma), tolerance = 1e-4)
  expect_equal(as.numeric(gpu), as.numeric(cpu), tolerance = 1e-4)
})

test_that("open water retains optical state after opaque ground reflection", {
  water = dielectric(refraction = 1.33, attenuation = c(.08, .03, .01))
  floor = xy_rect(
    z = -1,
    xwidth = 100,
    ywidth = 100,
    material = diffuse("white")
  )
  open = xy_rect(xwidth = 100, ywidth = 100, material = water) |>
    add_object(floor)
  closed = cube(
    z = -50,
    xwidth = 1000,
    ywidth = 1000,
    zwidth = 100,
    material = water
  ) |>
    add_object(floor)
  # The opaque floor closes every camera-visible water path; the box's other
  # faces lie far beyond it. A white sky makes return paths through the water visible.
  skyfile = tempfile(fileext = ".png")
  on.exit(unlink(skyfile))
  png::writePNG(array(1, c(2, 4, 3)), skyfile)
  sky = infinite_light(skyfile)
  open = add_infinite_light(open, sky)
  closed = add_infinite_light(closed, sky)
  gpu = render_open_dielectric(open, samples = 256)
  box = render_open_dielectric(closed, samples = 256)
  cpu = render_open_dielectric(open, "nee", samples = 256)
  means = function(x) apply(x[,, 1:3], 3, mean)
  expect_equal(means(gpu), means(box), tolerance = .035)
  expect_equal(means(gpu), means(cpu), tolerance = .04)
  expect_gt(mean(gpu), .2)
})

test_that("open sheets restore absorption after priority-selected closed glass exits", {
  scene = xy_rect(
    z = 1,
    xwidth = 100,
    ywidth = 100,
    material = dielectric(
      refraction = 1,
      attenuation = rep(.2, 3),
      priority = 2
    )
  ) |>
    add_object(cube(
      xwidth = 100,
      ywidth = 100,
      zwidth = 1,
      material = dielectric(
        refraction = 1,
        attenuation = rep(.6, 3),
        priority = 1
      )
    )) |>
    add_object(xy_rect(
      z = -1,
      xwidth = 100,
      ywidth = 100,
      material = light(intensity = 1)
    ))
  gpu = render_open_dielectric(scene, samples = 4)
  expect_equal(mean(gpu[,, 1:3]), exp(-.8), tolerance = 1e-4)
})

test_that("diffusion containment does not erase a higher priority open interface", {
  scene = xy_rect(
    z = 1,
    xwidth = 100,
    ywidth = 100,
    material = dielectric(
      refraction = 1,
      attenuation = rep(.2, 3),
      priority = 0
    )
  ) |>
    add_object(sphere(
      radius = .5,
      material = subsurface_diffusion(
        color = "white",
        radius = .1,
        priority = 2
      )
    )) |>
    add_object(xy_rect(
      z = -1,
      xwidth = 100,
      ywidth = 100,
      material = light(intensity = 1)
    ))
  gpu = render_open_dielectric(scene, samples = 4)
  expect_equal(mean(gpu[,, 1:3]), exp(-.4), tolerance = 1e-4)
  expect_equal(attr(gpu, "wavefront")$scattering_events, 0)
})

test_that("shared open mesh placements retain repeated material crossings", {
  skip_if_not_installed("rayvertex")
  mesh = rayvertex::xy_rect_mesh(scale = c(100, 100, 1))
  pane = raymesh_model(
    mesh,
    material = dielectric(refraction = 1, attenuation = rep(.2, 3))
  )
  # Enter the same material twice, then exit only one placement. Absorption
  # continues behind that exit until the emissive target.
  scene = create_instances(pane, z = c(1, .5)) |>
    add_object(create_instances(pane, z = 0, angle_y = 180)) |>
    add_object(xy_rect(
      z = -1,
      xwidth = 100,
      ywidth = 100,
      material = light(intensity = 1)
    ))
  gpu = render_open_dielectric(scene, samples = 4)
  expect_equal(mean(gpu[,, 1:3]), exp(-.4), tolerance = 1e-4)
  expect_equal(attr(gpu, "wavefront")$mesh_instances, 3)
})

test_that("tessellated open disks transmit and apply attenuation", {
  scene = disk(
    z = 1,
    radius = 5,
    angle = c(90, 0, 0),
    material = dielectric(refraction = 1, attenuation = rep(.3, 3))
  ) |>
    add_object(xy_rect(
      z = -1,
      xwidth = 10,
      ywidth = 10,
      material = light(intensity = 1)
    ))
  gpu = render_open_dielectric(scene, samples = 4)
  expect_equal(mean(gpu[,, 1:3]), exp(-.6), tolerance = 1e-4)
  expect_gt(attr(gpu, "wavefront")$triangles, 2)
})

test_that("an unpaired back face still produces total internal reflection", {
  # The camera is on the interior side, 60 degrees from the exit normal.
  # All rays must reflect to the black background rather than reach the emitter.
  scene = xy_rect(
    xwidth = 100,
    ywidth = 100,
    angle = c(0, 60, 0),
    flipped = TRUE,
    material = dielectric(refraction = 1.5)
  ) |>
    add_object(xy_rect(
      z = -10,
      xwidth = 1,
      ywidth = 100,
      material = light(intensity = 1)
    ))
  gpu = render_open_dielectric(scene, samples = 4)
  expect_equal(as.numeric(gpu[,, 1:3]), rep(0, 20 * 20 * 3), tolerance = 1e-6)
})

test_that("a first back-face hit attenuates the known camera-to-surface segment", {
  sigma = c(.1, .3, .7)
  scene = xy_rect(
    z = 1,
    xwidth = 10,
    ywidth = 10,
    flipped = TRUE,
    material = dielectric(refraction = 1, attenuation = sigma)
  ) |>
    add_object(xy_rect(
      z = -1,
      xwidth = 10,
      ywidth = 10,
      material = light(intensity = 1)
    ))
  gpu = render_open_dielectric(scene, samples = 4)
  # Camera z=4 to interface z=1 is three units; the next two units are air.
  expect_equal(apply(gpu[,, 1:3], 3, mean), exp(-3 * sigma), tolerance = 1e-4)
})

test_that("diffusion exits recover the surrounding open dielectric", {
  skyfile = tempfile(fileext = ".png")
  on.exit(unlink(skyfile))
  png::writePNG(array(1, c(2, 4, 3)), skyfile)
  water = dielectric(
    refraction = 1.33,
    attenuation = c(.6, .2, .05),
    priority = 2
  )
  body = sphere(
    z = -1,
    radius = .5,
    material = subsurface_diffusion(
      color = "white",
      radius = .1,
      refraction = 1.33,
      priority = 0
    )
  )
  open = xy_rect(xwidth = 1000, ywidth = 1000, material = water) |>
    add_object(body) |>
    add_infinite_light(infinite_light(skyfile))
  closed = cube(
    z = -500,
    xwidth = 1000,
    ywidth = 1000,
    zwidth = 1000,
    material = water
  ) |>
    add_object(body) |>
    add_infinite_light(infinite_light(skyfile))
  gpu = render_open_dielectric(open, samples = 128)
  reference = render_open_dielectric(closed, samples = 128)
  expect_gt(attr(gpu, "wavefront")$scattering_events, 0)
  expect_equal(
    apply(gpu[,, 1:3], 3, mean),
    apply(reference[,, 1:3], 3, mean),
    tolerance = .04
  )
})
