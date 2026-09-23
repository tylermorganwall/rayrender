test_that("media and clouds validate and retain their haze preference", {
  expect_true(homogeneous_medium()$haze)
  expect_false(homogeneous_medium(haze = FALSE)$haze)
  expect_false(grid_medium(array(1, c(2, 2, 2)), haze = FALSE)$haze)
  expect_equal(
    homogeneous_medium(haze_density_threshold = .2)$haze_density_threshold,
    .2
  )
  for (invalid in list(NA, Inf, 0, -1, TRUE, "0.2", c(.1, .2))) {
    expect_error(
      homogeneous_medium(haze_density_threshold = invalid),
      "haze_density_threshold"
    )
    expect_error(
      cloud(haze_density_threshold = invalid),
      "haze_density_threshold"
    )
  }
  for (invalid in list(NA, NULL, 0, "no", c(TRUE, FALSE))) {
    expect_error(homogeneous_medium(haze = invalid), "haze")
    expect_error(grid_medium(array(1, c(2, 2, 2)), haze = invalid), "haze")
  }
  medium = homogeneous_medium(haze = FALSE)
  expect_identical(unserialize(serialize(medium, NULL)), medium)
  skip_if_not_installed("ambient")
  puff = cloud(resolution = 24)
  expect_false(puff$shape_info[[1]]$medium$haze)
  expect_equal(puff$shape_info[[1]]$medium$haze_density_threshold, .05)
  expect_gt(puff$shape_info[[1]]$medium$sigma_s[1], 0)
  expect_error(cloud(haze = NA), "haze")
  puff = cloud(resolution = 24, haze = TRUE, haze_density_threshold = .2)
  expect_true(puff$shape_info[[1]]$medium$haze)
  expect_equal(puff$shape_info[[1]]$medium$haze_density_threshold, .2)
})

# Test each air interval analytically, rather than comparing noisy cloud images.
test_that("per-medium haze controls compose camera radiance, alpha, and nested intervals", {
  skip_if_not_installed("skymodelr")
  sky = sky_light(
    40.7,
    -74,
    as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York"),
    # Enable the global gate so these checks exercise each medium's controls.
    haze_in_volumes = TRUE,
    deferred_haze = FALSE,
    # These exact interval comparisons need no angular sampling noise.
    # Native filter tests separately check sampled means and subdivisions.
    haze_filter = FALSE,
    sampling_resolution = 16
  )
  description = tryCatch(prepare_infinite_light(sky), error = function(e) {
    if (
      grepl("dataset|download|file", conditionMessage(e), ignore.case = TRUE)
    ) {
      skip(conditionMessage(e))
    }
    stop(e)
  })
  settings = list(
    lookfrom = c(0, 200, 0),
    lookat = c(0, 200, 6000),
    fov = 0,
    ortho_dimensions = c(.002, .002),
    width = 4,
    height = 4,
    samples = 4,
    sample_method = "sobol",
    integrator_type = "nee",
    max_depth = 1,
    min_variance = 0,
    aperture = 0,
    denoise = FALSE,
    bloom = FALSE,
    clamp_value = Inf,
    tonemap = "raw",
    auto_exposure = FALSE,
    ambient_light = FALSE,
    preview = FALSE,
    interactive = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    parallel = FALSE
  )
  render = function(scene, light = sky, camera = settings) {
    set.seed(42)
    do.call(
      render_scene,
      c(list(scene = add_infinite_light(scene, light)), camera)
    )
  }
  panel = xy_rect(
    y = 200,
    z = 6000,
    xwidth = 20000,
    ywidth = 20000,
    material = light("white", intensity = 1),
    flipped = TRUE
  )
  shell = cube(y = 200, z = 3000, xwidth = 20000, ywidth = 20000, zwidth = 4000)
  excluded = set_medium(shell, homogeneous_medium(sigma_s = 0, haze = FALSE))
  enabled = set_medium(shell, homogeneous_medium(sigma_s = 0, haze = TRUE))
  outside_sky = sky
  outside_sky$haze_in_volumes = FALSE
  reference = render(add_object(panel, enabled), outside_sky)
  expect_equal(render(add_object(panel, excluded)), reference, tolerance = 1e-5)
  solid = set_scene_material(
    shell,
    subsurface(sigma_a = 0, sigma_s = 0, refraction = 1)
  )
  expect_equal(render(add_object(panel, solid)), reference, tolerance = 1e-5)
  expect_equal(
    render(add_object(panel, create_instances(excluded, x = 0))),
    reference,
    tolerance = 1e-5
  )
  # Missing fields preserve the default when loading older descriptions.
  legacy = homogeneous_medium(sigma_s = 0)
  legacy$haze = NULL
  expect_equal(
    render(add_object(panel, set_medium(shell, legacy))),
    render(add_object(panel, enabled)),
    tolerance = 1e-5
  )

  inner = set_medium(
    cube(y = 200, z = 3000, xwidth = 10000, ywidth = 10000, zwidth = 2000),
    homogeneous_medium(sigma_s = 0, haze = TRUE)
  )
  nested = add_object(add_object(panel, excluded), inner)
  q = query_prague_atmosphere(
    description,
    rbind(c(0, 200, 0), c(0, 200, 2000), c(0, 200, 5000)),
    matrix(rep(c(0, 0, 1), 3), ncol = 3, byrow = TRUE),
    c(1000, 2000, 1000)
  )
  expected = q$inscatter[1, ] +
    q$transmission[1, ] *
      (q$inscatter[2, ] +
        q$transmission[2, ] * (q$inscatter[3, ] + q$transmission[3, ]))
  for (deferred in c(FALSE, TRUE)) {
    sky$deferred_haze = deferred
    actual = render(nested)
    expect_equal(
      vapply(1:3, function(i) mean(actual[,, i]), numeric(1)),
      expected,
      tolerance = 3e-4
    )
  }
  # The global exclusion wins even over an explicitly enabled inner medium.
  expect_equal(render(nested, outside_sky), reference, tolerance = 3e-4)

  # The grid has no scattering candidates: its density ramp still excludes
  # exactly z = 3000..5000. This guards against classification only at events.
  ramp = grid_medium(
    array(c(0, 1), c(1, 1, 2)),
    sigma_s = 0,
    bounds = rbind(c(-10000, -10000, -2000), c(10000, 10000, 2000)),
    haze = TRUE,
    haze_density_threshold = .5
  )
  threshold_scene = add_object(panel, set_medium(shell, ramp))
  equivalent = add_object(
    panel,
    set_medium(
      cube(y = 200, z = 4000, xwidth = 20000, ywidth = 20000, zwidth = 2000),
      homogeneous_medium(sigma_s = 0, haze = FALSE)
    )
  )
  for (deferred in c(FALSE, TRUE)) {
    sky$deferred_haze = deferred
    expect_equal(render(threshold_scene), render(equivalent), tolerance = 3e-4)
    camera = settings
    camera$lookfrom = c(0, 200, 3500)
    expect_equal(
      render(threshold_scene, camera = camera),
      render(equivalent, camera = camera),
      tolerance = 3e-4
    )
    expect_equal(
      render(threshold_scene, outside_sky),
      reference,
      tolerance = 3e-4
    )
  }
  # Turning haze off overrides a positive cutoff, including empty ramp cells.
  ramp$haze = FALSE
  expect_equal(
    render(add_object(panel, set_medium(shell, ramp))),
    reference,
    tolerance = 3e-4
  )
  ramp$haze = TRUE
  ramp["haze_density_threshold"] = list(NULL)
  expect_equal(
    render(add_object(panel, set_medium(shell, ramp))),
    render(add_object(panel, enabled)),
    tolerance = 3e-4
  )
  ramp$haze_density_threshold = 2
  expect_equal(
    render(add_object(panel, set_medium(shell, ramp))),
    render(add_object(panel, enabled)),
    tolerance = 3e-4
  )

  # The camera and floor are below both boxes. Only direct-light connections
  # enter the ramp, so the same explicit half-box checks shadow selection.
  ramp$haze_density_threshold = .5
  floor = xz_rect(
    y = -10000,
    z = 3500,
    xwidth = 30000,
    zwidth = 30000,
    material = diffuse("white")
  )
  shadow_camera = settings
  shadow_camera$lookfrom = c(0, -9999, 3500)
  shadow_camera$lookat = c(0, -10000, 3500)
  shadow_camera$camera_up = c(0, 0, 1)
  shadow_camera$max_depth = 2
  shadow_camera$samples = 32
  shadow_reference = add_object(
    floor,
    set_medium(
      cube(y = 200, z = 4000, xwidth = 20000, ywidth = 20000, zwidth = 2000),
      homogeneous_medium(sigma_s = 0, haze = FALSE)
    )
  )
  for (deferred in c(FALSE, TRUE)) {
    sky$deferred_haze = deferred
    expect_equal(
      render(
        add_object(floor, set_medium(shell, ramp)),
        camera = shadow_camera
      ),
      render(shadow_reference, camera = shadow_camera),
      tolerance = 3e-4
    )
  }

  settings$lookat = c(0, 3800, 6000)
  settings$transparent_background = TRUE
  empty = sphere(x = 50000, radius = 1)
  for (origin in list(c(0, 200, 0), c(0, 200, 3000))) {
    settings$lookfrom = origin
    expect_equal(
      render(add_object(empty, excluded)),
      render(add_object(empty, enabled), outside_sky),
      tolerance = 3e-4
    )
  }
  # Scattering exercises direct-light connections as well as the camera path.
  scattering = homogeneous_medium(
    sigma_s = c(0.0002, 0.0004, 0.0006),
    g = 0.6,
    haze = FALSE
  )
  settings$lookfrom = c(0, 200, 0)
  settings$max_depth = 6
  settings$samples = 16
  for (deferred in c(FALSE, TRUE)) {
    sky$deferred_haze = outside_sky$deferred_haze = deferred
    tagged = render(add_object(empty, set_medium(shell, scattering)))
    scattering$haze = TRUE
    global = render(
      add_object(empty, set_medium(shell, scattering)),
      outside_sky
    )
    scattering$haze = FALSE
    expect_equal(tagged, global, tolerance = 1e-5)
  }
  # Native validation also protects manually edited or serialized descriptions.
  invalid = excluded
  invalid$shape_info[[1]]$medium$haze = NA
  expect_error(render(invalid), "haze")
})
