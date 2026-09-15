test_that("scene requirements select NEE without changing other integrators", {
  # Integrator selection must not require downloading sky data or rendering disks.
  local_mocked_bindings(prepare_scene_infinite_lights = function(lights) list())
  time = as.POSIXct("2026-06-21 12:00:00", tz = "UTC")
  sky = sky_light(0, 0, time, sun = FALSE, moon = FALSE)
  fixed_sky = sky_light(
    0,
    0,
    time,
    sun = FALSE,
    moon = FALSE,
    haze = FALSE,
    query_altitude = FALSE
  )
  image = withr::local_tempfile(fileext = ".png")
  png::writePNG(array(1, c(2, 4, 3)), image)
  surface = sphere()
  medium = set_medium(cube(), homogeneous_medium())
  unchanged = list(
    surface,
    sphere(material = diffuse(fog = TRUE)),
    add_infinite_light(surface, infinite_light(image)),
    add_infinite_light(surface, sky_light_image(0, 0, time)),
    add_infinite_light(surface, sun_light(0, 0, time)),
    add_infinite_light(surface, moon_light(0, 0, time))
  )
  requires_nee = list(
    medium,
    create_instances(medium, x = c(-2, 2)),
    add_infinite_light(surface, sky),
    add_infinite_light(surface, fixed_sky),
    add_infinite_light(medium, sky)
  )
  if (requireNamespace("ambient", quietly = TRUE)) {
    requires_nee = c(requires_nee, list(cloud(resolution = 24)))
  }
  for (integrator in c("nee", "rtiow", "basic")) {
    for (scene in unchanged) {
      info = prepare_scene_list(scene, integrator_type = integrator)$render_info
      expect_identical(
        info$integrator_type,
        match(integrator, c("nee", "rtiow", "basic"))
      )
    }
    for (scene in requires_nee) {
      info = prepare_scene_list(scene, integrator_type = integrator)$render_info
      expect_identical(info$integrator_type, 1L)
    }
  }
  for (scene in c(unchanged, requires_nee)) {
    expect_error(
      prepare_scene_list(scene, integrator_type = "typo"),
      "not recognized as valid"
    )
  }
})

test_that("automatically selected NEE recognizes emitting media as illumination", {
  scene = set_medium(cube(), homogeneous_medium(sigma_a = 1, emission = 1))
  for (integrator in c("nee", "rtiow", "basic")) {
    info = prepare_scene_list(
      scene,
      integrator_type = integrator,
      ambient_light = NULL
    )$render_info
    expect_false(info$ambient_light)
    expect_true(
      prepare_scene_list(
        scene,
        integrator_type = integrator,
        ambient_light = TRUE
      )$render_info$ambient_light
    )
  }
  expect_true(
    prepare_scene_list(sphere(), ambient_light = NULL)$render_info$ambient_light
  )
})

test_that("automatic NEE renders the same volume through still and animation APIs", {
  absorption = c(0.2, 0.4, 0.6)
  scene = set_medium(
    cube(),
    homogeneous_medium(sigma_a = absorption, sigma_s = 0)
  ) |>
    add_camera(camera(
      lookfrom = c(0, 0, 3),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(0.5, 0.5),
      aperture = 0
    ))
  args = list(
    scene = scene,
    width = 4,
    height = 4,
    samples = 1,
    ambient_light = FALSE,
    transparent_background = TRUE,
    min_variance = 0,
    tonemap = "raw",
    bloom = FALSE,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE
  )
  reference = do.call(render_scene, c(args, list(integrator_type = "nee")))
  expect_equal(
    as.numeric(reference[,, 4]),
    rep(1 - mean(exp(-absorption)), 16),
    tolerance = 2e-5
  )
  for (selection in list(
    list(),
    list(integrator_type = "rtiow"),
    list(integrator_type = "basic")
  )) {
    still = do.call(render_scene, c(args, selection))
    frames = do.call(render_scene, c(args, selection, list(mode = "animation")))
    animation = do.call(render_animation, c(args, selection))
    expect_equal(as.numeric(still), as.numeric(reference), tolerance = 1e-6)
    expect_equal(
      as.numeric(frames[[1]]),
      as.numeric(reference),
      tolerance = 1e-6
    )
    expect_equal(
      as.numeric(animation[[1]]),
      as.numeric(reference),
      tolerance = 1e-6
    )
  }
})
