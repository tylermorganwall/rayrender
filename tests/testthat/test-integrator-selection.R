test_that("automatic integration selects a usable Metal device or CPU NEE", {
  for (fun in list(
    render_scene,
    render_animation,
    prepare_scene_list,
    render_animation_camera,
    render_camera_batch
  )) {
    expect_identical(formals(fun)$integrator_type, "auto")
  }
  for (available in c(FALSE, TRUE)) {
    local_mocked_bindings(
      cpp_has_metal = function() available,
      .package = "rayrender"
    )
    info = prepare_scene_list(sphere())$render_info
    expect_identical(info$metal_wavefront, available)
    expect_identical(info$integrator_type, 1L)
    info = prepare_scene_list(sphere(), integrator_type = "auto")$render_info
    expect_identical(info$metal_wavefront, available)
  }
})

test_that("explicit integrators do not perform automatic device selection", {
  local_mocked_bindings(
    cpp_has_metal = function() stop("unexpected capability probe"),
    .package = "rayrender"
  )
  for (integrator in c("nee", "rtiow", "basic", "metal")) {
    info = prepare_scene_list(
      sphere(),
      integrator_type = integrator
    )$render_info
    expect_identical(info$metal_wavefront, integrator == "metal")
  }
})

test_that("PBRT render arguments leave device selection automatic", {
  path = withr::local_tempfile(fileext = ".pbrt")
  writeLines(
    c('WorldBegin', 'Shape "sphere" "float radius" [1]', 'WorldEnd'),
    path
  )
  expect_identical(read_pbrt(path)$render_args$integrator_type, "auto")
})

test_that("default still and animation renders match explicit device selection", {
  available = cpp_has_metal()
  selected = if (available) "metal" else "nee"
  args = list(
    scene = sphere(material = light("orange")),
    width = 12,
    height = 12,
    samples = 4,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    parallel = FALSE,
    bloom = FALSE
  )
  for (render in list(render_scene, render_animation)) {
    set.seed(71)
    automatic = do.call(render, args)
    set.seed(71)
    explicit = do.call(render, c(args, list(integrator_type = selected)))
    if (is.list(automatic)) {
      automatic = automatic[[1]]
      explicit = explicit[[1]]
    }
    expect_equal(as.numeric(automatic), as.numeric(explicit), tolerance = 0)
    expect_identical(isTRUE(attr(automatic, "wavefront")$used), available)
  }
})

test_that("scene requirements select NEE without changing other integrators", {
  # Integrator selection must not require downloading sky data or rendering disks.
  local_mocked_bindings(
    prepare_scene_infinite_lights = function(lights) list(),
    .package = "rayrender"
  )
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
    add_infinite_light(surface, infinite_light(image)),
    add_infinite_light(surface, sky_light_image(0, 0, time)),
    add_infinite_light(surface, sun_light(0, 0, time)),
    add_infinite_light(surface, moon_light(0, 0, time))
  )
  requires_nee = list(
    medium,
    set_medium(sphere(), homogeneous_medium(sigma_s = .01)),
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
  local_mocked_bindings(
    cpp_has_metal = function() FALSE,
    .package = "rayrender"
  )
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
