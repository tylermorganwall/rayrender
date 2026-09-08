infinite_light_fixture = function() {
  path = tempfile(fileext = ".png")
  png::writePNG(array(0.5, c(4, 8, 3)), path)
  path
}

test_that("image infinite lights validate their descriptions", {
  path = infinite_light_fixture()
  on.exit(unlink(path))
  light = infinite_light(path, intensity = 2, rotation = 45, name = "fill")
  expect_s3_class(light, "ray_infinite_light")
  expect_equal(light$intensity, 2)
  expect_equal(light$rotation, 45)
  expect_output(print(light), "Infinite light 'fill'")
  for (value in list(-1, NA_real_, Inf, c(1, 2), "1", NULL)) {
    expect_error(infinite_light(path, intensity = value), "intensity")
  }
  for (value in list(NA_real_, Inf, c(1, 2), "1", NULL)) {
    expect_error(infinite_light(path, rotation = value), "rotation")
  }
  for (value in list("", NA_character_, c("a", "b"), 1)) {
    expect_error(infinite_light(path, name = value), "name")
  }
  expect_error(infinite_light(tempdir()), "directory")
  expect_error(infinite_light(paste0(path, "missing")), "does not exist")
  expect_error(infinite_light(NA_character_), "filename")
})

test_that("scene lights can be added, replaced, queried and removed", {
  path = infinite_light_fixture()
  on.exit(unlink(path))
  scene = sphere() |>
    add_camera(camera(name = "main")) |>
    add_infinite_light(infinite_light(path, name = "key")) |>
    add_infinite_light(infinite_light(path, name = "fill", intensity = 0.25))
  expect_equal(names(list_infinite_lights(scene)), c("key", "fill"))
  expect_equal(get_infinite_light(scene, "fill")$intensity, 0.25)
  expect_equal(get_camera(scene)$name, "main")
  expect_equal(nrow(scene), 1)
  expect_error(
    add_infinite_light(scene, infinite_light(path, name = "key")),
    "already exists"
  )
  scene = add_infinite_light(
    scene,
    infinite_light(path, intensity = 3),
    name = "key",
    replace = TRUE
  )
  expect_equal(get_infinite_light(scene, "key")$intensity, 3)
  expect_error(get_infinite_light(scene, "missing"), "identify")
  expect_error(remove_infinite_light(scene, "missing"), "identify")
  scene = remove_infinite_light(scene, "key")
  expect_named(list_infinite_lights(scene), "fill")
  expect_length(list_infinite_lights(remove_infinite_light(scene, "fill")), 0)
})

test_that("scene operations and serialization preserve lights independently of geometry", {
  path = infinite_light_fixture()
  on.exit(unlink(path))
  a = sphere() |> add_infinite_light(infinite_light(path, name = "a"))
  b = cube() |> add_infinite_light(infinite_light(path, name = "b"))
  scene = add_object(a, b)
  expect_named(list_infinite_lights(scene), c("a", "b"))
  expect_error(add_object(a, a), "Duplicate infinite light")
  transformed = scene |>
    group_objects(translate = c(1, 2, 3), angle = c(20, 45, 0)) |>
    set_scene_material(diffuse(color = "red"))
  expect_identical(
    list_infinite_lights(transformed),
    list_infinite_lights(scene)
  )
  expect_identical(
    list_infinite_lights(unserialize(serialize(scene, NULL))),
    list_infinite_lights(scene)
  )
  expect_named(list_infinite_lights(add_object(sphere(), a)), "a")
})

test_that("preparation adds legacy images and disables implicit fallback", {
  path = infinite_light_fixture()
  on.exit(unlink(path))
  scene = sphere() |>
    add_infinite_light(infinite_light(path, intensity = 0, rotation = 20))
  info = prepare_scene_list(scene, ambient_light = NULL)$render_info
  expect_true(info$hasbackground)
  expect_false(info$ambient_light)
  expect_length(info$infinite_lights, 1)
  expect_equal(info$infinite_lights[[1]]$intensity, 0)
  info = prepare_scene_list(
    scene,
    environment_light = path,
    intensity_env = 2
  )$render_info
  expect_length(info$infinite_lights, 2)
  expect_equal(info$infinite_lights[[1]]$rotation, 20)
  expect_equal(info$infinite_lights[[2]]$intensity, 2)
  expect_true(
    prepare_scene_list(sphere(), ambient_light = NULL)$render_info$ambient_light
  )
  altered = scene
  attr(altered, "ray_infinite_lights")[[1]]$intensity = -1
  expect_error(prepare_scene_list(altered), "intensity")
  unlink(path)
  expect_error(prepare_scene_list(scene), "does not exist")
})

test_that("sky lights retain location, time and generation settings", {
  time = as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York")
  light = sky_light(
    40.7,
    -74,
    time,
    sky_args = list(hosek = FALSE, resolution = 64)
  )
  expect_s3_class(light, "ray_infinite_light")
  expect_identical(light$datetime, time)
  expect_identical(light$sky_args, list(hosek = FALSE, resolution = 64))
  scene = sphere() |> add_infinite_light(light)
  expect_identical(get_infinite_light(scene, "sky"), light)
  expect_output(print(light), "location: 40.7, -74")
  for (value in list(NA_real_, Inf, 91, c(0, 1))) {
    expect_error(sky_light(value, 0, time), "lat")
  }
  expect_error(sky_light(0, 181, time), "long")
  expect_error(sky_light(0, 0, "2026-06-21"), "POSIXct")
  expect_error(sky_light(0, 0, time, list(32)), "named list")
  expect_error(
    sky_light(0, 0, time, list(filename = "sky.exr")),
    "cannot override"
  )
  expect_error(sky_light(0, 0, time, list(lat = 1)), "cannot override")
  expect_error(sky_light(0, 0, time, intensity = -1), "intensity")
})

test_that("static sky generation is cached independently of intensity and rotation", {
  skip_if_not_installed("skymodelr")
  calls = 0
  seen = NULL
  testthat::local_mocked_bindings(
    generate_sky_latlong = function(
      lat,
      lon,
      datetime,
      filename,
      resolution,
      ...
    ) {
      calls <<- calls + 1
      seen <<- list(
        lat = lat,
        lon = lon,
        datetime = datetime,
        resolution = resolution
      )
      # This test checks preparation/caching; native image decoding is tested by renders.
      writeBin(as.raw(1), filename)
    },
    .package = "skymodelr"
  )
  time = as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York")
  # A unique timestamp separates this mocked cache from actual sky images.
  time = time + as.numeric(Sys.time()) / 1e6
  sky = sky_light(40.7, -74, time, sky_args = list(resolution = 32))
  a = prepare_infinite_light(sky)
  on.exit(unlink(a$filename))
  sky$intensity = 0.5
  sky$rotation = 45
  b = prepare_infinite_light(sky)
  expect_equal(calls, 1)
  expect_identical(a$filename, b$filename)
  expect_identical(seen$datetime, time)
  expect_equal(seen$lon, -74)
  expect_equal(b$rotation, 45)
  expect_equal(b$intensity, 0.5)
  expect_identical(b$type, "image")
})

test_that("all integrators sum scene image lights and retain transparent backgrounds", {
  skip_if_not_installed("libopenexr")
  directory = withr::local_tempdir()
  a = file.path(directory, "a.exr")
  b = file.path(directory, "b.exr")
  libopenexr::write_exr(
    a,
    matrix(0.25, 4, 8),
    matrix(0.125, 4, 8),
    matrix(0, 4, 8)
  )
  libopenexr::write_exr(
    b,
    matrix(0, 4, 8),
    matrix(0.125, 4, 8),
    matrix(0.5, 4, 8)
  )
  # The sphere is outside the camera view, so every pixel samples the background.
  scene = sphere(x = 100) |>
    add_infinite_light(infinite_light(a, name = "a")) |>
    add_infinite_light(infinite_light(b, name = "b", intensity = 0.5))
  render = function(scene, integrator, ...) {
    render_scene(
      scene,
      width = 4,
      height = 4,
      samples = 1,
      integrator_type = integrator,
      aperture = 0,
      interactive = FALSE,
      preview = FALSE,
      plot_scene = FALSE,
      parallel = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = "raw",
      ...
    )
  }
  for (integrator in c("nee", "rtiow", "basic")) {
    image = render(scene, integrator)
    expected = c(0.25, 0.1875, 0.25)
    for (channel in 1:3) {
      expect_equal(
        as.numeric(image[,, channel]),
        rep(expected[channel], 16),
        tolerance = 1e-6
      )
    }
    transparent = render(scene, integrator, transparent_background = TRUE)
    expect_equal(as.numeric(transparent[,, 4]), rep(0, 16))
  }
  corrupt = file.path(directory, "broken.exr")
  writeLines("not an EXR", corrupt)
  expect_error(
    render(sphere() |> add_infinite_light(infinite_light(corrupt)), "nee"),
    "load|Load|EXR|exr"
  )
})

test_that("environment rotation agrees across still, animation and camera batches", {
  skip_if_not_installed("libopenexr")
  directory = withr::local_tempdir()
  path = file.path(directory, "directional.exr")
  red = matrix(rep(seq(0, 1, length.out = 32), each = 16), 16, 32)
  libopenexr::write_exr(path, red, 1 - red, red * 0 + 0.1)
  scene = sphere(x = 100) |>
    add_infinite_light(infinite_light(path, rotation = 25))
  cam = camera(
    fov = 0,
    ortho_dimensions = c(2, 2),
    aperture = 0,
    name = "first"
  )
  scene = add_camera(scene, cam)
  args = list(
    scene = scene,
    width = 4,
    height = 4,
    samples = 1,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    parallel = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    rotate_env = 47,
    integrator_type = "nee",
    tonemap = "raw"
  )
  still = do.call(
    render_scene,
    c(args, list(mode = "image", interactive = FALSE))
  )
  frames = do.call(
    render_scene,
    c(args, list(mode = "animation", interactive = FALSE))
  )
  animation = do.call(render_animation, args)
  expect_equal(as.numeric(frames[[1]]), as.numeric(still), tolerance = 1e-6)
  expect_equal(as.numeric(animation[[1]]), as.numeric(still), tolerance = 1e-6)
  args$scene = add_camera(scene, cam, name = "second")
  batch = do.call(
    render_scene,
    c(args, list(camera = "all", mode = "image", interactive = FALSE))
  )
  expect_equal(as.numeric(batch$first), as.numeric(still), tolerance = 1e-6)
  expect_equal(as.numeric(batch$second), as.numeric(still), tolerance = 1e-6)
})
