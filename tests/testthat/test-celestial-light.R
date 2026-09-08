test_that("celestial descriptions validate and preserve observation settings", {
  time = as.POSIXct("2026-01-28 21:00:00", tz = "Pacific/Auckland")
  sun = sun_light(-36.9, 174.8, time)
  moon = moon_light(-36.9, 174.8, time, moon_args = list(earthshine = FALSE))
  expect_s3_class(sun, "ray_infinite_light")
  expect_equal(sun$type, "sun")
  expect_equal(moon$type, "moon")
  expect_identical(moon$datetime, time)
  expect_false(moon$moon_args$earthshine)
  expect_output(print(moon), "moon")
  scene = sphere() |> add_infinite_light(sun) |> add_infinite_light(moon)
  expect_named(list_infinite_lights(scene), c("sun", "moon"))
  expect_identical(
    get_infinite_light(unserialize(serialize(scene, NULL)), "moon"),
    moon
  )
  expect_equal(nrow(scene), 1)
  for (bad in list(0, -1, NA_real_, Inf, 15, 32.5, c(32, 32), "32")) {
    expect_error(sun_light(0, 0, time, resolution = bad), "resolution")
  }
  expect_error(sun_light(0, 0, time, sky_args = list(hosek = TRUE)), "Prague")
  expect_error(
    sun_light(0, 0, time, sky_args = list(render_mode = "all")),
    "Unsupported"
  )
  expect_error(
    moon_light(0, 0, time, sky_args = list(altitude = -1)),
    "altitude"
  )
  expect_error(
    moon_light(0, 0, time, moon_args = list(earthshine = NA)),
    "earthshine"
  )
  expect_error(
    moon_light(0, 0, time, moon_args = list(phase = 90)),
    "Unsupported"
  )
  expect_error(
    moon_light(0, 0, time, moon_args = list(moon_extinction_kV = Inf)),
    "moon_extinction"
  )
})

test_that("celestial preparation consumes public radiance and geometry unchanged", {
  skip_if_not_installed("skymodelr")
  skip_if_not(
    all(
      c("generate_sun_disk", "generate_moon_disk") %in%
        getNamespaceExports("skymodelr")
    ),
    "Installed skymodelr predates the disk exports"
  )
  calls = c(sun = 0, moon = 0)
  seen = NULL
  image = array(rep(c(3, 2, -0.1), each = 256), c(16, 16, 3))
  generator = function(body, ...) {
    calls[body] <<- calls[body] + 1
    seen <<- list(...)
    list(
      image = image,
      azimuth_deg = 90,
      elevation_deg = 30,
      angular_diameter_deg = 0.53,
      projection = "rectilinear"
    )
  }
  testthat::local_mocked_bindings(
    generate_sun_disk = function(...) generator("sun", ...),
    generate_moon_disk = function(...) generator("moon", ...),
    .package = "skymodelr"
  )
  time = as.POSIXct("2026-01-28", tz = "UTC") + as.numeric(Sys.time()) / 1e6
  for (body in c("sun", "moon")) {
    light = if (body == "sun") {
      sun_light(12, 34, time, resolution = 16)
    } else {
      moon_light(
        12,
        34,
        time,
        moon_args = list(earthshine = FALSE),
        resolution = 16
      )
    }
    light$sky_args = list(altitude = 500, hosek = FALSE)
    prepared = prepare_infinite_light(light)
    withr::defer(unlink(c(
      prepared$filename,
      sub("exr$", "rds", prepared$filename)
    )))
    expect_equal(prepared$type, "disk")
    expect_equal(prepared$direction, c(-sqrt(3) / 2, 0.5, 0), tolerance = 1e-12)
    expect_equal(prepared$angular_diameter, 0.53)
    expect_identical(seen$datetime, time)
    expect_equal(seen$lat, 12)
    expect_equal(seen$lon, 34)
    expect_equal(seen$altitude, 500)
    expect_equal(seen$resolution, 16)
    expect_null(seen$hosek)
    if (body == "moon") {
      expect_false(seen$earthshine)
    }
    pixels = rayimage::ray_read_image(prepared$filename)
    expect_equal(
      as.numeric(pixels[,, 1:3]),
      as.numeric(image),
      tolerance = 1e-6
    )
    expect_equal(as.numeric(pixels[,, 4]), rep(1, 256))
    light$intensity = 2
    light$rotation = 30
    second = prepare_infinite_light(light)
    expect_equal(unname(calls[body]), 1)
    expect_identical(prepared$filename, second$filename)
    expect_equal(second$intensity, 2)
    expect_equal(second$rotation, 30)
  }
})

test_that("invalid public disk output is rejected before texture preparation", {
  valid = list(
    image = array(1, c(16, 16, 3)),
    azimuth_deg = 0,
    elevation_deg = 30,
    angular_diameter_deg = 0.53,
    projection = "rectilinear"
  )
  for (update in list(
    list(image = array(1, c(16, 16, 4))),
    list(image = array(NA_real_, c(16, 16, 3))),
    list(projection = "latlong"),
    list(azimuth_deg = 360),
    list(elevation_deg = 91),
    list(angular_diameter_deg = 0),
    list(angular_diameter_deg = Inf)
  )) {
    expect_error(
      validate_skymodelr_disk(utils::modifyList(valid, update)),
      "invalid disk"
    )
  }
})

test_that("native disks preserve angular coverage and additive backgrounds in all integrators", {
  skip_if_not_installed("libopenexr")
  directory = withr::local_tempdir()
  file = file.path(directory, "disk.exr")
  libopenexr::write_exr(
    file,
    matrix(0.5, 8, 8),
    matrix(0.25, 8, 8),
    matrix(0.125, 8, 8)
  )
  scene = sphere(x = 100) |>
    add_infinite_light(celestial_disk(
      file,
      c(0, 0, 1),
      0.53,
      clip_horizon = FALSE
    ))
  settings = list(
    width = 4,
    height = 4,
    samples = 1,
    aperture = 0,
    fov = 0,
    ortho_dimensions = c(1, 1),
    lookfrom = c(0, 0, 0),
    lookat = c(0, 0, 1),
    preview = FALSE,
    interactive = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    parallel = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = "raw"
  )
  for (integrator in c("nee", "rtiow", "basic")) {
    render = function(scene, ...) {
      do.call(
        render_scene,
        c(
          list(scene = scene, integrator_type = integrator),
          settings,
          list(...)
        )
      )
    }
    image = render(scene)
    expect_equal(as.numeric(image[,, 1]), rep(0.5, 16), tolerance = 1e-6)
    image = render(scene, rotate_env = 1)
    expect_equal(as.numeric(image[,, 1:3]), rep(0, 48), tolerance = 1e-6)
    added = add_infinite_light(scene, infinite_light(file, name = "sky"))
    image = render(added)
    expect_equal(as.numeric(image[,, 1]), rep(1, 16), tolerance = 1e-6)
    transparent = render(added, transparent_background = TRUE)
    expect_equal(as.numeric(transparent[,, 4]), rep(0, 16))
  }
})
