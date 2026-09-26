test_that("disk lights validate and preserve scene metadata without files", {
  light = disk_light("#ff8000", 3, 12, c(0, 0, 1e300), name = "key")
  expect_s3_class(light, "ray_infinite_light")
  expect_equal(light$color, c(1, 128 / 255, 0))
  expect_equal(light$direction, c(0, 0, 1))
  expect_null(light$filename)
  expect_identical(prepare_infinite_light(light), light)
  scene = add_infinite_light(sphere(), light)
  expect_identical(get_infinite_light(scene, "key"), light)
  expect_output(print(light), "angular diameter: 12")
  expect_error(add_infinite_light(scene, light), "already exists")
  expect_length(list_infinite_lights(remove_infinite_light(scene, "key")), 0)
  for (diameter in list(0, -1, 180, Inf, NA_real_, c(1, 2))) {
    expect_error(disk_light(angular_diameter = diameter), "angular diameter")
  }
  for (direction in list(c(0, 0, 0), c(1, 2), c(0, NA, 1), c(0, Inf, 1))) {
    expect_error(disk_light(direction = direction), "direction")
  }
  for (color in list(
    numeric(),
    c(1, 0),
    c(1, 0, 0, 1),
    c(1, NA, 0),
    c(2, 0, 0)
  )) {
    expect_error(disk_light(color = color))
  }
  expect_error(disk_light(intensity = -1), "intensity")
  expect_error(disk_light(rotation = Inf), "rotation")
  expect_error(disk_light(name = ""), "name")
})

test_that("uniform disks render color, angular extent, rotation, and additive radiance", {
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
  light = disk_light(c(1, 0.5, 0.25), intensity = 2, direction = c(0, 0, 1))
  scene = add_infinite_light(sphere(x = 100), light)
  expected = rep(c(2, 1, 0.5), each = 16)
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
    expect_equal(as.numeric(render(scene)[,, 1:3]), expected, tolerance = 1e-6)
    expect_equal(as.numeric(render(scene, rotate_env = 1)[,, 1:3]), rep(0, 48))
    expect_equal(
      as.numeric(render(scene, rotate_env = 0.1)[,, 1:3]),
      expected,
      tolerance = 1e-6
    )
    added = add_infinite_light(scene, light, name = "fill")
    expect_equal(
      as.numeric(render(added)[,, 1:3]),
      2 * expected,
      tolerance = 1e-6
    )
    light$intensity = 0
    dark = add_infinite_light(scene, light, replace = TRUE)
    expect_equal(as.numeric(render(dark)[,, 1:3]), rep(0, 48))
    light$intensity = 2
  }
})
