test_that("attached point lights validate and survive scene composition", {
  light = spot_light(direction = c(0, -20, 0), intensity = 5)
  expect_equal(light$direction, c(0, -1, 0))
  expect_error(point_light(position = c(0, NA, 1)), "finite")
  expect_error(point_light(intensity = -1), "nonnegative")
  expect_error(spot_light(direction = c(0, 0, 0)), "zero")
  expect_error(spot_light(cone_angle = 10, falloff_angle = 11), "Require")
  scene = add_light(sphere(), light, name = "key")
  expect_equal(nrow(scene), 1)
  expect_equal(get_light(scene, "key")$intensity, 5)
  expect_error(add_light(scene, light, name = "key"), "already exists")
  scene = add_object(scene, cube(x = 4))
  expect_named(list_lights(scene), "key")
  expect_error(
    add_object(scene, add_light(sphere(), light, name = "key")),
    "Duplicate"
  )
  expect_length(list_lights(remove_light(scene, "key")), 0)
  scene = add_light(
    scene,
    point_light(intensity = 10),
    name = "key",
    replace = TRUE
  )
  expect_equal(get_light(scene, "key")$type, "point")
  prepared = prepare_scene_list(scene, integrator_type = "basic")
  expect_identical(prepared$render_info$integrator_type, 1L)
  expect_false(prepared$render_info$ambient_light)
  expect_length(prepared$render_info$point_lights, 1)
})

#' @param scene Scene to render.
#' @return Centered orthographic linear RGB mean.
#' @keywords internal
point_light_test_pixel = function(scene) {
  image = render_scene(
    scene,
    width = 4,
    height = 4,
    samples = 8,
    lookfrom = c(0, 0, 5),
    lookat = c(0, 0, 0),
    camera_up = c(0, 1, 0),
    fov = 0,
    ortho_dimensions = c(1e-4, 1e-4),
    aperture = 0,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = "raw",
    min_variance = 0,
    max_depth = 2
  )
  apply(image[,, 1:3, drop = FALSE], 3, mean)
}

test_that("native point lights match analytic irradiance and inverse square falloff", {
  skip_on_cran()
  surface = xy_rect(
    xwidth = 10,
    ywidth = 10,
    material = diffuse(color = "white")
  )
  near = point_light_test_pixel(add_light(
    surface,
    point_light(position = c(0, 0, 2), intensity = 4)
  ))
  far = point_light_test_pixel(add_light(
    surface,
    point_light(position = c(0, 0, 4), intensity = 4)
  ))
  expect_equal(near, rep(1 / pi, 3), tolerance = 2e-4)
  expect_equal(near / far, rep(4, 3), tolerance = 2e-4)
  two = add_light(
    add_light(
      surface,
      point_light(position = c(0, 0, 2), intensity = 1, name = "a")
    ),
    point_light(position = c(0, 0, 2), intensity = 3, name = "b")
  )
  expect_equal(point_light_test_pixel(two), near, tolerance = 2e-4)
})

test_that("point shadows respect blockers and RGB homogeneous extinction", {
  skip_on_cran()
  surface = xy_rect(
    xwidth = 10,
    ywidth = 10,
    material = diffuse(color = "white")
  )
  scene = add_light(surface, point_light(position = c(2, 0, 2), intensity = 8))
  reference = point_light_test_pixel(scene)
  expect_equal(reference, rep(sqrt(.5) / pi, 3), tolerance = 2e-4)
  block = cube(x = 1, z = 1, width = .6, material = diffuse("black"))
  expect_equal(
    point_light_test_pixel(add_object(scene, block)),
    rep(0, 3),
    tolerance = 1e-8
  )
  fog = set_medium(
    block,
    homogeneous_medium(sigma_a = c(.2, .5, 1), sigma_s = 0)
  )
  expect_equal(
    point_light_test_pixel(add_object(scene, fog)) / reference,
    exp(-c(.2, .5, 1) * .6 * sqrt(2)),
    tolerance = 3e-4
  )
  # A constant grid should agree with the homogeneous reference in expectation.
  grid = set_medium(
    block,
    grid_medium(
      array(1, c(2, 2, 2)),
      sigma_a = .5,
      sigma_s = 0,
      bounds = rbind(c(-2, -2, -2), c(2, 2, 2))
    )
  )
  set.seed(902)
  expect_equal(
    mean(point_light_test_pixel(add_object(scene, grid)) / reference),
    exp(-.5 * .6 * sqrt(2)),
    tolerance = .15
  )
})

test_that("spotlight cone axis and falloff affect rendered illumination", {
  skip_on_cran()
  surface = xy_rect(
    xwidth = 10,
    ywidth = 10,
    material = diffuse(color = "white")
  )
  beam = function(direction, cone_angle = 60, falloff_angle = 30) {
    point_light_test_pixel(add_light(
      surface,
      spot_light(
        position = c(0, 0, 2),
        direction = direction,
        cone_angle = cone_angle,
        falloff_angle = falloff_angle,
        intensity = 4
      )
    ))
  }
  expect_equal(beam(c(0, 0, -1)), rep(1 / pi, 3), tolerance = 2e-4)
  expect_equal(beam(c(0, 0, 1)), rep(0, 3), tolerance = 1e-8)
  cosine = (cos(pi / 3) + cos(pi / 6)) / 2
  expect_equal(
    beam(c(sqrt(1 - cosine^2), 0, -cosine)),
    rep(.5 / pi, 3),
    tolerance = 2e-4
  )
})
