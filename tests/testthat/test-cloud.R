test_that("clouds are centered volumes with local density and standard transforms", {
  skip_if_not_installed("ambient")
  original = cloud(resolution = 24)
  medium = original$shape_info[[1]]$medium
  expect_s3_class(original, "ray_scene")
  expect_equal(nrow(original), 1)
  expect_equal(c(original$x, original$y, original$z), c(0, 0, 0))
  expect_equal(medium$bounds, rbind(c(-50, -12.5, -37.5), c(50, 12.5, 37.5)))
  expect_equal(dim(medium$density), c(24, 8, 18))
  expect_true(all(is.finite(medium$density)))
  expect_true(all(medium$density >= 0 & medium$density <= 1))
  expect_gt(max(medium$density), 0)
  expect_true(all(medium$density[c(1, 24), , ] == 0))
  expect_true(all(medium$density[, c(1, 8), ] == 0))
  expect_true(all(medium$density[,, c(1, 18)] == 0))
  expect_false(original$shape_info[[1]]$medium_keep_surface)

  moved = cloud(
    x = 20,
    y = 30,
    z = -40,
    angle = c(10, 30, 20),
    order_rotation = c(2, 1, 3),
    scale = c(2, 1, 0.5),
    resolution = 24
  )
  expect_identical(moved$shape_info[[1]]$medium, medium)
  expect_equal(c(moved$x, moved$y, moved$z), c(20, 30, -40))
  expect_equal(moved$transforms[[1]]$angle[[1]], c(10, 30, 20))
  expect_equal(moved$transforms[[1]]$order_rotation[[1]], c(2, 1, 3))
  expect_equal(moved$transforms[[1]]$scale[[1]], c(2, 1, 0.5))
  expect_equal(
    cloud(scale = 2, resolution = 24)$transforms[[1]]$scale[[1]],
    rep(2, 3)
  )

  # Dimensions renormalize vertical extinction; transforms retain world rates.
  taller = cloud(height = 50, resolution = 24)$shape_info[[1]]$medium
  expect_equal(taller$sigma_s, medium$sigma_s / 2)
  expect_equal(taller$sigma_a, medium$sigma_a / 2)
  vacuum = cloud(optical_depth = 0, resolution = 24)$shape_info[[1]]$medium
  expect_equal(vacuum$sigma_s + vacuum$sigma_a, rep(0, 3))
  expect_error(prepare_scene_list(original, integrator_type = "basic"), "nee")
})

test_that("cloud shapes are repeatable without consuming the caller's RNG", {
  skip_if_not_installed("ambient")
  withr::local_seed(103)
  state = .Random.seed
  a = cloud(seed = 8, resolution = 24)
  expect_identical(.Random.seed, state)
  expect_identical(cloud(seed = 8, resolution = 24), a)
  b = cloud(seed = 9, resolution = 24)
  expect_false(identical(
    a$shape_info[[1]]$medium$density,
    b$shape_info[[1]]$medium$density
  ))
  bank = cloud(style = "stratus", seed = 8, resolution = 24)
  expect_gt(max(bank$shape_info[[1]]$medium$density), 0)
  expect_false(identical(
    a$shape_info[[1]]$medium$density,
    bank$shape_info[[1]]$medium$density
  ))
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(cloud(resolution = 24))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("cloud controls reject invalid positions, grids, and transforms", {
  skip_if_not_installed("ambient")
  for (name in c(
    "x",
    "y",
    "z",
    "width",
    "height",
    "depth",
    "seed",
    "resolution",
    "coverage",
    "detail",
    "optical_depth",
    "g"
  )) {
    for (value in list(NULL, NA_real_, Inf, TRUE, "1", c(1, 2))) {
      args = list(resolution = 24)
      args[name] = list(value)
      expect_error(do.call(cloud, args), name)
    }
  }
  invalid = list(
    width = 0,
    height = -1,
    depth = 0,
    seed = 0.5,
    resolution = 24.5,
    coverage = -0.1,
    detail = 1.1,
    optical_depth = -1,
    g = 1,
    angle = c(0, 1),
    order_rotation = c(1, 1, 2),
    scale = c(1, 0, 1)
  )
  for (name in names(invalid)) {
    args = list(resolution = 24)
    args[name] = invalid[name]
    expect_error(do.call(cloud, args), name)
  }
  expect_error(cloud(style = "smoke", resolution = 24), "arg")
})

test_that("cloud density follows object, group, and instance transformations in renders", {
  skip_if_not_installed("ambient")
  args = list(
    width = 4,
    height = 2,
    depth = 3,
    resolution = 24,
    optical_depth = 2
  )
  base = do.call(cloud, args)
  position = c(4, 3, 2)
  angle = c(10, 30, 20)
  scale = c(1.2, 0.8, 1.4)
  order = c(2, 1, 3)
  object = do.call(
    cloud,
    c(
      args,
      list(
        x = position[1],
        y = position[2],
        z = position[3],
        angle = angle,
        scale = scale,
        order_rotation = order
      )
    )
  )
  grouped = group_objects(
    base,
    translate = position,
    angle = angle,
    scale = scale,
    order_rotation = order
  )
  instanced = create_instances(
    base,
    x = position[1],
    y = position[2],
    z = position[3],
    angle_x = angle[1],
    angle_y = angle[2],
    angle_z = angle[3],
    scale_x = scale[1],
    scale_y = scale[2],
    scale_z = scale[3],
    order_rotation = order
  )
  render = function(scene) {
    set.seed(105)
    render_scene(
      scene,
      lookfrom = position + c(3, 4, -6),
      lookat = position,
      width = 16,
      height = 12,
      samples = 8,
      max_depth = 4,
      fov = 40,
      integrator_type = "nee",
      sample_method = "sobol",
      min_variance = 0,
      transparent_background = TRUE,
      aperture = 0,
      denoise = FALSE,
      bloom = FALSE,
      clamp_value = Inf,
      tonemap = "raw",
      auto_exposure = FALSE,
      ambient_light = TRUE,
      preview = FALSE,
      interactive = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      parallel = FALSE
    )
  }
  reference = render(object)
  expect_true(all(is.finite(reference)))
  expect_gt(max(reference[,, 4]), 0.1)
  expect_lt(min(reference[,, 4]), 0.01)
  expect_equal(render(grouped), reference, tolerance = 1e-4)
  expect_equal(render(instanced), reference, tolerance = 1e-4)
})
