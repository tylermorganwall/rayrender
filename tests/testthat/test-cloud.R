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
  expect_identical(
    prepare_scene_list(
      original,
      integrator_type = "basic"
    )$render_info$integrator_type,
    1L
  )
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
    "t",
    "animation_seed",
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
    animation_seed = -1,
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

test_that("cloud evolution is smooth, repeatable, and independent of placement", {
  skip_if_not_installed("ambient")
  for (style in c("cumulus", "stratus")) {
    initial = cloud(style = style, seed = 8, resolution = 32)
    at_zero = cloud(
      style = style,
      seed = 8,
      resolution = 32,
      t = 0,
      animation_seed = 101
    )
    expect_identical(at_zero, initial)
    middle = cloud(
      style = style,
      seed = 8,
      resolution = 32,
      t = 0.5,
      animation_seed = 101
    )
    next_frame = cloud(
      style = style,
      seed = 8,
      resolution = 32,
      t = 0.501,
      animation_seed = 101
    )
    later = cloud(
      style = style,
      seed = 8,
      resolution = 32,
      t = 1,
      animation_seed = 101
    )
    a = initial$shape_info[[1]]$medium$density
    b = middle$shape_info[[1]]$medium$density
    c = next_frame$shape_info[[1]]$medium$density
    d = later$shape_info[[1]]$medium$density
    expect_gt(mean(abs(b - a)), 1e-4)
    expect_gt(mean(abs(d - b)), 1e-4)
    expect_lt(mean(abs(c - b)), mean(abs(d - b)) / 20)
    near_zero = cloud(
      style = style,
      seed = 8,
      resolution = 32,
      t = 0.001,
      animation_seed = 101
    )$shape_info[[1]]$medium$density
    expect_lt(mean(abs(near_zero - a)), mean(abs(b - a)) / 20)

    # Revisiting a frame after later times does not depend on rendering order.
    expect_identical(
      cloud(
        style = style,
        seed = 8,
        resolution = 32,
        t = 0.5,
        animation_seed = 101
      ),
      middle
    )
    alternative = cloud(
      style = style,
      seed = 8,
      resolution = 32,
      t = 0.5,
      animation_seed = 102
    )
    expect_gt(mean(abs(alternative$shape_info[[1]]$medium$density - b)), 1e-4)
    moved = cloud(
      x = 20,
      y = 30,
      z = -40,
      angle = c(10, 30, 20),
      scale = 2,
      style = style,
      seed = 8,
      resolution = 32,
      t = 0.5,
      animation_seed = 101
    )
    expect_identical(
      moved$shape_info[[1]]$medium,
      middle$shape_info[[1]]$medium
    )
    expect_equal(c(middle$x, middle$y, middle$z), c(0, 0, 0))
    expect_identical(
      middle$shape_info[[1]]$medium$bounds,
      initial$shape_info[[1]]$medium$bounds
    )
    expect_identical(dim(b), dim(a))
    expect_true(all(is.finite(b) & b >= 0 & b <= 1))
    expect_true(all(b[c(1, dim(b)[1]), , ] == 0))
    expect_true(all(b[, c(1, dim(b)[2]), ] == 0))
    expect_true(all(b[,, c(1, dim(b)[3])] == 0))
  }
  # Broad features still evolve when fine detail is disabled.
  a = cloud(detail = 0, t = 0, resolution = 24)$shape_info[[1]]$medium$density
  b = cloud(detail = 0, t = 1, resolution = 24)$shape_info[[1]]$medium$density
  expect_gt(mean(abs(b - a)), 1e-4)
})

test_that("animated clouds preserve RNG state and support negative times and seed limits", {
  skip_if_not_installed("ambient")
  withr::local_seed(109)
  state = .Random.seed
  for (time in c(-1, 0.5, 1, 10)) {
    object = cloud(
      t = time,
      animation_seed = .Machine$integer.max - 1,
      resolution = 24
    )
    expect_identical(.Random.seed, state)
    expect_true(all(is.finite(object$shape_info[[1]]$medium$density)))
  }
  for (seed in c(0.5, -1, .Machine$integer.max)) {
    expect_error(
      cloud(t = 1, animation_seed = seed, resolution = 24),
      "animation_seed"
    )
  }
  rm(".Random.seed", envir = .GlobalEnv)
  invisible(cloud(t = 1, resolution = 24))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
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

test_that("cloud animation translates both Perlin fields without changing their seeds", {
  skip_if_not_installed("ambient")
  original = ambient::gen_perlin
  seen = list()
  local_mocked_bindings(
    gen_perlin = function(x, y, z, frequency, seed, ...) {
      seen[[length(seen) + 1L]] <<- list(
        coordinates = cbind(x, y, z),
        frequency = frequency,
        seed = seed
      )
      original(x, y, z, frequency = frequency, seed = seed, ...)
    },
    .package = "ambient"
  )
  invisible(cloud(resolution = 24, t = 0, animation_seed = 7))
  initial = seen
  seen = list()
  invisible(cloud(resolution = 24, t = 1, animation_seed = 7))
  expect_length(seen, length(initial))
  displacement = seen[[1]]$coordinates[1, ] - initial[[1]]$coordinates[1, ]
  expect_equal(sqrt(sum(displacement^2)), 0.1)
  for (i in seq_along(initial)) {
    expect_equal(seen[[i]]$seed, initial[[i]]$seed)
    expect_equal(seen[[i]]$frequency, initial[[i]]$frequency)
    expect_equal(
      seen[[i]]$coordinates,
      sweep(initial[[i]]$coordinates, 2, displacement, "+")
    )
  }
})
