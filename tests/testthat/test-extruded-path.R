sweep_test_line = function() rbind(c(0, 0, 0), c(0, 0, 10))
sweep_test_loop = function() {
  rbind(c(0, 0, 0), c(2, 0, 0), c(2, 2, 1), c(0, 2, 0), c(-1, 1, 2))
}

test_that("rayvertex sweep meshes pass through unchanged and empty sweeps stay empty", {
  args = list(
    points = sweep_test_line(),
    breaks = 6,
    smooth_normals = TRUE,
    width = 0,
    width_end = 1,
    twists = .25
  )
  expected = do.call(rayvertex::extruded_path_mesh, args)
  scene = do.call(extruded_path, args)
  expect_identical(scene$shape_info[[1]]$mesh_info[[1]], expected)
  expect_null(extruded_path(sweep_test_line(), width = 0, breaks = 5))
  expect_null(extruded_path(sweep_test_line(), u_min = .3, u_max = .3))
  expect_null(extruded_path(
    sweep_test_line(),
    width = 0,
    material_caps = diffuse("red")
  ))
})

test_that("wrapped intervals honor cap choices and share all material IDs", {
  for (caps in list(
    c(FALSE, FALSE),
    c(TRUE, FALSE),
    c(FALSE, TRUE),
    c(TRUE, TRUE)
  )) {
    object = extruded_path(
      sweep_test_loop(),
      closed = TRUE,
      breaks = 20,
      u_min = .8,
      u_max = 1.2,
      end_caps = caps
    )
    expect_equal(nrow(object), 1L)
    expect_length(object$shape_info[[1]]$mesh_info[[1]]$shapes, 2L + sum(caps))
    ids = vapply(object$shape_info, function(s) s$material_id, 0.0)
    expect_length(unique(ids), 1)
  }
  object = extruded_path(
    sweep_test_loop(),
    closed = TRUE,
    breaks = 20,
    u_min = .8,
    u_max = 1.2,
    material_caps = diffuse("red")
  )
  ids = vapply(object$shape_info, function(s) s$material_id, 0.0)
  expect_length(ids, 2)
  expect_false(ids[1] == ids[2])
  expect_length(object$shape_info[[1]]$mesh_info[[1]]$shapes, 2)
  expect_length(object$shape_info[[2]]$mesh_info[[1]]$shapes, 2)
})

test_that("rayvertex sweeps preserve renderer materials and transforms", {
  mat = diffuse("steelblue")
  scene = extruded_path(
    sweep_test_line(),
    breaks = 5,
    material = mat,
    x = 1,
    angle = c(0, 30, 0),
    scale = 2
  )
  expect_true(all(scene$x == 1))
  expect_equal(scene$material[[1]], mat[[1]])
  set.seed(1)
  image = render_scene(
    scene,
    width = 16,
    height = 16,
    samples = 2,
    lookfrom = c(4, 4, 25),
    lookat = c(1, 0, 10),
    fov = 60,
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_gt(diff(range(image)), 0)
})

test_that("sweep surfaces and caps form one absorbing SSS boundary", {
  scene = extruded_path(
    sweep_test_line(),
    breaks = 5,
    material = subsurface(sigma_a = .1, sigma_s = 0, refraction = 1)
  )
  set.seed(42)
  image = render_scene(
    scene,
    width = 4,
    height = 4,
    samples = 2048,
    sample_method = "random",
    lookfrom = c(0, 0, -2),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.25, .25),
    aperture = 0,
    ambient_light = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white",
    min_variance = 0,
    tonemap = "raw",
    denoise = FALSE,
    bloom = FALSE,
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_equal(mean(image[,, 1:3]), exp(-1), tolerance = .015)
})

test_that("smooth zero-width tips render with sparse normal indices", {
  scene = extruded_path(
    sweep_test_line(),
    width = 0,
    width_end = 2,
    breaks = 8,
    smooth_normals = TRUE
  )
  mesh = scene$shape_info[[1]]$mesh_info[[1]]
  expect_lt(
    length(unique(as.vector(mesh$shapes[[1]]$norm_indices))),
    nrow(mesh$normals[[1]])
  )
  set.seed(2)
  image = render_scene(
    scene,
    width = 8,
    height = 8,
    samples = 2,
    lookfrom = c(4, 2, 12),
    lookat = c(0, 0, 6),
    fov = 60,
    debug_channel = "normals",
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot = FALSE,
    denoise = FALSE,
    bloom = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_gt(diff(range(image)), 0)
})
