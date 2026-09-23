test_that("rayvertex lathes form absorbing SSS boundaries in rayrender", {
  mesh = rayvertex::lathe_mesh(
    rbind(c(0, 0), c(1, 0), c(1, 2), c(0, 2)),
    segments = 32
  )
  scene = raymesh_model(
    mesh,
    material = subsurface(sigma_a = .1, sigma_s = 0, refraction = 1)
  )
  expect_identical(scene$shape_info[[1]]$mesh_info[[1]], mesh)
  set.seed(42)
  image = render_scene(
    scene,
    width = 4,
    height = 4,
    samples = 2048,
    sample_method = "random",
    lookfrom = c(0, 1, -3),
    lookat = c(0, 1, 0),
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
    plot_scene = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_equal(mean(image[,, 1:3]), exp(-.2), tolerance = .015)
})
