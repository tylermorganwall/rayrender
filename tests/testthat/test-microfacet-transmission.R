test_that("rough microfacet glass transmits positive radiance", {
  # An emitter enclosed by rough glass isolates transmission from the direct
  # lighting strategy. The signed cosine product used to make this black.
  scene = sphere(
    radius = 1,
    material = microfacet(transmission = TRUE, eta = 1.333, roughness = .1)
  ) |>
    add_object(sphere(radius = .7, material = light(intensity = 1))) |>
    add_object(sphere(
      x = 20,
      radius = .1,
      material = subsurface(sigma_s = 1, sigma_a = .1)
    ))
  set.seed(42)
  image = render_scene(
    scene,
    width = 24,
    height = 24,
    samples = 32,
    lookfrom = c(0, 0, 4),
    lookat = c(0, 0, 0),
    fov = 35,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    tonemap = "raw",
    min_variance = 0,
    parallel = FALSE,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE
  )
  expect_true(all(is.finite(image)))
  expect_gt(mean(image[9:16, 9:16, 1:3]), .2)
})
