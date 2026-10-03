#' @param scene Scene to render.
#' @param samples Default `64`. Samples per pixel.
#' @param ambient Default `FALSE`. Enable the uniform environment.
#' @param depth Default `2`. Maximum scattering depth.
#' @return Mean linear RGB at a nearly point-sized orthographic film.
#' @keywords internal
matched_glass_pixel = function(
  scene,
  samples = 64,
  ambient = FALSE,
  depth = 2
) {
  set.seed(721)
  image = render_scene(
    scene,
    width = 8,
    height = 8,
    samples = samples,
    lookfrom = c(0, 0, 5),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(1e-4, 1e-4),
    aperture = 0,
    ambient_light = ambient,
    backgroundhigh = c(.5, .5, .5),
    backgroundlow = c(.5, .5, .5),
    integrator_type = "nee",
    parallel = FALSE,
    progress = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = "raw",
    min_variance = 0,
    max_depth = depth
  )
  apply(image[,, 1:3, drop = FALSE], 3, mean)
}

test_that("matched glass transmits point and spot light with analytic RGB absorption", {
  skip_on_cran()
  wall = xy_rect(xwidth = 10, ywidth = 10, material = diffuse("white"))
  sigma = c(.2, .5, 1)
  tint = c(.7, .8, .9)
  pane = cube(
    x = 1,
    z = 1,
    width = .6,
    material = dielectric(color = tint, refraction = 1, attenuation = sigma)
  )
  # Only the light connection crosses the pane; the camera looks along x = 0.
  for (lamp in list(
    point_light(position = c(2, 0, 2), intensity = 8),
    spot_light(position = c(2, 0, 2), direction = c(-1, 0, -1), intensity = 8)
  )) {
    scene = add_light(wall, lamp)
    reference = matched_glass_pixel(scene, samples = 8)
    expect_equal(reference, rep(sqrt(.5) / pi, 3), tolerance = 2e-4)
    expect_equal(
      matched_glass_pixel(add_object(scene, pane), samples = 8) / reference,
      tint * exp(-sigma * .6 * sqrt(2)),
      tolerance = 3e-4
    )
    # No approximate equality: genuinely refracting glass still blocks this
    # direct connection, even when the IOR difference is small.
    refracting = cube(
      x = 1,
      z = 1,
      width = .6,
      material = dielectric(refraction = 1.0001)
    )
    expect_equal(
      matched_glass_pixel(add_object(scene, refracting), samples = 8),
      rep(0, 3),
      tolerance = 1e-8
    )
  }
})

test_that("matched shadow connections use priority-selected absorption", {
  skip_on_cran()
  scene = add_light(
    xy_rect(xwidth = 10, ywidth = 10, material = diffuse("white")),
    point_light(position = c(2, 0, 2), intensity = 8)
  )
  reference = matched_glass_pixel(scene, samples = 8)
  outer_sigma = c(.1, .3, .5)
  inner_sigma = c(.8, .5, .2)
  outer = cube(
    x = 1,
    z = 1,
    width = .8,
    material = dielectric(
      refraction = 1,
      attenuation = outer_sigma,
      priority = 2
    )
  )
  inner = cube(
    x = 1,
    z = 1,
    width = .4,
    material = dielectric(
      refraction = 1,
      attenuation = inner_sigma,
      priority = 1
    )
  )
  expect_equal(
    matched_glass_pixel(
      add_object(add_object(scene, outer), inner),
      samples = 8
    ) /
      reference,
    exp(-.4 * sqrt(2) * (outer_sigma + inner_sigma)),
    tolerance = 3e-4
  )
  # A losing surface can have another IOR and black tint; neither is visible
  # while the outer material has higher priority.
  hidden = cube(
    x = 1,
    z = 1,
    width = .4,
    material = dielectric(
      color = "black",
      refraction = 1.5,
      attenuation = rep(100, 3),
      priority = 3
    )
  )
  expect_equal(
    matched_glass_pixel(
      add_object(add_object(scene, outer), hidden),
      samples = 8
    ) /
      reference,
    exp(-.8 * sqrt(2) * outer_sigma),
    tolerance = 3e-4
  )
})

test_that("the matched IOR is relative to the surrounding glass", {
  skip_on_cran()
  scene = add_light(
    xy_rect(xwidth = 10, ywidth = 10, material = diffuse("white")),
    point_light(position = c(2, 0, 2), intensity = 8)
  ) |>
    add_object(cube(
      width = 20,
      material = dielectric(refraction = 1.5, priority = 10)
    ))
  # Camera, receiver, and light are all inside the clear outer solid.
  reference = matched_glass_pixel(scene, samples = 8)
  expect_equal(reference, rep(sqrt(.5) / pi, 3), tolerance = 2e-4)
  pane = cube(
    x = 1,
    z = 1,
    width = .6,
    material = dielectric(refraction = 1.5, attenuation = c(.2, .5, 1))
  )
  expect_equal(
    matched_glass_pixel(add_object(scene, pane), samples = 8) / reference,
    exp(-c(.2, .5, 1) * .6 * sqrt(2)),
    tolerance = 3e-4
  )
  air_pane = cube(
    x = 1,
    z = 1,
    width = .6,
    material = dielectric(refraction = 1)
  )
  expect_equal(
    matched_glass_pixel(add_object(scene, air_pane), samples = 8),
    rep(0, 3),
    tolerance = 1e-8
  )
})

test_that("matched panes preserve area and environment MIS and scattering depth", {
  skip_on_cran()
  wall = xy_rect(xwidth = 10, ywidth = 10, material = diffuse("white"))
  area_scene = add_object(
    wall,
    sphere(x = 2, z = 2, radius = .5, material = light(intensity = 4))
  )
  # A broad slab crosses both primary and direct-light/BSDF continuation rays.
  # Zero absorption should leave illumination unchanged, not double it through
  # an erroneous specular flag at the emitter, or consume the scattering budget.
  pane = cube(
    z = 1,
    xwidth = 20,
    ywidth = 20,
    zwidth = .3,
    material = dielectric(refraction = 1)
  )
  for (ambient in c(FALSE, TRUE)) {
    scene = if (ambient) wall else area_scene
    reference = matched_glass_pixel(scene, ambient = ambient, samples = 256)
    through = matched_glass_pixel(
      add_object(scene, pane),
      ambient = ambient,
      samples = 256
    )
    expect_gt(min(reference), .01)
    expect_equal(through, reference, tolerance = .015)
    shallow = matched_glass_pixel(scene, ambient = ambient, depth = 1)
    expect_equal(
      matched_glass_pixel(
        add_object(scene, pane),
        ambient = ambient,
        depth = 1
      ),
      shallow,
      tolerance = .015
    )
  }
})
