render_primitive_instances = function(
  scene,
  backend = "metal",
  ambient = TRUE
) {
  set.seed(415)
  image = render_scene(
    scene,
    integrator_type = backend,
    width = 40,
    height = 40,
    samples = 256,
    lookfrom = c(0, 1, 6),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(4, 4),
    min_variance = 0,
    max_depth = 12,
    ambient_light = ambient,
    backgroundlow = if (ambient) "white" else "black",
    backgroundhigh = if (ambient) "white" else "black",
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = "raw"
  )
  if (backend == "metal") {
    report = attr(image, "wavefront")
    if (!isTRUE(report$used)) {
      expect_match(report$fallback, "not enabled|device|macOS 11")
      skip("Metal unavailable")
    }
    expect_equal(report$discarded_paths, 0)
  }
  image
}

test_that("primitive instances share faces and preserve mirrored surface mappings", {
  paint = diffuse(
    texture_checker("gold", "steelblue", texture_coordinates("uv", scale = 4)),
    bump_texture = outer(
      seq(0, 6, length.out = 16),
      seq(0, 6, length.out = 16),
      function(x, y) .5 + .2 * sin(x) * cos(y)
    ),
    bump_intensity = .02
  )
  shapes = list(
    cylinder(radius = .6, length = 1.5, material = paint),
    disk(radius = .7, angle = c(90, 0, 0), material = paint),
    cube(material = paint),
    xy_rect(material = paint)
  )
  for (shape in shapes) {
    scene = create_instances(
      shape,
      x = c(-.9, .9),
      scale_x = c(-.7, .7),
      scale_y = c(.8, 1.2),
      scale_z = .6,
      angle_y = c(-15, 20)
    )
    gpu = render_primitive_instances(scene)
    cpu = render_primitive_instances(scene, "nee")
    expect_equal(attr(gpu, "wavefront")$unique_meshes, 1)
    expect_equal(attr(gpu, "wavefront")$mesh_instances, 2)
    expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .035)
    expect_true(all(is.finite(gpu)))
  }
})

test_that("instanced primitive emitters retain placement-specific area sampling", {
  scene = cylinder(
    radius = .3,
    length = 1.1,
    material = light(intensity = 8)
  ) |>
    create_instances(x = c(-1.3, 1.3), y = c(1.1, 1.5), scale_x = c(-1, .7)) |>
    add_object(sphere(radius = .6, material = diffuse("grey60"))) |>
    add_object(xz_rect(
      y = -.6,
      xwidth = 8,
      zwidth = 8,
      material = diffuse("white")
    ))
  gpu = render_primitive_instances(scene, ambient = FALSE)
  cpu = render_primitive_instances(scene, "nee", ambient = FALSE)
  expect_equal(attr(gpu, "wavefront")$unique_meshes, 1)
  expect_equal(attr(gpu, "wavefront")$mesh_instances, 2)
  expect_gt(mean(gpu[,, 1:3]), .02)
  expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .05)
})

test_that("instanced primitive glass retains separate closed boundaries", {
  glass = dielectric(refraction = 1.3, attenuation = c(.1, .2, .3))
  for (shape in list(
    cylinder(radius = .5, length = 1.3, material = glass),
    cube(material = glass)
  )) {
    scene = create_instances(shape, x = c(-.7, .7), scale_x = c(-1, 1))
    gpu = render_primitive_instances(scene)
    cpu = render_primitive_instances(scene, "nee")
    expect_equal(attr(gpu, "wavefront")$mesh_instances, 2)
    expect_lt(mean(abs(gpu[,, 1:3] - cpu[,, 1:3])), .035)
  }
})
