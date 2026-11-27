#' @param scene Scene to render.
#' @param ... Overrides for the small deterministic render.
#' @return A raw image.
openpbr_test_render = function(scene, ...) {
  args = list(
    scene = scene,
    width = 8L,
    height = 8L,
    samples = 64L,
    sample_method = "random",
    lookfrom = c(0, 0, 4),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.5, .5),
    aperture = 0,
    ambient_light = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white",
    min_variance = 0,
    clamp_value = Inf,
    tonemap = "raw",
    denoise = FALSE,
    bloom = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE
  )
  set.seed(300926)
  do.call(render_scene, modifyList(args, list(...)))
}

test_that("OpenPBR descriptors validate the standard parameter groups", {
  mat = openpbr()
  expect_s3_class(mat, "ray_material")
  expect_identical(mat[[1]]$type, 11L)
  expect_identical(get_material_name(11L), "openpbr")
  expect_equal(mat[[1]]$openpbr$base_color, rep(.8, 3))
  expect_null(mat[[1]]$subsurface)
  expect_equal(openpbr(base_color = "red")[[1]]$openpbr$base_color, c(1, 0, 0))
  for (name in c(
    "base_weight",
    "base_metalness",
    "specular_roughness",
    "coat_weight",
    "thin_film_weight",
    "fuzz_roughness",
    "geometry_opacity"
  )) {
    for (value in list(
      -.1,
      1.1,
      NA_real_,
      Inf,
      numeric(),
      c(0, 1),
      matrix(.5),
      TRUE
    )) {
      expect_error(do.call(openpbr, setNames(list(value), name)), name)
    }
  }
  expect_error(openpbr(specular_ior = 0), "specular_ior")
  expect_error(openpbr(base_color = c(1, 2, 3)), "base_color")
  expect_error(openpbr(priority = .5), "priority")
  expect_error(openpbr(geometry_thin_walled = NA), "geometry_thin_walled")
  expect_error(openpbr(geometry_normal = c(0, 0, 0)), "geometry_normal")
  expect_s3_class(
    openpbr(subsurface_weight = 1, subsurface_radius = 0),
    "ray_material"
  )
  expect_error(
    openpbr(transmission_weight = 1, geometry_opacity = .5),
    "opacity"
  )
  expect_s3_class(
    openpbr(subsurface_weight = 1, subsurface_scatter_anisotropy = 1),
    "ray_material"
  )
  expect_error(openpbr(image_repeat = 0), "image_repeat")
})

test_that("OpenPBR owns and replaces solid interiors and selects NEE", {
  mat = openpbr(subsurface_weight = 1, priority = 2)
  ready = prepare_subsurface(sphere(material = mat))
  expect_equal(ready$shape_info[[1]]$medium$openpbr$priority, 2)
  expect_identical(prepare_subsurface(ready), ready)
  expect_null(
    prepare_subsurface(set_scene_material(ready, diffuse()))$shape_info[[
      1
    ]]$medium
  )
  expect_error(
    set_medium(sphere(material = mat), homogeneous_medium()),
    "conflicts"
  )
  expect_error(prepare_subsurface(xy_rect(material = mat)), "boundaries")
  thin = openpbr(subsurface_weight = 1, geometry_thin_walled = TRUE)
  expect_null(
    prepare_subsurface(xy_rect(material = thin))$shape_info[[1]]$medium
  )
  prep = prepare_scene_list(
    sphere(material = openpbr()),
    integrator_type = "basic"
  )
  expect_identical(prep$render_info$integrator_type, 1L)
})

test_that("opaque OpenPBR passes a white furnace and retains RGB", {
  white = openpbr_test_render(
    sphere(material = openpbr(base_color = "white")),
    samples = 256L
  )
  expect_true(all(is.finite(white)))
  expect_equal(mean(white[,, 1:3]), 1, tolerance = .025)
  red = openpbr_test_render(sphere(
    material = openpbr(base_color = c(.8, .1, .05))
  ))
  expect_gt(mean(red[,, 1]), 3 * mean(red[,, 2]))
  expect_null(attr(red, "path_warnings"))
})

test_that("OpenPBR transmission uses volume absorption and dielectric priorities", {
  glass = openpbr(
    transmission_weight = 1,
    transmission_color = c(.8, .5, .2),
    transmission_depth = 2,
    specular_roughness = 0,
    specular_ior = 1
  )
  image = openpbr_test_render(sphere(material = glass), samples = 128L)
  channels = apply(image[,, 1:3], 3, mean)
  expect_gt(channels[1], channels[2])
  expect_gt(channels[2], channels[3])
  expect_equal(channels, c(.8, .5, .2), tolerance = .035)
  expect_null(attr(image, "path_warnings"))
  black = sphere(
    material = openpbr(
      subsurface_weight = 1,
      subsurface_color = "black",
      subsurface_radius = .1,
      specular_ior = 1,
      priority = 1
    )
  )
  winner = sphere(
    radius = 1.2,
    material = dielectric(refraction = 1, priority = 0)
  )
  hidden = openpbr_test_render(add_object(winner, black), samples = 8L)
  expect_equal(mean(hidden[,, 1:3]), 1, tolerance = 1e-5)
  loser = set_scene_material(winner, dielectric(refraction = 1, priority = 2))
  visible = openpbr_test_render(add_object(loser, black), samples = 32L)
  expect_lt(mean(visible[,, 1:3]), .02)
  expect_null(attr(hidden, "path_warnings"))
  expect_null(attr(visible, "path_warnings"))
})

test_that("OpenPBR thin walls, opacity, emission and textures render", {
  thin = xy_rect(
    material = openpbr(
      geometry_thin_walled = TRUE,
      transmission_weight = 1,
      specular_roughness = 0,
      specular_ior = 1.5
    )
  )
  image = openpbr_test_render(thin, samples = 128L)
  expect_equal(mean(image[,, 1:3]), 1, tolerance = .035)
  transparent = xy_rect(
    material = openpbr(base_color = "black", geometry_opacity = 0)
  )
  expect_equal(
    mean(openpbr_test_render(transparent, samples = 8L)[,, 1:3]),
    1,
    tolerance = 1e-6
  )
  emission = openpbr_test_render(
    xy_rect(
      material = openpbr(
        emission_luminance = 2,
        emission_color = "white",
        base_color = "black"
      )
    ),
    samples = 8L,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  expect_equal(mean(emission[,, 1:3]), 2, tolerance = .01)
  map = array(.2, c(4, 4, 3))
  map[,, 1] = .9
  textured = sphere(
    material = openpbr(
      image_texture = map,
      roughness_texture = matrix(.5, 4, 4)
    )
  )
  image = openpbr_test_render(textured)
  expect_true(all(is.finite(image)))
  expect_gt(mean(image[,, 1]), mean(image[,, 2]))
  expect_null(attr(image, "path_warnings"))
})

test_that("OpenPBR instances and camera-inside rays use the volume integrator", {
  material = openpbr(transmission_weight = 1, specular_roughness = 0)
  copies = create_instances(sphere(material = material), x = 0)
  prepared = prepare_scene_list(copies, integrator_type = "rtiow")
  expect_identical(prepared$render_info$integrator_type, 1L)
  image = openpbr_test_render(copies, samples = 64L)
  expect_equal(mean(image[,, 1:3]), 1, tolerance = .035)
  inside = openpbr_test_render(
    sphere(material = material),
    lookfrom = c(0, 0, 0),
    lookat = c(0, 0, 1),
    fov = 10,
    samples = 64L
  )
  expect_equal(mean(inside[,, 1:3]), 1.5^2, tolerance = .05)
  expect_null(attr(inside, "path_warnings"))
  expect_null(attr(image, "path_warnings"))
})

test_that("OpenPBR area emission illuminates other surfaces through NEE", {
  scene = sphere(material = openpbr(base_color = "white")) |>
    add_object(sphere(
      x = 0,
      y = 3,
      z = 3,
      radius = .7,
      material = openpbr(base_color = "black", emission_luminance = 15)
    ))
  image = openpbr_test_render(
    scene,
    samples = 128L,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black"
  )
  expect_gt(mean(image[,, 1:3]), .1)
  expect_null(attr(image, "path_warnings"))
})

test_that("raw OpenPBR maps are independent of legacy texture remapping", {
  map = tempfile(fileext = ".png")
  on.exit(unlink(map))
  png::writePNG(matrix(c(.2, .4, .6, .8), 2), map)
  current = sphere(
    material = openpbr(base_color = c(.4, .2, .1), roughness_texture = map)
  )
  # An off-camera material loads the same file through the legacy remapping path.
  # Reordering scene rows must not mutate the OpenPBR material's authored map.
  legacy = sphere(x = 10000, material = microfacet(roughness_texture = map))
  forward = openpbr_test_render(add_object(current, legacy), samples = 32L)
  reverse = openpbr_test_render(add_object(legacy, current), samples = 32L)
  expect_equal(forward[,, 1:3], reverse[,, 1:3], tolerance = 1e-6)
})
