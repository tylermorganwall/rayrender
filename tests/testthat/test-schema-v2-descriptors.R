test_that("schema-v2 descriptors serialize deterministically", {
  descriptors = list(
    spectrum_constant(1),
    spectrum_rgb("#884422", role = "albedo"),
    spectrum_sampled(c(400, 500, 600), c(0.1, 0.5, 0.2), role = "unbounded"),
    spectrum_blackbody(6500),
    spectrum_named("stdillum-D65"),
    spectrum_cauchy_ior(1.5, B = 0.004),
    spectrum_sellmeier_ior(B = c(1, 0.5), C = c(0.01, 0.1)),
    texture_constant(c(0.2, 0.3, 0.4)),
    texture_image_color("albedo.png", role = "albedo"),
    texture_image_scalar("roughness.png"),
    texture_checker(0, 1),
    texture_mix(0, 1),
    texture_scale(texture_float(0.5), scale = 2),
    triplanar_mapping(scale = 2),
    object_mapping(),
    uniform_light_sampler(),
    path_integrator(max_depth = 8),
    volpath_integrator(max_depth = 8),
    random_walk_integrator(max_depth = 8),
    sobol_sampler(pixel_samples = 16),
    stratified_sampler(x_samples = 2, y_samples = 2),
    independent_sampler(pixel_samples = 16),
    perspective_camera(),
    orthographic_camera(scale = 2),
    realistic_camera(lens = "lens.dat"),
    cie1931_sensor(),
    rgb_sensor(),
    measured_sensor(
      spectrum_constant(1),
      spectrum_constant(1),
      spectrum_constant(1)
    ),
    rgb_film(width = 16, height = 16),
    spectral_options(),
    scene_validation("strict"),
    interface_material(),
    conductor(),
    dielectric_interface(),
    thin_dielectric(),
    coated_diffuse(reflectance = spectrum_constant(0.5)),
    coated_conductor(roughness = 0.1),
    diffuse_transmission(transmittance = spectrum_constant(0.5)),
    measured_material("material.bsdf"),
    mix_material(interface_material(), thin_dielectric()),
    henyey_greenstein_phase(0.1),
    homogeneous_medium(sigma_s = spectrum_constant(0.1)),
    optical_region("glass"),
    region_boundary("glass"),
    dielectric_region(id = "dielectric"),
    point_light(c(0, 1, 0)),
    spot_light(c(0, 1, 0), c(0, 0, 0)),
    distant_light(c(0, -1, 0)),
    uniform_infinite_light(),
    image_infinite_light("studio.hdr"),
    area_light(),
    named_texture("albedo"),
    named_material("surface")
  )

  for (descriptor in descriptors) {
    expect_identical(ray_schema_roundtrip(descriptor), descriptor)
  }
})

test_that("schema-v2 constructors reject invalid descriptors early", {
  expect_error(spectrum_constant(-1), "spectrum_constant\\(value\\)")
  expect_error(spectrum_sampled(c(500, 400), c(1, 1)), "strictly increasing")
  expect_error(texture_image("albedo.png"), "texture_image\\(role\\)")
  expect_error(point_light(c(0, 1)), "point_light\\(position\\)")
})

test_that("old scenes still construct through the legacy add_object path", {
  scene = add_object(sphere(), sphere(x = 1))

  expect_s3_class(scene, "ray_scene")
  expect_false(inherits(scene, "ray_scene_v2"))
  expect_equal(nrow(scene), 2)
})

test_that("schema-v2 scene decorators preserve rows and registries", {
  glass = optical_region("glass", eta = spectrum_constant(1.5))
  object = sphere(material = dielectric_interface()) |>
    with_light(area_light(scale = 4)) |>
    with_region_boundary(region_boundary(glass))

  expect_s3_class(object, "ray_scene_v2")
  expect_equal(nrow(object), 1)
  expect_s3_class(object$light[[1]], "ray_light")
  expect_equal(length(object$region_boundaries[[1]]), 1)

  scene = ray_scene() |>
    add_region(glass) |>
    add_named_texture("white", texture_constant(1)) |>
    add_named_material("interface", interface_material()) |>
    add_object(object) |>
    add_light(point_light(c(0, 4, 0))) |>
    set_environment(uniform_infinite_light())

  expect_s3_class(scene, "ray_scene_v2")
  expect_equal(nrow(scene), 1)
  expect_identical(names(attr(scene, "regions")), "glass")
  expect_identical(names(attr(scene, "named_textures")), "white")
  expect_identical(names(attr(scene, "named_materials")), "interface")
  expect_equal(length(attr(scene, "lights")), 1)
  expect_false(anyNA(scene$object_id))

  compiler_input = as_scene_compiler_input(
    scene,
    validation = scene_validation("strict")
  )
  expect_s3_class(compiler_input, "ray_scene_compiler_input")
  expect_equal(compiler_input$schema_version, 2)
  expect_equal(length(compiler_input$objects), 1)
  expect_s3_class(ray_schema_roundtrip(scene), "ray_scene_v2")
})

test_that("schema-v2 validation fails before rendering", {
  missing_region_scene = ray_scene_v2(
    sphere() |>
      with_region_boundary(region_boundary("missing"))
  )

  expect_error(
    as_scene_compiler_input(
      missing_region_scene,
      validation = scene_validation("strict")
    ),
    "object\\[1\\].region\\[missing\\]"
  )

  csg_light_scene = ray_scene_v2(
    csg_object(csg_sphere(), material = interface_material()) |>
      with_light(area_light())
  )

  expect_error(
    as_scene_compiler_input(
      csg_light_scene,
      validation = scene_validation("strict")
    ),
    "area light is not supported"
  )
})

test_that("legacy scene adaptation aggregates warnings once per scene", {
  legacy_scene = add_object(
    sphere(material = light()),
    sphere(x = 1, material = light())
  )

  expect_warning(
    {
      converted = legacy_scene_to_schema_v2(legacy_scene)
    },
    "legacy shorthand"
  )
  report = attr(
    suppressWarnings(legacy_scene_to_schema_v2(legacy_scene)),
    "conversion_report"
  )

  expect_s3_class(converted, "ray_scene_v2")
  expect_length(report$warnings, 1)
  expect_s3_class(converted$light[[1]], "ray_light")
  expect_s3_class(converted$material[[1]], "ray_material_v2")
})

test_that("render_scene exposes the schema-v2 spectral shell without rendering", {
  compiler_input = render_scene(
    sphere(),
    render_mode = "spectral",
    return_result = TRUE,
    width = 8,
    height = 8,
    samples = 4
  )

  expect_s3_class(compiler_input, "ray_scene_compiler_input")
  expect_s3_class(compiler_input$render_defaults$integrator, "ray_integrator")
  expect_s3_class(compiler_input$render_defaults$sampler, "ray_sampler")

  random_walk_input = render_scene(
    sphere(),
    render_mode = "spectral",
    integrator_type = "randomwalk",
    return_result = TRUE,
    width = 8,
    height = 8,
    samples = 4
  )

  expect_equal(random_walk_input$render_defaults$integrator$type, "random_walk")

  expect_error(
    render_scene(sphere(), render_mode = "spectral", return_result = FALSE),
    "validates schema-v2 input only"
  )
})
