test_that("schema-v2 descriptors serialize deterministically", {
  descriptors = list(
    spectrum_constant(1),
    spectrum_rgb("#884422", role = "albedo"),
    spectrum_rgb(c(0.2, 0.3, 0.4), role = "albedo", color_space = "DCI-P3"),
    spectrum_sampled(c(400, 500, 600), c(0.1, 0.5, 0.2), role = "unbounded"),
    spectrum_blackbody(6500),
    spectrum_named("stdillum-D65"),
    spectrum_cauchy_ior(1.5, B = 0.004),
    spectrum_sellmeier_ior(B = c(1, 0.5), C = c(0.01, 0.1)),
    texture_constant(c(0.2, 0.3, 0.4)),
    texture_image_color("albedo.png", role = "albedo"),
    texture_image_color(
      "rec2020-albedo.png",
      role = "albedo",
      color_space = "Rec.2020",
      encoding = "linear"
    ),
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
    rgb_sensor(color_space = "ACES2065-1"),
    measured_sensor(
      spectrum_constant(1),
      spectrum_constant(1),
      spectrum_constant(1)
    ),
    rgb_film(width = 16, height = 16),
    rgb_film(output_color_space = "Rec.2020"),
    spectral_options(),
    spectral_options(input_color_space = "DCI-P3"),
    scene_validation("strict"),
    interface_material(),
    conductor(),
    conductor(
      eta = spectrum_sampled(
        c(400, 500, 600),
        c(0.2, 0.4, 0.8),
        role = "unbounded"
      ),
      k = spectrum_sampled(c(400, 500, 600), c(2, 3, 4), role = "unbounded"),
      u_roughness = texture_float(0.1),
      v_roughness = texture_image_scalar("roughness.png"),
      remap_roughness = FALSE
    ),
    conductor(reflectance = spectrum_rgb("#ccaa55", role = "albedo")),
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
    image_infinite_light(
      "studio-aces.exr",
      color_space = "ACES2065-1",
      encoding = "linear"
    ),
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
  expect_error(
    spectrum_rgb(c(0.1, 0.2, 0.3), color_space = "AdobeRGB"),
    "must be one of"
  )
  expect_error(spectrum_sampled(c(500, 400), c(1, 1)), "strictly increasing")
  expect_error(texture_image("albedo.png"), "texture_image\\(role\\)")
  expect_error(
    conductor(
      reflectance = spectrum_constant(0.5),
      eta = spectrum_named("metal-Cu-eta")
    ),
    "cannot be combined"
  )
  expect_error(
    dielectric_region(sigma_a = spectrum_constant(0.1)),
    "nonzero absorption"
  )
  expect_error(point_light(c(0, 1)), "point_light\\(position\\)")
})

test_that("schema-v2 color spaces are explicit canonical names", {
  expect_identical(
    ray_supported_color_spaces(),
    c("sRGB", "DCI-P3", "Rec.2020", "ACES2065-1")
  )
  expect_equal(
    spectrum_rgb(c(0.1, 0.2, 0.3), color_space = "DCI-P3")$color_space,
    "DCI-P3"
  )
  expect_equal(
    rgb_film(output_color_space = "ACES2065-1")$output_color_space,
    "ACES2065-1"
  )
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
  expect_match(compiler_input$compiler_hash, "^[a-f0-9]{32}$")
  expect_gt(length(compiler_input$caches$materials), 0)
  expect_gt(length(compiler_input$caches$shapes), 0)
  expect_s3_class(ray_schema_roundtrip(scene), "ray_scene_v2")
})

test_that("compiler output caches and hashes are deterministic", {
  scene = ray_scene_v2(sphere(material = interface_material())) |>
    add_light(point_light(c(0, 4, 0))) |>
    set_environment(uniform_infinite_light())

  compiler_input1 = as_scene_compiler_input(
    scene,
    validation = scene_validation("strict")
  )
  compiler_input2 = as_scene_compiler_input(
    scene,
    validation = scene_validation("strict")
  )

  expect_identical(compiler_input1$compiler_hash, compiler_input2$compiler_hash)
  expect_match(compiler_input1$compiler_hash, "^[a-f0-9]{32}$")
  expect_gt(length(compiler_input1$caches$materials), 0)
  expect_gt(length(compiler_input1$caches$lights), 1)
  expect_gt(length(compiler_input1$caches$shapes), 0)
  expect_true(compiler_input1$compiled$objects[[1]]$normal_transform$finite)
})

test_that("compiler rejects legacy positional material data", {
  expect_error(
    as_scene_compiler_input(
      ray_scene_v2(sphere()),
      validation = scene_validation("none")
    ),
    "legacy positional material data"
  )
})

test_that("importers expose color policy and generated capabilities", {
  obj_file = r_obj()
  skip_if_not(file.exists(obj_file))

  imported = obj_model(
    obj_file,
    load_material = TRUE,
    load_textures = FALSE,
    load_normals = FALSE,
    calculate_consistent_normals = FALSE,
    material = interface_material()
  )
  compiler_input = as_scene_compiler_input(
    ray_scene_v2(imported),
    validation = scene_validation("strict")
  )
  caps = compiler_input$compiled$objects[[1]]$shape_capabilities
  policy = compiler_input$diagnostics$importers[[1]]$color_policy

  expect_false(caps$uv)
  expect_true(caps$generated_mapping)
  expect_true(caps$generated_normals)
  expect_equal(policy$base_color$role, "albedo")
  expect_equal(policy$emission$role, "illuminant")
  expect_equal(policy$scalar_maps$encoding, "linear")
  expect_false(policy$scalar_maps$gamma_decoded)
})

test_that("shared materials do not collapse region instances", {
  glass = optical_region("glass", eta = spectrum_constant(1.5))
  material = dielectric_interface()
  object1 = sphere(material = material) |>
    with_region_boundary(region_boundary(glass))
  object2 = sphere(x = 1, material = material) |>
    with_region_boundary(region_boundary(glass))
  compiler_input = ray_scene() |>
    add_region(glass) |>
    add_object(object1) |>
    add_object(object2) |>
    as_scene_compiler_input(validation = scene_validation("strict"))

  ids = vapply(
    compiler_input$compiled$objects,
    function(object) object$region_instances[[1]]$instance_id,
    character(1)
  )
  material_ids = vapply(
    compiler_input$compiled$objects,
    function(object) object$material_cache_id,
    character(1)
  )

  expect_length(unique(ids), 2)
  expect_identical(length(unique(material_ids)), 1L)
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
    "CSG area lights require sampling"
  )

  csg_textured_scene = ray_scene_v2(
    csg_object(
      csg_sphere(),
      material = coated_diffuse(
        reflectance = texture_image_color("albedo.png", role = "albedo")
      )
    )
  )

  expect_error(
    as_scene_compiler_input(
      csg_textured_scene,
      validation = scene_validation("strict")
    ),
    "CSG image textures require"
  )
})

test_that("sampled CSG area lights compile through render-mesh path", {
  csg_light_scene = ray_scene_v2(
    csg_object(csg_sphere(), material = interface_material()) |>
      with_light(area_light(sampling = "sampled"))
  )

  compiler_input = as_scene_compiler_input(
    csg_light_scene,
    validation = scene_validation("strict")
  )
  object = compiler_input$compiled$objects[[1]]

  expect_true(object$requires_mesh_area_light)
  expect_equal(object$csg_compilation$area_light, "render_mesh")
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

test_that("legacy metal adapts to one spectral compatibility warning", {
  legacy_scene = add_object(
    sphere(material = metal(color = "gold", fuzz = 0.2)),
    sphere(x = 1, material = metal())
  )

  expect_warning(
    {
      converted = legacy_scene_to_schema_v2(legacy_scene)
    },
    "compatibility conductor"
  )
  report = attr(
    suppressWarnings(legacy_scene_to_schema_v2(legacy_scene)),
    "conversion_report"
  )

  expect_length(report$warnings, 1)
  expect_equal(converted$material[[1]][[1]]$type, "compat_rgb_metal_conductor")
  expect_s3_class(
    converted$material[[1]][[1]]$params$reflectance,
    "ray_spectrum"
  )
  expect_equal(converted$material[[1]][[1]]$params$u_roughness$value, 0.2)
})

test_that("legacy glossy adapts to coated diffuse", {
  legacy_scene = add_object(
    sphere(
      material = glossy(color = "#336699", gloss = 0.5, reflectance = 0.04)
    )
  )

  expect_warning(
    {
      converted = legacy_scene_to_schema_v2(legacy_scene)
    },
    "coated diffuse"
  )
  report = attr(
    suppressWarnings(legacy_scene_to_schema_v2(legacy_scene)),
    "conversion_report"
  )

  expect_length(report$warnings, 1)
  expect_equal(converted$material[[1]][[1]]$type, "coated_diffuse")
  expect_s3_class(
    converted$material[[1]][[1]]$params$reflectance,
    "ray_spectrum"
  )
  expect_equal(converted$material[[1]][[1]]$params$u_roughness$value, 0.0625)
  expect_equal(converted$material[[1]][[1]]$params$v_roughness$value, 0.0625)
  expect_equal(converted$material[[1]][[1]]$params$eta$value, 1.5)
  expect_false(converted$material[[1]][[1]]$params$remap_roughness)
})

test_that("legacy dielectric adapts to optical region metadata", {
  legacy_scene = add_object(
    sphere(material = dielectric(refraction = 1.33, priority = -2))
  )

  expect_warning(
    {
      converted = legacy_scene_to_schema_v2(legacy_scene)
    },
    "dielectric\\(\\)"
  )
  report = attr(
    suppressWarnings(legacy_scene_to_schema_v2(legacy_scene)),
    "conversion_report"
  )
  regions = attr(converted, "regions")

  expect_length(report$warnings, 1)
  expect_equal(converted$material[[1]][[1]]$type, "dielectric_interface")
  expect_length(regions, 1)
  expect_equal(regions[[1]]$eta$value, 1.33)
  expect_equal(regions[[1]]$priority, -2L)
  expect_equal(
    converted$region_boundaries[[1]][[1]]$region,
    names(regions)[[1]]
  )
  expect_equal(converted$region_boundaries[[1]][[1]]$side, "negative_normal")
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
  expect_equal(
    compiler_input$render_defaults$spectral_policy$default_render_mode,
    "rgb_legacy"
  )
  expect_false(
    compiler_input$render_defaults$spectral_policy$runtime_rendering$still
  )
  expect_equal(
    compiler_input$render_defaults$runtime$runtime_bridge,
    "not_connected"
  )

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
    "spectral runtime renderer is not connected"
  )
})

test_that("spectral rollout policy documents unsupported runtime features", {
  capabilities = spectral_render_capabilities()

  expect_equal(capabilities$default_render_mode, "rgb_legacy")
  expect_equal(capabilities$default_policy$decision, "keep_rgb_legacy_default")
  expect_false(capabilities$compiler_input$direct_runtime_bridge)
  expect_false(capabilities$runtime_rendering$denoising)
  expect_equal(capabilities$denoiser$input_space, "output-linear RGB")
  expect_false(capabilities$adaptive_sampling$uses_packet_components)
  expect_false(capabilities$alpha$transmissive_dielectric_is_alpha)
})

test_that("spectral animation frame zero matches still compiler camera setup", {
  camera_motion = generate_camera_motion(
    positions = list(c(0, 1, -10), c(0, 1, -10)),
    lookats = list(c(0, 0, 0), c(0, 0, 0)),
    frames = 2,
    type = "linear",
    progress = FALSE
  )
  frame = list(
    lookfrom = c(
      camera_motion$x[[1]],
      camera_motion$y[[1]],
      camera_motion$z[[1]]
    ),
    lookat = c(
      camera_motion$dx[[1]],
      camera_motion$dy[[1]],
      camera_motion$dz[[1]]
    ),
    camera_up = c(
      camera_motion$upx[[1]],
      camera_motion$upy[[1]],
      camera_motion$upz[[1]]
    ),
    aperture = camera_motion$aperture[[1]],
    fov = camera_motion$fov[[1]],
    focal_distance = camera_motion$focal[[1]],
    ortho_dimensions = c(camera_motion$orthox[[1]], camera_motion$orthoy[[1]])
  )
  still_input = render_scene(
    sphere(),
    render_mode = "spectral",
    return_result = TRUE,
    width = 8,
    height = 8,
    samples = 4,
    fov = frame$fov,
    lookfrom = frame$lookfrom,
    lookat = frame$lookat,
    camera_up = frame$camera_up,
    aperture = frame$aperture,
    focal_distance = frame$focal_distance,
    ortho_dimensions = frame$ortho_dimensions,
    progress = FALSE
  )
  animation_input = render_animation(
    sphere(),
    camera_motion = camera_motion,
    render_mode = "spectral",
    return_result = TRUE,
    start_frame = 1,
    end_frame = 1,
    width = 8,
    height = 8,
    samples = 4,
    progress = FALSE
  )

  expect_equal(
    animation_input$render_defaults$camera,
    still_input$render_defaults$camera
  )
  expect_equal(animation_input$render_defaults$frame_index, 1)
  expect_equal(
    animation_input$render_defaults$runtime$request$context,
    "render_animation"
  )
  expect_error(
    render_animation(
      sphere(),
      camera_motion = camera_motion,
      render_mode = "spectral",
      return_result = FALSE,
      progress = FALSE
    ),
    "spectral runtime renderer is not connected"
  )
})
