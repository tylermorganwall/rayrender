#' @param text PBRT text to write.
#' @return Temporary PBRT filename.
#' @keywords internal
pbrt_test_file = function(text) {
  file = tempfile(fileext = ".pbrt")
  writeLines(text, file)
  file
}

test_that("PBRT execution options accept a name and an untyped value", {
  file = pbrt_test_file(c(
    'Option "wavefront" true',
    'Option "disable-pixel-jitter" false',
    'WorldBegin Shape "sphere"'
  ))
  commands = pbrt_parse_file(file)
  expect_identical(commands[[1]]$args, c("wavefront", "true"))
  expect_identical(commands[[2]]$args, c("disable-pixel-jitter", "false"))
  imported = suppressWarnings(read_pbrt(file, strict = FALSE))
  expect_equal(nrow(imported$scene), 1)
  expect_equal(sum(imported$diagnostics$directive == "Option"), 2)
  expect_error(
    pbrt_parse_file(pbrt_test_file('Option "wavefront"')),
    "Expected 2 operands"
  )
})

test_that("PBRT parameters can continue across nested Include boundaries", {
  directory = tempfile("pbrt-continuation-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  root = file.path(directory, "scene.pbrt")
  writeLines(
    c(
      'WorldBegin Include "nested.pbrt" "float radius" 2',
      'Include "nested.pbrt" "float radius" 3'
    ),
    root
  )
  writeLines('Include "sphere.pbrt"', file.path(directory, "nested.pbrt"))
  writeLines('Shape "sphere"', file.path(directory, "sphere.pbrt"))
  imported = read_pbrt(root)
  radii = vapply(
    imported$scene$shape_info,
    function(shape) shape$shape_properties$radius,
    numeric(1)
  )
  expect_equal(radii, c(2, 3))
  expect_length(imported$source_files, 3)
  writeLines(
    'Shape "sphere" "float radius" 1',
    file.path(directory, "sphere.pbrt")
  )
  expect_error(read_pbrt(root), "Duplicate parameter: radius")
  writeLines('Translate 0 1 0', file.path(directory, "sphere.pbrt"))
  expect_error(read_pbrt(root), "does not end in a parameterized directive")
})

test_that("PBRT singleton environment axes remain safe to export as EXR", {
  context = new.env(parent = emptyenv())
  context$asset_dir = tempfile("pbrt-single-pixel-")
  context$assets = character()
  on.exit(unlink(context$asset_dir, recursive = TRUE))
  for (size in list(c(1, 1), c(1, 5), c(5, 1))) {
    image = array(rep(c(.2, .4, .6), each = prod(size)), c(size, 3))
    file = pbrt_write_asset(image, ".exr", context)
    output = rayimage::ray_read_image(file)
    expect_equal(dim(output)[1:2], pmax(size, 2))
    for (channel in 1:3) {
      expect_equal(
        as.numeric(output[,, channel]),
        rep(c(.2, .4, .6)[channel], prod(pmax(size, 2))),
        tolerance = 1e-3
      )
    }
  }
})

test_that("PBRT permissive previews report unresolved and repeated names", {
  file = pbrt_test_file(c(
    'WorldBegin NamedMaterial "later" Shape "sphere"',
    'Texture "color" "spectrum" "constant" "rgb value" [.2 .3 .4]',
    'Material "diffuse" "texture reflectance" "color" Shape "sphere"',
    'Texture "color" "spectrum" "constant" "rgb value" [.7 .6 .5]',
    'Material "diffuse" "texture reflectance" "color" Shape "sphere"'
  ))
  expect_error(read_pbrt(file), "Undefined material")
  result = suppressWarnings(read_pbrt(file, strict = FALSE))
  colors = lapply(result$scene$material, function(material) {
    material$properties[[1]][1:3]
  })
  expect_equal(colors, list(rep(.5, 3), c(.2, .3, .4), c(.7, .6, .5)))
  expect_true(any(grepl("Undefined material", result$diagnostics$message)))
  expect_true(any(grepl("Duplicate named texture", result$diagnostics$message)))
  duplicate = pbrt_test_file(c(
    'WorldBegin Texture "x" "float" "constant" "float value" 1',
    'Texture "x" "float" "constant" "float value" 2'
  ))
  expect_error(read_pbrt(duplicate), "Duplicate named texture")
})

test_that("PBRT identical repeated parameters are harmless", {
  result = read_pbrt(pbrt_test_file(
    'WorldBegin Shape "sphere" "float radius" 2 "float radius" 2'
  ))
  expect_equal(result$scene$shape_info[[1]]$shape_properties$radius, 2)
  expect_error(
    read_pbrt(pbrt_test_file(
      'WorldBegin Shape "sphere" "float radius" 2 "float radius" 3'
    )),
    "Duplicate parameter"
  )
})

test_that("PBRT point and spot lights retain world transforms and cone parameters", {
  imported = read_pbrt(pbrt_test_file(c(
    'WorldBegin',
    'Translate 1 2 3',
    'LightSource "point" "point3 from" [1 0 0] "rgb I" [2 4 6] "float scale" 2',
    'LightSource "spot" "point3 from" [0 0 0] "point3 to" [0 1 0] "float coneangle" 40 "float conedeltaangle" 10',
    'Shape "sphere"'
  )))
  lights = list_lights(imported$scene)
  expect_equal(nrow(imported$scene), 1)
  expect_length(lights, 2)
  expect_equal(lights[[1]]$position, c(2, 2, 3))
  expect_equal(lights[[1]]$color * lights[[1]]$intensity, c(4, 8, 12))
  expect_equal(lights[[2]]$direction, c(0, 1, 0))
  expect_equal(lights[[2]]$cone_angle, 40)
  expect_equal(lights[[2]]$falloff_angle, 10)
  expect_equal(nrow(imported$diagnostics), 0)
})

test_that("PBRT coated conductors and hair use existing material models", {
  imported = suppressWarnings(read_pbrt(
    pbrt_test_file(c(
      'WorldBegin',
      'Material "coatedconductor" "rgb conductor.reflectance" [.8 .5 .2] "float interface.eta" 1.4',
      '"float conductor.roughness" .0625 "float interface.roughness" .0016 Shape "sphere"',
      'Material "hair" "rgb sigma_a" [.1 .2 .3] "float beta_m" .4 Shape "sphere"',
      'Material "hair" "rgb color" [.2 .3 .4] Shape "sphere"'
    )),
    strict = FALSE
  ))
  coat = imported$scene$material[[1]]$openpbr
  expect_equal(coat$base_metalness, 1)
  expect_equal(coat$base_color, c(.8, .5, .2))
  expect_equal(coat$coat_ior, 1.4)
  expect_equal(coat$specular_roughness, .5)
  expect_equal(coat$coat_roughness, .2)
  expect_equal(imported$scene$material[[2]]$properties[[1]][1:3], c(.1, .2, .3))
  expect_equal(imported$scene$material[[3]], hair(color = c(.2, .3, .4))[[1]])
  expect_equal(hair(sigma_a = .2)[[1]]$properties[[1]][1:3], rep(.2, 3))
  expect_error(hair(sigma_a = c(1, NA, 2)), "finite")
})

test_that("PBRT Disney parameters produce OpenPBR lobes and preserve GGX widths", {
  file = pbrt_test_file(c(
    'WorldBegin',
    'MakeNamedMaterial "default" "string type" "disney" NamedMaterial "default" Shape "sphere"',
    'Material "disney" "rgb color" [.25 .5 .75] "float metallic" .2',
    '"float eta" 1.4 "float roughness" .3 "float anisotropic" .7',
    '"float speculartint" .6 "float sheen" .4 "float sheentint" .8',
    '"float clearcoat" .9 "float clearcoatgloss" .2 "float spectrans" .35 Shape "sphere"',
    'Material "disney" "rgb color" [0 0 0] "float roughness" 0',
    '"float speculartint" 1 "float sheen" 1 "float sheentint" 1 Shape "sphere"'
  ))
  expect_error(read_pbrt(file), "Disney material converted to OpenPBR")
  imported = suppressWarnings(read_pbrt(file, strict = FALSE))
  expect_true(all(vapply(
    imported$scene$material,
    function(x) x$type == 11L,
    logical(1)
  )))
  defaults = imported$scene$material[[1]]$openpbr
  expect_equal(defaults$base_color, rep(.5, 3))
  expect_equal(defaults$specular_roughness, .5)
  expect_equal(defaults$specular_ior, 1.5)
  expect_equal(
    defaults$base_metalness +
      defaults$coat_weight +
      defaults$fuzz_weight +
      defaults$subsurface_weight +
      defaults$transmission_weight,
    0
  )
  expect_null(imported$scene$material[[1]]$subsurface)
  mapped = imported$scene$material[[2]]$openpbr
  expect_equal(mapped$base_color, c(.25, .5, .75))
  expect_equal(mapped$base_metalness, .2)
  expect_equal(mapped$base_diffuse_roughness, .3)
  expect_equal(mapped$specular_ior, 1.4)
  # Check the resulting GGX widths, rather than only restating a parameter map.
  alpha_x = mapped$specular_roughness^2 *
    sqrt(2 / (1 + (1 - mapped$specular_roughness_anisotropy)^2))
  alpha_y = (1 - mapped$specular_roughness_anisotropy) * alpha_x
  expect_equal(c(alpha_x, alpha_y), .3^2 * c(1 / sqrt(.37), sqrt(.37)))
  expect_equal(mapped$coat_weight, .225)
  expect_equal(mapped$coat_ior, 1.5)
  expect_equal(mapped$coat_roughness^2, .0802)
  expect_equal(mapped$transmission_weight, .35)
  expect_equal(mapped$transmission_color^2, c(.25, .5, .75))
  expect_equal(mapped$transmission_depth, 0)
  expect_equal(mapped$fuzz_weight, .4 * .8 * .65)
  expect_lt(mapped$specular_color[1], mapped$specular_color[3])
  expect_true(all(mapped$specular_color <= 1))
  black = imported$scene$material[[3]]$openpbr
  expect_equal(black$specular_color, rep(1, 3))
  expect_equal(black$fuzz_color, rep(1, 3))
  expect_equal(black$specular_roughness^2, .001)
  expect_false(any(grepl(
    "Unsupported material|Unsupported parameter",
    imported$diagnostics$message
  )))
})

test_that("Disney specular tint is inactive at both metalness endpoints", {
  for (metallic in c(0, 1)) {
    file = pbrt_test_file(c(
      'WorldBegin Material "disney" "rgb color" [.9 .1 .01]',
      paste(
        '"float metallic"',
        metallic,
        '"float speculartint" 1 Shape "sphere"'
      )
    ))
    imported = suppressWarnings(read_pbrt(file, strict = FALSE))
    mapped = imported$scene$material[[1]]$openpbr
    expect_equal(mapped$base_metalness, metallic)
    expect_equal(mapped$base_color, c(.9, .1, .01))
    expect_equal(mapped$specular_color, rep(1, 3))
    expect_false(any(grepl("tint exceeds one", imported$diagnostics$message)))
  }
})

test_that("Disney partial metals preserve diffuse color and normal Fresnel", {
  color = c(.13, .37, .81)
  f0 = ((1.6 - 1) / (1.6 + 1))^2
  for (m in c(0, .001, .2, .5, .9, .999, 1)) {
    file = pbrt_test_file(c(
      'WorldBegin Material "disney" "rgb color" [.13 .37 .81]',
      paste('"float metallic"', m, '"float eta" 1.6 Shape "sphere"')
    ))
    p = suppressWarnings(read_pbrt(file, strict = FALSE))$scene$material[[
      1
    ]]$openpbr
    base = p$base_weight * p$base_color
    converted_f0 = ((p$specular_ior - 1) / (p$specular_ior + 1))^2
    expect_equal((1 - p$base_metalness) * base, (1 - m) * color)
    expect_equal(
      (1 - p$base_metalness) * converted_f0 + p$base_metalness * base,
      (1 - m^2) * f0 + m^2 * color
    )
    expect_equal(p$base_diffuse_roughness, 0)
  }
})

test_that("Disney fitted tint approaches metal endpoints continuously", {
  for (m in c(1e-8, 1 - 1e-8)) {
    file = pbrt_test_file(c(
      'WorldBegin Material "disney" "rgb color" [.9 .1 .01]',
      paste('"float metallic"', m, '"float speculartint" 1 Shape "sphere"')
    ))
    p = suppressWarnings(read_pbrt(file, strict = FALSE))$scene$material[[
      1
    ]]$openpbr
    expect_lt(max(abs(p$specular_color - 1)), 1e-7)
  }
})

test_that("Disney fitted sheen preserves tint ratios without clipping its color", {
  file = pbrt_test_file(c(
    'WorldBegin Material "disney" "rgb color" [.9 .1 .01]',
    '"float sheen" 1 "float sheentint" 1 Shape "sphere"'
  ))
  imported = suppressWarnings(read_pbrt(file, strict = FALSE))
  p = imported$scene$material[[1]]$openpbr
  expect_equal(p$fuzz_color / p$fuzz_color[1], c(1, 1 / 9, 1 / 90))
  expect_equal(max(p$fuzz_color), 1)
  expect_lt(p$fuzz_weight, .5)
  expect_false(any(grepl("tint exceeds one", imported$diagnostics$message)))
})

test_that("Disney fitted coat is monotone and retains calibrated endpoints", {
  file = pbrt_test_file(c(
    'WorldBegin',
    vapply(
      c(0, .5, 1),
      function(gloss) {
        paste(
          'Material "disney" "float clearcoat" 1 "float clearcoatgloss"',
          gloss,
          'Shape "sphere"'
        )
      },
      character(1)
    )
  ))
  materials = suppressWarnings(read_pbrt(file, strict = FALSE))$scene$material
  r = vapply(materials, function(x) x$openpbr$coat_roughness, numeric(1))
  expect_true(all(diff(r) < 0))
  expect_equal(r[c(1, 3)], .936 * c(.1, .001)^.43)
  expect_equal(materials[[1]]$openpbr$coat_weight, .0605)
  expect_equal(materials[[1]]$openpbr$coat_darkening, 0)
})

test_that("Disney surface fits leave physical volume and thin interfaces alone", {
  for (parameter in c(
    '"float spectrans" .2',
    '"float scatterdistance" .1',
    '"bool thin" true',
    '"float eta" .8'
  )) {
    file = pbrt_test_file(c(
      'WorldBegin Material "disney" "float metallic" .3',
      if (!grepl("eta", parameter)) '"float eta" 1.4',
      '"float clearcoat" 1 "float sheen" .5',
      parameter,
      'Shape "sphere"'
    ))
    p = suppressWarnings(read_pbrt(file, strict = FALSE))$scene$material[[
      1
    ]]$openpbr
    expect_equal(p$specular_ior, if (grepl('eta', parameter)) .8 else 1.4)
    expect_equal(p$base_metalness, .3)
    expect_equal(p$coat_weight, .25)
  }
})

test_that("Disney surface fits diagnose saturation for extreme Fresnel and sheen", {
  file = pbrt_test_file(c(
    'WorldBegin Material "disney" "float eta" 100 "float metallic" .9 Shape "sphere"',
    'Material "disney" "rgb color" [0 0 1] "float sheen" 1 "float sheentint" 1 Shape "sphere"'
  ))
  imported = suppressWarnings(read_pbrt(file, strict = FALSE))
  expect_true(any(grepl("Fresnel exceeds", imported$diagnostics$message)))
  expect_true(any(grepl("sheen weight exceeds", imported$diagnostics$message)))
  expect_equal(imported$scene$material[[2]]$openpbr$fuzz_weight, 1)
})

test_that("Disney solid scattering and thin diffuse transmission remain distinct", {
  imported = suppressWarnings(read_pbrt(
    pbrt_test_file(c(
      'WorldBegin Material "disney" "rgb color" [.2 .4 .6]',
      '"rgb scatterdistance" [.1 .02 0] Shape "sphere"',
      'Material "disney" "bool thin" true "float difftrans" .8 "float flatness" .7',
      '"rgb color" [.2 .4 .6] "float spectrans" .3 Shape "sphere"',
      'Material "disney" "float scatterdistance" 0 Shape "sphere"'
    )),
    strict = FALSE
  ))
  solid = imported$scene$material[[1]]
  thin = imported$scene$material[[2]]
  expect_equal(solid$openpbr$subsurface_weight, 1)
  # Physical extinction MFPs, not the old diffusion-profile distances.
  expect_equal(
    solid$openpbr$subsurface_radius * solid$openpbr$subsurface_radius_scale,
    c(.04016832, .004955904, 0)
  )
  expect_gt(solid$openpbr$subsurface_color[1], .2)
  expect_equal(solid$openpbr$subsurface_scatter_anisotropy, 0)
  expect_equal(solid$openpbr$specular_ior, 1.5)
  expect_equal(solid$subsurface$openpbr, solid$openpbr)
  expect_true(thin$openpbr$geometry_thin_walled)
  expect_equal(thin$openpbr$subsurface_weight, .8)
  expect_equal(thin$openpbr$transmission_color, c(.2, .4, .6))
  expect_null(thin$subsurface)
  expect_null(imported$scene$material[[3]]$subsurface)
  expect_true(any(grepl("flatness", imported$diagnostics$message)))
  expect_true(any(grepl(
    "internal-reflection IOR 1.4",
    imported$diagnostics$message
  )))
  expect_true(any(grepl("max_depth", imported$diagnostics$message)))
})

test_that("Disney SSS conversion reproduces fitted volume albedo and RGB lengths", {
  # Independent values from Hyperion 2018, section 4.4.2, for A=(.08,.5,.8).
  mapped = pbrt_disney_subsurface(c(.08, .5, .8), c(.1, .2, .3))
  expect_equal(
    mapped$radius * mapped$radius_scale,
    c(.059699791872, .041055, .02897856),
    tolerance = 1e-9
  )
  # Decode using the actual native OpenPBR equation, so this catches passing
  # volume alpha directly or forgetting to re-encode the authored color.
  color = mapped$color
  s = 4.09712 +
    4.20863 * color -
    sqrt(9.59217 + 41.6808 * color + 17.7126 * color^2)
  alpha = 1 - s^2
  expect_equal(
    alpha,
    c(.5609270688, .9729143354, .9983758191),
    tolerance = 1e-5
  )
  scaled = pbrt_disney_subsurface(c(.08, .5, .8), c(1, 2, 3))
  expect_equal(scaled$color, mapped$color)
  expect_equal(scaled$radius, 10 * mapped$radius)
  expect_equal(scaled$radius_scale, mapped$radius_scale)

  endpoints = pbrt_disney_subsurface(c(0, 1, .5), c(.1, .1, 0))
  expect_equal(endpoints$color[1:2], c(0, 1))
  expect_equal(endpoints$radius_scale[3], 0)
  zero = pbrt_disney_subsurface(rep(.5, 3), rep(0, 3))
  expect_equal(zero$radius, 0)
  expect_equal(zero$radius_scale, rep(0, 3))
  ramp = vapply(
    seq(0, 1, length.out = 1001),
    function(a) {
      pbrt_disney_subsurface(rep(a, 3), rep(1, 3))$color[1]
    },
    numeric(1)
  )
  expect_true(all(is.finite(ramp) & ramp >= 0 & ramp <= 1))
  expect_true(all(diff(ramp) >= 0))
})

test_that("Disney texture imports retain base color, raw roughness, and bump amplitude", {
  directory = tempfile("disney-textures-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  png::writePNG(
    array(rep(c(.2, .4, .6), each = 4), c(2, 2, 3)),
    file.path(directory, "color.png")
  )
  png::writePNG(matrix(c(.2, .8, .2, .8), 2), file.path(directory, "rough.png"))
  file = file.path(directory, "scene.pbrt")
  writeLines(
    c(
      'WorldBegin Texture "color" "spectrum" "imagemap" "string filename" "color.png"',
      '"string encoding" "linear" "float uscale" 2 "float vscale" 3',
      'Texture "rough" "float" "imagemap" "string filename" "rough.png" "string encoding" "linear"',
      'Texture "bump" "float" "checkerboard" "float tex1" 0 "float tex2" .02 "float uscale" 2 "float vscale" 2',
      'Texture "metal" "float" "constant" "float value" .7',
      'Material "disney" "texture color" "color" "texture roughness" "rough"',
      '"texture bumpmap" "bump" "texture metallic" "metal" Shape "sphere"'
    ),
    file
  )
  imported = suppressWarnings(read_pbrt(
    file,
    strict = FALSE,
    asset_dir = file.path(directory, "assets")
  ))
  material = imported$scene$material[[1]]
  expect_equal(material$openpbr$base_weight, .79)
  expect_equal(material$openpbr$base_metalness, .49 / .79)
  expect_equal(material$image_repeat[[1]], c(2, 3))
  expect_true(file.exists(material$image))
  roughness = png::readPNG(material$roughness_texture)
  expect_equal(as.numeric(roughness), c(.2, .8, .2, .8), tolerance = 1 / 255)
  expect_equal(material$bump_intensity, .02, tolerance = 1e-4)
  expect_true(file.exists(material$bump_texture))
  expect_true(any(grepl(
    "roughness-map UV repeat",
    imported$diagnostics$message
  )))
  expect_no_error(process_scene(imported$scene))
})

test_that("Disney conversion validates ranges and diagnoses roughness clipping", {
  for (parameter in c(
    '"float metallic" -1',
    '"float clearcoatgloss" 2',
    '"float eta" 0',
    '"float roughness" -1',
    '"rgb scatterdistance" [1 -1 1]'
  )) {
    file = pbrt_test_file(paste(
      'WorldBegin Material "disney"',
      parameter,
      'Shape "sphere"'
    ))
    expect_error(read_pbrt(file, strict = FALSE), "Disney")
  }
  imported = suppressWarnings(read_pbrt(
    pbrt_test_file(
      'WorldBegin Material "disney" "float roughness" 1 "float anisotropic" 1 Shape "sphere"'
    ),
    strict = FALSE
  ))
  expect_equal(imported$scene$material[[1]]$openpbr$specular_roughness, 1)
  expect_true(any(grepl("clamped to one", imported$diagnostics$message)))
})

test_that("PBRT animation retains scoped endpoints and transform times", {
  imported = read_pbrt(pbrt_test_file(c(
    'TransformTimes 2 3',
    'Camera "perspective" "float shutteropen" 2 "float shutterclose" 3',
    'WorldBegin',
    'ObjectBegin "ball" Shape "sphere" ObjectEnd',
    'AttributeBegin ActiveTransform EndTime Translate 2 0 0',
    'CoordinateSystem "moving" ObjectInstance "ball" AttributeEnd',
    'Shape "sphere"',
    'CoordSysTransform "moving" Shape "sphere"'
  )))
  animation = imported$scene$animation_info[[1]]
  expect_equal(animation$start_transform_animation[[1]], diag(4))
  expect_equal(animation$end_transform_animation[[1]][1, 4], 2)
  expect_equal(c(animation$start_time, animation$end_time), c(2, 3))
  expect_true(all(is.na(imported$scene$animation_info[[
    2
  ]]$start_transform_animation[[1]])))
  expect_equal(
    imported$scene$animation_info[[3]]$end_transform_animation[[1]][1, 4],
    2
  )
  expect_equal(imported$render_args$shutteropen, 2)
})

test_that("motion bounds include a shutter outside the default time interval", {
  skip_on_cran()
  imported = read_pbrt(pbrt_test_file(c(
    'TransformTimes 2 3',
    'LookAt 0 0 5 0 0 0 0 1 0 Camera "orthographic" "float shutteropen" 3 "float shutterclose" 3',
    'WorldBegin AreaLightSource "diffuse" "rgb L" [1 1 1]',
    'Material "diffuse" "rgb reflectance" [0 0 0]',
    'ActiveTransform StartTime Translate 4 0 0 Shape "sphere" "float radius" .2'
  )))
  args = imported$render_args
  args$width = 8
  args$height = 8
  args$samples = 1
  args$ortho_dimensions = c(.01, .01)
  args$preview = args$progress = args$parallel = args$plot_scene = FALSE
  args$denoise = args$bloom = FALSE
  args$tonemap = "raw"
  image = do.call(render_scene, c(list(scene = imported$scene), args))
  expect_equal(mean(image[,, 1:3]), 1, tolerance = 1e-5)
})

test_that("PBRT moving cameras become one blurred exposure", {
  imported = read_pbrt(pbrt_test_file(c(
    'LookAt 0 0 5 0 0 0 0 1 0',
    'ActiveTransform EndTime Translate 1 0 0 Camera "perspective"',
    'WorldBegin Shape "sphere"'
  )))
  camera = imported$render_args$camera
  expect_s3_class(camera, "ray_camera")
  expect_equal(nrow(camera$motion), 2)
  expect_true(camera$camera_motion_blur)
  expect_equal(camera$shutter_speed, 1)
  expect_equal(imported$render_args$mode, "image")
  expect_equal(abs(diff(camera$motion$x)), 1)
})

test_that("PBRT grid definitions keep their world placement independently of boundaries", {
  imported = read_pbrt(pbrt_test_file(c(
    'WorldBegin Translate 2 0 0',
    'MakeNamedMedium "fog" "string type" "uniformgrid" "integer nx" 2 "integer ny" 1 "integer nz" 1',
    '"float density" [.2 .8] "rgb sigma_a" [.1 .2 .3] "float temperatureoffset" 100 "float temperaturescale" 2',
    'Translate 3 0 0 Material "interface" MediumInterface "fog" "" Shape "sphere"'
  )))
  medium = imported$scene$shape_info[[1]]$medium
  expect_equal(as.vector(medium$density), c(.2, .8))
  expect_equal(medium$medium_transform[1, 4], -3)
  expect_equal(medium$temperature_offset, 100)
  expect_equal(medium$temperature_scale, 2)
})

test_that("PBRT displacement preserves signed bump amplitude", {
  imported = read_pbrt(pbrt_test_file(c(
    'WorldBegin',
    'Texture "check" "float" "checkerboard" "float tex1" 0 "float tex2" 1 "float uscale" 2 "float vscale" 2',
    'Texture "height" "float" "scale" "texture tex" "check" "float scale" -.05',
    'Material "diffuse" "texture displacement" "height"',
    'Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0] "integer indices" [0 1 2]'
  )))
  material = imported$scene$material[[1]]
  expect_equal(material$bump_intensity, .05, tolerance = 1e-6)
  expect_true(file.exists(material$bump_texture))
  expect_equal(range(png::readPNG(material$bump_texture)), c(0, 1))
})

test_that("PBRT PLY displacement moves the mesh along its normals", {
  file = tempfile(fileext = ".ply")
  writeLines(
    c(
      "ply",
      "format ascii 1.0",
      "element vertex 3",
      "property float x",
      "property float y",
      "property float z",
      "element face 1",
      "property list uchar int vertex_indices",
      "end_header",
      "0 0 0",
      "1 0 0",
      "0 1 0",
      "3 0 1 2"
    ),
    file
  )
  input = pbrt_test_file(c(
    'WorldBegin',
    'Texture "raise" "float" "constant" "float value" .25',
    sprintf(
      'Shape "plymesh" "string filename" "%s" "texture displacement" "raise" "float edgelength" 2',
      file
    )
  ))
  imported = suppressWarnings(read_pbrt(input, strict = FALSE))
  mesh = imported$scene$shape_info[[1]]$mesh_info[[1]]
  expect_equal(mesh$vertices[[1]][, 3], rep(.25, 3))
  expect_true(any(grepl("PLY displacement", imported$diagnostics$message)))
})

test_that("PBRT realistic cameras retain lens files and physical parameters", {
  lens = tempfile(fileext = ".dat")
  writeLines(c("# test lens", "20 5 1.5 12", "-20 30 1 12"), lens)
  imported = read_pbrt(pbrt_test_file(c(
    'Film "rgb" "float diagonal" 36',
    sprintf(
      'Camera "realistic" "string lensfile" "%s" "float aperturediameter" 8 "float focusdistance" 4',
      lens
    ),
    'WorldBegin Shape "sphere"'
  )))
  expect_equal(
    imported$render_args$camera_description_file,
    normalizePath(lens)
  )
  expect_equal(imported$render_args$film_size, 36)
  expect_equal(imported$render_args$aperture, 8)
  expect_equal(imported$render_args$focal_distance, 4)
  expect_equal(imported$render_args$fov, -1)
})

test_that("PBRT NanoVDB imports permit density-only files", {
  skip_on_cran()
  file = normalizePath(test_path("fixtures", "volumes", "tiles.nvdb"))
  imported = suppressWarnings(read_pbrt(
    pbrt_test_file(c(
      'WorldBegin',
      sprintf(
        'MakeNamedMedium "fog" "string type" "nanovdb" "string filename" "%s" "rgb sigma_s" [0 0 0]',
        file
      ),
      'Material "interface" MediumInterface "fog" "" Shape "sphere"'
    )),
    strict = FALSE
  ))
  medium = imported$scene$shape_info[[1]]$medium
  expect_equal(medium$filename, file)
  expect_true(medium$temperature_optional)
  expect_no_error(render_scene(
    imported$scene,
    width = 8,
    height = 8,
    samples = 1,
    preview = FALSE,
    progress = FALSE,
    parallel = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE
  ))
})

test_that("a still exposure uses the camera motion endpoint", {
  skip_on_cran()
  imported = read_pbrt(pbrt_test_file(c(
    'LookAt 0 0 5 0 0 0 0 1 0',
    'ActiveTransform EndTime Translate 2 0 0 Camera "orthographic"',
    'WorldBegin AreaLightSource "diffuse" "rgb L" [1 1 1]',
    'Material "diffuse" "rgb reflectance" [0 0 0] Shape "sphere" "float radius" .2'
  )))
  render = function(blur, groups = c(1L, 1L), enabled = NULL) {
    args = imported$render_args
    args$camera$camera_motion_blur = blur
    args$camera$motion$camera_motion_blur_group = groups
    if (!is.null(enabled)) {
      args$camera$motion$camera_motion_blur = rep(enabled, 2)
    }
    args = c(
      args,
      list(
        width = 32,
        height = 16,
        samples = 64,
        preview = FALSE,
        progress = FALSE,
        parallel = FALSE,
        plot_scene = FALSE,
        denoise = FALSE,
        bloom = FALSE,
        min_variance = 0,
        tonemap = "raw"
      )
    )
    args = args[!duplicated(names(args), fromLast = TRUE)]
    set.seed(4)
    do.call(render_scene, c(list(scene = imported$scene), args))[,, 1]
  }
  still = render(FALSE)
  blurred = render(TRUE)
  expect_gt(sum(colSums(blurred) > .01), sum(colSums(still) > .01))
  expect_true(all(is.finite(blurred)))
  expect_equal(render(TRUE, groups = c(1L, 2L)), still)
  expect_equal(render(TRUE, enabled = FALSE), still)
})

test_that("PBRT parsing handles comments, escapes, arrays and multiple directives", {
  file = pbrt_test_file(c(
    '# a comment',
    'WorldBegin Material "diffuse" "rgb reflectance" [.2 3e-1 +.4]',
    'MakeNamedMaterial "a#b\\\"c" "string type" "diffuse"',
    'NamedMaterial "a#b\\\"c" Shape "sphere" "float radius" [2]'
  ))
  commands = pbrt_parse_file(file)
  expect_length(commands, 5)
  expect_equal(commands[[2]]$params$reflectance$value, c(.2, .3, .4))
  expect_equal(commands[[3]]$args, 'a#b"c')
  scene = read_pbrt(file)
  expect_equal(scene$scene$shape_info[[1]]$shape_properties$radius, 2)
  expect_equal(nrow(scene$diagnostics), 0)
})

test_that("PBRT invalid syntax has a useful file and line", {
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin\nShape "sphere" "float radius" [1')),
    ':2: Shape: Unclosed'
  )
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin\nShape "sphere')),
    ':2: malformed quoted'
  )
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin Shape "sphere" "float radius" Inf')),
    'nonfinite'
  )
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin Shape "sphere" "float radius" [1 2]')),
    'radius requires 1'
  )
  expect_error(
    read_pbrt(pbrt_test_file(
      'WorldBegin Shape "sphere" "float radius" 1 "float radius" 2'
    )),
    'Duplicate parameter'
  )
  expect_error(
    read_pbrt(pbrt_test_file(
      'WorldBegin Shape "trianglemesh" "point3 P" [1 2]'
    )),
    'tuple length'
  )
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin AttributeEnd')),
    'Unmatched'
  )
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin AttributeBegin')),
    'Unclosed'
  )
  expect_error(
    read_pbrt(pbrt_test_file('Shape "sphere"')),
    'must follow WorldBegin'
  )
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin NamedMaterial "missing"')),
    'Undefined material'
  )
  expect_error(
    read_pbrt(pbrt_test_file(
      'WorldBegin Material "diffuse" "texture reflectance" "missing"'
    )),
    'Undefined texture'
  )
  expect_error(read_pbrt(pbrt_test_file(character())), 'Empty PBRT file')
})

test_that("PBRT Includes use the main file directory, cache parsing and detect cycles", {
  directory = tempfile()
  dir.create(file.path(directory, "parts"), recursive = TRUE)
  root = file.path(directory, "scene.pbrt")
  writeLines(
    'WorldBegin Include "parts/one.pbrt" Include "parts/one.pbrt"',
    root
  )
  writeLines(
    'Translate 1 0 0 Include "two.pbrt.gz"',
    file.path(directory, "parts/one.pbrt")
  )
  connection = gzfile(file.path(directory, "two.pbrt.gz"), "wt")
  writeLines('Shape "sphere"', connection)
  close(connection)
  result = read_pbrt(root)
  expect_length(result$source_files, 3)
  expect_equal(nrow(result$scene), 2)
  expect_equal(result$scene$transforms[[2]]$group_transform[[1]][1, 4], 2)
  writeLines('Include "scene.pbrt"', file.path(directory, "parts/one.pbrt"))
  expect_error(read_pbrt(root), 'Cyclic PBRT Include')
  writeLines('WorldBegin Include "missing.pbrt"', root)
  expect_error(read_pbrt(root, strict = FALSE), 'Asset not found')
})

test_that("PBRT transforms compose on the right and restore scopes", {
  scene = read_pbrt(pbrt_test_file(c(
    'WorldBegin Translate 1 2 3',
    'AttributeBegin Scale 2 3 4 Material "diffuse" "rgb reflectance" [1 0 0]',
    'ReverseOrientation Shape "sphere" AttributeEnd Shape "sphere"',
    'CoordinateSystem "saved" Identity CoordSysTransform "saved" Shape "sphere"',
    'Transform [1 0 0 0  .25 1 0 0  0 0 1 0  7 8 9 1] Shape "sphere"'
  )))$scene
  expect_equal(diag(scene$transforms[[1]]$group_transform[[1]]), c(2, 3, 4, 1))
  expect_equal(scene$transforms[[1]]$group_transform[[1]][1:3, 4], c(1, 2, 3))
  expect_true(scene$shape_info[[1]]$flipped)
  expect_false(scene$shape_info[[2]]$flipped)
  expect_equal(scene$material[[2]]$properties[[1]][1:3], rep(.5, 3))
  expect_equal(scene$transforms[[2]], scene$transforms[[3]])
  expect_equal(scene$transforms[[4]]$group_transform[[1]][1, 2], .25)
  expect_equal(scene$transforms[[4]]$group_transform[[1]][1:3, 4], c(7, 8, 9))
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin Scale 0 1 1 Shape "sphere"')),
    'nonsingular'
  )
  expect_error(
    read_pbrt(pbrt_test_file('LookAt 0 0 0 0 0 1 0 0 1 WorldBegin')),
    'parallel'
  )
})

test_that("PBRT camera transform and portrait FOV preserve framing", {
  result = read_pbrt(pbrt_test_file(c(
    'LookAt 2 3 -5 2 3 0 0 1 0 Camera "perspective" "float fov" 60',
    'Film "rgb" "integer xresolution" 300 "integer yresolution" 600 "string filename" "do-not-write.exr"',
    'Sampler "halton" "integer pixelsamples" 12 WorldBegin Shape "sphere"'
  )))
  expect_equal(result$render_args$lookfrom, c(2, 3, -5))
  expect_equal(result$render_args$lookat, c(2, 3, -4))
  expect_equal(result$render_args$camera_up, c(0, 1, 0))
  expect_equal(result$render_args$fov, 360 / pi * atan(tan(pi / 6) * 2))
  expect_equal(result$render_args$samples, 12)
  expect_null(result$render_args$filename)
  orthographic = read_pbrt(pbrt_test_file(
    'Camera "orthographic" Film "rgb" "integer xresolution" 600 "integer yresolution" 300 WorldBegin Shape "sphere"'
  ))
  expect_equal(orthographic$render_args$fov, 0)
  expect_equal(orthographic$render_args$ortho_dimensions, c(4, 2))
})

test_that("PBRT object definitions are instanced and do not appear on their own", {
  result = read_pbrt(pbrt_test_file(c(
    'WorldBegin ObjectBegin "pair" Shape "sphere" Translate 2 0 0 Shape "sphere" ObjectEnd',
    'Translate 4 0 0 ObjectInstance "pair" Translate 0 3 0 ObjectInstance "pair"'
  )))
  expect_equal(result$scene$shape, rep("instance", 2))
  original = result$scene$shape_info[[1]]$shape_properties$original_scene[[1]]
  expect_equal(nrow(original), 2)
  expect_equal(
    result$scene$transforms[[2]]$group_transform[[1]][1:3, 4],
    c(4, 3, 0)
  )
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin ObjectInstance "missing"')),
    'Undefined object'
  )
})

test_that("PBRT geometry keeps zero-based topology, normals, UVs and axes", {
  result = read_pbrt(pbrt_test_file(c(
    'WorldBegin Shape "trianglemesh" "point3 P" [0 0 0  1 0 0  0 1 0]',
    '"normal3 N" [0 0 1 0 0 1 0 0 1] "point2 uv" [0 0 1 0 0 1]',
    'Shape "bilinearmesh" "point3 P" [0 0 0 1 0 0 0 1 0 1 1 0]',
    'Shape "cylinder" "float zmin" 2 "float zmax" 6',
    'Shape "disk" "float height" 3'
  )))
  mesh = result$scene$shape_info[[1]]$mesh_info[[1]]
  expect_equal(unclass(mesh$shapes)[[1]]$indices, matrix(0:2, 1))
  expect_equal(
    mesh$vertices[[1]],
    matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0), 3, byrow = TRUE)
  )
  quad = result$scene$shape_info[[2]]$mesh_info[[1]]
  expect_equal(unclass(quad$shapes)[[1]]$indices, rbind(c(0, 1, 3), c(0, 3, 2)))
  expect_false(result$scene$shape_info[[3]]$shape_properties$has_cap)
  expect_equal(
    result$scene$transforms[[3]]$group_transform[[1]][1:3, 4],
    c(0, 0, 4)
  )
  expect_equal(
    result$scene$transforms[[4]]$group_transform[[1]][1:3, 2],
    c(0, 0, 1)
  )
  expect_error(
    read_pbrt(pbrt_test_file(
      'WorldBegin Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0] "integer indices" [0 1 3]'
    )),
    'indices'
  )
})

test_that("PBRT roughness is converted through GGX alpha", {
  result = read_pbrt(pbrt_test_file(c(
    'WorldBegin',
    'Material "dielectric" "float eta" 1.4 Shape "sphere"',
    'Material "dielectric" "float roughness" .0625 Shape "sphere"',
    'Material "conductor" "rgb reflectance" [.4 .5 .6] "bool remaproughness" false "float roughness" .25 Shape "sphere"'
  )))
  expect_equal(result$scene$material[[1]]$type, 3L)
  expect_equal(result$scene$material[[1]]$properties[[1]][4], 1.4)
  expected = microfacet(transmission = TRUE, eta = 1.5, roughness = .5)
  expect_equal(result$scene$material[[2]], expected[[1]])
  expect_equal(result$scene$material[[3]]$type, 6L)
})

test_that("PBRT named textures and inherited parameters are scoped", {
  result = read_pbrt(pbrt_test_file(c(
    'WorldBegin',
    'Texture "red" "spectrum" "constant" "rgb value" [1 0 0]',
    'MakeNamedMaterial "red" "string type" "diffuse" "texture reflectance" "red"',
    'AttributeBegin Attribute "shape" "float radius" 3 NamedMaterial "red" Shape "sphere" AttributeEnd',
    'Shape "sphere"'
  )))
  expect_equal(result$scene$material[[1]]$properties[[1]][1:3], c(1, 0, 0))
  expect_equal(result$scene$shape_info[[1]]$shape_properties$radius, 3)
  expect_equal(result$scene$shape_info[[2]]$shape_properties$radius, 1)
})

test_that("constant scalar texture references work for material parameters", {
  result = read_pbrt(pbrt_test_file(c(
    'WorldBegin Texture "rough" "float" "constant" "float value" .0625',
    'Material "dielectric" "texture roughness" "rough" Shape "sphere"'
  )))
  expect_equal(
    result$scene$material[[1]],
    microfacet(transmission = TRUE, eta = 1.5, roughness = .5)[[1]]
  )
  expect_error(
    read_pbrt(
      pbrt_test_file(
        'WorldBegin Material "dielectric" "texture roughness" "undefined" Shape "sphere"'
      ),
      strict = FALSE
    ),
    'Undefined texture'
  )
})

test_that("PBRT media attach to bounded interfaces and preserve named definitions", {
  result = read_pbrt(pbrt_test_file(c(
    'WorldBegin',
    'MakeNamedMedium "fog" "string type" "homogeneous" "rgb sigma_a" [.1 .2 .3] "rgb sigma_s" [1 2 3]',
    'MediumInterface "fog" "" Material "interface" Shape "sphere"'
  )))
  medium = result$scene$shape_info[[1]]$medium
  expect_equal(medium$sigma_a, c(.1, .2, .3))
  expect_false(result$scene$shape_info[[1]]$medium_keep_surface)
  expect_error(
    read_pbrt(pbrt_test_file('WorldBegin MediumInterface "undefined" ""')),
    'Undefined medium'
  )
})

test_that("PBRT unsupported features fail strictly and are inspectable permissively", {
  file = pbrt_test_file(
    'WorldBegin Shape "mystery" Shape "sphere" "float unknown" 3'
  )
  expect_error(read_pbrt(file), 'Unsupported shape')
  expect_warning(result <- read_pbrt(file, strict = FALSE), '2 diagnostic')
  expect_equal(nrow(result$scene), 1)
  expect_equal(nrow(result$diagnostics), 2)
  expect_match(result$diagnostics$message[2], 'unknown')
  expect_warning(
    result <- read_pbrt(
      pbrt_test_file('WorldBegin Material "coateddiffuse" Shape "sphere"'),
      strict = FALSE
    ),
    'OpenPBR'
  )
  expect_equal(result$scene$material[[1]]$type, 11L)
})

test_that("PBRT textures preserve linear color values and resolve image paths", {
  directory = tempfile()
  dir.create(directory)
  png::writePNG(array(.5, c(4, 4, 3)), file.path(directory, "color#map.png"))
  file = file.path(directory, "scene.pbrt")
  writeLines(
    c(
      'WorldBegin',
      'Texture "image" "spectrum" "imagemap" "string filename" "color#map.png" "string encoding" "linear"',
      'Texture "checks" "spectrum" "checkerboard" "rgb tex1" [.2 .4 .6] "rgb tex2" [.8 .6 .4] "float uscale" 2 "float vscale" 2',
      'Material "diffuse" "texture reflectance" "image"',
      'Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0]',
      'Material "diffuse" "texture reflectance" "checks"',
      'Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0]'
    ),
    file
  )
  x = read_pbrt(file, asset_dir = file.path(directory, "assets"))
  expect_true(all(file.exists(x$assets)))
  image = rayimage::ray_read_image(x$scene$material[[1]]$image)
  expect_equal(as.numeric(image[1, 1, 1]), 128 / 255, tolerance = .001)
  checker = rayimage::ray_read_image(x$scene$material[[2]]$image)
  expect_equal(
    sort(unique(as.numeric(checker[,, 1]))),
    c(.2, .8),
    tolerance = .001
  )
  expect_no_error(process_scene(x$scene))
})

test_that("PBRT image-map declarations accept offsets and retain independent repeat", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  image = array(seq_len(48) / 50, c(4, 4, 3))
  png::writePNG(image, file.path(directory, "map.png"))
  source_image = pbrt_read_image(
    file.path(directory, "map.png"),
    source_linear = TRUE
  )
  file = file.path(directory, "scene.pbrt")
  for (offset in list(c(0, 0), c(.25, -.25), c(-.25, .25), c(.125, .125))) {
    writeLines(
      c(
        'WorldBegin',
        sprintf(
          paste(
            'Texture "image" "spectrum" "imagemap" "string filename" "map.png"',
            '"string encoding" "linear" "float uscale" 2 "float vscale" 3',
            '"float udelta" %s "float vdelta" %s'
          ),
          offset[1],
          offset[2]
        ),
        'Material "diffuse" "texture reflectance" "image"',
        'Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0]'
      ),
      file
    )
    expect_no_warning(
      result <- read_pbrt(file, asset_dir = file.path(directory, "assets"))
    )
    expect_equal(nrow(result$diagnostics), 0L)
    material = result$scene$material[[1]]
    expect_equal(material$image_repeat[[1]], c(2, 3))
    expect_equal(material$image_offset[[1]], offset)
    baked = rayimage::ray_read_image(material$image, convert_to_array = TRUE)
    expect_equal(
      baked,
      source_image,
      tolerance = .001,
      ignore_attr = TRUE
    )
  }
})

test_that("PBRT equal-area environments are resampled with their transforms", {
  image = array(.5, c(64, 64, 3))
  constant = pbrt_environment_image(image, diag(4))
  expect_equal(dim(constant), c(64, 128, 3))
  expect_equal(range(constant), c(.5, .5))
  image[,, 1] = matrix(rep((seq_len(64) - .5) / 64, each = 64), 64)
  image[,, 2] = t(image[,, 1])
  result = pbrt_environment_image(image, diag(4))
  expect_gt(result[32, 96, 1], .95) # world +X -> PBRT square (+1, .5)
  expect_gt(result[1, 32, 2], .95) # world +Y -> PBRT square (.5, +1)
  expect_equal(result[32, 1, 1:2], c(.5, .5), tolerance = .1) # world +Z -> center
  rotation = diag(4)
  rotation[1:3, 1:3] = matrix(c(1, 0, 0, 0, 0, 1, 0, -1, 0), 3)
  rotated = pbrt_environment_image(image, rotation)
  expect_equal(rotated[64, 32, 1:2], c(.5, .5), tolerance = .1)
  directory = tempfile()
  dir.create(directory)
  rayimage::ray_write_image(
    image,
    file.path(directory, 'environment.exr'),
    write_linear = TRUE
  )
  file = file.path(directory, 'scene.pbrt')
  writeLines(
    'WorldBegin LightSource "infinite" "string filename" "environment.exr" Shape "sphere"',
    file
  )
  x = read_pbrt(file)
  expect_length(attr(x$scene, "ray_infinite_lights"), 1)
  expect_true(file.exists(x$assets[1]))
})

test_that("PBRT PLY and cubic curve assets convert without external tools", {
  directory = tempfile()
  dir.create(directory)
  writeLines(
    c(
      'ply',
      'format ascii 1.0',
      'element vertex 3',
      'property float x',
      'property float y',
      'property float z',
      'element face 1',
      'property list uchar int vertex_indices',
      'end_header',
      '0 0 0',
      '1 0 0',
      '0 1 0',
      '3 0 1 2'
    ),
    file.path(directory, 'triangle.ply')
  )
  file = file.path(directory, 'scene.pbrt')
  writeLines(
    c(
      'WorldBegin Shape "plymesh" "string filename" "triangle.ply"',
      'Shape "curve" "point3 P" [0 0 0 0 1 0 1 2 0 1 3 0 2 4 0 2 5 0 3 6 0] "float width" .1'
    ),
    file
  )
  x = read_pbrt(file)
  expect_equal(x$scene$shape, c('ply', 'curve', 'curve'))
  expect_equal(
    x$scene$shape_info[[1]]$fileinfo,
    normalizePath(file.path(directory, 'triangle.ply'))
  )
})

test_that("a converted PBRT fixture renders geometry and area emission", {
  skip_on_cran()
  withr::local_options(cores = 2L)
  result = read_pbrt(pbrt_test_file(c(
    'LookAt 0 0 -5 0 0 0 0 1 0 Camera "perspective" "float fov" 30',
    'Film "rgb" "integer xresolution" 24 "integer yresolution" 24',
    'Sampler "independent" "integer pixelsamples" 4 WorldBegin',
    'Material "diffuse" "rgb reflectance" [0 0 0] AreaLightSource "diffuse" "rgb L" [2 1 .5] Shape "sphere"'
  )))
  image = do.call(
    render_scene,
    c(
      list(scene = result$scene),
      result$render_args,
      list(
        preview = FALSE,
        interactive = FALSE,
        plot_scene = FALSE,
        progress = FALSE,
        denoise = FALSE,
        bloom = FALSE
      )
    )
  )
  expect_true(all(is.finite(image)))
  expect_gt(max(image), .5)
  expect_gt(mean(image[10:15, 10:15, 1]), mean(image[1:3, 1:3, 1]))
})

test_that("PBRT row buffers preserve scene and diagnostic order across chunks", {
  count = 600L
  file = pbrt_test_file(c(
    "WorldBegin",
    as.vector(rbind(
      sprintf('Shape "sphere" "float radius" [%d]', seq_len(count)),
      rep("UnsupportedDirective", count)
    ))
  ))
  imported = suppressWarnings(read_pbrt(file, strict = FALSE))
  expect_equal(nrow(imported$scene), count)
  expect_equal(
    vapply(
      imported$scene$shape_info,
      function(x) x$shape_properties$radius,
      numeric(1)
    ),
    seq_len(count)
  )
  expect_equal(imported$diagnostics$line, seq.int(3L, 2L * count + 1L, by = 2L))
  expect_true(all(imported$diagnostics$directive == "UnsupportedDirective"))
  expect_s3_class(imported$scene, "ray_scene")
  expect_s3_class(imported$scene$shape_info, "ray_shape_info")
  expect_s3_class(imported$scene$material, "ray_material")
})

test_that("PBRT definitions have independent buffers and can be reused", {
  file = pbrt_test_file(c(
    "WorldBegin",
    rep('Shape "sphere" "float radius" [1]', 270),
    'ObjectBegin "many"',
    rep('Shape "sphere" "float radius" [2]', 300),
    'ObjectEnd',
    'ObjectBegin "empty"',
    'ObjectEnd',
    'ObjectInstance "empty"',
    'ObjectInstance "many"',
    'Translate 3 0 0',
    'ObjectInstance "many"',
    'Shape "sphere" "float radius" [4]'
  ))
  scene = read_pbrt(file)$scene
  expect_equal(nrow(scene), 273L)
  expect_equal(scene$shape[c(271L, 272L)], c("instance", "instance"))
  children = lapply(scene$shape_info[c(271L, 272L)], function(x) {
    x$shape_properties$original_scene[[1]]
  })
  expect_identical(children[[1]], children[[2]])
  expect_equal(nrow(children[[1]]), 300L)
  expect_true(all(
    vapply(
      children[[1]]$shape_info,
      function(x) x$shape_properties$radius,
      numeric(1)
    ) ==
      2
  ))
  expect_equal(scene$shape_info[[273L]]$shape_properties$radius, 4)
  expect_true(all(
    vapply(
      scene$shape_info[seq_len(270L)],
      function(x) x$shape_properties$radius,
      numeric(1)
    ) ==
      1
  ))
})

test_that("PBRT grammar preserves large numeric arrays and following source lines", {
  values = seq(-1, 1, length.out = 10001)
  text = c(
    'WorldBegin',
    'Shape "trianglemesh" "float payload" [',
    sprintf('%.17g', values),
    ']',
    'Translate 1 2 3'
  )
  commands = pbrt_parse_file(pbrt_test_file(text))
  # The reference is R's conversion of the text actually stored in the file;
  # R_strtod and the platform strtod can differ by one ULP after formatting.
  expect_identical(
    commands[[2]]$params$payload$value,
    as.numeric(sprintf("%.17g", values))
  )
  expect_equal(commands[[3]]$line, length(text))
  expect_identical(commands[[3]]$args, c('1', '2', '3'))
  for (value in c('1e', '1foo', '@array:1')) {
    expect_error(pbrt_parse_file(pbrt_test_file(paste0(
      'Shape "sphere" "float radius" [',
      value,
      ']'
    ))))
  }
  empty = pbrt_parse_file(pbrt_test_file('Shape "sphere" "string names" []'))
  expect_identical(empty[[1]]$params$names$value, character())
})

test_that("PBRT gzip grammar errors identify the original file", {
  path = tempfile(fileext = ".pbrt.gz")
  output = gzfile(path, "wt")
  writeLines(c("# comment", "stray_token"), output)
  close(output)
  error = tryCatch(pbrt_parse_file(path), error = identity)
  expect_s3_class(error, "error")
  expect_match(conditionMessage(error), paste0(path, ":2:"), fixed = TRUE)
})

test_that("PBRT uniform grids retain spatial emission scale", {
  scene = read_pbrt(pbrt_test_file(c(
    'WorldBegin',
    'MakeNamedMedium "glow" "string type" "uniformgrid"',
    '"integer nx" 2 "integer ny" 1 "integer nz" 1',
    '"float density" [1 1] "rgb Le" [1 2 3] "float Lescale" [0 2]',
    'Material "interface" MediumInterface "glow" "" Shape "sphere"'
  )))
  medium = scene$scene$shape_info[[1]]$medium
  expect_equal(as.vector(medium$emission_scale_grid), c(0, 2))
  expect_equal(unname(dim(medium$emission_scale_grid)), c(2, 1, 1))
  expect_equal(medium$emission, c(1, 2, 3))
})

test_that("PBRT PFM images preserve orientation, channels, scale and byte order", {
  expected = array(seq_len(18) / 20, c(3, 2, 3))
  for (endian in c("little", "big")) {
    path = tempfile(fileext = ".pfm")
    output = file(path, "wb")
    writeLines(c("PF", "2 3", if (endian == "little") "-2" else "2"), output)
    writeBin(
      as.vector(aperm(expected[3:1, , , drop = FALSE], c(3, 2, 1))) / 2,
      output,
      size = 4L,
      endian = endian
    )
    close(output)
    expect_equal(pbrt_read_image(path), expected, tolerance = 1e-7)
    expect_equal(
      pbrt_read_image(path, source_linear = FALSE),
      ifelse(
        expected <= .04045,
        expected / 12.92,
        ((expected + .055) / 1.055)^2.4
      ),
      tolerance = 1e-7
    )
  }
  path = tempfile(fileext = ".pfm")
  output = file(path, "wb")
  writeLines(c("Pf", "2 1", "-1"), output)
  writeBin(c(.25, .5), output, size = 4L, endian = "little")
  close(output)
  expect_equal(pbrt_read_image(path), array(rep(c(.25, .5), 3), c(1, 2, 3)))
  writeLines(c("PF", "2 3", "-1"), path)
  expect_error(pbrt_read_image(path), "Invalid PFM dimensions|Truncated")
  writeLines(c("PF", "2 3", "0"), path)
  expect_error(pbrt_read_image(path), "Invalid PFM dimensions")
})

test_that("large PBRT light lists preserve names, order, and independent values", {
  count = 513L
  imported = read_pbrt(pbrt_test_file(c(
    "WorldBegin",
    sprintf('LightSource "point" "point3 from" [%d 0 0]', seq_len(count)),
    'LightSource "spot" "point3 from" [0 2 0] "point3 to" [0 0 0]',
    'Shape "sphere"'
  )))
  lights = list_lights(imported$scene)
  expect_length(lights, count + 1L)
  expect_identical(names(lights), paste0("pbrt-point-", seq_len(count + 1L)))
  expect_equal(
    vapply(lights[seq_len(count)], function(x) x$position[1], numeric(1)),
    setNames(seq_len(count), names(lights)[seq_len(count)])
  )
  expect_identical(lights[[count + 1L]]$type, "spot")
})

test_that("PBRT retains independent native offsets for color, bump, and roughness", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  png::writePNG(
    array(seq_len(48) / 50, c(4, 4, 3)),
    file.path(directory, "map.png")
  )
  file = file.path(directory, "scene.pbrt")
  writeLines(
    c(
      'WorldBegin',
      'Texture "color" "spectrum" "imagemap" "string filename" "map.png" "float udelta" .1 "float vdelta" .2',
      'Texture "bump" "float" "imagemap" "string filename" "map.png" "float udelta" -.3 "float vdelta" .4',
      'Texture "rough" "float" "imagemap" "string filename" "map.png" "float udelta" .5 "float vdelta" -.6',
      'Material "disney" "texture color" "color" "texture roughness" "rough" "texture bumpmap" "bump"',
      'Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0]'
    ),
    file
  )
  result = suppressWarnings(read_pbrt(
    file,
    strict = FALSE,
    asset_dir = file.path(directory, "assets")
  ))
  material = result$scene$material[[1]]
  expect_equal(material$image_offset[[1]], c(.1, .2))
  expect_equal(material$texture_offsets$bump, c(-.3, .4))
  expect_equal(material$texture_offsets$roughness, c(.5, -.6))
  expect_false(any(grepl("udelta|vdelta", result$diagnostics$message)))
})
