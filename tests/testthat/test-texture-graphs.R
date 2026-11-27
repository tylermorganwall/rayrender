#' @param scene Test scene.
#' @param x Default `0`. Camera x translation.
#' @param debug Default `"color"`. Renderer diagnostic mode.
#' @return Small deterministic linear image.
texture_test_render = function(scene, x = 0, debug = "color") {
  set.seed(52)
  image = render_scene(
    scene,
    width = 12,
    height = 12,
    samples = 16,
    sample_method = "random",
    lookfrom = c(x, 0, 4),
    lookat = c(x, 0, 0),
    fov = 0,
    ortho_dimensions = c(1, 1),
    ambient_light = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white",
    min_variance = 0,
    clamp_value = Inf,
    tonemap = "raw",
    denoise = FALSE,
    bloom = FALSE,
    debug_channel = debug,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE
  )
  image[,, 1:3, drop = FALSE]
}

test_that("texture graphs compose, validate types, and survive serialization", {
  noise = texture_noise(scale = 12, seed = 9)
  color = texture_mix("red", "blue", noise)
  rough = texture_mix(.02, .3, noise)
  expect_identical(color$type, "color")
  expect_identical(rough$type, "scalar")
  expect_identical(unserialize(serialize(color, NULL)), color)
  expect_equal(texture_mix(.1, .5, .25)$value, .2)
  expect_equal(texture_scale(.2, 3)$value, .6)
  expect_error(texture_mix(0, 1, "red"), "scalar")
  expect_error(texture_direction_mix(0, 1, c(0, 0, 0)), "nonzero")
  expect_error(texture_noise(octaves = 0), "octave")
  expect_error(texture_coordinates(offset = c(0, NA, 0)), "coordinate")
  expect_error(texture_image("does-not-exist.png"), "existing")
  expect_error(microfacet(roughness = color), "scalar")
  expect_error(microfacet(roughness = texture_constant("red")), "scalar")
  expect_identical(
    texture_coordinates("uv", scale = c(2, 3), offset = c(.2, .3)),
    texture_coordinates("uv", scale = c(2, 3, 1), offset = c(.2, .3, 0))
  )
  expect_error(diffuse(color, checkercolor = "red"), "Choose")
  expect_error(
    openpbr(base_color = color, image_texture = "other.png"),
    "Choose"
  )
  expect_error(openpbr(coat_weight = noise), "not supported")
  expect_identical(diffuse(color)[[1]]$texture_graphs$color, color)
  expect_identical(
    microfacet(roughness = rough)[[1]]$texture_graphs$roughness,
    rough
  )
  expect_identical(
    openpbr(base_color = color, specular_roughness = rough)[[1]]$texture_graphs,
    list(color = color, roughness = rough)
  )
})

test_that("PBRT directionmix preserves its transform, nesting, and roughness remap", {
  path = tempfile(fileext = ".pbrt")
  on.exit(unlink(path))
  writeLines(
    c(
      'WorldBegin AttributeBegin Rotate 90 0 0 1',
      'Texture "a" "float" "constant" "float value" .0005',
      'Texture "mix" "float" "directionmix" "texture tex1" "a" "float tex2" .005 "vector3 dir" [0 2 0]',
      'AttributeEnd',
      'Texture "scaled" "float" "scale" "texture tex" "mix" "float scale" 2',
      'Material "dielectric" "texture roughness" "scaled"',
      'Shape "sphere"',
      'Translate 3 0 0 Material "dielectric" "texture roughness" "mix" "bool remaproughness" false',
      'Shape "sphere"'
    ),
    path
  )
  imported = read_pbrt(path)
  expect_equal(nrow(imported$diagnostics), 0)
  graphs = imported$scene$material[[1]]$texture_graphs
  expect_true(graphs$roughness_is_alpha)
  expect_equal(graphs$roughness$exponent, .5)
  expect_equal(graphs$roughness$child$factor$value, 2)
  expect_equal(
    graphs$roughness$child$child$weight$direction,
    c(-1, 0, 0),
    tolerance = 1e-7
  )
  expect_equal(
    imported$scene$material[[2]]$texture_graphs$roughness$exponent,
    1
  )
})

test_that("graph images decode before filtering and preserve independent UV mappings", {
  file = tempfile(fileext = ".png")
  on.exit(unlink(file))
  png::writePNG(array(.5, c(2, 2, 3)), file)
  color = texture_test_render(xy_rect(material = diffuse(texture_image(file))))
  scalar = texture_test_render(xy_rect(
    material = diffuse(texture_image(file, type = "scalar"))
  ))
  byte = 128 / 255
  expect_equal(mean(color), ((byte + .055) / 1.055)^2.4, tolerance = .002)
  expect_equal(mean(scalar), byte, tolerance = .002)
  black_white = array(0, c(2, 2, 3))
  black_white[, 2, ] = 1
  png::writePNG(black_white, file)
  mix = texture_image(
    file,
    coordinates = texture_coordinates("uv", scale = 0, offset = c(.5, .5, 0))
  )
  expect_equal(
    mean(texture_test_render(xy_rect(material = diffuse(mix)))),
    .5,
    tolerance = .002
  )
})

test_that("world and object patterns render correctly on translated geometry and instances", {
  object = texture_gradient(
    0,
    1,
    coordinates = texture_coordinates("object"),
    axis = "x"
  )
  world = texture_gradient(
    0,
    1,
    coordinates = texture_coordinates("world"),
    axis = "x"
  )
  base = xy_rect(material = diffuse(object))
  moved = xy_rect(x = 2, material = diffuse(object))
  expect_equal(
    texture_test_render(base),
    texture_test_render(moved, x = 2),
    tolerance = 1e-6
  )
  expect_equal(
    mean(texture_test_render(xy_rect(x = 2, material = diffuse(world)), x = 2)),
    1,
    tolerance = 1e-6
  )
  instanced = create_instances(base, x = 2)
  expect_equal(
    texture_test_render(base),
    texture_test_render(instanced, x = 2),
    tolerance = 1e-6
  )
})

test_that("roughness graphs agree with numeric inputs in rendered microfacet and OpenPBR", {
  # A nonconstant graph with equal leaves exercises native graph evaluation.
  roughness = texture_direction_mix(.3, .3)
  for (constructor in list(
    function(r) {
      microfacet(roughness = r, eta = c(.2, .9, 1.1), kappa = c(3.9, 2.4, 2.1))
    },
    function(r) openpbr(specular_roughness = r, base_metalness = 1)
  )) {
    a = texture_test_render(sphere(material = constructor(.3)), debug = "none")
    b = texture_test_render(
      sphere(material = constructor(roughness)),
      debug = "none"
    )
    expect_equal(b, a, tolerance = 1e-5)
    expect_true(all(is.finite(b)))
  }
})
