test_that("local texture processing preserves input columns and material IDs", {
  texture = array(seq(.1, .9, length.out = 12), c(2, 2, 3))
  alpha = matrix(c(.2, .4, .6, .8), 2, 2)
  bump = matrix(c(.1, .3, .5, .7), 2, 2)
  material = microfacet(
    image_texture = texture,
    alpha_texture = alpha,
    bump_texture = bump,
    roughness_texture = texture
  )
  first = sphere(material = material)
  second = sphere(x = 2, material = material)
  first$shape_info[[1]]$material_id = 41L
  second$shape_info[[1]]$material_id = 41L
  scene = do.call(vctrs::vec_rbind, list(first, second))
  before = serialize(scene, NULL)
  processed = process_scene(scene)
  expect_identical(serialize(scene, NULL), before)
  expect_s3_class(processed$scene$material, "ray_material")
  expect_s3_class(processed$scene$shape_info, "ray_shape_info")
  paths = unlist(lapply(processed$scene$material, function(m) {
    c(m$image, m$alphaimage, m$bump_texture, m$roughness_texture)
  }))
  on.exit(unlink(paths), add = TRUE)
  expect_true(all(file.exists(paths)))
  expect_true(all(
    vapply(processed$scene$shape_info, function(x) x$material_id, integer(1)) ==
      0L
  ))
  # Both rows' input descriptors can be reused after processing independently.
  expect_identical(scene$material[[1]], scene$material[[2]])
  expect_equal(processed$typevec, rep(material[[1]]$type, 2))
})

test_that("nested instances do not shift subsequent bump texture assignments", {
  withr::local_options(list(cores = 1))
  bump = outer(seq(0, 1, length.out = 8), seq(0, 1, length.out = 8))
  background = sphere(x = -4)
  nested = create_instances(create_instances(background))
  render_bump = function(scene) {
    set.seed(930)
    render_scene(
      scene,
      width = 24,
      height = 24,
      samples = 1,
      lookfrom = c(0, 0, 5),
      lookat = c(0, 0, 0),
      fov = 35,
      parallel = FALSE,
      preview = FALSE,
      interactive = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      debug_channel = "bump",
      tonemap = "raw"
    )
  }
  # The first case previously lost its bump; the second dereferenced the null
  # buffer assigned to the preceding, off-camera object's bump texture.
  for (tail in list(
    sphere(material = diffuse(bump_texture = bump)),
    add_object(
      sphere(x = -4, material = diffuse(bump_texture = bump)),
      sphere()
    )
  )) {
    reference = render_bump(add_object(background, tail))
    actual = render_bump(add_object(nested, tail))
    expect_true(all(is.finite(actual)))
    # Instance bounds can change the last float bit of camera/BVH arithmetic.
    expect_equal(actual, reference, tolerance = 1e-6)
  }
})
