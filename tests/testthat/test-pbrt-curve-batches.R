curve_batch_fixture = function(material = diffuse('orange')) {
  lapply(seq_len(20L), function(i) {
    x = (i - 10.5) * .12
    bezier_curve(
      p1 = rbind(
        c(x, 0, 0),
        c(x - .2, .45, .15),
        c(x + .2, .95, -.1),
        c(x, 1.4, 0)
      ),
      width = .06 + i * .001,
      width_end = .025,
      u_min = .02,
      u_max = .98,
      split_depth = i %% 3,
      type = c('flat', 'cylinder', 'ribbon')[(i - 1L) %% 3L + 1L],
      normal = c(0, 0, 1),
      normal_end = c(.1, 0, 1),
      material = material
    )
  })
}

test_that('curve batches preserve every descriptor and common row field', {
  rows = curve_batch_fixture()
  packed = pbrt_batch_curve_rows(rows)
  expect_length(packed, 1L)
  values = packed[[1]]$shape_info[[1]]$shape_properties$curve_data
  expected = vapply(
    rows,
    function(x) {
      as.numeric(unlist(x$shape_info[[1]]$shape_properties, use.names = FALSE))
    },
    numeric(24)
  )
  expect_identical(values, expected)
  common = packed[[1]]
  common$shape_info[[1]]$shape_properties = list()
  for (row in rows) {
    row$shape_info[[1]]$shape_properties = list()
    expect_identical(common, row)
  }
  chains = list(vctrs::list_unchop(rows[1:10]), vctrs::list_unchop(rows[11:20]))
  expect_identical(pbrt_batch_curve_rows(chains), packed)
  expect_identical(pbrt_batch_curve_rows(packed), packed)
  expect_identical(pbrt_batch_curve_rows(rows[1:15]), rows[1:15])
  shifted = lapply(rows, group_objects, translate = c(1, 0, 0))
  expect_length(pbrt_batch_curve_rows(c(rows, shifted)), 2L)
  # A multi-row fragment with different common fields must not inherit its first row.
  heterogeneous = vctrs::list_unchop(c(rows[1:10], shifted[11:20]))
  expect_identical(
    pbrt_batch_curve_rows(list(heterogeneous)),
    list(heterogeneous)
  )
})

test_that('curve batching leaves emitters and medium-related rows separate', {
  for (material in list(
    light(),
    openpbr(emission_luminance = 2),
    diffuse(fog = TRUE),
    subsurface()
  )) {
    rows = curve_batch_fixture(material)
    expect_identical(pbrt_batch_curve_rows(rows), rows)
  }
  rows = curve_batch_fixture()
  for (i in seq_along(rows)) {
    rows[[i]]$shape_info[[1]]$medium = homogeneous_medium()
    rows[[i]]$shape_info[[1]]$medium_keep_surface = TRUE
  }
  expect_identical(pbrt_batch_curve_rows(rows), rows)
})

test_that('packed curves preserve raw transport and animation behind priority glass', {
  old_options = options(cores = 1L)
  on.exit(options(old_options))
  settings = list(
    width = 40L,
    height = 32L,
    samples = 8L,
    lookfrom = c(0, .8, 5),
    lookat = c(0, .7, 0),
    fov = 35,
    aperture = 0,
    sample_method = 'random',
    min_variance = 0,
    max_depth = 12L,
    parallel = FALSE,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = 'raw'
  )
  for (variant in c(
    'diffuse',
    'openpbr',
    'translucent',
    'hair',
    'glass',
    'animated'
  )) {
    material = switch(
      variant,
      openpbr = openpbr(base_color = 'orange', specular_roughness = .4),
      translucent = translucent(
        reflectance = c(.3, .15, .03),
        transmittance = c(.4, .2, .05)
      ),
      hair = hair(color = 'orange'),
      diffuse('orange')
    )
    rows = curve_batch_fixture(material)
    if (variant == 'animated') {
      rows = lapply(
        rows,
        animate_objects,
        end_position = c(.2, 0, 0),
        end_angle = c(0, 10, 0)
      )
    }
    rows = lapply(rows, group_objects, angle = c(0, 10, 0), scale = 1.1)
    packed = pbrt_batch_curve_rows(rows)
    expect_equal(length(packed), 1L, info = variant)
    before = vctrs::list_unchop(rows)
    after = packed[[1]]
    if (variant == 'glass') {
      glass = cube(
        y = .75,
        xwidth = 3.5,
        ywidth = 2,
        zwidth = 1,
        material = dielectric(refraction = 1.3, priority = 5)
      )
      before = add_object(before, glass)
      after = add_object(after, glass)
    }
    illumination = disk_light(
      direction = c(-1, 2, 3),
      angular_diameter = 45,
      intensity = 20
    )
    before = add_infinite_light(before, illumination)
    after = add_infinite_light(after, illumination)
    set.seed(20261004)
    reference = do.call(render_scene, c(list(scene = before), settings))
    set.seed(20261004)
    actual = do.call(render_scene, c(list(scene = after), settings))
    expect_true(all(is.finite(actual)), info = variant)
    expect_true(max(actual[,, 1:3]) > 0, info = variant)
    # Compare every image channel, excluding nondeterministic build-time attributes.
    expect_equal(actual[,,], reference[,,], tolerance = 1e-6, info = variant)
  }
})

test_that('native packed curve input validates its schema', {
  row = pbrt_batch_curve_rows(curve_batch_fixture())[[1]]
  settings = list(
    width = 4L,
    height = 4L,
    samples = 1L,
    parallel = FALSE,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE
  )
  data = row$shape_info[[1]]$shape_properties$curve_data
  invalid = list(data[-1, ], data)
  invalid[[2]][1, 1] = NA_real_
  for (index in c(13L, 15L, 17L, 18L)) {
    value = data
    value[index, 1] = -1
    invalid[[length(invalid) + 1L]] = value
  }
  value = data
  value[18, 1] = 3
  value[19:24, 1] = 0
  invalid[[length(invalid) + 1L]] = value
  for (value in invalid) {
    bad = row
    bad$shape_info[[1]]$shape_properties$curve_data = value
    expect_error(
      do.call(render_scene, c(list(scene = bad), settings)),
      '[Pp]acked'
    )
  }
  bad = set_scene_material(row, light())
  expect_error(
    do.call(render_scene, c(list(scene = bad), settings)),
    'importance-sampled'
  )
})
