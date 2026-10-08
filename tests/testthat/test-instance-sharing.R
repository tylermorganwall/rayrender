test_that("shared prototypes retain independent transforms, lights and media", {
  options_before = options(cores = 1L)
  on.exit(options(options_before))
  texture = array(1, c(8, 8, 4))
  texture[,, 1] = diag(8)
  texture[,, 2] = 0.2
  texture[1:3, , 4] = 0
  prototypes = list(
    sphere(radius = .4, material = diffuse('steelblue')),
    sphere(
      radius = .4,
      material = diffuse(
        image_texture = texture,
        bump_texture = diag(8),
        bump_intensity = .05
      )
    ),
    sphere(radius = .25, material = light('orange', intensity = 3)),
    create_instances(
      sphere(radius = .16, material = diffuse('tomato')),
      y = c(-.2, .2)
    ),
    set_medium(
      sphere(radius = .4),
      homogeneous_medium(sigma_a = c(.2, .5, .1), sigma_s = .3)
    )
  )
  for (prototype in prototypes) {
    instance = create_instances(prototype)
    shared = add_object(
      group_objects(instance, translate = c(-.65, 0, 0)),
      group_objects(instance, translate = c(.65, 0, 0))
    )
    shared$transforms[[2]]$group_transform[[1]][1, 1] = -1.2
    shared$transforms[[2]]$group_transform[[1]][1, 2] = .15
    packed = instance
    packed$shape_info[[1]]$shape_properties$instance_transforms =
      vapply(
        shared$transforms,
        function(x) as.vector(x$group_transform[[1]]),
        numeric(16)
      )
    separate = shared
    for (i in 1:2) {
      separate$shape_info[[i]]$shape_properties$original_scene =
        unserialize(serialize(
          shared$shape_info[[i]]$shape_properties$original_scene,
          NULL
        ))
    }
    ground = xz_rect(
      y = -.45,
      xwidth = 6,
      zwidth = 6,
      material = diffuse('white')
    )
    settings = list(
      width = 24,
      height = 16,
      samples = 4,
      lookfrom = c(0, 1, 4),
      lookat = c(0, 0, 0),
      fov = 35,
      sample_method = 'random',
      min_variance = 0,
      max_depth = 8,
      parallel = FALSE,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = 'raw'
    )
    shared = add_infinite_light(
      add_object(shared, ground),
      disk_light(direction = c(-1, 2, 3), intensity = 2)
    )
    separate = add_infinite_light(
      add_object(separate, ground),
      disk_light(direction = c(-1, 2, 3), intensity = 2)
    )
    packed = add_infinite_light(
      add_object(packed, ground),
      disk_light(direction = c(-1, 2, 3), intensity = 2)
    )
    set.seed(812)
    actual = do.call(render_scene, c(list(scene = shared), settings))
    set.seed(812)
    expected = do.call(render_scene, c(list(scene = separate), settings))
    expect_true(all(is.finite(actual)))
    expect_gt(max(actual), 0)
    expect_equal(actual, expected, tolerance = 1e-6)
    set.seed(812)
    compact = do.call(render_scene, c(list(scene = packed), settings))
    expect_equal(compact, expected, tolerance = 1e-6)
  }
})

test_that("large static PBRT populations pack affine matrices and retain animated uses", {
  file = tempfile(fileext = '.pbrt')
  on.exit(unlink(file))
  placements = sprintf(
    'AttributeBegin Translate %g 0 0 ObjectInstance "part" AttributeEnd',
    1:80
  )
  writeLines(
    c(
      'WorldBegin ObjectBegin "part" Shape "sphere" ObjectEnd',
      placements,
      'ActiveTransform EndTime Translate 1 0 0 ObjectInstance "part"'
    ),
    file
  )
  scene = read_pbrt(file)$scene
  expect_equal(nrow(scene), 2)
  matrices = scene$shape_info[[1]]$shape_properties$instance_transforms
  expect_equal(dim(matrices), c(16L, 80L))
  expect_equal(matrices[13, ], 1:80)
  expect_equal(scene$animation_info[[2]]$end_transform_animation[[1]][1, 4], 1)
  expect_null(scene$shape_info[[2]]$shape_properties$instance_transforms)
})
