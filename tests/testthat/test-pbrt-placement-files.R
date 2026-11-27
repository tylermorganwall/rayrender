test_that('repeated restoring placement files share a nested prototype', {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  part = file.path(root, 'part.pbrt')
  path = file.path(root, 'scene.pbrt')
  placements = sprintf(
    'AttributeBegin Translate %g %g 0 ObjectInstance "ball" AttributeEnd',
    rep(seq(-.875, .875, length.out = 8), 8),
    rep(seq(-.875, .875, length.out = 8), each = 8)
  )
  writeLines(placements, part)
  writeLines(
    c(
      'WorldBegin ObjectBegin "ball" Shape "sphere" "float radius" [.1] "float alpha" [.7] ObjectEnd',
      'AttributeBegin Translate -2 0 0 Include "part.pbrt" AttributeEnd',
      'Include "part.pbrt"',
      'AttributeBegin Translate 2 0 0 Import "part.pbrt" AttributeEnd'
    ),
    path
  )
  candidate = read_pbrt(path)$scene
  control_env = new.env(parent = environment(read_pbrt))
  control_env$pbrt_share_placement_file = function(...) FALSE
  control_env$pbrt_execute_file = pbrt_execute_file
  environment(control_env$pbrt_execute_file) = control_env
  control_env$read_pbrt = read_pbrt
  environment(control_env$read_pbrt) = control_env
  control_env$pbrt_execute = pbrt_execute
  environment(control_env$pbrt_execute) = control_env
  control_env$pbrt_import_file = pbrt_import_file
  environment(control_env$pbrt_import_file) = control_env
  control = control_env$read_pbrt(path)$scene
  expect_equal(nrow(candidate), 3L)
  expect_equal(nrow(control), 1L)
  expect_equal(
    ncol(candidate$shape_info[[1]]$shape_properties$instance_transforms),
    64L
  )
  expect_identical(
    candidate$shape_info[[2]]$shape_properties$original_scene,
    candidate$shape_info[[3]]$shape_properties$original_scene
  )
  old = options(cores = 1L)
  on.exit(options(old), add = TRUE)
  render = function(scene) {
    scene = add_infinite_light(
      scene,
      disk_light(direction = c(-1, 2, 3), intensity = 20000)
    )
    set.seed(443)
    render_scene(
      scene,
      width = 32,
      height = 16,
      samples = 8,
      max_depth = 6,
      lookfrom = c(0, 0, 10),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(7, 3.5),
      sample_method = 'random',
      min_variance = 0,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      parallel = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = 'raw'
    )
  }
  expect_equal(render(candidate), render(control), tolerance = 1e-6)
  # A non-restoring translation changes the state after every Include, so it
  # must retain ordinary execution even after the file has been seen before.
  writeLines(c('Translate 1 0 0', placements), part)
  expect_equal(read_pbrt(path)$scene, control_env$read_pbrt(path)$scene)
  writeLines(c('Identity', placements), part)
  expect_equal(read_pbrt(path)$scene, control_env$read_pbrt(path)$scene)
  writeLines(
    c('Material "diffuse" "rgb reflectance" [.2 .3 .4]', placements),
    part
  )
  expect_equal(read_pbrt(path)$scene, control_env$read_pbrt(path)$scene)
  writeLines(c('ActiveTransform EndTime Translate 1 0 0', placements), part)
  expect_equal(read_pbrt(path)$scene, control_env$read_pbrt(path)$scene)
  writeLines(placements, part)
  writeLines(
    c(
      'WorldBegin Include "part.pbrt" Import "part.pbrt"',
      'ObjectBegin "ball" Shape "sphere" ObjectEnd'
    ),
    path
  )
  expect_equal(read_pbrt(path)$scene, control_env$read_pbrt(path)$scene)
  writeLines(
    c(
      'WorldBegin ObjectBegin "ball" ObjectEnd',
      'Include "part.pbrt" Import "part.pbrt"'
    ),
    path
  )
  expect_equal(read_pbrt(path)$scene, control_env$read_pbrt(path)$scene)
  writeLines(c('AttributeBegin', placements), part)
  expect_error(read_pbrt(path), 'Unclosed|Unmatched')
})
