test_that('zero alpha preserves hidden area illumination but not ordinary geometry', {
  path = tempfile(fileext = '.pbrt')
  on.exit(unlink(path))
  floor = 'Shape "trianglemesh" "point3 P" [-2 0 -2 2 0 -2 2 0 2 -2 0 2] "integer indices" [0 2 1 0 3 2]'
  emitter = 'Shape "trianglemesh" "point3 P" [-1 2 -1 1 2 -1 1 2 1 -1 2 1] "integer indices" [0 1 2 0 2 3]'
  prefix = c(
    'WorldBegin Material "diffuse" "rgb reflectance" [.7 .7 .7]',
    floor,
    'AreaLightSource "diffuse" "rgb L" [20 20 20]'
  )
  writeLines(c(prefix, paste(emitter, '"float alpha" [0]')), path)
  imported = suppressWarnings(read_pbrt(path, strict = FALSE))
  scene = imported$scene
  expect_equal(nrow(scene), 2)
  expect_equal(scene$material[[2]]$properties[[1]][4], 1)
  expect_equal(scene$material[[2]]$alpha_value, 1)
  expect_true(any(grepl(
    'Constant-zero-alpha area emitter',
    imported$diagnostics$message
  )))
  writeLines(
    c(
      'WorldBegin Texture "hidden" "float" "constant" "float value" [0]',
      prefix[-1],
      paste(emitter, '"texture alpha" "hidden"')
    ),
    path
  )
  named = suppressWarnings(read_pbrt(path, strict = FALSE))$scene
  expect_equal(named$material[[2]]$properties[[1]][4], 1)
  writeLines(c(prefix, paste(emitter, '"float alpha" [-.1]')), path)
  expect_equal(nrow(suppressWarnings(read_pbrt(path, strict = FALSE))$scene), 1)
  old = options(cores = 1L)
  on.exit(options(old), add = TRUE)
  render = function(x, below = FALSE) {
    set.seed(482)
    image = render_scene(
      x,
      width = 24,
      height = 24,
      samples = 16,
      lookfrom = if (below) c(0, 0, 0) else c(0, 4, 4),
      lookat = if (below) c(0, 2, 0) else c(0, 0, 0),
      camera_up = if (below) c(0, 0, 1) else c(0, 1, 0),
      fov = 45,
      backgroundhigh = 'black',
      backgroundlow = 'black',
      ambient_light = FALSE,
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
    image[,, 1:3, drop = FALSE]
  }
  expect_gt(mean(render(scene)), .1)
  expect_lt(max(render(scene[1, , drop = FALSE])), 1e-12)
  hidden = scene[2, , drop = FALSE]
  expect_lt(max(render(hidden, below = TRUE)), 1e-12)
  hidden$material[[1]]$properties[[1]][4] = 0
  # See the emitter's emitting side directly for the visible-light control.
  visible = render(hidden, below = TRUE)
  expect_gt(max(visible), 1)
})
