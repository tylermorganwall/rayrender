test_that('static restoring matrix blocks pack without changing placement semantics', {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  scene_file = file.path(root, 'scene.pbrt')
  part = file.path(root, 'part.pbrt')
  matrices = lapply(seq_len(130), function(i) {
    m = diag(4)
    m[1, 1] = if (i %% 2) -1 else 1
    m[1, 2] = .13
    m[1:3, 4] = c((i %% 13) * .3, (i %/% 13) * .3, 0)
    m
  })
  blocks = vapply(
    seq_along(matrices),
    function(i) {
      sprintf(
        'AttributeBegin ConcatTransform [%s] ObjectInstance "%s" AttributeEnd',
        paste(as.vector(matrices[[i]]), collapse = ' '),
        if (i <= 128) 'ball' else 'small'
      )
    },
    character(1)
  )
  prefix = c(
    'WorldBegin',
    'ObjectBegin "ball" Shape "sphere" "float radius" [.1] ObjectEnd',
    'ObjectBegin "small" Shape "sphere" "float radius" [.05] ObjectEnd'
  )
  writeLines(blocks, part)
  writeLines(
    c(
      prefix,
      'Scale 1.2 .8 1 Include "part.pbrt"',
      'Shape "sphere" "float radius" [.12]'
    ),
    scene_file
  )
  control_env = new.env(parent = environment(read_pbrt))
  control_env$pbrt_pack_placement_blocks = function(...) FALSE
  for (name in c(
    'read_pbrt',
    'pbrt_execute_file',
    'pbrt_execute',
    'pbrt_import_file',
    'pbrt_share_placement_file'
  )) {
    fn = get(name, envir = environment(read_pbrt))
    environment(fn) = control_env
    assign(name, fn, envir = control_env)
  }
  candidate = read_pbrt(scene_file)$scene
  control = control_env$read_pbrt(scene_file)$scene
  expect_equal(candidate, control)
  old = options(cores = 1L)
  on.exit(options(old), add = TRUE)
  render = function(scene) {
    scene = add_infinite_light(
      scene,
      disk_light(direction = c(-1, 2, 3), intensity = 20000)
    )
    set.seed(614)
    render_scene(
      scene,
      width = 32,
      height = 32,
      samples = 4,
      lookfrom = c(2, 1.5, 10),
      lookat = c(2, 1.5, 0),
      fov = 0,
      ortho_dimensions = c(5, 5),
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
  candidate_image = render(candidate)
  expect_gt(max(candidate_image[,, 1:3]), .1)
  expect_identical(candidate_image, render(control))
  # Cached repeated includes still share their nested prototype.
  writeLines(
    c(
      prefix,
      'Include "part.pbrt"',
      'Translate 5 0 0 Import "part.pbrt"',
      'Translate 5 0 0 Include "part.pbrt"'
    ),
    scene_file
  )
  expect_equal(
    read_pbrt(scene_file)$scene,
    control_env$read_pbrt(scene_file)$scene
  )
  # Forward definitions, animation, and non-restoring files use the executor.
  writeLines(c('WorldBegin Include "part.pbrt"', prefix[-1]), scene_file)
  expect_equal(
    read_pbrt(scene_file)$scene,
    control_env$read_pbrt(scene_file)$scene
  )
  writeLines(
    c(
      prefix,
      'ActiveTransform EndTime Translate 1 0 0 ActiveTransform All',
      'Include "part.pbrt"'
    ),
    scene_file
  )
  expect_equal(
    read_pbrt(scene_file)$scene,
    control_env$read_pbrt(scene_file)$scene
  )
  writeLines(c(prefix, 'Include "part.pbrt"'), scene_file)
  writeLines(c(blocks, 'Translate 1 0 0'), part)
  expect_equal(
    read_pbrt(scene_file)$scene,
    control_env$read_pbrt(scene_file)$scene
  )
  invalid = matrices[[1]]
  invalid[1, ] = 0
  broken = blocks
  broken[70] = sprintf(
    'AttributeBegin ConcatTransform [%s] ObjectInstance "ball" AttributeEnd',
    paste(as.vector(invalid), collapse = ' ')
  )
  writeLines(broken, part)
  expect_error(read_pbrt(scene_file), 'part.pbrt.*70.*nonsingular')
  expect_error(control_env$read_pbrt(scene_file), 'part.pbrt.*70.*nonsingular')
})
