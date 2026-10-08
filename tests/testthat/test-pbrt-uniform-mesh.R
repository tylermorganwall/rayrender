test_that('uniform PBRT mesh coverage matches the public constructor exactly', {
  context = new.env(parent = emptyenv())
  context$mesh_material = rayvertex::material_list()
  vertices = matrix(c(0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0), 4, 3, byrow = TRUE)
  indices = matrix(c(0L, 1L, 2L, 1L, 3L, 2L), 2, 3, byrow = TRUE)
  normals = matrix(rep(c(0, 0, 1), 4), 4, 3, byrow = TRUE)
  uv = vertices[, 1:2]
  saved = NULL
  for (has_normals in c(TRUE, FALSE)) {
    for (has_uv in c(TRUE, FALSE)) {
      n = if (has_normals) normals else NULL
      t = if (has_uv) uv else NULL
      candidate = pbrt_uniform_mesh(vertices, indices, n, t, context)
      control = rayvertex::construct_mesh(
        vertices,
        indices,
        n,
        if (has_normals) indices else NULL,
        t,
        if (has_uv) indices else NULL,
        material = context$mesh_material
      )
      expect_identical(candidate, control)
      if (is.null(saved)) saved = candidate
    }
  }
  expect_identical(saved$normals[[1]], normals)
  expect_identical(saved$texcoords[[1]], uv)
  expect_equal(nrow(context$mesh_template$vertices[[1]]), 0)
})

test_that('uniform construction preserves subdivision, patches, and raw output', {
  path = tempfile(fileext = '.pbrt')
  on.exit(unlink(path))
  control_env = new.env(parent = environment(read_pbrt))
  control_env$pbrt_uniform_mesh = function(
    vertices,
    indices,
    normals,
    uv,
    context
  ) {
    rayvertex::construct_mesh(
      vertices,
      indices,
      normals,
      if (!is.null(normals)) indices else NULL,
      uv,
      if (!is.null(uv)) indices else NULL,
      material = context$mesh_material
    )
  }
  for (name in c(
    'read_pbrt',
    'pbrt_execute_file',
    'pbrt_execute',
    'pbrt_shape',
    'pbrt_mesh'
  )) {
    fn = get(name, envir = environment(read_pbrt))
    environment(fn) = control_env
    assign(name, fn, control_env)
  }
  points = '"point3 P" [0 0 0 1 0 0 0 1 0 1 1 0]'
  shapes = c(
    paste('Shape "bilinearmesh"', points, '"integer indices" [0 1 2 3]'),
    paste(
      'Shape "loopsubdiv"',
      points,
      '"integer indices" [0 1 2 1 3 2] "integer levels" [1]'
    ),
    paste(
      'Shape "loopsubdiv"',
      points,
      '"integer indices" [0 1 2 1 3 2] "integer levels" [2]'
    )
  )
  for (shape in shapes) {
    writeLines(
      c('WorldBegin Material "diffuse" "rgb reflectance" [.8 .3 .1]', shape),
      path
    )
    candidate = suppressWarnings(read_pbrt(path, strict = FALSE)$scene)
    control = suppressWarnings(
      control_env$read_pbrt(path, strict = FALSE)$scene
    )
    expect_identical(candidate, control)
  }
  old = options(cores = 1L)
  on.exit(options(old), add = TRUE)
  render = function(scene) {
    scene = add_infinite_light(
      scene,
      disk_light(direction = c(-1, 2, 3), intensity = 20000)
    )
    set.seed(318)
    render_scene(
      scene,
      width = 24,
      height = 24,
      samples = 4,
      lookfrom = c(.5, .5, 4),
      lookat = c(.5, .5, 0),
      fov = 0,
      ortho_dimensions = c(1.4, 1.4),
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
  im = render(candidate)
  expect_gt(max(im[,, 1:3]), .1)
  expect_identical(im, render(control))
})
