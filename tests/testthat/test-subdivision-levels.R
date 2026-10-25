test_that('mesh subdivision counts are validated before loading geometry', {
  for (constructor in list(obj_model, ply_model, mesh3d_model, raymesh_model)) {
    expect_identical(formals(constructor)$subdivision_levels, 0)
    for (invalid in list(
      -1,
      .5,
      NA_real_,
      Inf,
      NaN,
      numeric(),
      c(0, 1),
      '1',
      TRUE
    )) {
      expect_error(
        do.call(constructor, list(subdivision_levels = invalid)),
        'subdivision_levels must be a single nonnegative integer'
      )
    }
  }
  expect_identical(validate_subdivision_levels(0L), 0)
})

test_that('all mesh loaders perform exactly zero, one or two subdivision passes', {
  skip_if_not_installed('rayvertex')
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  vertices = rbind(c(-1, -1, 0), c(1, -1, 0), c(0, 1, 0))
  indices = matrix(0:2, 1, 3)
  mesh = rayvertex::construct_mesh(vertices, indices)
  rgl_mesh = structure(
    list(vb = rbind(t(vertices), 1), it = t(indices + 1L)),
    class = c('mesh3d', 'shape3d')
  )
  obj = file.path(directory, 'triangle.obj')
  ply = file.path(directory, 'triangle.ply')
  writeLines(c('v -1 -1 0', 'v 1 -1 0', 'v 0 1 0', 'f 1 2 3'), obj)
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
      '-1 -1 0',
      '1 -1 0',
      '0 1 0',
      '3 0 1 2'
    ),
    ply
  )
  constructors = list(
    obj = function(...) obj_model(obj, load_material = FALSE, ...),
    ply = function(...) ply_model(ply, ...),
    mesh3d = function(...) mesh3d_model(rgl_mesh, ...),
    raymesh = function(...) raymesh_model(mesh, ...)
  )
  render = function(scene, backend = 'nee') {
    set.seed(417)
    suppressWarnings(render_scene(
      scene,
      integrator_type = backend,
      width = 32,
      height = 32,
      samples = 4,
      lookfrom = c(0, 0, 3),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(2.5, 2.5),
      aperture = 0,
      max_depth = 2,
      ambient_light = FALSE,
      parallel = FALSE,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      min_variance = 0
    ))
  }
  for (constructor in constructors) {
    expect_identical(
      constructor(material = light()),
      constructor(subdivision_levels = 0, material = light())
    )
  }
  for (level in 0:2) {
    reference_mesh = rayvertex::subdivide_mesh(mesh, subdivision_levels = level)
    expect_equal(nrow(reference_mesh$shapes[[1]]$indices), 4^level)
    reference = render(raymesh_model(reference_mesh, material = light()))
    for (name in names(constructors)) {
      scene = constructors[[name]](
        subdivision_levels = level,
        material = light()
      )
      cpu = render(scene)
      expect_equal(
        as.numeric(cpu),
        as.numeric(reference),
        tolerance = 1e-6,
        info = paste(name, 'level', level)
      )
      expect_true(all(is.finite(cpu)))
    }
  }
  # The displacement importer previously rounded a single refinement up to two.
  pbrt = file.path(directory, 'displaced.pbrt')
  writeLines(
    c(
      'WorldBegin',
      sprintf('Shape "plymesh" "string filename" "%s"', ply),
      '"float displacement" .25 "float edgelength" 2'
    ),
    pbrt
  )
  displaced = suppressWarnings(read_pbrt(pbrt, strict = FALSE))
  displaced_mesh = displaced$scene$shape_info[[1]]$mesh_info[[1]]
  expect_equal(nrow(displaced_mesh$shapes[[1]]$indices), 4)
  expect_equal(displaced_mesh$vertices[[1]][, 3], rep(.25, 6))
  # Count the actual triangles handed to Metal, rather than just the R settings.
  for (level in 0:2) {
    for (constructor in constructors) {
      gpu = render(
        constructor(subdivision_levels = level, material = light()),
        'metal'
      )
      if (!isTRUE(attr(gpu, 'wavefront')$used)) {
        expect_match(
          attr(gpu, 'wavefront')$fallback,
          '^(Metal was not enabled|this Metal device does not support|Metal wavefront requires macOS 11)'
        )
        skip('Metal unavailable; CPU subdivision comparisons completed')
      }
      expect_equal(attr(gpu, 'wavefront')$triangles, 4^level)
    }
  }
})

test_that('PBRT preserves exact subdivision counts and keeps refined meshes out of lazy storage', {
  path = tempfile(fileext = '.pbrt')
  assets = tempfile()
  on.exit(unlink(c(path, assets), recursive = TRUE))
  for (level in 0:2) {
    writeLines(
      c(
        'WorldBegin',
        paste0('Shape "loopsubdiv" "integer levels" [', level, ']'),
        '"point3 P" [-1 -1 0 1 -1 0 0 1 0] "integer indices" [0 1 2]'
      ),
      path
    )
    imported = suppressWarnings(read_pbrt(
      path,
      mesh_storage = 'ply',
      asset_dir = assets,
      strict = FALSE
    ))
    expect_equal(
      imported$scene$shape_info[[1]]$shape_properties$subdivision_levels,
      level
    )
    expect_identical(imported$scene$shape, if (level == 0) 'ply' else 'raymesh')
  }
})
