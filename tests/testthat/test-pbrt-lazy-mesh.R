test_that('lazy PBRT meshes retain smooth shading, UVs and transforms', {
  path = tempfile(fileext = '.pbrt')
  assets = tempfile()
  on.exit(unlink(c(path, assets), recursive = TRUE))
  old = options(cores = 1L)
  on.exit(options(old), add = TRUE)
  shape = paste(
    'Shape "trianglemesh"',
    '"point3 P" [-1 -1 0 1 -1 0 1 1 0 -1 1 0 8 8 8]',
    '"integer indices" [0 1 2 0 2 3]',
    '"normal N" [-.2 0 1 .2 0 1 .2 .2 1 -.2 .2 1 0 0 1]',
    '"point2 uv" [0 0 1 0 1 1 0 1 0 0]'
  )
  render = function(scene) {
    set.seed(827)
    scene = add_infinite_light(
      scene,
      disk_light(direction = c(-1, 2, 3), intensity = 20)
    )
    render_scene(
      scene,
      width = 32,
      height = 32,
      samples = 8,
      lookfrom = c(0, 0, 5),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(3, 3),
      sample_method = 'random',
      min_variance = 0,
      parallel = FALSE,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = 'raw'
    )
  }
  for (transform in c('', 'Scale -1 1.2 .8')) {
    writeLines(c('WorldBegin', transform, shape), path)
    memory = read_pbrt(path)$scene
    lazy = read_pbrt(path, asset_dir = assets, mesh_storage = 'ply')
    expect_identical(lazy$scene$shape, 'ply')
    expect_length(lazy$assets, 1)
    expect_true(all(file.exists(lazy$assets)))
    expect_identical(lazy$scene$transforms, memory$transforms)
    expect_identical(lazy$scene$animation_info, memory$animation_info)
    # A structured bump and alpha pattern exercises UV lookup, not just color.
    material = diffuse(
      image_texture = array(rep(c(.8, .3, .1), each = 64), c(8, 8, 3)),
      alpha_texture = matrix(rep(c(.3, 1), 32), 8, 8),
      bump_texture = diag(8),
      bump_intensity = .02
    )
    memory$material = material
    lazy$scene$material = material
    a = render(memory)
    b = render(lazy$scene)
    expect_gt(max(a[,, 1:3]), 0)
    expect_equal(as.vector(a), as.vector(b), tolerance = 1e-6)
  }
})

test_that('lazy mesh assets are owned by the import and cleaned after failure', {
  path = tempfile(fileext = '.pbrt')
  assets = tempfile()
  on.exit(unlink(c(path, assets), recursive = TRUE))
  shape = 'Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0] "integer indices" [0 1 2]'
  writeLines(
    c(
      'WorldBegin ObjectBegin "leaf"',
      rep(shape, 32),
      'ObjectEnd ObjectInstance "leaf" ObjectInstance "leaf"'
    ),
    path
  )
  result = read_pbrt(path, asset_dir = assets, mesh_storage = 'ply')
  expect_length(result$assets, 1)
  retained = result$assets
  # Finalizing an object has already written its cache before this bad reference.
  writeLines(
    c(
      'WorldBegin ObjectBegin "leaf"',
      shape,
      'ObjectEnd ObjectInstance "undefined"'
    ),
    path
  )
  expect_error(
    read_pbrt(path, asset_dir = assets, mesh_storage = 'ply'),
    'Undefined'
  )
  expect_setequal(
    list.files(assets, recursive = TRUE, full.names = TRUE),
    retained
  )
  expect_error(read_pbrt(path, mesh_storage = 'unknown'), 'arg')
})

test_that('lazy storage leaves transport boundaries and emitters in memory', {
  mesh = rayvertex::construct_mesh(
    vertices = rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0)),
    indices = matrix(0:2, 1, 3)
  )
  cache = new.env(parent = emptyenv())
  cache$directory = tempfile()
  on.exit(unlink(cache$directory, recursive = TRUE))
  for (material in list(
    subsurface(),
    dielectric(),
    light(),
    openpbr(emission_luminance = 1)
  )) {
    row = raymesh_model(mesh, material = material)
    expect_identical(pbrt_store_mesh_rows(list(row), cache), list(row))
  }
  row = set_medium(raymesh_model(mesh), homogeneous_medium(sigma_s = .1))
  expect_identical(pbrt_store_mesh_rows(list(row), cache), list(row))
  expect_false(dir.exists(cache$directory))
})
