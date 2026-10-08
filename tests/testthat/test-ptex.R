#' @param scene Scene under test.
#' @param debug Default `"color"`. Debug channel or path tracing.
#' @return Deterministic linear RGB pixels.
#' @keywords internal
#' @noRd
ptex_test_render = function(scene, debug = "color") {
  withr::local_options(list(cores = 2L))
  set.seed(827)
  scene = add_infinite_light(
    scene,
    disk_light(direction = c(-1, 2, 3), intensity = 2)
  )
  image = render_scene(
    scene,
    width = 16,
    height = 16,
    samples = 4,
    lookfrom = c(0, 0, 5),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(1.5, 1.5),
    sample_method = "random",
    min_variance = 0,
    parallel = TRUE,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = "raw",
    debug_channel = debug
  )
  image[,, 1:3, drop = FALSE]
}

test_that("Ptex descriptors validate and survive serialization without opening textures", {
  file = tempfile(fileext = ".ptx")
  file.create(file)
  on.exit(unlink(file))
  graph = texture_ptex(file, encoding = "gamma 2.2")
  expect_identical(unserialize(serialize(graph, NULL)), graph)
  expect_equal(graph$gamma, 2.2)
  expect_equal(graph$filter, 5L)
  expect_error(texture_ptex(file, encoding = "gamma"), "encoding")
  expect_error(texture_ptex(file, encoding = "gamma -1"), "encoding")
  expect_error(texture_ptex(file, encoding = NA), "encoding")
  expect_error(texture_ptex("missing.ptx"), "existing")
})

test_that("Ptex source faces survive memory, PLY, batching and instancing", {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  file = file.path(root, "faces.ptx")
  colors = list(c(.8, .1, .2), c(.15, .7, .3))
  faces = lapply(colors, function(color) {
    array(rep(color, each = 64), c(8, 8, 3))
  })
  ptex::ptex_write(file, faces)
  path = file.path(root, "scene.pbrt")
  shape = paste(
    'Shape "trianglemesh" "point3 P" [-1 -1 0 1 -1 0 1 1 0 -1 1 0]',
    '"integer indices" [0 1 2 0 2 3] "point2 uv" [0 0 1 0 1 1 0 1]'
  )
  for (face in 0:1) {
    writeLines(
      c(
        'WorldBegin Texture "paint" "spectrum" "ptex" "string filename" "faces.ptx" "string encoding" "linear"',
        'ObjectBegin "panel" Material "diffuse" "texture reflectance" "paint"',
        paste(shape, '"integer faceIndices" [', face, face, ']'),
        # Two same-material meshes batch, but must retain the per-source IDs.
        paste(shape, '"integer faceIndices" [', face, face, ']'),
        'ObjectEnd ObjectInstance "panel"'
      ),
      path
    )
    memory = read_pbrt(path, mesh_storage = "memory")
    lazy = read_pbrt(
      path,
      mesh_storage = "ply",
      asset_dir = file.path(root, paste0("assets", face))
    )
    a = ptex_test_render(memory$scene)
    b = ptex_test_render(lazy$scene)
    expect_equal(as.vector(a), as.vector(b), tolerance = 1e-6)
    expect_equal(as.numeric(a[8, 8, ]), colors[[face + 1L]], tolerance = 1e-5)
    expect_false(any(grepl("faceIndices|ptex", memory$diagnostics$message)))
  }
  writeLines(c("WorldBegin", paste(shape, '"integer faceIndices" [0]')), path)
  expect_error(read_pbrt(path), "one nonnegative integer")
  writeLines(
    c("WorldBegin", paste(shape, '"integer faceIndices" [-1 0]')),
    path
  )
  expect_error(read_pbrt(path), "one nonnegative integer")
})

test_that("Ptex imports match PBRT encoding, scalar averaging and scale", {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  color = c(.502, .302, .702)
  file = file.path(root, "color.ptx")
  ptex::ptex_write(file, list(array(rep(color, each = 16), c(4, 4, 3))))
  path = file.path(root, "scene.pbrt")
  for (type in c("spectrum", "float")) {
    writeLines(
      c(
        sprintf(
          'WorldBegin Texture "paint" "%s" "ptex" "string filename" "color.ptx" "float scale" .7',
          type
        ),
        'Material "diffuse" "texture reflectance" "paint"',
        'Shape "bilinearmesh" "point3 P" [-1 -1 0 1 -1 0 -1 1 0 1 1 0]'
      ),
      path
    )
    scene = read_pbrt(path)$scene
    expected = (floor(color * 255 + .5) / 255)^2.2 * .7
    if (type == "float") {
      expected = rep(mean(expected), 3)
    }
    expect_equal(
      as.numeric(ptex_test_render(scene)[8, 8, ]),
      expected,
      tolerance = 1e-5
    )
  }
})

test_that("Ptex preserves face-local UVs without an image v flip", {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  size = 32L
  face = array(0, c(size, size, 3))
  face[,, 1] = (row(face[,, 1]) - .5) / size
  face[,, 2] = (col(face[,, 2]) - .5) / size
  face[,, 3] = .25
  ptex::ptex_write(file.path(root, "uv.ptx"), list(face))
  path = file.path(root, "uv.pbrt")
  writeLines(
    c(
      'WorldBegin Texture "paint" "spectrum" "ptex" "string filename" "uv.ptx" "string encoding" "linear"',
      'Material "diffuse" "texture reflectance" "paint"',
      'Shape "bilinearmesh" "point3 P" [-1 -1 0 1 -1 0 -1 1 0 1 1 0]'
    ),
    path
  )
  scene = read_pbrt(path)$scene
  color = ptex_test_render(scene)
  uv = ptex_test_render(scene, "uv")
  # Both reconstruction and interpolation must preserve an affine ramp away
  # from chart edges, on both halves of the triangulated default quad chart.
  expect_equal(
    as.vector(color[4:12, 4:12, 1:2]),
    as.vector(uv[4:12, 4:12, 1:2]),
    tolerance = 1e-4
  )
  expect_equal(range(color[,, 3]), c(.25, .25), tolerance = 1e-6)
})

test_that("ray differentials filter Ptex detail at a distance", {
  file = tempfile(fileext = ".ptx")
  on.exit(unlink(file))
  face = array(0, c(64, 64, 1))
  face[,, 1] = (row(face[,, 1]) + col(face[,, 1])) %% 2
  ptex::ptex_write(file, list(face))
  indices = rbind(c(0, 1, 2), c(0, 2, 3))
  mesh = rayvertex::construct_mesh(
    vertices = rbind(c(-1, -1, 0), c(1, -1, 0), c(1, 1, 0), c(-1, 1, 0)),
    indices = indices,
    tex_indices = indices,
    texcoords = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  filtered = ptex_test_render(raymesh_model(
    mesh,
    material = diffuse(texture_ptex(file))
  ))
  point = ptex_test_render(raymesh_model(
    mesh,
    material = diffuse(texture_ptex(file, filter = "point"))
  ))
  expect_lt(max(abs(filtered[4:12, 4:12, ] - .5)), .02)
  expect_gt(diff(range(point)), .8)
})

test_that("PLY polygon triangulation retains the source Ptex face", {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  file = file.path(root, "faces.ptx")
  ptex::ptex_write(file, list(array(.1, c(4, 4, 1)), array(.7, c(4, 4, 1))))
  ply = file.path(root, "quad.ply")
  writeLines(
    c(
      'ply',
      'format ascii 1.0',
      'element vertex 4',
      'property float x',
      'property float y',
      'property float z',
      'property float u',
      'property float v',
      'element face 1',
      'property list uchar int vertex_indices',
      'property int face_indices',
      'end_header',
      '-1 -1 0 0 0',
      '1 -1 0 1 0',
      '1 1 0 1 1',
      '-1 1 0 0 1',
      '4 0 1 2 3 1'
    ),
    ply
  )
  scene = ply_model(ply, material = diffuse(texture_ptex(file)))
  image = ptex_test_render(scene)
  expect_equal(range(image), c(.7, .7), tolerance = 1e-6)
  expect_error(
    ptex_test_render(ply_model(
      ply,
      subdivision_levels = 2,
      material = diffuse(texture_ptex(file))
    )),
    "Subdivision of Ptex"
  )
})

test_that("Ptex colors reach both diffuse-transmission lobes", {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  color = c(.2, .4, .6)
  ptex::ptex_write(
    file.path(root, "leaf.ptx"),
    list(array(rep(color, each = 16), c(4, 4, 3)))
  )
  path = file.path(root, "leaf.pbrt")
  writeLines(
    c(
      'WorldBegin Texture "leaf" "spectrum" "ptex" "string filename" "leaf.ptx" "string encoding" "linear"',
      'Material "diffusetransmission" "texture reflectance" "leaf" "texture transmittance" "leaf" "float scale" .5',
      'Shape "bilinearmesh" "point3 P" [-1 -1 0 1 -1 0 -1 1 0 1 1 0]'
    ),
    path
  )
  imported = read_pbrt(path)
  expect_equal(
    as.numeric(ptex_test_render(imported$scene)[8, 8, ]),
    color,
    tolerance = 1e-6
  )
  expect_length(imported$diagnostics$message, 0)
})

test_that("Ptex bump maps affect shading equally through memory and lazy PLY", {
  root = tempfile()
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  height = array(rep(seq(0, .3, length.out = 32), 32), c(32, 32, 1))
  ptex::ptex_write(file.path(root, "height.ptx"), list(height))
  path = file.path(root, "bump.pbrt")
  writeLines(
    c(
      'WorldBegin Texture "height" "float" "ptex" "string filename" "height.ptx" "string encoding" "linear"',
      'Material "diffuse" "rgb reflectance" [.6 .6 .6] "texture displacement" "height"',
      'Shape "trianglemesh" "point3 P" [-1 -1 0 1 -1 0 1 1 0 -1 1 0]',
      '"integer indices" [0 1 2 0 2 3] "point2 uv" [0 0 1 0 1 1 0 1] "integer faceIndices" [0 0]'
    ),
    path
  )
  memory = read_pbrt(path)$scene
  lazy = read_pbrt(
    path,
    mesh_storage = "ply",
    asset_dir = file.path(root, "assets")
  )$scene
  a = ptex_test_render(memory, "none")
  b = ptex_test_render(lazy, "none")
  flat = memory
  flat$material[[1]]$texture_graphs$bump = NULL
  c = ptex_test_render(flat, "none")
  expect_true(all(is.finite(a)))
  expect_equal(as.vector(a), as.vector(b), tolerance = 1e-6)
  expect_gt(mean(abs(a - c)) / mean(abs(c)), .0001)
})

test_that("Ptex failure warnings are aggregated after worker completion", {
  file = tempfile(fileext = ".ptx")
  on.exit(unlink(file))
  ptex::ptex_write(file, list(array(.5, c(4, 4, 1))))
  mesh = rayvertex::construct_mesh(
    vertices = rbind(c(-1, -1, 0), c(1, -1, 0), c(1, 1, 0), c(-1, 1, 0)),
    indices = rbind(c(0, 1, 2), c(0, 2, 3))
  )
  mesh$shapes[[1]]$ptex_face_indices = c(10L, 10L)
  scene = raymesh_model(mesh, material = diffuse(texture_ptex(file)))
  expect_warning(
    image <- ptex_test_render(scene, "none"),
    "Ptex:.*failed and returned zero"
  )
  expect_true(all(is.finite(image)))
  expect_equal(max(image), 0)
  unlink(file)
  expect_error(ptex_test_render(scene), "Cannot open Ptex")
})
