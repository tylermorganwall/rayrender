test_that('PBRT meshes share defaults without changing indices or material overrides', {
  path = tempfile(fileext = '.pbrt')
  on.exit(unlink(path))
  shape = 'Shape "trianglemesh" "point3 P" [0 0 0 1 0 0 0 1 0] "integer indices" [0 1 2] "normal N" [0 0 1 0 0 1 0 0 1] "point2 uv" [0 0 1 0 0 1]'
  writeLines(
    c(
      'WorldBegin Material "diffuse" "rgb reflectance" [.8 .1 .2]',
      shape,
      'Translate 2 0 0',
      shape
    ),
    path
  )
  scene = read_pbrt(path)$scene
  a = scene$shape_info[[1]]$mesh_info[[1]]
  b = scene$shape_info[[2]]$mesh_info[[1]]
  expect_type(a$shapes[[1]]$indices, 'integer')
  expect_equal(as.vector(a$shapes[[1]]$indices), 0:2)
  expect_identical(a$materials, b$materials)
  original = b$materials
  a$materials[[1]][[1]]$diffuse = c(0, 0, 0)
  expect_identical(b$materials, original)
  expect_equal(scene$material[[1]], scene$material[[2]])
})
