test_that("PBRT masks retain independent UVs, cache images, and stay shape-local", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  png::writePNG(
    array(rep(c(.2, .4, .6), each = 4), c(2, 2, 3)),
    file.path(directory, "mask.png")
  )
  file = file.path(directory, "scene.pbrt")
  triangle = 'Shape "trianglemesh" "integer indices" [0 1 2] "point3 P" [0 0 0 1 0 0 0 1 0]'
  writeLines(
    c(
      'WorldBegin',
      'Texture "mask" "float" "imagemap" "string filename" "mask.png" "string encoding" "linear" "float uscale" 2 "float vscale" 1.8 "float udelta" .2 "float vdelta" -.12',
      'MakeNamedMaterial "shared" "string type" "diffuse" NamedMaterial "shared"',
      paste(triangle, '"texture alpha" "mask"'),
      paste(triangle, '"texture alpha" "mask"'),
      triangle,
      paste(triangle, '"float alpha" .37'),
      paste(triangle, '"float alpha" 0')
    ),
    file
  )
  imported = read_pbrt(file)
  expect_equal(nrow(imported$diagnostics), 0)
  expect_equal(nrow(imported$scene), 4)
  materials = imported$scene$material
  expect_identical(materials[[1]]$alphaimage, materials[[2]]$alphaimage)
  expect_equal(materials[[1]]$alpha_repeat, c(2, 1.8))
  expect_equal(materials[[1]]$texture_offsets$alpha, c(.2, -.12))
  expect_equal(materials[[1]]$image_repeat[[1]], c(1, 1))
  expect_equal(materials[[3]]$alphaimage, "")
  expect_null(materials[[3]]$alpha_value)
  expect_equal(materials[[4]]$alpha_value, .37)
  expect_equal(
    png::readPNG(materials[[1]]$alphaimage)[,, 4],
    matrix(.4, 2, 2),
    tolerance = 1 / 255
  )
})

test_that("PBRT float alpha maps decode encoding, scale, and embedded coverage", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  rgba = array(0, c(2, 2, 4))
  rgba[,, 4] = matrix(c(.2, .4, .6, .8), 2, 2)
  png::writePNG(rgba, file.path(directory, "rgba.png"))
  png::writePNG(array(.5, c(2, 2, 3)), file.path(directory, "gray.png"))
  file = file.path(directory, "scene.pbrt")
  triangle = 'Shape "trianglemesh" "integer indices" [0 1 2] "point3 P" [0 0 0 1 0 0 0 1 0]'
  writeLines(
    c(
      'WorldBegin',
      'Texture "rgba" "float" "imagemap" "string filename" "rgba.png" "float scale" .5',
      'Texture "gray" "float" "imagemap" "string filename" "gray.png"',
      'Texture "linear" "float" "imagemap" "string filename" "gray.png" "string encoding" "linear"',
      paste(triangle, '"texture alpha" "rgba"'),
      paste(triangle, '"texture alpha" "gray"'),
      paste(triangle, '"texture alpha" "linear"')
    ),
    file
  )
  imported = read_pbrt(file)
  coverage = vapply(
    imported$scene$material,
    function(mat) png::readPNG(mat$alphaimage)[1, 1, 4],
    numeric(1)
  )
  expect_lt(max(abs(coverage - c(.1, .214, .5))), .006)
  expect_lte(
    max(abs(
      png::readPNG(imported$scene$material[[1]]$alphaimage)[,, 4] -
        rgba[,, 4] * .5
    )),
    1 / 255
  )
})

test_that("a PBRT cutout transmits camera and shadow rays with mapped UVs", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  mask = array(1, c(16, 16, 3))
  mask[5:12, 5:12, ] = 0
  png::writePNG(mask, file.path(directory, "mask.png"))
  file = file.path(directory, "scene.pbrt")
  writeLines(
    c(
      'LookAt 0 0 -3 0 0 0 0 1 0 Camera "orthographic" "float screenwindow" [-1 1 -1 1]',
      'Film "rgb" "integer xresolution" 32 "integer yresolution" 32',
      'Integrator "path" "integer maxdepth" 2 WorldBegin',
      'LightSource "point" "point3 from" [0 0 -2] "rgb I" [40 40 40]',
      'Texture "mask" "float" "imagemap" "string filename" "mask.png" "string encoding" "linear"',
      'Material "diffuse" "rgb reflectance" [0 0 0]',
      'Shape "trianglemesh" "integer indices" [0 2 1 0 3 2] "point3 P" [-1 -1 0 1 -1 0 1 1 0 -1 1 0] "point2 uv" [0 0 1 0 1 1 0 1] "texture alpha" "mask"',
      'Material "diffuse" "rgb reflectance" [1 1 1]',
      'Shape "trianglemesh" "integer indices" [0 2 1 0 3 2] "point3 P" [-2 -2 1 2 -2 1 2 2 1 -2 2 1]'
    ),
    file
  )
  imported = read_pbrt(file)
  args = imported$render_args
  args$samples = 8L
  args$preview = args$progress = args$plot_scene = args$denoise = args$parallel = FALSE
  args$tonemap = "raw"
  open = do.call(render_scene, c(list(scene = imported$scene), args))
  expect_gt(mean(open[15:18, 15:18, 1:3]), .7)
  expect_lt(mean(open[3:6, 3:6, 1:3]), .01)
  # A transparent first hit must not discard opaque surfaces behind it inside
  # the same object instance. Exercise both camera and shadow continuation.
  text = readLines(file)
  first_material = which(startsWith(text, 'Material'))[1]
  writeLines(
    c(
      text[seq_len(first_material - 1L)],
      'ObjectBegin "layers"',
      text[first_material:length(text)],
      'ObjectEnd ObjectInstance "layers"'
    ),
    file
  )
  instanced = read_pbrt(file)
  image = do.call(render_scene, c(list(scene = instanced$scene), args))
  expect_gt(mean(image[15:18, 15:18, 1:3]), .7)
  expect_lt(mean(image[3:6, 3:6, 1:3]), .01)
  imported$scene$material[[1]]$texture_offsets$alpha = c(.5, 0)
  shifted = do.call(render_scene, c(list(scene = imported$scene), args))
  expect_lt(mean(shifted[15:18, 15:18, 1:3]), .01)
})
