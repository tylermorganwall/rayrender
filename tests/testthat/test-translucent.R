test_that("translucent sheets validate colors and retain separate textures", {
  x = translucent(.25, c(.6, .5, .3), image_offset = c(.2, -.1))[[1]]
  expect_equal(x$type, get_material_enum("translucent"))
  expect_equal(x$properties[[1]], rep(.25, 3))
  expect_equal(x$transmittance, c(.6, .5, .3))
  expect_equal(x$transmission_offset, c(.2, -.1))
  expect_error(translucent(.8, .8), "must not exceed")
  expect_error(translucent(numeric(), .5), "finite scalar")
  expect_error(translucent(.2, .5, image_repeat = Inf), "finite")
})

test_that("PBRT ISO, coated diffuse defaults, and diffuse transmission are retained", {
  file = tempfile(fileext = ".pbrt")
  on.exit(unlink(file))
  writeLines(
    c(
      'Film "rgb" "float iso" 200',
      'WorldBegin',
      'Material "coateddiffuse" Shape "sphere"',
      'Material "coateddiffuse" "float roughness" .1 Shape "sphere"',
      'Material "diffusetransmission" "rgb transmittance" [.6 .5 .3] Shape "sphere"'
    ),
    file
  )
  x = suppressWarnings(read_pbrt(file, strict = FALSE))
  expect_equal(x$render_args$iso, 200)
  expect_equal(x$scene$material[[1]]$openpbr$coat_roughness, 0)
  expect_equal(x$scene$material[[2]]$openpbr$coat_roughness, .1^.25)
  expect_equal(x$scene$material[[3]]$transmittance, c(.6, .5, .3))
  expect_equal(x$scene$material[[3]]$properties[[1]], rep(.25, 3))
  expect_false(any(grepl(
    "iso|Unsupported material|transmittance",
    x$diagnostics$message
  )))
  writeLines('Film "rgb" "float iso" -1 WorldBegin Shape "sphere"', file)
  expect_error(read_pbrt(file), "iso must be")
})

test_that("PBRT sheet textures keep independent UVs and radiometric scaling", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  png::writePNG(
    array(rep(c(.1, .2, .3), each = 4), c(2, 2, 3)),
    file.path(directory, "color.png")
  )
  file = file.path(directory, "sheet.pbrt")
  writeLines(
    c(
      'WorldBegin',
      'Texture "r" "spectrum" "imagemap" "string filename" "color.png" "string encoding" "linear" "float uscale" 2 "float vscale" 3 "float udelta" .1',
      'Texture "t" "spectrum" "imagemap" "string filename" "color.png" "string encoding" "linear" "float uscale" 4 "float vscale" 5 "float vdelta" .2',
      'Material "diffusetransmission" "texture reflectance" "r" "texture transmittance" "t" "float scale" 2',
      'Shape "trianglemesh" "integer indices" [0 1 2] "point3 P" [0 0 0 1 0 0 0 1 0] "point2 uv" [0 0 1 0 0 1]'
    ),
    file
  )
  x = read_pbrt(file)$scene$material[[1]]
  expect_equal(x$image_repeat[[1]], c(2, 3))
  expect_equal(x$image_offset[[1]], c(.1, 0))
  expect_equal(x$transmission_repeat, c(4, 5))
  expect_equal(x$transmission_offset, c(0, .2))
  a = rayimage::ray_read_image(x$image)
  b = rayimage::ray_read_image(x$transmission_texture)
  expect_equal(as.numeric(a[1, 1, 1:3]), c(.2, .4, .6), tolerance = .01)
  expect_equal(a, b)
  writeLines(
    c(
      'WorldBegin',
      'Texture "t" "spectrum" "imagemap" "string filename" "color.png" "string encoding" "linear"',
      'Material "diffusetransmission" "texture transmittance" "t" Shape "sphere"'
    ),
    file
  )
  expect_error(read_pbrt(file), "native UV convention")
})

test_that("a backlit sheet transmits diffuse light with next-event estimation", {
  file = tempfile(fileext = ".pbrt")
  on.exit(unlink(file))
  writeLines(
    c(
      'LookAt 0 0 -3 0 0 0 0 1 0 Camera "perspective" "float fov" 20',
      'Film "rgb" "integer xresolution" 16 "integer yresolution" 16',
      'Integrator "path" "integer maxdepth" 2 WorldBegin',
      'LightSource "point" "point3 from" [0 0 2] "rgb I" [4 4 4]',
      'Material "diffusetransmission" "rgb reflectance" [0 0 0] "rgb transmittance" [.8 .4 .2]',
      'Shape "trianglemesh" "integer indices" [0 1 2 0 2 3]',
      '"point3 P" [-2 -2 0 2 -2 0 2 2 0 -2 2 0]'
    ),
    file
  )
  x = read_pbrt(file)
  args = x$render_args
  args$samples = 16L
  args$preview = args$progress = args$plot_scene = args$denoise = args$parallel = FALSE
  args$tonemap = "raw"
  set.seed(1)
  image = do.call(render_scene, c(list(scene = x$scene), args))[,, 1:3]
  expect_gt(mean(image[,, 1]), .15)
  expect_equal(mean(image[,, 1]) / mean(image[,, 2]), 2, tolerance = 1e-4)
  expect_equal(mean(image[,, 2]) / mean(image[,, 3]), 2, tolerance = 1e-4)
  args$iso = 200
  set.seed(1)
  exposed = do.call(render_scene, c(list(scene = x$scene), args))[,, 1:3]
  expect_equal(exposed, 2 * image, tolerance = 1e-6)
  # Array images follow the existing PNG color-texture decoding convention.
  texture = array(rep(c(.8, .4, .2)^(1 / 2.2), each = 4), c(2, 2, 3))
  x$scene$material[[1]] = translucent(
    reflectance = 0,
    transmittance = 0,
    transmission_texture = texture
  )[[1]]
  args$iso = 100
  set.seed(1)
  textured = do.call(render_scene, c(list(scene = x$scene), args))[,, 1:3]
  expect_lt(max(abs(textured - image)), .003)
  x$scene$material[[1]] = diffuse(color = "black")[[1]]
  set.seed(1)
  blocked = do.call(render_scene, c(list(scene = x$scene), args))[,, 1:3]
  expect_equal(max(blocked), 0)
})
