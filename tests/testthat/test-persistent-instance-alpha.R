test_that("instance preparation retains usable alpha PNGs", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  file = file.path(directory, "alpha.png")
  for (coverage in list(matrix(c(0, .25, .75, 1), 2), matrix(1, 2, 2))) {
    mask = array(1, c(2, 2, 4))
    mask[,, 4] = coverage
    png::writePNG(mask, file)
    prepared = process_scene(sphere(material = diffuse(alpha_texture = file)))
    expect_identical(prepared$scene$material[[1]]$alphaimage, file)
    expect_equal(png::readPNG(file)[,, 4], coverage, tolerance = 1 / 255)
  }
  # Preserve the legacy convention for an opaque RGBA grayscale mask.
  mask[,, 1:3] = .25
  png::writePNG(mask, file)
  prepared = process_scene(sphere(material = diffuse(alpha_texture = file)))
  converted = prepared$scene$material[[1]]$alphaimage
  on.exit(unlink(converted), add = TRUE)
  expect_false(identical(converted, file))
  expect_equal(
    png::readPNG(converted)[,, 4],
    matrix(.25, 2, 2),
    tolerance = 1 / 255
  )
  expect_equal(png::readPNG(file)[,, 4], matrix(1, 2, 2))

  color_file = file.path(directory, "color.png")
  mask[,, 4] = .5
  png::writePNG(mask, color_file)
  before = tools::md5sum(file)
  prepared = process_scene(sphere(
    material = diffuse(
      image_texture = color_file,
      alpha_texture = file
    )
  ))
  expect_identical(prepared$scene$material[[1]]$alphaimage, color_file)
  expect_identical(tools::md5sum(file), before)
})

test_that("serialized nested PBRT instances render with persistent alpha assets", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  mask = array(.5, c(2, 2, 3))
  png::writePNG(mask, file.path(directory, "mask.png"))
  file = file.path(directory, "scene.pbrt")
  writeLines(
    c(
      'WorldBegin',
      'Texture "mask" "float" "imagemap" "string filename" "mask.png" "string encoding" "linear"',
      'ObjectBegin "leaf"',
      'Shape "trianglemesh" "integer indices" [0 1 2] "point3 P" [-1 -1 0 1 -1 0 0 1 0] "texture alpha" "mask"',
      'ObjectEnd',
      'ObjectBegin "branch" ObjectInstance "leaf" ObjectEnd',
      'ObjectInstance "branch"'
    ),
    file
  )
  imported = read_pbrt(file, asset_dir = file.path(directory, "assets"))
  branch = imported$scene$shape_info[[1]]$shape_properties$original_scene[[1]]
  leaf = branch$shape_info[[1]]$shape_properties$original_scene[[1]]
  alpha = leaf$material[[1]]$alphaimage
  expect_true(alpha %in% imported$assets)
  expect_true(file.exists(alpha))
  # Simulate deletion of any session-generated mask before another R process
  # loads the scene. The persistent files declared by the importer remain.
  unlink(setdiff(alpha, imported$assets))
  cache = file.path(directory, "scene.rds")
  saveRDS(imported, cache)
  script = file.path(directory, "render.R")
  writeLines(
    c(
      'library(rayrender)',
      'scene = readRDS(commandArgs(TRUE)[1])$scene',
      'image = render_scene(scene, lookfrom = c(0,0,3), lookat = c(0,0,0),',
      '  width = 8, height = 8, samples = 1, parallel = FALSE,',
      '  preview = FALSE, plot = FALSE, progress = FALSE, denoise = FALSE)',
      'stopifnot(all(is.finite(image)))'
    ),
    script
  )
  output = system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script), shQuote(cache)),
    stdout = TRUE,
    stderr = TRUE,
    env = paste0(
      "R_LIBS=",
      shQuote(paste(.libPaths(), collapse = .Platform$path.sep))
    )
  )
  expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
})
