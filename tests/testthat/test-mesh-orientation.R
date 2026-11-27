test_that("PBRT ReverseOrientation points flat mesh emitters toward the intended side", {
  directory = tempfile()
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  file = file.path(directory, "lamp.pbrt")
  render = function(reverse, mirrored) {
    writeLines(
      c(
        'LookAt 0 0 -3 0 0 0 0 1 0 Camera "perspective" "float fov" 25',
        'Film "rgb" "integer xresolution" 16 "integer yresolution" 16',
        'WorldBegin Material "diffuse" "rgb reflectance" [0 0 0]',
        'AreaLightSource "diffuse" "rgb L" [2 2 2]',
        if (reverse) 'ReverseOrientation' else '',
        if (mirrored) 'Scale -1 1 1' else '',
        'Shape "trianglemesh" "integer indices" [0 1 2 0 2 3]',
        '"point3 P" [-1 -1 0 1 -1 0 1 1 0 -1 1 0]'
      ),
      file
    )
    imported = read_pbrt(file)
    args = imported$render_args
    args$samples = 1
    args$parallel = args$denoise = args$preview = args$progress = args$plot_scene = FALSE
    args$tonemap = "raw"
    set.seed(1)
    image = do.call(render_scene, c(list(scene = imported$scene), args))
    image[,, 1:3, drop = FALSE]
  }
  for (mirrored in c(FALSE, TRUE)) {
    expect_equal(max(render(FALSE, mirrored)), 0)
    expect_gt(mean(render(TRUE, mirrored)), 1)
  }
})
