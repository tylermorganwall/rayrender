compare_image = function(path1, path2) {
  image1 = png::readPNG(path1)
  image2 = png::readPNG(path2)
  return(identical(image1, image2))
}


test_that("rayshader tests", {
  testthat::expect_no_error({
    volcano %>%
      rayshader::sphere_shade() %>%
      rayshader::plot_3d(volcano,zscale = 2)
    temp_hq = tempfile(fileext = ".png")
    set.seed(1)
    rayshader::render_highquality(filename=temp_hq, samples = 100)
    rgl::close3d()
  })
  if(!isTRUE(as.logical(Sys.getenv("IS_GHA", "false")))) {
    announce_snapshot_file(temp_hq, name = "rayshader.png")
    expect_snapshot_file(temp_hq, name = "rayshader.png", compare = compare_image)
  }
})
