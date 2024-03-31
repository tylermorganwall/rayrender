compare_image = function(path1, path2) {
  image1 = png::readPNG(path1)
  image2 = png::readPNG(path2)
  return(identical(image1, image2))
}


test_that("rayshader tests", {
  testthat::skip_if_not_installed("rayshader")
  set.seed(1)
  
  volcano %>%
    rayshader::sphere_shade() %>%
    rayshader::plot_3d(volcano,zscale = 2)
  temp_hq = tempfile(fileext = ".png")
  rayshader::render_highquality(filename=temp_hq, samples = 100)
  rgl::close3d()
  expect_snapshot_file(temp_hq, compare = compare_image)
})
