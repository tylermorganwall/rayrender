test_that("all texture-capable materials validate and retain UV offsets", {
  constructors = list(
    diffuse,
    metal,
    dielectric,
    microfacet,
    light,
    glossy,
    openpbr
  )
  for (material in constructors) {
    expect_equal(material()[[1]]$image_offset[[1]], c(0, 0))
    expect_equal(
      material(image_offset = c(-.2, .35))[[1]]$image_offset[[1]],
      c(-.2, .35)
    )
    for (invalid in list(0, c(1, 2, 3), c(NA, 0), c(Inf, 0), c("a", "b"))) {
      expect_error(material(image_offset = invalid), "two finite numeric")
    }
  }
  expect_equal(
    suppressWarnings(lambertian(image_offset = c(.1, .2)))[[1]]$image_offset[[
      1
    ]],
    c(.1, .2)
  )
})

test_that("native scene construction passes material offsets into image lookup", {
  pixels = array(0, c(8, 8, 3))
  pixels[, 1:4, 1] = 1
  pixels[, 5:8, 3] = 1
  filename = tempfile(fileext = ".png")
  on.exit(unlink(filename))
  png::writePNG(pixels, filename)
  render = function(offset) {
    set.seed(17)
    render_scene(
      xy_rect(
        material = light(image_texture = filename, image_offset = offset)
      ),
      width = 16,
      height = 16,
      samples = 1,
      lookfrom = c(0, 0, 3),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(2, 2),
      ambient_light = FALSE,
      backgroundhigh = "black",
      backgroundlow = "black",
      denoise = FALSE,
      tonemap = "raw",
      parallel = FALSE,
      plot_scene = FALSE,
      progress = FALSE,
      preview = FALSE
    )
  }
  base = render(c(0, 0))
  expect_equal(render(c(2, -3)), base)
  expect_gt(max(abs(render(c(.5, 0)) - base)), .1)
})
