test_that("bump arrays preserve signed sub-byte heights and singleton axes", {
  for (dimensions in list(c(1L, 1L), c(1L, 4L), c(3L, 1L))) {
    bump = matrix(
      seq(-.0004, .0006, length.out = prod(dimensions)),
      dimensions[1],
      dimensions[2]
    )
    processed = process_scene(sphere(material = diffuse(bump_texture = bump)))
    filename = processed$scene$material[[1]]$bump_texture
    height = rayimage::ray_read_image(filename)[,, 1]
    expect_lt(abs(min(height) - min(bump)), 1e-9)
    expect_lt(abs(max(height) - max(bump)), 1e-9)
    unlink(filename)
  }
})
