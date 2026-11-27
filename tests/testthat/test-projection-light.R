test_that("projection descriptors decode once and validate image and camera inputs", {
  pixels = array(.5, c(2, 4, 3))
  lamp = projection_light(pixels)
  expect_equal(lamp$image, pixels)
  expect_equal(
    projection_light(pixels, encoding = "srgb")$image,
    array(((.5 + .055) / 1.055)^2.4, dim(pixels)),
    tolerance = 1e-6
  )
  expect_error(projection_light(pixels, direction = c(0, 0, 0)), "nonzero")
  expect_error(projection_light(pixels, up = c(0, 0, 1)), "nonparallel")
  expect_error(projection_light(pixels, fov = 180), "fov")
  expect_error(projection_light(pixels, intensity = -1), "nonnegative")
  expect_error(projection_light(array(-1, c(1, 1, 3))), "nonnegative")
  expect_error(projection_light(array(NA_real_, c(1, 1, 3))), "finite")
  expect_error(projection_light("nonexistent-projection.png"), "existing")
  file = tempfile(fileext = ".png")
  png::writePNG(pixels, file)
  decoded = projection_light(file)
  unlink(file)
  scene = add_light(xy_rect(), decoded, name = "slide")
  expect_equal(get_light(scene, "slide")$image, decoded$image)
  expect_equal(
    list_lights(add_object(scene, sphere()))$slide$type,
    "projection"
  )
  expect_length(list_lights(remove_light(scene, "slide")), 0)
})

#' @param light One light or a list of lights to test against a white diffuse screen.
#' @param backend Integrator name.
#' @param at Default `c(0, 0, 0)`. Screen sample position.
#' @param objects Default `NULL`. Additional occluding or absorbing objects.
#' @param samples Default `8`. Samples per pixel for light-selection convergence.
#' @return Mean linear RGB of a very small orthographic patch.
#' @keywords internal
projection_test_pixel = function(
  light,
  backend,
  at = c(0, 0, 0),
  objects = NULL,
  samples = 8
) {
  scene = xy_rect(xwidth = 100, ywidth = 100, material = diffuse("white"))
  lights = if (inherits(light, "ray_light")) list(light) else light
  for (item in lights) {
    scene = add_light(scene, item)
  }
  if (!is.null(objects)) {
    scene = add_object(scene, objects)
  }
  set.seed(128)
  image = suppressWarnings(render_scene(
    scene,
    integrator_type = backend,
    width = 4,
    height = 4,
    samples = samples,
    sample_method = "sobol",
    lookfrom = at + c(0, 0, .1),
    lookat = at,
    fov = 0,
    ortho_dimensions = c(1e-5, 1e-5),
    aperture = 0,
    max_depth = 2,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    parallel = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    tonemap = "raw",
    min_variance = 0
  ))
  if (backend == "metal" && !isTRUE(attr(image, "wavefront")$used)) {
    reason = attr(image, "wavefront")$fallback
    if (grepl("not enabled|device|macOS 11", reason)) {
      skip("Metal unavailable")
    }
    stop("Unexpected projection fallback: ", reason)
  }
  apply(image[,, 1:3, drop = FALSE], 3, mean)
}

for (backend in c("nee", "metal")) {
  test_that(
    paste(
      backend,
      "projection has correct orientation, clipping, and attenuation"
    ),
    {
      skip_on_cran()
      slide = array(0, c(2, 2, 3))
      slide[1, 1, ] = c(1, 0, 0)
      slide[1, 2, ] = c(0, 1, 0)
      slide[2, 1, ] = c(0, 0, 1)
      slide[2, 2, ] = 1
      lamp = projection_light(
        slide,
        position = c(0, 0, 2),
        fov = 90,
        intensity = 4
      )
      points = rbind(
        c(-1, 1, 0),
        c(1, 1, 0),
        c(-1, -1, 0),
        c(1, -1, 0),
        c(0, 0, 0)
      )
      colors = rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), rep(1, 3), rep(.5, 3))
      for (i in seq_len(nrow(points))) {
        distance = sqrt(sum((points[i, ] - lamp$position)^2))
        expected = colors[i, ] * 4 * (2 / distance) / (pi * distance^2)
        expect_equal(
          projection_test_pixel(lamp, backend, points[i, ]),
          expected,
          tolerance = 2e-4
        )
      }
      expect_equal(
        projection_test_pixel(lamp, backend, c(2.1, 0, 0)),
        rep(0, 3),
        tolerance = 1e-8
      )
      turned = projection_light(
        slide,
        position = c(0, 0, 2),
        direction = c(0, 0, 1),
        fov = 90
      )
      expect_equal(
        projection_test_pixel(turned, backend),
        rep(0, 3),
        tolerance = 1e-8
      )
      wide = projection_light(
        array(1, c(2, 4, 3)),
        position = c(0, 0, 2),
        fov = 90,
        intensity = 4
      )
      expect_gt(min(projection_test_pixel(wide, backend, c(3, 0, 0))), 0)
      expect_equal(
        projection_test_pixel(wide, backend, c(0, 2.1, 0)),
        rep(0, 3),
        tolerance = 1e-8
      )
      rotated = lamp
      rotated$up = c(1, 0, 0)
      expect_equal(
        projection_test_pixel(rotated, backend, c(1, 1, 0)),
        colors[1, ] * 8 / (pi * 6^1.5),
        tolerance = 2e-4
      )
      farther = lamp
      farther$position = c(0, 0, 4)
      reference = projection_test_pixel(lamp, backend)
      # Eight emitters exercise Metal's finite-light BVH rather than its
      # small-scene uniform proposal. Use more samples for the light-selection
      # mixture (which can also select the black environment on Metal).
      lamps = lapply(seq_len(8), function(i) {
        item = lamp
        item$name = paste0("slide-", i)
        item
      })
      expect_equal(
        projection_test_pixel(lamps, backend, samples = 256),
        8 * reference,
        tolerance = .005
      )
      dark = lamp
      dark$image[] = 0
      dark$name = "dark"
      expect_equal(
        projection_test_pixel(c(lamps, list(dark)), backend, samples = 256),
        8 * reference,
        tolerance = .005
      )
      expect_equal(
        reference / projection_test_pixel(farther, backend),
        rep(4, 3),
        tolerance = 2e-4
      )
      blocker = cube(z = 1, width = .5, material = diffuse("black"))
      expect_equal(
        projection_test_pixel(lamp, backend, objects = blocker),
        rep(0, 3),
        tolerance = 1e-8
      )
      medium = set_medium(
        blocker,
        homogeneous_medium(sigma_a = c(.2, .5, 1), sigma_s = 0)
      )
      expect_equal(
        projection_test_pixel(lamp, backend, objects = medium) / reference,
        exp(-c(.2, .5, 1) * .5),
        tolerance = 2e-4
      )
      cosine = (cos(pi / 3) + cos(pi / 6)) / 2
      spot = spot_light(
        position = c(0, 0, 2),
        direction = c(sqrt(1 - cosine^2), 0, -cosine),
        cone_angle = 60,
        falloff_angle = 30,
        intensity = 4
      )
      expect_equal(
        projection_test_pixel(spot, backend),
        rep(.5 / pi, 3),
        tolerance = 2e-4
      )
    }
  )
}

test_that("PBRT projection import preserves aspect, transform, and v3/v4 image conventions", {
  directory = tempfile("projection-pbrt-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  png::writePNG(array(1, c(2, 4, 3)), file.path(directory, "slide.png"))
  file = file.path(directory, "scene.pbrt")
  for (legacy in c(FALSE, TRUE)) {
    writeLines(
      c(
        'WorldBegin',
        'AttributeBegin',
        'Translate 1 2 3',
        'Scale -2 3 1',
        sprintf(
          'LightSource "projection" "string %s" "slide.png" "float scale" 2 "float fov" 60',
          if (legacy) "mapname" else "filename"
        ),
        'AttributeEnd',
        'Shape "sphere"'
      ),
      file
    )
    imported = read_pbrt(file, strict = TRUE)
    light = list_lights(imported$scene)[[1]]
    expect_equal(light$type, "projection")
    expect_equal(light$position, c(1, 2, 3))
    expect_equal(dim(light$image), c(2L, 4L, 3L))
    expect_equal(light$intensity, 2)
    expect_equal(
      light$projection_basis,
      diag(c(-.5, if (legacy) -1 / 3 else 1 / 3, 1))
    )
  }
})
