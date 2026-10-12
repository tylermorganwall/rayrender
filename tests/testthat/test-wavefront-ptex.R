#' @param scene Scene to compare without lighting variance.
#' @param backend Default `"metal"`. Integrator under test.
#' @param settings Default `list()`. Additional render settings.
#' @return Linear image with backend diagnostics.
#' @keywords internal
#' @noRd
wavefront_ptex_render = function(scene, backend = "metal", settings = list()) {
  withr::local_options(list(cores = 6L))
  set.seed(825)
  # A narrow frontal disk gives near-constant unit irradiance / pi.
  scene = add_infinite_light(
    scene,
    disk_light(
      direction = c(0, 0, 1),
      angular_diameter = 1,
      intensity = 1 / sin(pi / 360)^2
    )
  )
  image = do.call(
    render_scene,
    modifyList(
      list(
        scene = scene,
        integrator_type = backend,
        width = 32,
        height = 32,
        samples = 32,
        lookfrom = c(0, 0, 5),
        lookat = c(0, 0, 0),
        fov = 0,
        ortho_dimensions = c(1.5, 1.5),
        sample_method = "sobol",
        min_variance = 0,
        preview = FALSE,
        progress = FALSE,
        plot_scene = FALSE,
        denoise = FALSE,
        bloom = FALSE,
        tonemap = "raw",
        max_depth = 2,
        ambient_light = FALSE,
        backgroundhigh = "black",
        backgroundlow = "black"
      ),
      settings
    )
  )
  if (backend == "metal" && !isTRUE(attr(image, "wavefront")$used)) {
    expect_match(
      attr(image, "wavefront")$fallback,
      "not enabled|device|macOS 11"
    )
    skip("Metal unavailable")
  }
  image
}

#' @param face Default `0L`. Zero-based source face ID.
#' @return A quad split into two triangles sharing a Ptex face.
#' @keywords internal
#' @noRd
wavefront_ptex_panel = function(face = 0L) {
  indices = rbind(c(0, 1, 2), c(0, 2, 3))
  mesh = rayvertex::construct_mesh(
    vertices = rbind(c(-1, -1, 0), c(1, -1, 0), c(1, 1, 0), c(-1, 1, 0)),
    indices = indices,
    tex_indices = indices,
    texcoords = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  )
  mesh$shapes[[1]]$ptex_face_indices = rep(face, 2)
  mesh
}

test_that("Metal Ptex keeps source faces, encoding and compositional scalar values", {
  skip_if_not_installed("rayvertex")
  file = tempfile(fileext = ".ptx")
  on.exit(unlink(file))
  colors = list(c(.17, .43, .72), c(.81, .24, .12))
  ptexr::ptex_write(
    file,
    lapply(colors, function(x) array(rep(x, each = 64), c(8, 8, 3)))
  )
  for (face in 0:1) {
    for (encoding in c("linear", "srgb", "gamma 2.2")) {
      tex = texture_ptex(file, encoding = encoding)
      scene = raymesh_model(
        wavefront_ptex_panel(face),
        material = diffuse(tex)
      )
      gpu = wavefront_ptex_render(scene)
      cpu = wavefront_ptex_render(scene, "nee")
      expect_equal(as.vector(gpu), as.vector(cpu), tolerance = .001)
      expect_equal(attr(gpu, "wavefront")$discarded_paths, 0)
    }
  }
  tex = texture_scale(
    texture_ptex(file, type = "scalar", encoding = "gamma 2.2"),
    .7
  )
  scene = raymesh_model(
    wavefront_ptex_panel(1L),
    material = diffuse(tex)
  )
  expect_equal(
    as.vector(wavefront_ptex_render(scene)),
    as.vector(wavefront_ptex_render(scene, "nee")),
    tolerance = .001
  )
  point = raymesh_model(
    wavefront_ptex_panel(1L),
    material = diffuse(texture_ptex(file, filter = "point"))
  )
  expect_equal(
    as.vector(wavefront_ptex_render(point)),
    as.vector(wavefront_ptex_render(point, "nee")),
    tolerance = .001
  )
})

test_that("Metal Ptex preserves UV orientation, minification and atlas limits", {
  skip_if_not_installed("rayvertex")
  file = tempfile(fileext = ".ptx")
  on.exit(unlink(file))
  n = 64L
  pixels = array(0, c(n, n, 3))
  pixels[,, 1] = (row(pixels[,, 1]) - .5) / n
  pixels[,, 2] = (col(pixels[,, 2]) - .5) / n
  pixels[,, 3] = .25
  ptexr::ptex_write(file, list(pixels))
  mesh = wavefront_ptex_panel()
  scene = raymesh_model(
    mesh,
    material = diffuse(texture_ptex(file))
  )
  gpu = wavefront_ptex_render(scene)
  cpu = wavefront_ptex_render(scene, "nee")
  expect_lt(max(abs(gpu[4:28, 4:28, 1:3] - cpu[4:28, 4:28, 1:3])), .03)
  withr::local_options(list(rayrender.metal_ptex_memory = 4096))
  expect_warning(
    small <- wavefront_ptex_render(scene),
    "reduced its maximum face resolution"
  )
  expect_lt(max(abs(small[8:24, 8:24, 1:3] - cpu[8:24, 8:24, 1:3])), .025)
  withr::local_options(list(rayrender.metal_ptex_memory = 256 * 1024^2))
  pixels[,,] = (row(pixels[,, 1]) + col(pixels[,, 1])) %% 2
  unlink(file)
  ptexr::ptex_write(file, list(pixels))
  scene = raymesh_model(
    mesh,
    material = diffuse(texture_ptex(file))
  )
  filtered = wavefront_ptex_render(scene)
  expect_lt(max(abs(filtered[4:28, 4:28, 1:3] - .5)), .03)
})

test_that("Metal Ptex follows triangle charts and moving instances", {
  skip_if_not_installed("rayvertex")
  file = tempfile(fileext = ".ptx")
  on.exit(unlink(file))
  ptexr::ptex_write(file, list(array(.375, c(4, 4, 1))), mesh_type = "triangle")
  mesh = rayvertex::construct_mesh(
    vertices = rbind(c(-2, -2, 0), c(4, -2, 0), c(-2, 4, 0)),
    indices = matrix(0:2, nrow = 1),
    tex_indices = matrix(0:2, nrow = 1),
    texcoords = rbind(c(0, 0), c(1, 0), c(0, 1))
  )
  mesh$shapes[[1]]$ptex_face_indices = 0L
  panel = raymesh_model(
    mesh,
    material = diffuse(texture_ptex(file))
  )
  scene = create_instances(panel, x = 0) |>
    animate_objects(end_position = c(.2, .1, 0), end_angle = c(0, 5, 0))
  gpu = wavefront_ptex_render(scene, settings = list(shutter_speed = 1))
  cpu = wavefront_ptex_render(scene, "nee", list(shutter_speed = 1))
  expect_equal(as.vector(gpu), as.vector(cpu), tolerance = .001)
  expect_gt(attr(gpu, "wavefront")$mesh_instances, 0)
})

test_that("Metal Ptex rejects an atlas that cannot fit without corrupting the scene", {
  skip_if_not_installed("rayvertex")
  # Check device availability before exercising an intentional fallback.
  wavefront_ptex_render(raymesh_model(
    wavefront_ptex_panel(),
    material = diffuse()
  ))
  file = tempfile(fileext = ".ptx")
  on.exit(unlink(file))
  ptexr::ptex_write(file, rep(list(array(.4, c(1, 1, 1))), 16))
  withr::local_options(list(cores = 6L, rayrender.metal_ptex_memory = 1024))
  scene = raymesh_model(
    wavefront_ptex_panel(),
    material = diffuse(texture_ptex(file))
  )
  expect_warning(
    image <- render_scene(
      scene,
      integrator_type = "metal",
      width = 8,
      height = 8,
      samples = 1,
      lookfrom = c(0, 0, 5),
      lookat = c(0, 0, 0),
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      denoise = FALSE,
      bloom = FALSE
    ),
    "Ptex face metadata exceeds the atlas memory budget"
  )
  expect_false(attr(image, "wavefront")$used)
  expect_true(all(is.finite(image)))
})
