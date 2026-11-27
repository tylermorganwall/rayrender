#' @param scene Scene to render.
#' @param ... Render arguments overriding the small deterministic test defaults.
#' @return A raw rendered image.
diffusion_test_render = function(scene, ...) {
  args = list(
    scene = scene,
    width = 6L,
    height = 6L,
    samples = 128L,
    sample_method = "random",
    lookfrom = c(0, 0, 4),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.5, .5),
    aperture = 0,
    ambient_light = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white",
    min_variance = 0,
    clamp_value = Inf,
    tonemap = "raw",
    denoise = FALSE,
    bloom = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE
  )
  set.seed(192709)
  do.call(render_scene, modifyList(args, list(...)))
}

test_that("diffusion descriptors validate and retain automatic medium ownership", {
  material = subsurface_diffusion(
    color = c(.3, .7, 1),
    radius = c(.1, .2, .3),
    scale = 2,
    priority = 2
  )
  expect_s3_class(material, "ray_material")
  medium = material[[1]]$subsurface
  expect_identical(medium$subsurface$method, "diffusion")
  expect_equal(medium$subsurface$color, c(.3, .7, 1))
  expect_equal(medium$subsurface$radius, c(.2, .4, .6))
  expect_equal(medium$sigma_a + medium$sigma_s, rep(0, 3))
  expect_false(medium$haze)
  expect_equal(material[[1]]$properties[[1]][8], 2)
  for (name in c("radius", "scale", "refraction")) {
    for (value in list(0, -1, Inf, NA_real_, numeric(), "white")) {
      args = setNames(list(value), name)
      expect_error(do.call(subsurface_diffusion, args))
    }
  }
  expect_error(subsurface_diffusion(color = c(1, 2, 3)), "color")
  expect_error(subsurface_diffusion(roughness = .00001), "roughness")
  expect_error(subsurface_diffusion(priority = .5), "priority")
  ready = prepare_subsurface(sphere(material = material))
  expect_identical(ready$shape_info[[1]]$medium$subsurface$method, "diffusion")
  expect_identical(prepare_subsurface(ready), ready)
  expect_error(
    set_medium(sphere(material = material), homogeneous_medium()),
    "conflicts"
  )
})

test_that("diffusion normalizes a thick planar furnace and keeps RGB reflectance", {
  for (ior in c(1, 1.333)) {
    image = diffusion_test_render(
      cube(
        xwidth = 100,
        ywidth = 100,
        zwidth = 100,
        z = -50,
        material = subsurface_diffusion(
          color = c(.25, .6, 1),
          radius = .05,
          refraction = ior
        )
      ),
      samples = 1024L
    )
    fresnel = ((ior - 1) / (ior + 1))^2
    expected = fresnel + (1 - fresnel) * c(.25, .6, 1)
    expect_equal(apply(image[,, 1:3], 3, mean), expected, tolerance = .035)
    expect_null(attr(image, "path_warnings"))
  }
})

test_that("winning glass hides diffusion and a lower-priority glass cannot hide it", {
  glass = sphere(
    radius = 1.2,
    material = dielectric(refraction = 1, priority = 0)
  )
  milk = sphere(
    material = subsurface_diffusion(
      color = "black",
      radius = .1,
      refraction = 1,
      priority = 1
    )
  )
  hidden = diffusion_test_render(add_object(glass, milk), samples = 8L)
  expect_equal(as.numeric(hidden[,, 1:3]), rep(1, 108), tolerance = 1e-6)
  glass = set_scene_material(glass, dielectric(refraction = 1, priority = 2))
  visible = diffusion_test_render(add_object(glass, milk), samples = 8L)
  expect_lt(max(visible[,, 1:3]), 1e-6)
})

test_that("overlapped liquid and a geometrically clipped liquid agree at glass contact", {
  glass = cube(
    z = .6,
    zwidth = 1,
    xwidth = 5,
    ywidth = 5,
    material = dielectric(refraction = 1.5, priority = 0)
  )
  material = subsurface_diffusion(
    color = c(.6, .8, 1),
    radius = .1,
    refraction = 1.333,
    priority = 1
  )
  # Glass ends at z=.1; one milk box extends into it, the other ends at that face.
  overlap = add_object(
    glass,
    cube(z = -.35, zwidth = 1.5, xwidth = 4, ywidth = 4, material = material)
  )
  a = diffusion_test_render(overlap, samples = 1024L)
  expect_null(attr(a, "path_warnings"))
  for (extension in c(0, 1e-7, 1e-5)) {
    clipped = add_object(
      glass,
      cube(
        z = -.5 + extension / 2,
        zwidth = 1.2 + extension,
        xwidth = 4,
        ywidth = 4,
        material = material
      )
    )
    b = diffusion_test_render(clipped, samples = 1024L)
    expect_equal(
      apply(a[,, 1:3], 3, mean),
      apply(b[,, 1:3], 3, mean),
      tolerance = .035
    )
    expect_null(attr(b, "path_warnings"))
  }
})

test_that("diffusion has opaque coverage, finite rough boundaries and explicit camera-inside handling", {
  for (roughness in c(0, .3)) {
    body = sphere(
      material = subsurface_diffusion(radius = .04, roughness = roughness)
    )
    image = diffusion_test_render(
      body,
      samples = 16L,
      transparent_background = TRUE
    )
    expect_true(all(is.finite(image)))
    expect_equal(as.numeric(image[,, 4]), rep(1, 36))
    expect_null(attr(image, "path_warnings"))
  }
})

test_that("active diffusion interiors terminate rays and warn once without aborting", {
  withr::local_envvar(RAYRENDER_DEBUG_PATHS = "false")
  withr::local_options(warn = 2, cores = 2L)
  body = sphere(material = subsurface_diffusion(radius = .04, priority = 1))
  for (parallel in c(FALSE, TRUE)) {
    output = capture.output(
      image <- diffusion_test_render(
        body,
        samples = 4L,
        lookfrom = c(0, 0, 0),
        lookat = c(0, 0, 1),
        parallel = parallel
      )
    )
    expect_equal(sum(grepl("Warning: rays starting inside", output)), 1L)
    expect_true(all(is.finite(image)))
    expect_equal(as.numeric(image[,, 1:3]), rep(0, 6 * 6 * 3))
    diagnostics = attr(image, "path_warnings")
    expect_equal(diagnostics$terminated_paths, 6 * 6 * 4)
    expect_equal(unname(diagnostics$counts["diffusion_interior"]), 6 * 6 * 4)
    expect_length(diagnostics$examples$diffusion_interior, 8L)
    expect_match(
      diagnostics$examples$diffusion_interior[1],
      "stage=normalized_diffusion"
    )
  }
  # Ordinary outside views remain usable after a camera-inside render.
  outside = diffusion_test_render(body, samples = 4L)
  expect_true(all(is.finite(outside)))
  expect_null(attr(outside, "path_warnings"))
  expect_gt(max(outside[,, 1:3]), 0)

  # Geometric containment alone is insufficient: higher-priority glass wins.
  covered = add_object(
    body,
    sphere(radius = 1.2, material = dielectric(refraction = 1, priority = 0))
  )
  output = capture.output(
    image <- diffusion_test_render(
      covered,
      samples = 4L,
      lookfrom = c(0, 0, 0),
      lookat = c(0, 0, 1)
    )
  )
  expect_false(any(grepl("Warning: rays starting inside", output)))
  expect_null(attr(image, "path_warnings"))
  expect_gt(max(image[,, 1:3]), 0)
})

test_that("instanced diffusion stays with its object placement", {
  body = sphere(radius = .3, material = subsurface_diffusion(radius = .03))
  copies = create_instances(body, x = c(-.4, .4), scale_y = 1.3) |>
    create_instances(y = c(-.5, .5))
  image = diffusion_test_render(
    copies,
    samples = 16L,
    ortho_dimensions = c(2, 2)
  )
  expect_true(all(is.finite(image)))
  expect_gt(mean(image[,, 1:3]), .5)
  expect_null(attr(image, "path_warnings"))
})


test_that("unresolvably distant profile samples miss without nonfinite rays", {
  image = diffusion_test_render(
    sphere(
      material = subsurface_diffusion(
        radius = 1e300,
        refraction = 1
      )
    ),
    samples = 8L
  )
  expect_true(all(is.finite(image)))
  expect_lt(max(image[,, 1:3]), 1e-6)
  expect_null(attr(image, "path_warnings"))
})
