test_that("BVH selectors are validated and deterministic renders agree", {
  scene = generate_ground(depth = -1) |>
    add_object(sphere(-.7, material = diffuse('tomato'))) |>
    add_object(sphere(.8, .2, .1, radius = .7, material = diffuse('steelblue')))
  render = function(method) {
    set.seed(814)
    render_scene(
      scene,
      bvh_type = method,
      width = 32,
      height = 32,
      samples = 4,
      lookfrom = c(0, 2, 7),
      lookat = c(0, 0, 0),
      aperture = 0,
      mode = 'image',
      denoise = FALSE,
      plot_scene = FALSE,
      progress = FALSE
    )
  }
  expect_error(render('not-a-builder'), 'arg')
  reference = render('sah')
  expect_equal(unname(render('hlbvh')), unname(reference), tolerance = 1e-6)
  expect_equal(unname(render('equal')), unname(reference), tolerance = 1e-6)
  # Device access can be denied by a test sandbox even in a Metal-enabled build.
  gpu = tryCatch(render('metal'), error = identity)
  if (inherits(gpu, 'error')) {
    expect_match(
      conditionMessage(gpu),
      'Metal.*(unavailable|device)|No Metal device'
    )
  } else {
    expect_equal(unname(gpu), unname(reference), tolerance = 1e-6)
  }
})

test_that("HLBVH retains priority exclusion and interior medium initialization", {
  scene = cube(
    xwidth = 3,
    ywidth = 3,
    zwidth = 3,
    material = dielectric(refraction = 1, priority = 0)
  ) |>
    add_object(sphere(
      material = subsurface(
        sigma_a = 100,
        sigma_s = 100,
        refraction = 1.8,
        priority = 1
      )
    ))
  for (method in c('sah', 'hlbvh', 'metal')) {
    for (case in list(
      list(scene = scene, origin = c(0, 0, 4)),
      list(scene = scene, origin = c(0, 0, 0))
    )) {
      image = tryCatch(
        render_scene(
          case$scene,
          bvh_type = method,
          width = 8,
          height = 8,
          samples = 4,
          lookfrom = case$origin,
          lookat = c(0, 0, -2),
          aperture = 0,
          fov = 0,
          ortho_dimensions = c(.5, .5),
          ambient_light = TRUE,
          backgroundhigh = 'white',
          backgroundlow = 'white',
          mode = 'image',
          denoise = FALSE,
          plot_scene = FALSE,
          progress = FALSE,
          max_depth = 4
        ),
        error = identity
      )
      if (inherits(image, 'error')) {
        expect_identical(method, 'metal')
        expect_match(
          conditionMessage(image),
          'Metal.*(unavailable|device)|No Metal device'
        )
      } else {
        expect_equal(
          as.numeric(image[,, 1:3]),
          rep(1, 8 * 8 * 3),
          tolerance = 1e-6
        )
      }
    }
  }
})
