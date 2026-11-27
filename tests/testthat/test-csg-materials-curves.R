test_that('CSG components inherit NULL and reject incompatible mixtures', {
  inherited = csg_combine(csg_sphere(), csg_box())
  expect_null(inherited$object1$material)
  expect_length(
    rayrender:::prepare_csg_materials(inherited, dielectric())$materials,
    0
  )
  colored = csg_combine(csg_sphere(material = diffuse('red')), csg_box())
  palette = rayrender:::prepare_csg_materials(colored, diffuse('blue'))
  expect_length(palette$materials, 1)
  expect_equal(palette$object$object1$surface_material, 1)
  expect_null(palette$object$object2$surface_material)
  expect_error(csg_object(colored, material = dielectric()), 'dielectric')
  for (mat in list(
    dielectric(),
    subsurface(),
    subsurface_diffusion(),
    microfacet(transmission = TRUE),
    openpbr(transmission_weight = .5),
    openpbr(subsurface_weight = .5)
  )) {
    expect_error(
      csg_object(csg_sphere(material = mat)),
      'dielectric|subsurface'
    )
  }
  expect_error(
    csg_object(csg_sphere(material = diffuse(bump_texture = matrix(0, 2, 2)))),
    'bump'
  )
  expect_error(
    csg_combine(csg_sphere(), csg_box(), 'mix', radius = 2),
    'radius'
  )
  expect_error(
    csg_combine(csg_sphere(), csg_box(), blend_materials = NA),
    'blend_materials'
  )
  expect_error(csg_combine(csg_sphere(), csg_box(), 'typo'), 'arg')
})

test_that('3D CSG curves are closed capped solids with controlled taper', {
  points = rbind(c(0, 0, 0), c(0, 1, 0), c(1, 1, 1), c(1, 0, 1))
  curve = csg_curve(points, radius = c(.2, .05), segments = 8)
  expect_length(curve$shapes, 8)
  expect_equal(curve$shapes[[1]]$start, points[1, ])
  expect_equal(curve$shapes[[8]]$end, points[4, ])
  expect_equal(curve$shapes[[1]]$radius, .2)
  expect_equal(curve$shapes[[8]]$upper_radius, .05)
  expect_null(curve$material)
  expect_silent(csg_object(curve, material = dielectric()))
  stationary = csg_curve(matrix(0, 4, 3), radius = c(.2, .1), segments = 4)
  expect_true(all(vapply(
    stationary$shapes,
    function(s) s$csg_type == 2,
    logical(1)
  )))
  expect_error(csg_curve(points[, 1:2]), 'points')
  expect_error(csg_curve(points, radius = 0), 'radius')
  expect_error(csg_curve(points, segments = 1.5), 'segments')
  expect_error(csg_curve(points, material = subsurface()), 'subsurface')
})

# Constant environment tests measure the mixture without a denoiser or exposure.
# The output is linear radiance, so averaging inherited red and explicit blue
# Lambertian surfaces has an analytical reference at equal-distance joins.
render_csg_material_test = function(
  shape,
  backend,
  material = diffuse('red'),
  ...
) {
  set.seed(39)
  image = render_scene(
    csg_object(shape, material = material),
    width = 24,
    height = 24,
    samples = 512,
    sample_method = 'sobol',
    lookfrom = c(0, 0, 4),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.1, .1),
    integrator_type = backend,
    ambient_light = TRUE,
    backgroundlow = 'white',
    backgroundhigh = 'white',
    parallel = TRUE,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    min_variance = 0,
    denoise = FALSE,
    bloom = FALSE,
    ...
  )
  if (backend == 'metal' && !isTRUE(attr(image, 'wavefront')$used)) {
    reason = attr(image, 'wavefront')$fallback
    if (grepl('not enabled|device|macOS 11', reason)) {
      skip('Metal unavailable')
    }
    stop('Unexpected fallback: ', reason)
  }
  if (backend == 'metal') {
    expect_equal(attr(image, 'wavefront')$csg_failed_queries, 0)
  }
  expect_true(all(is.finite(image)))
  apply(image[,, 1:3], 3, mean)
}

for (backend in c('nee', 'metal')) {
  local({
    backend = backend
    test_that(
      paste(
        backend,
        'resolves inherited, hard, and smooth component materials'
      ),
      {
        red = csg_sphere()
        blue = csg_sphere(material = diffuse('blue'))
        expect_equal(
          render_csg_material_test(csg_combine(red, blue), backend),
          c(1, 0, 0),
          tolerance = .03
        )
        expect_equal(
          render_csg_material_test(
            csg_combine(red, blue, 'blend', radius = .4),
            backend
          ),
          c(.5, 0, .5),
          tolerance = .035
        )
        expect_equal(
          render_csg_material_test(
            csg_combine(
              red,
              blue,
              'blend',
              radius = .4,
              blend_materials = FALSE
            ),
            backend
          ),
          c(1, 0, 0),
          tolerance = .03
        )
        expect_equal(
          render_csg_material_test(
            csg_combine(red, blue, 'mix', radius = .25),
            backend
          ),
          c(.75, 0, .25),
          tolerance = .035
        )
        # The cutter defines the visible concave face, not the enclosing object.
        cut = csg_combine(
          red,
          csg_sphere(z = 1, radius = .65, material = diffuse('blue')),
          'subtract'
        )
        expect_gt(render_csg_material_test(cut, backend)[3], .5)
        # A component material must also survive the coordinate modifier stack.
        transformed = csg_translate(
          csg_scale(csg_group(list(blue)), scale = 1.5),
          z = -.2
        )
        expect_equal(
          render_csg_material_test(transformed, backend),
          c(0, 0, 1),
          tolerance = .03
        )
        # Tapered 3D curve, viewed at its end cap.
        tube = csg_curve(
          rbind(c(0, 0, 0), c(0, 0, -1), c(0, 1, -2), c(1, 1, -3)),
          radius = c(.25, .1),
          material = diffuse('blue')
        )
        expect_equal(
          render_csg_material_test(tube, backend),
          c(0, 0, 1),
          tolerance = .04
        )
      }
    )
  })
}
