test_that('Metal applies consistent mesh normals to glossy reflection', {
  mesh = tempfile(fileext = '.obj')
  on.exit(unlink(mesh))
  # An intentionally coarse smooth face makes the view-dependent correction
  # observable without relying on a dense production asset or narrow highlights.
  writeLines(
    c(
      'v -4 -4 0',
      'v 4 -4 0',
      'v 0 4 0',
      'vn -0.6 0 0.8',
      'vn 0.6 0 0.8',
      'vn 0 0.6 0.8',
      'f 1//1 2//2 3//3'
    ),
    mesh
  )
  images = list()
  for (instanced in c(FALSE, TRUE)) {
    for (consistent in c(FALSE, TRUE)) {
      scene = obj_model(
        mesh,
        load_material = FALSE,
        calculate_consistent_normals = consistent,
        material = glossy(color = c(.6, .7, .8), gloss = .5)
      )
      if (instanced) {
        scene = create_instances(
          scene,
          scale_x = .8,
          scale_y = 1.4,
          scale_z = .9
        )
      }
      scene = add_infinite_light(
        scene,
        disk_light(
          direction = c(.8, 0, 1),
          angular_diameter = 30,
          intensity = 8
        )
      )
      for (backend in c('nee', 'metal')) {
        set.seed(519)
        img = render_scene(
          scene,
          width = 24,
          height = 24,
          samples = 512,
          sample_method = 'sobol',
          lookfrom = c(2, 0, 4),
          lookat = c(0, 0, 0),
          fov = 0,
          ortho_dimensions = c(.8, .8),
          integrator_type = backend,
          max_depth = 2,
          min_variance = 0,
          denoise = FALSE,
          tonemap = 'raw',
          bloom = FALSE,
          progress = FALSE,
          preview = FALSE,
          plot_scene = FALSE
        )
        if (backend == 'metal' && !isTRUE(attr(img, 'wavefront')$used)) {
          expect_match(
            attr(img, 'wavefront')$fallback,
            'not enabled|device|macOS 11'
          )
          skip('Metal is unavailable')
        }
        if (backend == 'metal' && instanced) {
          expect_equal(attr(img, 'wavefront')$mesh_instances, 1)
        }
        images[[paste(consistent, backend)]] = mean(img[,, 1:3])
      }
      expect_equal(
        images[[paste(consistent, 'metal')]],
        images[[paste(consistent, 'nee')]],
        tolerance = .015
      )
    }
    expect_gt(abs(images[['TRUE nee']] - images[['FALSE nee']]), .02)
  }
})
