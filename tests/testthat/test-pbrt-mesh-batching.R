test_that('PBRT mesh batching preserves shading and leaves boundaries separate', {
  texture = array(1, c(8, 8, 4))
  texture[,, 1] = diag(8)
  texture[,, 2] = .3
  texture[1:2, , 4] = .2
  material = diffuse(
    image_texture = texture,
    bump_texture = diag(8),
    bump_intensity = .03
  )
  rows = lapply(0:31, function(i) {
    x = (i %% 8) * .25 - .875
    y = (i %/% 8) * .25 - .375
    mesh = rayvertex::construct_mesh(
      vertices = rbind(
        c(x - .1, y - .1, 0),
        c(x + .1, y - .1, 0),
        c(x + .1, y + .1, 0),
        c(x - .1, y + .1, 0)
      ),
      indices = rbind(c(0L, 1L, 2L), c(0L, 2L, 3L)),
      normals = matrix(c(0, 0, 1), 4, 3, byrow = TRUE),
      norm_indices = rbind(c(0L, 1L, 2L), c(0L, 2L, 3L)),
      texcoords = rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)),
      tex_indices = rbind(c(0L, 1L, 2L), c(0L, 2L, 3L))
    )
    raymesh_model(mesh, material = material)
  })
  combined = pbrt_batch_mesh_rows(rows)
  expect_length(combined, 1L)
  expect_length(combined[[1]]$shape_info[[1]]$mesh_info[[1]]$shapes, 32L)
  separate = vctrs::list_unchop(rows)
  old = options(cores = 1L)
  on.exit(options(old), add = TRUE)
  render = function(scene) {
    set.seed(619)
    scene = add_infinite_light(
      scene,
      disk_light(direction = c(-1, 2, 3), intensity = 20000)
    )
    render_scene(
      scene,
      width = 32,
      height = 16,
      samples = 8,
      max_depth = 6,
      lookfrom = c(0, 0, 5),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(2.4, 1.2),
      sample_method = 'random',
      min_variance = 0,
      preview = FALSE,
      progress = FALSE,
      plot_scene = FALSE,
      parallel = FALSE,
      denoise = FALSE,
      bloom = FALSE,
      tonemap = 'raw'
    )
  }
  expect_equal(render(combined[[1]]), render(separate), tolerance = 1e-6)
  for (replacement in list(dielectric(), subsurface(), light(), openpbr())) {
    altered = lapply(rows, function(x) {
      x$material = replacement
      x
    })
    expect_identical(pbrt_batch_mesh_rows(altered), altered)
  }
  altered = lapply(rows, function(x) {
    set_medium(x, homogeneous_medium(sigma_s = .1))
  })
  expect_identical(pbrt_batch_mesh_rows(altered), altered)
  for (setting in c(
    'displacement_texture',
    'subdivision_levels',
    'recalculate_normals'
  )) {
    altered = lapply(rows, function(x) {
      x$shape_info[[1]]$shape_properties[[setting]] = switch(
        setting,
        displacement_texture = 'unused.png',
        subdivision_levels = 2,
        recalculate_normals = TRUE
      )
      x
    })
    expect_identical(pbrt_batch_mesh_rows(altered), altered)
  }
  # Successive compactions must keep the shape-count bound, not just row count.
  twice = pbrt_batch_mesh_rows(rep(combined, 32))
  expect_true(all(vapply(
    twice,
    function(x) length(x$shape_info[[1]]$mesh_info[[1]]$shapes) <= 256L,
    logical(1)
  )))
  # A mixed complete/missing normal group must not disable consistency tables
  # on the previously complete meshes.
  altered = rows
  for (i in seq(2, 32, 2)) {
    altered[[i]]$shape_info[[1]]$mesh_info[[1]]$shapes[[
      1
    ]]$has_vertex_normals = c(FALSE, FALSE)
  }
  expect_identical(pbrt_batch_mesh_rows(altered), altered)
})
