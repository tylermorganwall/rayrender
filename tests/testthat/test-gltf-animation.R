gltf_fixture = function(name = "rigid.gltf") {
  testthat::test_path("..", "..", "inst", "animation", "fixtures", name)
}

test_that("glTF, external buffers, and GLB decode the same owned source arrays", {
  source = read_gltf_native(gltf_fixture())
  expect_identical(
    source,
    read_gltf_native(gltf_fixture("rigid-external.gltf"))
  )
  expect_identical(source, read_gltf_native(gltf_fixture("rigid.glb")))
  attributes = source$meshes[[1]]$primitives[[1]]$attributes
  # Sparse replacements are tightly packed even though base vertices are strided.
  expect_equal(
    attributes$POSITION,
    rbind(c(0, 0, 0), c(2, 0, 0), c(0, 1, 0), c(0, 0, 2))
  )
  expect_equal(attributes$TEXCOORD_0, rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1)))
  expect_equal(
    attributes$NORMAL,
    rbind(rep(-sqrt(1 / 3), 3), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)),
    tolerance = 1e-7
  )
  expect_equal(
    source$meshes[[1]]$primitives[[1]]$indices,
    c(1, 3, 2, 1, 2, 4, 1, 4, 3, 2, 3, 4)
  )
  # The parser preserves future deformation source data, but the public rigid
  # importer diagnoses it before a misleading static asset can be rendered.
  expect_length(
    read_gltf_native(gltf_fixture("skin.gltf"))$skins[[1]]$joints,
    2
  )
  expect_length(
    read_gltf_native(gltf_fixture("morph.gltf"))$meshes[[1]]$primitives[[
      1
    ]]$targets,
    1
  )
})

test_that("imported clips preserve hierarchy, clocks, and authored cubic signs", {
  asset = read_gltf(gltf_fixture())
  clips = list_animations(asset)
  expect_equal(clips$id, c("clip-0", "clip-1", "clip-2"))
  expect_equal(clips$duration, c(3, 2, 2))
  expect_error(
    animation_model_instance(asset, animation = clip_player("motion")),
    "Ambiguous"
  )
  instance = animation_model_instance(
    asset,
    animation = clip_player("clip-0", loop = "repeat")
  )
  poses = animation_model_pose(instance, c(1, 0, 4, -2, 1))
  expect_identical(poses[[1]]$global, poses[[3]]$global)
  expect_identical(poses[[1]]$global, poses[[4]]$global)
  expect_identical(poses[[1]]$global, poses[[5]]$global)
  expect_equal(
    poses[[1]]$global[["node-0"]][1:3, 4],
    c(2 / 3, 0, 0),
    tolerance = 1e-6
  )
  expect_equal(
    poses[[1]]$global[["node-1"]][1:3, 4],
    c(2 / 3, 2, 0),
    tolerance = 1e-6
  )
  expect_equal(
    poses[[1]]$local[["node-1"]][1:2, 1:2],
    matrix(c(0, 1, -1, 0), 2),
    tolerance = 1e-6
  )
  cubic = animation_model_pose(
    animation_model_instance(asset, animation = clip_player("clip-2")),
    1
  )[[1]]
  expected = matrix(c(-sqrt(.5), -sqrt(.5), sqrt(.5), -sqrt(.5)), 2)
  expect_equal(cubic$local[["node-1"]][1:2, 1:2], expected, tolerance = 1e-6)
  step = animation_model_pose(
    animation_model_instance(asset, animation = clip_player("clip-1")),
    c(.5, 1)
  )
  expect_equal(
    diag(step[[1]]$local[["node-0"]])[1:3],
    rep(1.25, 3),
    tolerance = 1e-6
  )
  expect_equal(step[[1]]$local[["node-1"]][1:3, 4], c(0, 0, 0))
  expect_equal(step[[2]]$local[["node-1"]][1:3, 4], c(0, 0, 1))
  matrix_asset = read_gltf(gltf_fixture(), scene = 2)
  matrix_pose = animation_model_pose(
    animation_model_instance(matrix_asset),
    0
  )[[1]]$world[[
    "node-2"
  ]]
  expect_equal(diag(matrix_pose), c(-1, 1, 1, 1))
  expect_equal(matrix_pose[1:3, 4], c(-3, 0, 0))
  expect_null(matrix_asset$nodes[[1]]$scene)
  expect_error(read_gltf(gltf_fixture(), scene = 8), "out of range")
  file = tempfile(fileext = ".gltf")
  file.copy(gltf_fixture(), file)
  independent = read_gltf(file)
  unlink(file)
  saved = unserialize(serialize(
    animation_model_instance(independent, animation = clip_player("clip-0")),
    NULL
  ))
  expect_equal(animation_model_pose(saved, 1)[[1]]$global, poses[[1]]$global)
})

test_that("invalid glTF storage and unsupported animation fail explicitly", {
  for (name in c(
    "invalid-offset.gltf",
    "invalid-overflow.gltf",
    "invalid-index.gltf",
    "invalid-sparse.gltf"
  )) {
    expect_error(
      read_gltf(gltf_fixture(name)),
      "glTF.*(exceeds|range|validation)"
    )
  }
  expect_error(
    read_gltf(gltf_fixture("invalid-time.gltf")),
    "times must increase"
  )
  expect_error(
    read_gltf(gltf_fixture("invalid-extension.gltf")),
    "Required extensions"
  )
  expect_error(
    read_gltf(gltf_fixture("invalid-matrix-animation.gltf")),
    "authored matrix"
  )
  expect_error(
    read_gltf(gltf_fixture("invalid-cycle.gltf")),
    "parse failed|validation"
  )
  expect_error(
    read_gltf(gltf_fixture("invalid-duplicate-channel.gltf")),
    "Duplicate animation"
  )
  expect_error(
    read_gltf(gltf_fixture("network.gltf")),
    "network loading is disabled"
  )
  expect_error(
    read_gltf(gltf_fixture("skin.gltf")),
    "Skin requires contiguous paired"
  )
})

test_that("two imported instances render independent paths and phases", {
  asset = read_gltf(
    gltf_fixture("rigid.glb"),
    material = light("orange", intensity = 1)
  )
  path = transform_track(
    c(0, .4, 1.5, 3),
    position = rbind(c(-3, 0, 0), c(-2, 0, .5), c(-3, .5, 0), c(-3, 0, 0))
  )
  first = animation_model_instance(
    asset,
    path,
    clip_player("clip-0", loop = "repeat")
  )
  second = animation_model_instance(
    asset,
    transform_track(
      c(0, 1, 2, 3),
      position = rbind(c(2, 0, 0), c(2, .3, 1), c(1, 0, 0), c(2, 0, 0))
    ),
    clip_player("clip-0", offset = .5, speed = -1, loop = "repeat")
  )
  scene = unserialize(serialize(add_object(first, second), NULL))
  settings = list(
    width = 48,
    height = 36,
    samples = 16,
    camera = camera(lookfrom = c(1, 5, 13), lookat = c(1, 1, 0)),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    bloom = FALSE,
    denoise = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    tonemap = "raw",
    sample_method = "sobol",
    min_variance = 0,
    parallel = FALSE
  )
  times = c(.25, 1, -.5, .25)
  cpu = do.call(
    render_animation,
    c(list(scene = scene, times = times, integrator_type = "nee"), settings)
  )
  expect_identical(as.numeric(cpu[[1]]), as.numeric(cpu[[4]]))
  expect_gt(mean(abs(cpu[[1]] - cpu[[2]])), .01)
  for (i in seq_along(times)) {
    expected = add_object(
      flatten_model_pose(first, times[i]),
      flatten_model_pose(second, times[i])
    )
    reference = do.call(
      render_scene,
      c(
        list(scene = expected, time = times[i], integrator_type = "nee"),
        settings
      )
    )
    expect_lt(mean(abs(as.numeric(cpu[[i]]) - as.numeric(reference))), .015)
  }
  gpu = suppressWarnings(do.call(
    render_animation,
    c(list(scene = scene, times = times, integrator_type = "metal"), settings)
  ))
  if (!isTRUE(attr(gpu[[1]], "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  for (i in seq_along(times)) {
    report = attr(gpu[[i]], "wavefront")
    expect_true(report$used)
    expect_equal(report$unique_meshes, 1)
    expect_equal(report$mesh_instances, 4)
    expect_equal(report$triangles, 4)
    expect_identical(report$scene_reused, i > 1L)
    expect_lt(mean(abs(as.numeric(cpu[[i]]) - as.numeric(gpu[[i]]))), .025)
  }
})
