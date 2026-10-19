test_that("shared assets retain rest properties and deterministic hierarchy poses", {
  asset = rigid_model_fixture()
  clips = list_animations(asset)
  expect_identical(clips$id, "flap")
  expect_equal(clips$duration, 2)
  expect_setequal(clips$targets[[1]], c("body", "left", "right"))
  path = transform_track(
    c(0, .7, 2, 3),
    position = rbind(c(-1, 0, 0), c(0, 1, 0), c(1, .2, 1), c(2, 0, 0)),
    position_interpolation = "cubic"
  )
  first = animation_model_instance(
    asset,
    path,
    clip_player("flap", loop = "repeat")
  )
  second = animation_model_instance(
    asset,
    path,
    clip_player("Flight", offset = .5, speed = -2, loop = "pingpong")
  )
  expect_identical(
    first$animation_info[[1]]$model_instance$asset,
    second$animation_info[[1]]$model_instance$asset
  )
  expect_error(asset$nodes <- list(), "locked")
  times = c(8, -1, .125, 3, 0, 1e9, .125)
  poses = animation_model_pose(first, times)
  expect_identical(poses[[3]]$local, poses[[7]]$local)
  expect_identical(poses[[3]]$global, poses[[7]]$global)
  expect_identical(poses[[3]]$world, poses[[7]]$world)
  for (pose in poses) {
    expect_equal(pose$local$left[1:3, 4], c(-.2, 0, 0))
    expect_equal(
      pose$global$left,
      pose$global$body %*% pose$local$left,
      tolerance = 1e-6
    )
  }
  expect_false(identical(
    animation_model_pose(first, .4)[[1]]$local,
    animation_model_pose(second, .4)[[1]]$local
  ))
  frozen = animation_model_instance(
    asset,
    path,
    clip_player("flap", offset = .2, speed = 0)
  )
  flags = animation_model_pose(frozen, c(0, 1, 2))
  expect_true(flags[[2]]$placement_changed)
  expect_false(flags[[2]]$pose_changed)
  expect_identical(flags[[1]]$global, flags[[3]]$global)
  legacy = animate_objects(
    animation_model_instance(asset),
    end_position = c(1, 0, 0)
  )
  expect_error(animation_model_pose(legacy, .5), "legacy matrix attachments")
})

test_that("one clip clock spans all independently scheduled node channels", {
  asset = rigid_model_fixture()
  for (loop in c("hold", "repeat", "pingpong")) {
    for (speed in c(-2, 0, 1.5)) {
      player = clip_player(
        "flap",
        start_time = .3,
        offset = .7,
        speed = speed,
        loop = loop
      )
      instance = animation_model_instance(asset, animation = player)
      times = c(-3, 0, .3, 1.1, 8, .3)
      x = speed * (times - .3) + .7
      expected = 2 +
        switch(
          loop,
          hold = pmax(0, pmin(2, x)),
          "repeat" = x %% 2,
          pingpong = 2 - abs(x %% 4 - 2)
        )
      poses = animation_model_pose(instance, times)
      expect_equal(vapply(poses, `[[`, numeric(1), "clip_time"), expected)
      expect_identical(poses[[3]]$global, poses[[6]]$global)
    }
  }
})

test_that("composition remaps instance identities and preserves cameras and serialization", {
  asset = rigid_model_fixture()
  first = animation_model_instance(
    asset,
    animation = clip_player("flap", loop = "repeat")
  )
  combined = add_object(first, group_objects(first, translate = c(3, 0, 0)))
  bindings = lapply(combined$animation_info, `[[`, "model_instance")
  expect_false(identical(bindings[[1]]$id, bindings[[2]]$id))
  expect_identical(bindings[[1]]$asset, bindings[[2]]$asset)
  bound = rbind(first, first)
  expect_false(identical(
    bound$animation_info[[1]]$model_instance$id,
    bound$animation_info[[2]]$model_instance$id
  ))
  nested = add_object(create_instances(first), create_instances(first))
  nested_ids = model_instance_ids(nested)
  expect_length(nested_ids, 2)
  expect_equal(anyDuplicated(nested_ids), 0L)
  check_nodes = function(scene, id) {
    for (i in seq_len(nrow(scene))) {
      node = scene$animation_info[[i]]$model_node
      if (!is.null(node)) {
        expect_identical(node$instance$id, id)
      }
      prototype = scene$shape_info[[i]]$shape_properties$original_scene
      if (!is.null(prototype)) check_nodes(prototype[[1]], id)
    }
  }
  for (i in 1:2) {
    check_nodes(combined[i, ], bindings[[i]]$id)
  }
  scene = add_camera(combined, camera(lookfrom = c(0, 2, 6)), name = "view")
  scene = add_object(scene, sphere(y = -1))
  expect_named(ray_scene_cameras(scene), "view")
  loaded = unserialize(serialize(scene, NULL))
  expect_identical(
    loaded$animation_info[[1]]$model_instance$asset,
    loaded$animation_info[[2]]$model_instance$asset
  )
  for (i in 1:2) {
    expect_equal(
      animation_model_pose(loaded[i, ], c(2, .3)),
      animation_model_pose(scene[i, ], c(2, .3))
    )
  }
  p = animation_model_pose(combined[2, ], .3)[[1]]
  q = animation_model_pose(combined[1, ], .3)[[1]]
  expect_equal(p$world$body[1, 4] - q$world$body[1, 4], 3)
  expect_error(
    render_scene(create_instances(first), preview = FALSE),
    "explicit scene time"
  )
})

test_that("root motion removes only selected locomotion and preserves bob and articulation", {
  asset = rigid_model_fixture()
  player = clip_player(
    "flap",
    loop = "repeat",
    root_motion = "in_place",
    locomotion_root = "body",
    remove_yaw = TRUE
  )
  instance = animation_model_instance(asset, animation = player)
  poses = animation_model_pose(instance, c(1, 5, 0, 1))
  expect_equal(poses[[1]]$local$body[1:3, 4], c(0, .6, 0), tolerance = 1e-6)
  expect_equal(poses[[1]]$local$body[1:3, 1:3], diag(3), tolerance = 1e-6)
  expect_identical(poses[[1]]$global, poses[[4]]$global)
  expect_false(identical(poses[[1]]$local$left, poses[[3]]$local$left))
  player$translation_axes = "x"
  selected = animation_model_pose(
    animation_model_instance(asset, animation = player),
    1
  )[[1]]
  expect_equal(selected$local$body[1:3, 4], c(0, .6, .2), tolerance = 1e-6)
  expect_error(
    clip_player(root_motion = "in_place"),
    "explicit locomotion root"
  )
  expect_error(
    animation_model_instance(
      asset,
      animation = clip_player(
        "flap",
        root_motion = "in_place",
        locomotion_root = "left"
      )
    ),
    "top-level"
  )
  expect_error(
    animate_keyframes(sphere(), transform_track(), player = player),
    "bare transform|hierarchy"
  )
  expect_error(clip_player(root_motion = "extract"), "not implemented")
  expect_error(clip_player(root_motion = NA_character_), "Unsupported")
  expect_error(
    clip_player(root_motion = c("preserve", "in_place")),
    "Unsupported"
  )
})

test_that("invalid graphs, targets, and ambiguous names fail before rendering", {
  expect_error(
    animation_model_asset(list(
      animation_model_node("a", sphere(), "b"),
      animation_model_node("b", parent = "a")
    )),
    "cycle"
  )
  expect_error(
    animation_model_asset(list(animation_model_node("a", sphere(), "missing"))),
    "Unknown model node"
  )
  expect_error(
    animation_model_asset(list(
      animation_model_node("a"),
      animation_model_node("a")
    )),
    "unique"
  )
  expect_error(
    animation_model_asset(
      list(animation_model_node("a", sphere())),
      list(animation_model_clip("clip", list(missing = transform_track())))
    ),
    "Unknown model node"
  )
  expect_error(
    animation_model_clip(
      "clip",
      list(a = transform_track(), a = transform_track())
    ),
    "unique"
  )
  asset = animation_model_asset(
    list(animation_model_node("body", sphere())),
    list(
      animation_model_clip("a", list(), "same"),
      animation_model_clip("b", list(), "same")
    )
  )
  expect_error(
    animation_model_instance(asset, animation = clip_player("same")),
    "Ambiguous"
  )
  expect_equal(
    animation_model_pose(
      animation_model_instance(asset, animation = clip_player("b")),
      10
    )[[
      1
    ]]$clip_time,
    0
  )
  expect_error(
    animation_model_asset(list(animation_model_node(
      "body",
      animate_objects(sphere())
    ))),
    "static"
  )
})

test_that("hierarchy snapshots and loop exposures render through CPU and Metal", {
  asset = rigid_model_fixture(emissive = TRUE)
  first = animation_model_instance(
    asset,
    transform_track(
      c(0, 1, 2, 3),
      position = rbind(
        c(-1.3, -.5, 0),
        c(-1, 0, .2),
        c(-1.2, .2, -.3),
        c(-1.3, -.5, 0)
      )
    ),
    clip_player(
      "flap",
      loop = "repeat",
      root_motion = "in_place",
      locomotion_root = "body",
      remove_yaw = TRUE
    )
  )
  second = animation_model_instance(
    asset,
    transform_track(
      c(0, .8, 2, 3),
      position = rbind(
        c(1.2, -.5, .3),
        c(1.3, .1, -.2),
        c(1.2, -.1, .4),
        c(1.2, -.5, .3)
      )
    ),
    clip_player("flap", offset = .5, speed = -1.5, loop = "pingpong")
  )
  scene = add_object(first, second)
  settings = list(
    width = 48,
    height = 36,
    samples = 16,
    camera = camera(
      lookfrom = c(0, 2.5, 7),
      lookat = c(0, .3, 0),
      fov = 40,
      aperture = 0
    ),
    sample_method = "sobol",
    min_variance = 0,
    parallel = FALSE,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    tonemap = "raw"
  )
  times = c(.25, 2, -.25, 1, .25)
  set.seed(714)
  cpu = do.call(
    render_animation,
    c(list(scene = scene, times = times, integrator_type = "nee"), settings)
  )
  expect_identical(as.numeric(cpu[[1]]), as.numeric(cpu[[5]]))
  expect_gt(mean(abs(cpu[[1]] - cpu[[4]])), .005)
  loaded = unserialize(serialize(scene, NULL))
  set.seed(714)
  loaded_cpu = do.call(
    render_animation,
    c(list(scene = loaded, times = times, integrator_type = "nee"), settings)
  )
  for (i in seq_along(times)) {
    expect_identical(as.numeric(loaded_cpu[[i]]), as.numeric(cpu[[i]]))
  }
  for (i in seq_along(times)) {
    expected_scene = add_object(
      flatten_model_pose(first, times[i]),
      flatten_model_pose(second, times[i])
    )
    expected = do.call(
      render_scene,
      c(
        list(scene = expected_scene, time = times[i], integrator_type = "nee"),
        settings
      )
    )
    expect_lt(mean(abs(as.numeric(expected) - as.numeric(cpu[[i]]))), .01)
  }
  gpu = suppressWarnings(do.call(
    render_animation,
    c(list(scene = loaded, times = times, integrator_type = "metal"), settings)
  ))
  if (!isTRUE(attr(gpu[[1]], "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  for (i in seq_along(times)) {
    expect_true(attr(gpu[[i]], "wavefront")$used)
    # Both models use the same two wing meshes: four placements, 24 source
    # triangles total. Duplicating source geometry would produce four meshes.
    expect_equal(attr(gpu[[i]], "wavefront")$unique_meshes, 2)
    expect_equal(attr(gpu[[i]], "wavefront")$mesh_instances, 4)
    expect_equal(attr(gpu[[i]], "wavefront")$triangles, 24)
    expect_identical(attr(gpu[[i]], "wavefront")$scene_reused, i > 1L)
    expect_lt(mean(abs(as.numeric(cpu[[i]]) - as.numeric(gpu[[i]]))), .02)
  }
  for (backend in c("nee", "metal")) {
    blurred = do.call(
      render_scene,
      c(
        list(
          scene = scene,
          time = 2,
          shutter_offsets = c(-.4, .4),
          integrator_type = backend
        ),
        settings
      )
    )
    expect_true(all(is.finite(blurred)))
    expect_gt(mean(blurred), .005)
    if (backend == "nee") {
      reference = blurred
    } else {
      expect_lt(mean(abs(as.numeric(reference) - as.numeric(blurred))), .02)
    }
  }
})

test_that("copied model players are owned by the outer instance after reload", {
  asset = rigid_model_fixture(emissive = TRUE)
  first = animation_model_instance(
    asset,
    animation = clip_player("flap", speed = 0)
  )
  second = group_objects(first, translate = c(2, 0, 0))
  # Ordinary R replacement must not leave stale descendant clocks. This also
  # exercises equal prototype list identities under different outer players.
  second$animation_info[[1]]$model_instance$player = clip_player(
    "flap",
    offset = .75,
    speed = -1,
    loop = "repeat"
  )
  scene = unserialize(serialize(add_object(first, second), NULL))
  settings = list(
    time = .2,
    width = 40,
    height = 30,
    samples = 16,
    lookfrom = c(1, 2, 7),
    lookat = c(1, .3, 0),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    sample_method = "sobol",
    min_variance = 0,
    parallel = FALSE,
    tonemap = "raw",
    integrator_type = "nee"
  )
  expected = add_object(
    flatten_model_pose(first, .2),
    flatten_model_pose(second, .2)
  )
  reference = do.call(render_scene, c(list(scene = expected), settings))
  set.seed(507)
  result = do.call(render_scene, c(list(scene = scene), settings))
  expect_lt(mean(abs(as.numeric(reference) - as.numeric(result))), .01)
  single_reference = do.call(
    render_scene,
    c(list(scene = flatten_model_pose(second, .2)), settings)
  )
  single = do.call(render_scene, c(list(scene = second), settings))
  expect_lt(mean(abs(as.numeric(single_reference) - as.numeric(single))), .01)
  expect_false(identical(
    animation_model_pose(scene[1, ], .2)[[1]]$local,
    animation_model_pose(scene[2, ], .2)[[1]]$local
  ))
  set.seed(507)
  repeated = do.call(
    render_scene,
    c(list(scene = create_instances(scene)), settings)
  )
  expect_identical(as.numeric(result), as.numeric(repeated))
  tinted = set_scene_material(second, light(color = "red", intensity = 1))
  tinted_reference = set_scene_material(
    flatten_model_pose(second, .2),
    light(color = "red", intensity = 1)
  )
  tinted_image = do.call(render_scene, c(list(scene = tinted), settings))
  expected_image = do.call(
    render_scene,
    c(list(scene = tinted_reference), settings)
  )
  expect_lt(
    mean(abs(as.numeric(tinted_image) - as.numeric(expected_image))),
    .01
  )
  # The override removes orange's green emission, independently of how much
  # screen area the bird occupies. A whole-image difference would conflate it
  # with framing and antialiasing coverage.
  expect_gt(max(single[,, 2]), .1)
  expect_lt(max(tinted_image[,, 2]), 1e-6)
  expect_gt(max(tinted_image[,, 1]), .5)
  expect_identical(
    first$animation_info[[1]]$model_instance$asset,
    tinted$animation_info[[1]]$model_instance$asset
  )
})

test_that('Metal snapshots preserve STEP sides next to nonrepresentable repeat periods', {
  times = readBin(
    writeBin(c(1 / 24, 25 / 24, 49 / 24), raw(), size = 4),
    numeric(),
    n = 3,
    size = 4
  )
  duration = diff(range(times))
  asset = animation_model_asset(
    list(animation_model_node(
      'part',
      sphere(material = light('orange', intensity = 1))
    )),
    list(animation_model_clip(
      'loop',
      list(
        part = transform_track(
          times,
          position = rbind(c(0, 0, 0), c(10, 0, 0), c(100, 0, 0)),
          position_interpolation = 'step'
        )
      )
    ))
  )
  scene = animation_model_instance(
    asset,
    animation = clip_player('loop', loop = 'repeat')
  )
  settings = list(
    scene = scene,
    times = c(2, duration + 1e-4, 0, 2),
    width = 48,
    height = 24,
    samples = 32,
    sample_method = "sobol",
    min_variance = 0,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = 'black',
    backgroundlow = 'black',
    tonemap = 'raw',
    camera = camera(
      lookfrom = c(5, 0, 8),
      lookat = c(5, 0, 0),
      fov = 0,
      ortho_dimensions = c(14, 7)
    )
  )
  cpu = do.call(render_animation, c(settings, list(integrator_type = 'nee')))
  gpu = do.call(render_animation, c(settings, list(integrator_type = 'metal')))
  if (!isTRUE(attr(gpu[[1]], 'wavefront')$used)) {
    skip('Metal runtime unavailable')
  }
  for (i in seq_along(gpu)) {
    expect_true(attr(gpu[[i]], 'wavefront')$used)
    expect_lt(mean(abs(as.numeric(cpu[[i]]) - as.numeric(gpu[[i]]))), .002)
  }
  expect_gt(mean(abs(as.numeric(cpu[[1]]) - as.numeric(cpu[[2]]))), .005)
  expect_identical(as.numeric(gpu[[1]]), as.numeric(gpu[[4]]))
})
