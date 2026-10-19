deformation_fixture = function(name, extension = "gltf") {
  path = system.file(
    "animation",
    "fixtures",
    paste0("deform-", name, ".", extension),
    package = "rayrender"
  )
  if (!nzchar(path)) {
    path = file.path(
      "..",
      "..",
      "inst",
      "animation",
      "fixtures",
      paste0("deform-", name, ".", extension)
    )
  }
  path
}

test_that("native skin and morph positions use the correct spaces and all influence sets", {
  for (extension in c("gltf", "glb")) {
    for (kind in c("skin", "morph", "combined")) {
      asset = read_gltf(deformation_fixture(kind, extension))
      player = clip_player("clip-0")
      result = animation_model_surface(asset, player, 1)[[1]][[1]]
      expected = rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
      if (kind != "morph") {
        expected[, 1] = expected[, 1] - 1.3
      }
      if (kind != "skin") {
        expected[, 2] = expected[, 2] + .75
        expected[, 3] = expected[, 3] + .375
      }
      expect_equal(result$positions, expected, tolerance = 1e-6)
      expect_equal(
        result$normals,
        matrix(rep(c(0, 0, 1), 3), ncol = 3, byrow = TRUE)
      )
      rest = animation_model_surface(asset, clip_player(), 0)[[1]][[
        1
      ]]$positions
      expect_equal(rest[, 1], c(0, 1, 0), tolerance = 1e-6)
      expect_equal(
        rest[, 2],
        c(0, 0, 1) + if (kind != "skin") .5 else 0,
        tolerance = 1e-6
      )
      instance = animation_model_instance(
        asset,
        transform_track(0, position = c(10, 0, 0)),
        player
      )
      matrix = animation_model_pose(instance, 1)[[1]]$world[["node-0"]]
      world = t(matrix %*% t(cbind(result$positions, 1)))[, 1:3]
      expect_equal(
        world[, 1],
        c(0, 1, 0) + if (kind == "morph") 15 else 13.7,
        tolerance = 1e-6
      )
    }
  }
})

test_that("morph tracks keep defaults, independent schedules, cubic tangents and STEP boundaries", {
  asset = read_gltf(deformation_fixture("morph"))
  expect_equal(list_animations(asset)$duration, c(2, 2, 1))
  expect_equal(
    animation_model_pose(animation_model_instance(asset), 4)[[1]]$weights[[1]],
    c(.5, 0)
  )
  cubic = animation_model_surface(asset, clip_player("clip-1"), 1)[[1]][[
    1
  ]]$positions
  expect_equal(cubic[1, ], c(0, .6, .15), tolerance = 1e-6)
  step = clip_player("clip-2")
  before = animation_model_surface(asset, step, 1 - 1e-8)[[1]][[1]]$positions
  after = animation_model_surface(asset, step, 1)[[1]][[1]]$positions
  expect_equal(before[1, ], c(0, 0, 0))
  expect_equal(after[1, ], c(0, 1, .5))
})

test_that("deformed copies have independent clocks and repeated backward seeks do not drift", {
  asset = unserialize(serialize(
    read_gltf(deformation_fixture("combined")),
    NULL
  ))
  first = clip_player("clip-0", loop = "repeat")
  second = clip_player("clip-0", loop = "repeat", offset = .7, speed = -1)
  reference = animation_model_surface(asset, first, .5)
  for (time in c(5, -1, 0, .5, 2, 1e6, .5)) {
    animation_model_surface(asset, second, time)
    expect_identical(animation_model_surface(asset, first, .5), reference)
  }
  expect_false(identical(animation_model_surface(asset, second, .5), reference))
  expect_identical(
    animation_model_surface(asset, first, 0),
    animation_model_surface(asset, first, 2)
  )
})

# Independent analytic reference: the two skin joints contribute .3/.7, while
# mesh-node motion cancels in the palette. The two morph deltas are unit Y/Z.
deformation_reference_scene = function(instances, time, material) {
  parts = lapply(instances, function(instance) {
    binding = instance$animation_info[[1]]$model_instance
    player = binding$player
    local = (player$speed * (time - player$start_time) + player$offset) %% 2
    weight = max(0, min(1, local - .25))
    positions = rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
    positions = sweep(positions, 2, c(3 + .7 * local, weight, .5 * weight), '+')
    placement = animation_evaluate_transform(
      instance$animation_info[[1]]$keyframes,
      time
    )[[1]]
    positions = t(placement %*% t(cbind(positions, 1)))[, 1:3]
    mesh = structure(
      list(
        vb = t(positions),
        it = matrix(1:3, 3),
        normals = matrix(rep(c(0, 0, 1), 3), 3)
      ),
      class = 'mesh3d'
    )
    mesh3d_model(mesh, material = material, override_material = TRUE)
  })
  Reduce(add_object, parts)
}

test_that("CPU and Metal snapshots update mesh geometry and preserve independent phases", {
  material = light('orange', intensity = 1)
  asset = read_gltf(deformation_fixture('combined'), material = material)
  first = animation_model_instance(
    asset,
    transform_track(
      c(0, .4, 1, 1.6, 2),
      position = rbind(
        c(-1.5, 0, 0),
        c(-1.3, .1, 0),
        c(-1.5, .3, 0),
        c(-1.7, .1, 0),
        c(-1.5, 0, 0)
      ),
      position_interpolation = 'cubic'
    ),
    clip_player('clip-0', loop = 'repeat')
  )
  second = animation_model_instance(
    asset,
    transform_track(
      c(0, .5, 1.2, 1.7, 2),
      position = rbind(
        c(1.3, 0, 0),
        c(1.1, .2, 0),
        c(1.5, -.1, 0),
        c(1.2, .2, 0),
        c(1.3, 0, 0)
      ),
      position_interpolation = 'cubic'
    ),
    clip_player('clip-0', offset = .7, speed = -1, loop = 'repeat')
  )
  settings = list(
    width = 64,
    height = 48,
    samples = 16,
    camera = camera(
      lookfrom = c(3.8, 1.5, 8),
      lookat = c(3.8, 1.2, 0),
      fov = 45,
      aperture = 0
    ),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    bloom = FALSE,
    denoise = FALSE,
    ambient_light = FALSE,
    backgroundhigh = 'black',
    backgroundlow = 'black',
    tonemap = 'raw',
    sample_method = 'sobol',
    parallel = FALSE
  )
  scene = unserialize(serialize(add_object(first, second), NULL))
  times = c(.25, 1, 1.75, .25)
  cpu = do.call(
    render_animation,
    c(list(scene = scene, times = times, integrator_type = 'nee'), settings)
  )
  expect_identical(as.numeric(cpu[[1]]), as.numeric(cpu[[4]]))
  expect_gt(mean(abs(cpu[[1]] - cpu[[2]])), .002)
  for (i in seq_along(times)) {
    expected = deformation_reference_scene(
      list(first, second),
      times[i],
      material
    )
    reference = do.call(
      render_scene,
      c(
        list(scene = expected, time = times[i], integrator_type = 'nee'),
        settings
      )
    )
    expect_lt(mean(abs(as.numeric(cpu[[i]]) - as.numeric(reference))), .005)
  }
  expect_message(
    do.call(
      render_scene,
      c(
        list(
          scene = scene,
          time = 0,
          shutter_offsets = c(0, .1),
          integrator_type = 'nee'
        ),
        settings
      )
    ),
    'Internal model poses are sampled at frame time'
  )
  gpu = do.call(
    render_animation,
    c(list(scene = scene, times = times, integrator_type = 'metal'), settings)
  )
  if (!isTRUE(attr(gpu[[1]], 'wavefront')$used)) {
    skip('Metal runtime unavailable')
  }
  for (i in seq_along(times)) {
    expect_true(attr(gpu[[i]], 'wavefront')$used)
    report = attr(gpu[[i]], 'wavefront')
    expect_identical(report$scene_reused, i > 1L)
    if (i > 1L) {
      expect_equal(report$acceleration_builds, 0)
      expect_gt(report$acceleration_refits, 0)
      expect_gt(report$geometry_bytes_updated, 0)
    }
    expect_lt(mean(abs(as.numeric(cpu[[i]]) - as.numeric(gpu[[i]]))), .025)
  }
})

test_that("placement-only animation reuses a frozen deformed mesh", {
  asset = read_gltf(deformation_fixture('combined'), material = light('white'))
  instance = animation_model_instance(
    asset,
    transform_track(
      c(0, 1, 2),
      position = rbind(c(-2, 0, 0), c(-1, 0, 0), c(-2, 0, 0))
    ),
    clip_player('clip-0', speed = 0, offset = .5)
  )
  frames = render_animation(
    instance,
    times = c(0, .5, 1),
    width = 24,
    height = 24,
    samples = 4,
    camera = camera(lookfrom = c(3, 2, 7), lookat = c(3, 1, 0)),
    integrator_type = 'metal',
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE
  )
  if (!isTRUE(attr(frames[[1]], 'wavefront')$used)) {
    skip('Metal runtime unavailable')
  }
  expect_false(attr(frames[[1]], 'wavefront')$scene_reused)
  expect_true(attr(frames[[2]], 'wavefront')$scene_reused)
  expect_true(attr(frames[[3]], 'wavefront')$scene_reused)
})

test_that("morph normals and emitter geometry update with changing triangle area", {
  material = light('white', intensity = 4)
  asset = read_gltf(deformation_fixture('bend'), material = material)
  instance = animation_model_instance(
    asset,
    transform_track(0, position = c(0, 3, 0), angle = c(90, 0, 0)),
    clip_player('clip-0')
  )
  floor = xz_rect(x = 3.5, xwidth = 5, zwidth = 5, material = diffuse('white'))
  settings = list(
    width = 32,
    height = 32,
    samples = 128,
    camera = camera(
      lookfrom = c(6, 5, 7),
      lookat = c(3.5, 0, 0),
      fov = 40,
      aperture = 0
    ),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = 'black',
    backgroundlow = 'black',
    tonemap = 'raw',
    sample_method = 'sobol',
    parallel = FALSE
  )
  cpu = do.call(
    render_animation,
    c(
      list(
        scene = add_object(instance, floor),
        times = c(0, 1, 2, 0),
        integrator_type = 'nee'
      ),
      settings
    )
  )
  expect_identical(as.numeric(cpu[[1]]), as.numeric(cpu[[4]]))
  expect_gt(mean(abs(cpu[[1]] - cpu[[3]])), .001)
  for (time in c(0, 1, 2)) {
    weight = time / 2
    surface = animation_model_surface(asset, clip_player('clip-0'), time)[[1]][[
      1
    ]]
    normal = c(-weight, 0, 1) / sqrt(1 + weight^2)
    expect_equal(
      surface$normals,
      matrix(rep(normal, 3), ncol = 3, byrow = TRUE),
      tolerance = 1e-6
    )
    reference = structure(
      list(
        vb = t(rbind(c(0, 0, 0), c(1, 0, weight), c(0, 1, 0))),
        it = matrix(1:3, 3),
        normals = matrix(rep(normal, 3), 3)
      ),
      class = 'mesh3d'
    )
    reference = mesh3d_model(
      reference,
      material = material,
      override_material = TRUE
    ) |>
      group_objects(translate = c(3, 3, 0), angle = c(90, 0, 0)) |>
      add_object(floor)
    image = do.call(
      render_scene,
      c(list(scene = reference, time = time, integrator_type = 'nee'), settings)
    )
    expect_lt(mean(abs(as.numeric(cpu[[time + 1]]) - as.numeric(image))), .02)
  }
})

test_that('missing binds default to identity and malformed influences fail', {
  asset = read_gltf(deformation_fixture('combined'))
  # Deliberately mutable native-input copies exercise validation below the API.
  asset = list2env(as.list.environment(asset), parent = emptyenv())
  asset$nodes[[1]]$deformation$skin['inverse_bind'] = list(NULL)
  positions = animation_model_surface(asset, clip_player(), 0)[[1]][[
    1
  ]]$positions
  expect_equal(positions[, 2], c(0, 0, 1) + 1.2, tolerance = 1e-6)
  for (kind in c('negative', 'zero', 'joint', 'nonfinite', 'unpaired')) {
    broken = list2env(as.list.environment(asset), parent = emptyenv())
    attrs = broken$nodes[[1]]$deformation$primitives[[1]]$attributes
    if (kind == 'negative') {
      attrs$WEIGHTS_0[1, 1] = -1
    }
    if (kind == 'zero') {
      attrs$WEIGHTS_0[,] = 0
      attrs$WEIGHTS_1[,] = 0
    }
    if (kind == 'joint') {
      attrs$JOINTS_0[1, 1] = 99
    }
    if (kind == 'nonfinite') {
      attrs$WEIGHTS_0[1, 1] = Inf
    }
    if (kind == 'unpaired') {
      attrs$WEIGHTS_1 = NULL
    }
    broken$nodes[[1]]$deformation$primitives[[1]]$attributes = attrs
    expect_error(
      animation_model_surface(broken, clip_player(), 0),
      'weight|attribute|WEIGHTS'
    )
  }
})

test_that('opt-in phase timing reports exact-pose reuse without changing images', {
  previous = options(rayrender.benchmark_timing = TRUE)
  on.exit(options(previous))
  asset = read_gltf(deformation_fixture('combined'), material = light('white'))
  instance = animation_model_instance(asset, animation = clip_player('clip-0'))
  frames = render_animation(
    instance,
    times = c(0, 1, 1, .25),
    width = 16,
    height = 16,
    samples = 2,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    camera = camera(lookfrom = c(3, 2, 7), lookat = c(3, 1, 0))
  )
  timing = lapply(frames, attr, 'animation_timing')
  expect_true(all(vapply(
    timing,
    function(x) {
      x$compile_seconds > 0 &&
        x$pose_seconds >= 0 &&
        x$refit_seconds >= 0
    },
    logical(1)
  )))
  expect_equal(
    vapply(timing, function(x) x$meshes_updated, numeric(1)),
    c(0, 1, 0, 1)
  )
  expect_identical(as.numeric(frames[[2]]), as.numeric(frames[[3]]))
})


test_that("posed normal recalculation survives serialization and reverse seeking", {
  for (invalid in list(NA, 1, "yes", logical(), c(TRUE, FALSE))) {
    expect_error(
      read_gltf(deformation_fixture("bend"), recalculate_normals = invalid),
      "recalculate_normals"
    )
  }
  asset = read_gltf(deformation_fixture("bend"), recalculate_normals = TRUE)
  # Deliberately remove the authored normal deltas: geometry alone defines the
  # expected slope, including after serialization and nonmonotonic seeks.
  nodes = asset$nodes
  nodes[[1]]$deformation$primitives[[1]]$targets[[1]]$NORMAL = NULL
  asset = animation_model_asset(nodes, asset$clips)
  asset = unserialize(serialize(asset, NULL))
  player = clip_player("clip-0")
  for (time in c(2, 0, 1, .25, 2)) {
    surface = animation_model_surface(asset, player, time)[[1]][[1]]
    expected = c(-time / 2, 0, 1) / sqrt(1 + (time / 2)^2)
    expect_equal(
      surface$normals,
      matrix(rep(expected, 3), ncol = 3, byrow = TRUE),
      tolerance = 1e-6
    )
  }
  # Exercise the dynamic-mesh binding and frame publication, including repeats.
  scene = animation_model_instance(asset, animation = player) |>
    add_object(sphere(
      x = 3.5,
      y = 2,
      z = 3,
      radius = .3,
      material = light(intensity = 10)
    ))
  frames = render_animation(
    scene,
    times = c(2, 0, 2),
    width = 16,
    height = 16,
    samples = 4,
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    parallel = FALSE,
    integrator_type = "nee",
    camera = camera(lookfrom = c(3.5, .5, 4), lookat = c(3.5, .5, 0))
  )
  expect_identical(as.numeric(frames[[1]]), as.numeric(frames[[3]]))
  expect_gt(mean(abs(frames[[1]] - frames[[2]])), .001)
})

test_that("frozen internal poses retain placement and unrelated object blur", {
  material = light("orange", intensity = 1)
  asset = animation_read_gltf(deformation_fixture("combined"), material)
  path = animation_transform_track(
    c(0, .5, 1, 2),
    position = rbind(c(-2, 0, 0), c(-.5, 0, 0), c(-2, .5, 0), c(-2, 0, 0)),
    position_interpolation = "cubic"
  )
  first = animation_model_instance(
    asset,
    path,
    animation_clip_player("clip-0", loop = "repeat")
  )
  second = animation_model_instance(
    asset,
    animation_transform_track(
      c(0, 1, 2),
      position = rbind(c(0, 1, 0), c(1, .7, 0), c(0, 1, 0))
    ),
    animation_clip_player("clip-0", speed = -1, offset = .7, loop = "repeat")
  )
  flyer = animation_keyframes(
    sphere(radius = .2, material = light("cyan", intensity = 1)),
    animation_transform_track(
      c(0, .5, 1, 2),
      position = rbind(c(0, -1, 0), c(5, -1, 0), c(1, -1, 0), c(0, -1, 0))
    )
  )
  scene = unserialize(serialize(
    add_object(add_object(first, second), flyer),
    NULL
  ))
  settings = list(
    width = 48,
    height = 40,
    samples = 32,
    camera = camera(lookfrom = c(3, 1, 10), lookat = c(3, 1, 0), fov = 45),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    tonemap = "raw",
    parallel = FALSE,
    sample_method = "sobol",
    min_variance = 0
  )
  times = c(.5, 1.8, -.1, .5)
  withr::local_options(warn = 2)
  for (backend in c("nee", "metal")) {
    # Deliberately exclude frame time from the shutter: pose time must still
    # be the requested frame, while placement follows the later exposure.
    offsets = c(.125, .625)
    messages = character()
    set.seed(617)
    frames = withCallingHandlers(
      do.call(
        render_animation,
        c(
          list(
            scene = scene,
            times = times,
            shutter_offsets = offsets,
            integrator_type = backend
          ),
          settings
        )
      ),
      message = function(m) {
        messages <<- c(messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      }
    )
    expect_equal(sum(grepl("Internal model poses", messages)), 1)
    if (backend == "metal" && !isTRUE(attr(frames[[1]], "wavefront")$used)) {
      skip("Metal runtime unavailable")
    }
    expect_identical(as.numeric(frames[[1]]), as.numeric(frames[[4]]))
    for (i in seq_along(times)) {
      # Bake the analytic combined-skin/morph fixture at frame time, then
      # reattach only the authored outer path. No native pose evaluator is
      # used to construct this reference geometry.
      expected = lapply(list(first, second), function(instance) {
        placement = instance$animation_info[[1]]$keyframes
        instance$animation_info[[1]]$keyframes = animation_transform_track()
        static = deformation_reference_scene(list(instance), times[i], material)
        animation_keyframes(static, placement)
      })
      expected = add_object(Reduce(add_object, expected), flyer)
      set.seed(617)
      reference = do.call(
        render_animation,
        c(
          list(
            scene = expected,
            times = times[i],
            shutter_offsets = offsets,
            integrator_type = backend
          ),
          settings
        )
      )
      expect_lt(
        mean(abs(as.numeric(frames[[i]]) - as.numeric(reference[[1]]))),
        .003
      )
      if (backend == "metal") expect_true(attr(frames[[i]], "wavefront")$used)
    }
    still = do.call(
      render_scene,
      c(
        list(scene = scene, time = times[1], integrator_type = backend),
        settings
      )
    )
    expect_gt(mean(abs(as.numeric(frames[[1]]) - as.numeric(still))), .005)
    set.seed(617)
    shifted = do.call(
      render_scene,
      c(
        list(
          scene = scene,
          time = .5,
          shutter_offsets = c(.25, .25),
          integrator_type = backend
        ),
        settings
      )
    )
    set.seed(617)
    exact = do.call(
      render_scene,
      c(list(scene = scene, time = .75, integrator_type = backend), settings)
    )
    expect_equal(as.numeric(shifted), as.numeric(exact), tolerance = 1e-6)
  }
})

test_that("skin-only, morph-only, and rigid models freeze their full internal hierarchy", {
  settings = list(
    width = 32,
    height = 32,
    samples = 16,
    camera = camera(lookfrom = c(4, 1, 10), lookat = c(4, 1, 0), fov = 40),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    tonemap = "raw",
    parallel = FALSE,
    sample_method = "sobol",
    min_variance = 0
  )
  for (kind in c("skin", "morph", "rigid")) {
    asset = if (kind == "rigid") {
      animation_model_asset(
        list(animation_model_node(
          "body",
          sphere(
            radius = .5,
            material = light("white", intensity = 1)
          )
        )),
        list(animation_model_clip(
          "clip-0",
          list(
            body = animation_transform_track(
              c(0, 1, 2),
              position = rbind(c(2, 0, 0), c(5, 1, 0), c(2, 0, 0))
            )
          )
        ))
      )
    } else {
      animation_read_gltf(
        deformation_fixture(kind),
        light("white", intensity = 1)
      )
    }
    instance = animation_model_instance(
      asset,
      animation = animation_clip_player("clip-0", loop = "repeat")
    )
    for (backend in c("nee", "metal")) {
      messages = character()
      set.seed(617)
      result = withCallingHandlers(
        do.call(
          render_animation,
          c(
            list(
              scene = instance,
              times = c(.5, 1.5, .5),
              shutter_offsets = c(-.4, .7),
              integrator_type = backend
            ),
            settings
          )
        ),
        message = function(m) {
          messages <<- c(messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      )
      expect_equal(sum(grepl("Internal model poses", messages)), 1)
      if (backend == "metal" && !isTRUE(attr(result[[1]], "wavefront")$used)) {
        skip("Metal runtime unavailable")
      }
      for (i in 1:2) {
        set.seed(617)
        reference = do.call(
          render_animation,
          c(
            list(
              scene = instance,
              times = c(.5, 1.5)[i],
              integrator_type = backend
            ),
            settings
          )
        )
        expect_equal(
          as.numeric(result[[i]]),
          as.numeric(reference[[1]]),
          tolerance = 1e-6
        )
      }
      expect_identical(as.numeric(result[[1]]), as.numeric(result[[3]]))
    }
  }
})
