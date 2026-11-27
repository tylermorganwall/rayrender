test_that("native full-pose blur integrates independent node, skin and morph clips", {
  withr::local_options(rayrender.benchmark_timing = TRUE, cores = 2L)
  settings = list(
    width = 40,
    height = 32,
    samples = 64,
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
    sample_method = "sobol",
    min_variance = 0,
    integrator_type = "nee"
  )
  for (kind in c("skin", "morph", "combined", "rigid")) {
    asset = if (kind == "rigid") {
      animation_model_asset(
        list(animation_model_node(
          "body",
          sphere(radius = .5, material = light("white", intensity = 1))
        )),
        list(animation_model_clip(
          "move",
          list(
            body = animation_transform_track(
              c(0, 1, 2),
              position = rbind(c(2, 0, 0), c(5, 1, 0), c(2, 0, 0))
            )
          )
        ))
      )
    } else {
      read_gltf(
        system.file(
          "animation",
          "fixtures",
          paste0("deform-", kind, ".glb"),
          package = "rayrender"
        ),
        material = light("white", intensity = 1)
      )
    }
    clip = if (kind == "rigid") "move" else "clip-0"
    first = animation_model_instance(
      asset,
      animation_transform_track(
        c(0, 1, 2),
        position = rbind(c(-1, 0, 0), c(-.3, .1, 0), c(-1, 0, 0)),
        position_interpolation = "cubic"
      ),
      animation_clip_player(clip, loop = "repeat")
    )
    second = animation_model_instance(
      asset,
      animation_transform_track(
        c(0, .7, 2),
        position = rbind(c(1, 0, 0), c(.2, -.3, 0), c(1, 0, 0))
      ),
      animation_clip_player(clip, speed = -1, offset = .7, loop = "repeat")
    )
    scene = add_object(first, second)
    set.seed(621)
    frames = do.call(
      render_animation,
      c(
        list(
          scene = scene,
          times = c(.5, 1.7, -.1, .5),
          shutter_offsets = c(-.7, .7)
        ),
        settings
      )
    )
    expect_identical(as.numeric(frames[[1]]), as.numeric(frames[[4]]))
    timing = attr(frames[[1]], "animation_timing")
    expect_equal(timing$shutter_samples, settings$samples)
    expect_length(timing$shutter_times, settings$samples)
    expect_true(all(timing$shutter_times > -.2 & timing$shutter_times < 1.2))
    expect_equal(
      length(unique(floor((timing$shutter_times + .2) / 1.4 * 64))),
      64
    )
    expect_gt(timing$shutter_refit_seconds, 0)
    # Independently rendered exact-time images at the recorded shutter times.
    # Different spatial samples introduce edge noise, bounded below in mean RGB.
    reference_settings = settings
    reference_settings$samples = 8
    reference = do.call(
      render_animation,
      c(
        list(
          scene = scene,
          times = timing$shutter_times,
          shutter_offsets = c(0, 0)
        ),
        reference_settings
      )
    )
    averaged = Reduce(`+`, lapply(reference, as.numeric)) / length(reference)
    expect_true(
      mean(abs(as.numeric(frames[[1]]) - averaged)) < .015,
      info = kind
    )
    set.seed(621)
    frozen = suppressMessages(do.call(
      render_animation,
      c(
        list(
          scene = scene,
          times = .5,
          shutter_offsets = c(-.7, .7),
          animation_motion_blur = "placement"
        ),
        settings
      )
    ))[[1]]
    expect_true(
      mean(abs(as.numeric(frames[[1]]) - as.numeric(frozen))) > .001,
      info = kind
    )
    # A high adaptive threshold must not stop temporal integration early.
    adaptive_settings = settings
    adaptive_settings$min_variance = .5
    set.seed(621)
    adaptive = do.call(
      render_animation,
      c(
        list(
          scene = scene,
          times = .5,
          shutter_offsets = c(-.7, .7)
        ),
        adaptive_settings
      )
    )[[1]]
    expect_identical(as.numeric(adaptive), as.numeric(frames[[1]]))
  }
})

test_that("point exposures and ordinary scenes retain their sampling behavior", {
  withr::local_options(rayrender.benchmark_timing = TRUE)
  asset = animation_model_asset(list(animation_model_node("body", sphere())))
  scene = animation_model_instance(asset)
  settings = list(
    scene = scene,
    time = 1,
    shutter_offsets = c(.125, .125),
    width = 16,
    height = 16,
    samples = 4,
    integrator_type = "nee",
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    bloom = FALSE,
    denoise = FALSE,
    parallel = FALSE
  )
  set.seed(20)
  full = do.call(
    render_scene,
    c(settings, list(animation_motion_blur = "full"))
  )
  set.seed(20)
  placement = do.call(
    render_scene,
    c(settings, list(animation_motion_blur = "placement"))
  )
  expect_identical(as.numeric(full), as.numeric(placement))
  expect_equal(attr(full, "animation_timing")$shutter_samples, 0)
  settings$scene = sphere()
  settings$shutter_offsets = c(-.25, .25)
  set.seed(20)
  ordinary = do.call(render_scene, settings)
  expect_equal(attr(ordinary, "animation_timing")$shutter_samples, 0)
})

test_that("full-pose still, timeline and camera batches share camera and lighting time", {
  withr::local_options(rayrender.benchmark_timing = TRUE, cores = 2L)
  asset = animation_model_asset(
    list(
      animation_model_node(
        "body",
        sphere(radius = .5, material = diffuse("red"))
      ),
      animation_model_node(
        "lamp",
        sphere(y = 2, radius = .3, material = light("white", intensity = 4)),
        parent = "body"
      )
    ),
    list(animation_model_clip(
      "move",
      list(
        body = animation_transform_track(
          c(0, 1, 2),
          position = rbind(c(-1, 0, 0), c(1, 0, 0), c(-1, 0, 0))
        )
      )
    ))
  )
  scene = add_object(
    animation_model_instance(
      asset,
      animation = animation_clip_player("move", loop = "repeat")
    ),
    xz_rect(y = -.5, xwidth = 12, zwidth = 12, material = diffuse("white"))
  )
  motion = generate_camera_motion(
    rbind(c(-2, 3, 7), c(2, 3, 7)),
    lookats = rbind(c(0, .5, 0), c(0, .5, 0)),
    frames = 3
  )
  motion$time = c(0, 1, 2)
  scene = add_camera(scene, camera(motion = motion), name = "first")
  settings = list(
    scene = scene,
    time = .75,
    shutter_offsets = c(-.4, .4),
    width = 32,
    height = 32,
    samples = 64,
    integrator_type = "nee",
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    bloom = FALSE,
    denoise = FALSE,
    tonemap = "raw",
    min_variance = 0
  )
  set.seed(328)
  still = do.call(render_scene, settings)
  set.seed(328)
  timeline = do.call(render_scene, c(settings, list(timeline = c(0, 2))))
  expect_identical(as.numeric(still), as.numeric(timeline))
  timing = attr(still, "animation_timing")
  expect_equal(timing$shutter_samples, 64)
  reference_settings = settings
  reference_settings$time = NULL
  reference_settings$shutter_offsets = c(0, 0)
  reference_settings$samples = 8
  reference = do.call(
    render_animation,
    c(
      reference_settings,
      list(times = timing$shutter_times)
    )
  )
  averaged = Reduce(`+`, lapply(reference, as.numeric)) / length(reference)
  expect_lt(mean(abs(as.numeric(still) - averaged)), .03)
  expect_gt(mean(still[,, 1:3]), .001)
  settings$scene = add_camera(scene, camera(motion = motion), name = "second")
  batches = do.call(render_scene, c(settings, list(camera = "all")))
  expect_length(batches, 2)
  expect_true(all(vapply(
    batches,
    function(x) attr(x, "animation_timing")$shutter_samples == 64,
    logical(1)
  )))
  # A binary-exact exposure must not depend on the batch's first frame anchor.
  anchor_settings = settings
  anchor_settings$time = NULL
  anchor_settings$shutter_offsets = c(-1, 1) / 16
  set.seed(329)
  solo = do.call(render_animation, c(anchor_settings, list(times = .75)))
  set.seed(329)
  sequence = do.call(
    render_animation,
    c(anchor_settings, list(times = c(0, .75)))
  )
  expect_identical(as.numeric(solo[[1]]), as.numeric(sequence[[2]]))
  expect_identical(
    attr(solo[[1]], "animation_timing")$shutter_times,
    attr(sequence[[2]], "animation_timing")$shutter_times
  )
  settings$denoise = TRUE
  denoised = do.call(render_scene, settings)
  expect_true(all(is.finite(denoised)))
})

test_that("Metal full-pose exposures refit resident geometry and match CPU integration", {
  withr::local_options(rayrender.benchmark_timing = TRUE, cores = 2L)
  asset = read_gltf(
    system.file(
      "animation",
      "fixtures",
      "deform-combined.glb",
      package = "rayrender"
    ),
    material = light("white", intensity = 1)
  )
  scene = add_object(
    animation_model_instance(
      asset,
      animation = animation_clip_player("clip-0", loop = "repeat")
    ),
    animation_model_instance(
      asset,
      animation_transform_track(
        c(0, 1, 2),
        position = rbind(c(-2, 0, 0), c(-1, 0, 0), c(-2, 0, 0))
      ),
      animation_clip_player("clip-0", speed = -1, offset = .7, loop = "repeat")
    )
  )
  settings = list(
    scene = scene,
    times = c(.5, 1.5, .5),
    shutter_offsets = c(-.4, .4),
    width = 40,
    height = 32,
    samples = 64,
    camera = camera(lookfrom = c(3, 1, 10), lookat = c(3, 1, 0), fov = 45),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    bloom = FALSE,
    denoise = FALSE,
    tonemap = "raw",
    min_variance = 0
  )
  set.seed(500)
  gpu = withCallingHandlers(
    do.call(render_animation, c(settings, list(integrator_type = "metal"))),
    warning = function(w) {
      if (
        grepl(
          "^Metal wavefront: this Metal device does not support ray tracing",
          conditionMessage(w)
        )
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )
  if (!isTRUE(attr(gpu[[1]], "wavefront")$used)) {
    skip(
      "Metal runtime unavailable; full-pose device acceptance remains pending"
    )
  }
  set.seed(500)
  cpu = do.call(render_animation, c(settings, list(integrator_type = "nee")))
  for (i in seq_along(gpu)) {
    expect_lt(mean(abs(as.numeric(gpu[[i]]) - as.numeric(cpu[[i]]))), .015)
    expect_equal(
      attr(gpu[[i]], "animation_timing")$shutter_times,
      attr(cpu[[i]], "animation_timing")$shutter_times
    )
    info = attr(gpu[[i]], "wavefront")
    expect_equal(info$completed_samples, 64)
    expect_gt(info$acceleration_refits, 0)
    expect_gt(info$geometry_bytes_updated, 0)
    if (i > 1) expect_equal(info$acceleration_builds, 0)
  }
  expect_identical(as.numeric(gpu[[1]]), as.numeric(gpu[[3]]))
})
