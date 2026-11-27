test_that("mixed deforming emitters refit selectively and keep light sampling current", {
  asset = animation_read_gltf(
    system.file(
      "animation",
      "fixtures",
      "deform-combined.glb",
      package = "rayrender"
    ),
    material = light("white", intensity = 1),
    recalculate_normals = TRUE
  )
  copies = lapply(seq_len(10), function(i) {
    animation_model_instance(
      asset,
      animation_transform_track(
        0,
        position = c(-7 + (i - 1) %% 5 * 1.3, (i - 1) %/% 5, 0)
      ),
      animation_clip_player(
        "clip-0",
        loop = "repeat",
        speed = if (i == 1) 0 else 1
      )
    )
  })
  scene = Reduce(add_object, copies) |>
    add_object(xz_rect(
      y = -.1,
      xwidth = 30,
      zwidth = 30,
      material = diffuse("grey")
    )) |>
    add_object(sphere(x = 5, radius = .3, material = diffuse("red")))
  settings = list(
    scene = scene,
    width = 48,
    height = 48,
    samples = 32,
    camera = camera(lookfrom = c(0, 3, 9), lookat = c(0, 1, 0), fov = 45),
    preview = FALSE,
    progress = FALSE,
    plot_scene = FALSE,
    denoise = FALSE,
    bloom = FALSE,
    min_variance = 0,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    tonemap = "raw",
    sample_method = "sobol",
    parallel = FALSE,
    integrator_type = "metal"
  )
  times = c(0, .75, 1.5, .75, .75)
  set.seed(75)
  frames = do.call(render_animation, c(settings, list(times = times)))
  if (!isTRUE(attr(frames[[1]], "wavefront")$used)) {
    skip("Metal runtime unavailable")
  }
  expect_gt(attr(frames[[1]], "wavefront")$acceleration_builds, 0)
  for (i in 2:5) {
    report = attr(frames[[i]], "wavefront")
    expect_true(report$used)
    expect_true(report$scene_reused)
    expect_equal(report$acceleration_builds, 0)
    expect_equal(report$acceleration_refits, if (i == 5) 0 else 10)
    expect_equal(report$geometry_bytes_updated, if (i == 5) 0 else 9 * 160)
    set.seed(75)
    fresh = do.call(
      render_scene,
      c(settings, list(time = times[[i]], interactive = FALSE))
    )
    expect_true(attr(fresh, "wavefront")$used)
    # Refits keep the first tree partition; a fresh build can choose another
    # partition and traversal order, so compare radiance with sampling tolerance.
    expect_lt(mean(abs(as.numeric(frames[[i]]) - as.numeric(fresh))), .025)
  }
  expect_identical(as.numeric(frames[[4]]), as.numeric(frames[[5]]))
})
