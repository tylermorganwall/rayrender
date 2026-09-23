test_that("subsurface descriptors validate both parameter modes", {
  s = subsurface()
  expect_s3_class(s, "ray_material")
  m = s[[1]]$subsurface
  expect_s3_class(m, "ray_medium")
  expect_identical(m$sigma_a, rep(0, 3))
  expect_identical(m$sigma_s, rep(1, 3))
  expect_false(m$haze)
  expect_equal(
    subsurface(color = "black", radius = 2)[[1]]$subsurface$sigma_a,
    rep(.5, 3)
  )
  m = subsurface(color = c(.2, .4, .7), radius = c(1, 2, 3), scale = 2)[[
    1
  ]]$subsurface
  expect_equal(m$sigma_a + m$sigma_s, 1 / c(2, 4, 6))
  expect_equal(
    subsurface(sigma_a = 0, sigma_s = c(0, 1, 2))[[1]]$subsurface$sigma_s,
    c(0, 1, 2)
  )
  for (arg in c("color", "radius", "scale")) {
    args = list(sigma_a = 0, sigma_s = 1)
    args[[arg]] = if (arg == "color") "white" else 1
    expect_error(do.call(subsurface, args), "explicit")
  }
  expect_error(subsurface(sigma_a = 1), "both")
  expect_error(subsurface(sigma_s = 1), "both")
  for (arg in c("radius", "scale", "refraction")) {
    for (value in list(0, -1, NA_real_, Inf, "red", numeric())) {
      args = list()
      args[[arg]] = value
      expect_error(do.call(subsurface, args))
    }
  }
  for (arg in c("sigma_a", "sigma_s")) {
    for (value in list(-1, NA_real_, Inf, "red", c(1, 2), matrix(1, 1, 3))) {
      args = list(sigma_a = 0, sigma_s = 1)
      args[[arg]] = value
      expect_error(do.call(subsurface, args))
    }
  }
  expect_error(subsurface(g = 1), "g")
  expect_error(subsurface(g = -1), "g")
  expect_error(subsurface(roughness = 1.1), "roughness")
  expect_error(subsurface(color = c(1, 2, 3)), "color")
  expect_error(subsurface(color = c("red", "green")), "color")
  expect_error(subsurface(method = "fake"), "arg")
  expect_equal(subsurface(priority = 7)[[1]]$properties[[1]][8], 7)
  for (value in list(-1, .5, NA_real_, Inf, 2^31, "high", c(0, 1))) {
    expect_error(subsurface(priority = value), "priority")
  }
})

test_that("auto interiors preserve ownership through preparation and replacement", {
  prepare = rayrender:::prepare_subsurface
  object = sphere(material = subsurface())
  original = serialize(object, NULL)
  ready = prepare(object)
  expect_identical(serialize(object, NULL), original)
  expect_true(ready$shape_info[[1]]$medium_keep_surface)
  expect_identical(ready$shape_info[[1]]$medium_owner, "subsurface")
  expect_identical(prepare(ready), ready)
  expect_identical(unserialize(serialize(ready, NULL)), ready)
  expect_null(set_scene_material(ready, diffuse())$shape_info[[1]]$medium)
  explicit = set_medium(sphere(), homogeneous_medium())
  expect_error(
    prepare(set_scene_material(explicit, subsurface())),
    "explicit medium"
  )
  expect_error(set_medium(object, homogeneous_medium()), "conflicts")
  expect_identical(
    set_scene_material(explicit, diffuse("red"))$shape_info[[1]]$medium,
    explicit$shape_info[[1]]$medium
  )
  grouped = prepare(group_objects(object, scale = 2))
  expect_identical(grouped$shape_info[[1]]$medium, ready$shape_info[[1]]$medium)
  instanced = create_instances(object, x = c(-2, 2)) |>
    create_instances(y = c(-3, 3))
  expect_true(rayrender:::scene_medium_features(prepare(instanced))$attached)
  expect_false(
    rayrender:::scene_medium_features(prepare(set_scene_material(
      instanced,
      diffuse()
    )))$attached
  )
  replaced = set_scene_material(
    create_instances(sphere()),
    subsurface(sigma_a = .2, sigma_s = 2)
  )
  expect_true(rayrender:::scene_medium_features(prepare(replaced))$attached)
  expect_error(prepare(xy_rect(material = subsurface())), "closed")
  ready = rayrender:::prepare_scene_list(
    object,
    integrator_type = "basic",
    denoise = FALSE
  )
  expect_equal(ready$render_info$integrator_type, 1L)
})

sss_test_render = function(scene, samples = 64L, ...) {
  args = list(
    scene = scene,
    width = 6L,
    height = 6L,
    samples = samples,
    sample_method = "random",
    lookfrom = c(0, 0, 3),
    lookat = c(0, 0, 0),
    fov = 0,
    ortho_dimensions = c(.5, .5),
    aperture = 0,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
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
  do.call(render_scene, modifyList(args, list(...)))
}

test_that("dielectric priority controls overlap extinction from both directions", {
  glass_absorption = c(.1, .2, .3)
  milk_absorption = c(.2, .3, .4)
  for (milk_wins in c(FALSE, TRUE)) {
    glass_priority = if (milk_wins) 2 else 0
    scene = cube(
      z = .5,
      zwidth = 2,
      xwidth = 3,
      ywidth = 3,
      material = dielectric(
        refraction = 1,
        attenuation = glass_absorption,
        priority = glass_priority
      )
    ) |>
      add_object(cube(
        z = -.5,
        zwidth = 2,
        xwidth = 2,
        ywidth = 2,
        material = subsurface(
          sigma_a = milk_absorption,
          sigma_s = 0,
          refraction = 1,
          priority = 1
        )
      ))
    glass_distance = if (milk_wins) 1 else 2
    milk_distance = if (milk_wins) 2 else 1
    expected = exp(
      -glass_absorption * glass_distance - milk_absorption * milk_distance
    )
    for (side in c(-1, 1)) {
      set.seed(20260923)
      image = sss_test_render(
        scene,
        samples = 2048L,
        lookfrom = c(0, 0, 4 * side),
        ambient_light = TRUE,
        backgroundhigh = "white",
        backgroundlow = "white"
      )
      expect_equal(apply(image[,, 1:3], 3, mean), expected, tolerance = .015)
    }
    # Camera starts inside both overlapping solids.
    set.seed(20260924)
    image = sss_test_render(
      scene,
      samples = 2048L,
      lookfrom = c(0, 0, 0),
      lookat = c(0, 0, -1),
      ambient_light = TRUE,
      backgroundhigh = "white",
      backgroundlow = "white"
    )
    expected_inside = if (milk_wins) {
      exp(-1.5 * milk_absorption)
    } else {
      exp(-.5 * glass_absorption - milk_absorption)
    }
    expect_equal(
      apply(image[,, 1:3], 3, mean),
      expected_inside,
      tolerance = .015
    )
  }
})

test_that("a higher-priority glass excludes a hidden scattering body entirely", {
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
  image = sss_test_render(
    scene,
    samples = 8L,
    max_depth = 3L,
    ambient_light = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  expect_equal(as.numeric(image[,, 1:3]), rep(1, 6 * 6 * 3), tolerance = 1e-6)
})

test_that("overlapping SSS solids can exit out of entry order", {
  scene = cube(
    z = .5,
    zwidth = 2,
    xwidth = 3,
    ywidth = 3,
    material = subsurface(
      sigma_a = .1,
      sigma_s = 0,
      refraction = 1,
      priority = 0
    )
  ) |>
    add_object(cube(
      z = -.5,
      zwidth = 2,
      xwidth = 2,
      ywidth = 2,
      material = subsurface(
        sigma_a = .4,
        sigma_s = 0,
        refraction = 1,
        priority = 1
      )
    ))
  for (side in c(-1, 1)) {
    set.seed(4192)
    image = sss_test_render(
      scene,
      samples = 2048L,
      lookfrom = c(0, 0, side * 4),
      ambient_light = TRUE,
      backgroundhigh = "white",
      backgroundlow = "white"
    )
    expect_equal(mean(image[,, 1:3]), exp(-.6), tolerance = .012)
  }
})

test_that("priority-resolved glass and SSS conserve a uniform white furnace", {
  for (glass_priority in c(0, 2)) {
    scene = sphere(
      z = .4,
      radius = .8,
      material = dielectric(refraction = 1.5, priority = glass_priority)
    ) |>
      add_object(sphere(
        material = subsurface(
          sigma_a = 0,
          sigma_s = 3,
          refraction = 1.333,
          priority = 1
        )
      ))
    set.seed(7123)
    image = sss_test_render(
      scene,
      samples = 32L,
      max_depth = 128L,
      ambient_light = TRUE,
      backgroundhigh = "white",
      backgroundlow = "white"
    )
    expect_equal(as.numeric(image[,, 1:3]), rep(1, 108), tolerance = 2e-5)
  }
})

test_that("vacuum and absorbing matched SSS slabs obey Beer-Lambert", {
  for (coefficients in list(c(0, 0, 0), c(0, .5, 1.5))) {
    scene = cube(
      material = subsurface(sigma_a = coefficients, sigma_s = 0, refraction = 1)
    ) |>
      add_object(xy_rect(
        z = -2,
        xwidth = 20,
        ywidth = 20,
        material = light(intensity = 1)
      ))
    set.seed(313)
    image = sss_test_render(scene, samples = 256L)
    expect_equal(
      apply(image[,, 1:3], 3, mean),
      exp(-coefficients),
      tolerance = .035
    )
    expect_true(all(is.finite(image)))
  }
})

test_that("independent mixture and scalar analog slab references agree", {
  rgb = sss_reference_slab(.2, 2, n = 4000L)
  analog = sss_analog_slab(.2, 2, n = 8000L)
  for (j in seq_len(2L)) {
    expect_lt(
      abs(rgb$mean[c(1L, 4L)[j]] - analog$mean[j + 1L]),
      5 * sqrt(rgb$se[c(1L, 4L)[j]]^2 + analog$se[j + 1L]^2) + .005
    )
  }
  vacuum = sss_reference_slab(0, 0, n = 2L)
  expect_equal(vacuum$mean, c(0, 0, 0, 1, 1, 1))
  conservative = sss_analog_slab(0, 2, n = 1000L)
  expect_equal(sum(conservative$mean[2:3]), 1)
})

test_that("internal SSS depth is independent of the ordinary surface budget", {
  material = subsurface(sigma_a = .1, sigma_s = 5, refraction = 1)
  scene = sphere(material = material) |>
    add_object(xy_rect(
      z = -3,
      xwidth = 40,
      ywidth = 40,
      material = light(intensity = 1)
    ))
  set.seed(229)
  shallow = sss_test_render(scene, max_depth = 1L)
  set.seed(229)
  deep = sss_test_render(scene, max_depth = 50L)
  expect_identical(shallow, deep)
})

test_that("a conservative smooth SSS body preserves a uniform white furnace", {
  image = sss_test_render(
    sphere(material = subsurface(sigma_a = 0, sigma_s = 8)),
    samples = 16L,
    backgroundhigh = "white",
    backgroundlow = "white",
    ambient_light = TRUE
  )
  expect_equal(as.numeric(image[,, 1:3]), rep(1, 108), tolerance = 2e-5)
})

test_that("SSS placement, nested cavities, embedded surfaces and motion use scene geometry", {
  m = subsurface(sigma_a = .2, sigma_s = 2, method = "guided")
  emitter = xy_rect(
    z = -3,
    xwidth = 20,
    ywidth = 20,
    material = light(intensity = 1)
  )
  bodies = list(
    sphere(material = m),
    cube(xwidth = 2, ywidth = 2, zwidth = .05, material = m),
    group_objects(sphere(material = m), scale = c(-1, 1, 1)),
    create_instances(sphere(material = m), scale_x = 1.4, scale_y = .7),
    sphere(material = m) |>
      add_object(sphere(radius = .3, material = diffuse("black"))),
    sphere(material = m) |>
      add_object(set_medium(
        sphere(radius = .4),
        homogeneous_medium(sigma_s = 0)
      )),
    sphere(material = m) |>
      add_object(sphere(
        radius = .4,
        material = subsurface(sigma_a = 0, sigma_s = 0, refraction = 1)
      )),
    animate_objects(sphere(material = m), end_position = c(.1, 0, 0)),
    set_medium(sphere(radius = 2), homogeneous_medium(sigma_s = .1)) |>
      add_object(sphere(material = m))
  )
  for (body in bodies) {
    set.seed(82)
    image = sss_test_render(add_object(body, emitter), samples = 8L)
    expect_true(all(is.finite(image)))
  }
  inside = sss_test_render(
    add_object(sphere(material = m), emitter),
    samples = 8L,
    lookfrom = c(0, 0, 0),
    lookat = c(0, 0, -1)
  )
  expect_true(all(is.finite(inside)))
  expect_gt(mean(inside[,, 1:3]), 0)
  blocked = cube(
    material = subsurface(sigma_a = 0, sigma_s = 0, refraction = 1)
  ) |>
    add_object(xy_rect(material = diffuse("black"))) |>
    add_object(emitter)
  blocked_image = sss_test_render(blocked, samples = 4L)
  expect_lt(mean(blocked_image[,, 1:3]), 1e-6)
  alpha = sss_test_render(
    add_object(bodies[[1]], emitter),
    samples = 4L,
    transparent_background = TRUE
  )
  expect_equal(as.numeric(alpha[,, 4]), rep(1, 36))
})

test_that("SSS raw transport preserves scale and random-sampler reproducibility", {
  for (method in c("random_walk", "guided")) {
    images = lapply(c(1, .1, 10), function(k) {
      scene = sphere(
        radius = k,
        material = subsurface(
          sigma_a = .2 / k,
          sigma_s = 2 / k,
          method = method
        )
      ) |>
        add_object(xy_rect(
          z = -3 * k,
          xwidth = 20 * k,
          ywidth = 20 * k,
          material = light(intensity = 1)
        ))
      set.seed(739)
      sss_test_render(
        scene,
        samples = 32L,
        lookfrom = c(0, 0, 3 * k),
        ortho_dimensions = c(.5, .5) * k
      )
    })
    for (j in 2:3) {
      expect_equal(
        as.numeric(images[[j]]),
        as.numeric(images[[1]]),
        tolerance = 2e-4
      )
    }
  }
  scene = sphere(material = subsurface(sigma_a = .1, sigma_s = 2))
  set.seed(834)
  a = sss_test_render(
    scene,
    samples = 8L,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  set.seed(834)
  b = sss_test_render(
    scene,
    samples = 8L,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  expect_identical(a, b)
  set.seed(835)
  c = sss_test_render(
    scene,
    samples = 8L,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  expect_false(identical(a, c))
})

test_that("camera rays on a mesh face classify both sides consistently", {
  file = tempfile(fileext = ".obj")
  on.exit(unlink(file))
  writeLines(
    c(
      "v -1 0 -1",
      "v 1 0 -1",
      "v 1 0 1",
      "v -1 0 1",
      "v -1 2.715 -1",
      "v 1 2.715 -1",
      "v 1 2.715 1",
      "v -1 2.715 1",
      "f 1 2 3",
      "f 1 3 4",
      "f 5 7 6",
      "f 5 8 7",
      "f 1 5 6",
      "f 1 6 2",
      "f 2 6 7",
      "f 2 7 3",
      "f 3 7 8",
      "f 3 8 4",
      "f 4 8 5",
      "f 4 5 1"
    ),
    file
  )
  scene = obj_model(
    file,
    load_material = FALSE,
    material = subsurface(
      sigma_a = 0,
      sigma_s = .5,
      refraction = 1
    )
  )
  for (origin in list(c(.2, 2.715, .1), c(-.75, 2.715, .67))) {
    for (side in c(-1, 1)) {
      image = sss_test_render(
        scene,
        samples = 4L,
        width = 3L,
        height = 3L,
        fov = 20,
        lookfrom = origin,
        lookat = origin + c(.317, side, .129),
        camera_up = c(0, 0, 1),
        ambient_light = TRUE,
        backgroundhigh = "white",
        backgroundlow = "white"
      )
      expect_equal(as.numeric(image[,, 1:3]), rep(1, 27), tolerance = 2e-5)
    }
  }
})

test_that("closed concave meshes, disconnected placements and open-boundary errors are supported", {
  m = subsurface(sigma_a = .2, sigma_s = 2, method = "guided")
  outline = rbind(c(-1, -1), c(1, -1), c(1, 0), c(0, 0), c(0, 1), c(-1, 1))
  concave = extruded_polygon(
    outline,
    plane = "xy",
    bottom = -.4,
    top = .4,
    material = m
  )
  emitter = xy_rect(
    z = -3,
    xwidth = 20,
    ywidth = 20,
    material = light(intensity = 1)
  )
  expect_true(all(is.finite(sss_test_render(
    add_object(concave, emitter),
    samples = 8L
  ))))
  file = tempfile(fileext = ".obj")
  on.exit(unlink(file))
  writeLines(
    c(
      "v -.5 -1 -.5",
      "v -.5 1 -.5",
      "v -.5 1 .5",
      "v .5 -1 -.5",
      "v .5 1 -.5",
      "v .5 1 .5",
      "f 1 3 2",
      "f 4 5 6",
      "f 1 2 5",
      "f 1 5 4",
      "f 2 3 6",
      "f 2 6 5",
      "f 3 1 4",
      "f 3 4 6"
    ),
    file
  )
  wedge = obj_model(file, material = m)
  expect_true(all(is.finite(sss_test_render(
    add_object(wedge, emitter),
    samples = 8L
  ))))
  copies = create_instances(sphere(radius = .3, material = m), x = c(-.35, .35))
  expect_true(all(is.finite(sss_test_render(
    add_object(copies, emitter),
    samples = 8L,
    ortho_dimensions = c(1.5, 1)
  ))))
  writeLines(c("v -1 -1 0", "v 1 -1 0", "v 0 1 0", "f 1 2 3"), file)
  expect_error(
    sss_test_render(obj_model(file, material = m), samples = 1L),
    "watertight|closed"
  )
  # Two disconnected closed tetrahedra within one mesh share one material.
  vertices = .3 * rbind(c(-1, -1, -1), c(1, -1, -1), c(0, 1, -1), c(0, 0, 1))
  left = right = vertices
  left[, 1] = left[, 1] - .5
  right[, 1] = right[, 1] + .5
  faces = rbind(c(1, 3, 2), c(1, 2, 4), c(2, 3, 4), c(3, 1, 4))
  writeLines(
    c(
      apply(rbind(left, right), 1, function(x) {
        paste("v", paste(x, collapse = " "))
      }),
      apply(rbind(faces, faces + 4), 1, function(x) {
        paste("f", paste(x, collapse = " "))
      })
    ),
    file
  )
  withr::local_envvar(RAYRENDER_VOLUME_STATS = "true")
  disconnected = sss_test_render(
    add_object(obj_model(file, material = m), emitter),
    samples = 8L,
    ortho_dimensions = c(2, 1)
  )
  expect_true(all(is.finite(disconnected)))
  expect_gt(attr(disconnected, "volume_statistics")$subsurface_events, 0)
})

test_that("guide fallback and internal-event statistics describe actual camera paths", {
  withr::local_envvar(RAYRENDER_VOLUME_STATS = "true")
  scene = sphere(
    material = subsurface(sigma_a = .2, sigma_s = 3, method = "guided")
  )
  inside = sss_test_render(
    scene,
    samples = 16L,
    lookfrom = c(0, 0, 0),
    lookat = c(0, 0, -1),
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  stats = attr(inside, "volume_statistics")
  expect_gt(stats$guide_fallback, 0)
  expect_equal(stats$guide_eligible, 0)
  expect_gte(stats$total_subsurface_events, stats$subsurface_events)
  expect_gte(
    stats$subsurface_events_p99_upper,
    stats$subsurface_events_p95_upper
  )
  expect_lte(stats$subsurface_events_p99_upper, 2 * stats$max_subsurface_events)
  scene = sphere(
    material = subsurface(sigma_a = 0, sigma_s = 3, method = "guided")
  )
  outside = sss_test_render(
    scene,
    samples = 8L,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  expect_equal(attr(outside, "volume_statistics")$guide_eligible, 0)
  expect_gt(attr(outside, "volume_statistics")$guide_fallback, 0)
})

test_that("SSS preserves the random sampler across thread counts", {
  withr::local_options(cores = 2L)
  scene = sphere(
    material = subsurface(sigma_a = .2, sigma_s = 3, method = "guided")
  )
  set.seed(909)
  serial = sss_test_render(
    scene,
    samples = 16L,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  set.seed(909)
  threaded = sss_test_render(
    scene,
    samples = 16L,
    parallel = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  expect_equal(as.numeric(serial), as.numeric(threaded), tolerance = 0)
})

test_that("anisotropic guidance falls back before sampling without changing HG transport", {
  withr::local_envvar(RAYRENDER_VOLUME_STATS = "true")
  images = lapply(c("random_walk", "guided"), function(method) {
    set.seed(671)
    sss_test_render(
      sphere(
        material = subsurface(
          sigma_a = .2,
          sigma_s = 3,
          g = .8,
          method = method
        )
      ),
      samples = 16L,
      backgroundhigh = "white",
      backgroundlow = "white"
    )
  })
  expect_identical(as.numeric(images[[1]]), as.numeric(images[[2]]))
  stats = attr(images[[2]], "volume_statistics")
  expect_equal(stats$guide_eligible, 0)
  expect_gt(stats$guide_fallback, 0)
})

test_that("subsurface uses the common preparation for stills, camera previews and animation", {
  scene = cube(
    material = subsurface(sigma_a = 0, sigma_s = 0, refraction = 1)
  ) |>
    add_object(xy_rect(
      z = -2,
      xwidth = 20,
      ywidth = 20,
      material = light(intensity = 1)
    )) |>
    add_camera(camera(
      lookfrom = c(0, 0, 3),
      lookat = c(0, 0, 0),
      fov = 0,
      ortho_dimensions = c(.5, .5),
      aperture = 0
    ))
  args = list(
    scene = scene,
    width = 4,
    height = 4,
    samples = 2,
    ambient_light = FALSE,
    min_variance = 0,
    tonemap = "raw",
    bloom = FALSE,
    denoise = FALSE,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE,
    integrator_type = "basic"
  )
  still = do.call(render_scene, args)
  animated = do.call(render_animation, args)
  camera_animation = do.call(render_scene, c(args, list(mode = "animation")))
  expect_equal(as.numeric(still), as.numeric(animated[[1]]), tolerance = 1e-6)
  expect_equal(
    as.numeric(still),
    as.numeric(camera_animation[[1]]),
    tolerance = 1e-6
  )
  expect_gt(mean(still), .99)
})

test_that("SSS preview, AO and denoising keep their established auxiliary policies", {
  scene = sphere(material = subsurface(sigma_a = .1, sigma_s = 2))
  args = list(
    scene = scene,
    width = 4,
    height = 4,
    samples = 2,
    preview = FALSE,
    plot_scene = FALSE,
    parallel = FALSE,
    progress = FALSE,
    denoise = FALSE,
    bloom = FALSE
  )
  expect_true(all(is.finite(do.call(render_preview, args))))
  local({
    grDevices::pdf(tempfile(fileext = ".pdf"))
    on.exit(grDevices::dev.off())
    expect_true(all(is.finite(render_ao(
      scene,
      width = 4,
      height = 4,
      samples = 1,
      parallel = FALSE,
      progress = FALSE
    ))))
  })
  set.seed(93)
  image = sss_test_render(
    scene,
    samples = 4L,
    denoise = TRUE,
    backgroundhigh = "white",
    backgroundlow = "white"
  )
  expect_true(all(is.finite(image)))
})
