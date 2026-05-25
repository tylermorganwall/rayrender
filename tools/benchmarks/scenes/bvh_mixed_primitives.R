build_mixed_primitives_scene = function(settings) {
  object_count = as.integer(settings$object_count %||% 180)
  set.seed(as.integer(settings$seed) + 101)

  alpha_texture = matrix(rep(c(1, 0.25, 0.75, 1), length.out = 32 * 32), 32, 32)
  alpha_material = rayrender::diffuse(
    color = "white",
    alpha_texture = alpha_texture,
    sigma = 25
  )

  scene = rayrender::generate_ground(
    depth = -0.55,
    spheresize = 1000,
    material = rayrender::diffuse(color = "grey25", checkercolor = "grey55")
  )

  objects = vector("list", object_count)
  for (i in seq_len(object_count)) {
    x = runif(1, -4.8, 4.8)
    z = runif(1, -3.6, 4.8)
    selector = i %% 6
    material = rayrender::diffuse(
      color = grDevices::hsv((i %% 23) / 23, 0.55, 0.72),
      sigma = 35
    )
    objects[[i]] = switch(
      as.character(selector),
      "0" = rayrender::sphere(
        x = x,
        y = -0.25,
        z = z,
        radius = runif(1, 0.14, 0.32),
        material = material
      ),
      "1" = rayrender::cube(
        x = x,
        y = -0.25,
        z = z,
        width = runif(1, 0.22, 0.42),
        material = material,
        angle = c(0, runif(1, 0, 180), 0)
      ),
      "2" = rayrender::cylinder(
        x = x,
        y = -0.2,
        z = z,
        radius = runif(1, 0.08, 0.18),
        length = runif(1, 0.35, 0.75),
        material = material,
        angle = c(0, runif(1, 0, 180), runif(1, -25, 25))
      ),
      "3" = rayrender::disk(
        x = x,
        y = -0.05,
        z = z,
        radius = runif(1, 0.18, 0.34),
        material = if (i %% 18 == 3) alpha_material else material,
        angle = c(90, runif(1, 0, 180), 0)
      ),
      "4" = rayrender::triangle(
        v1 = c(x - 0.22, -0.42, z),
        v2 = c(x + 0.22, -0.42, z + 0.1),
        v3 = c(x, runif(1, -0.02, 0.28), z - 0.22),
        material = material
      ),
      rayrender::xy_rect(
        x = x,
        y = runif(1, -0.1, 0.45),
        z = z,
        xwidth = runif(1, 0.25, 0.6),
        ywidth = runif(1, 0.25, 0.6),
        material = if (i %% 17 == 5) alpha_material else material,
        angle = c(0, runif(1, 0, 180), 0)
      )
    )
  }

  scene = rayrender::add_object(scene, do.call(rbind, objects))
  scene = rayrender::add_object(
    scene,
    rayrender::xz_rect(
      x = 0,
      y = 5.5,
      z = 0.5,
      xwidth = 5,
      zwidth = 5,
      flipped = TRUE,
      material = rayrender::light(intensity = 16, importance_sample = TRUE)
    )
  )
  scene
}

render_mixed_primitives_scene = function(scene, settings) {
  rayrender::render_scene(
    scene,
    width = settings$width,
    height = settings$height,
    samples = settings$samples,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    parallel = settings$threads > 1,
    ambient_light = FALSE,
    lookfrom = c(6, 3.5, -8),
    lookat = c(0, -0.1, 0.3),
    aperture = 0,
    fov = 38,
    clamp_value = 8,
    sample_method = "random",
    integrator_type = "nee",
    new_page = FALSE
  )
}

run_benchmark = function(settings) {
  if (isTRUE(settings$time_build)) {
    elapsed = system.time({
      scene = build_mixed_primitives_scene(settings)
      image = render_mixed_primitives_scene(scene, settings)
    })[["elapsed"]]
  } else {
    scene = build_mixed_primitives_scene(settings)
    elapsed = system.time({
      image = render_mixed_primitives_scene(scene, settings)
    })[["elapsed"]]
  }
  attr(image, "render_seconds") = as.numeric(elapsed)
  image
}

`%||%` = function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    y
  } else {
    x
  }
}
