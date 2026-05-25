build_many_spheres_scene = function(settings) {
  sphere_count = as.integer(settings$sphere_count %||% 360)
  rng_seed = as.integer(settings$seed)
  set.seed(rng_seed)

  scene = rayrender::generate_ground(
    depth = -0.45,
    spheresize = 1000,
    material = rayrender::diffuse(color = "grey35", checkercolor = "grey70")
  )

  grid_side = ceiling(sqrt(sphere_count))
  x_values = seq(-5.5, 5.5, length.out = grid_side)
  z_values = seq(-5.5, 5.5, length.out = grid_side)
  positions = expand.grid(x = x_values, z = z_values)
  positions = positions[seq_len(sphere_count), , drop = FALSE]

  objects = vector("list", sphere_count)
  for (i in seq_len(sphere_count)) {
    jitter = runif(2, -0.12, 0.12)
    radius = runif(1, 0.08, 0.2)
    color = grDevices::hsv(
      h = (i %% 29) / 29,
      s = 0.45 + 0.35 * runif(1),
      v = 0.55 + 0.35 * runif(1)
    )
    material = if (i %% 11 == 0) {
      rayrender::metal(color = color, fuzz = 0.12)
    } else {
      rayrender::diffuse(color = color, sigma = 20)
    }
    objects[[i]] = rayrender::sphere(
      x = positions$x[[i]] + jitter[[1]],
      y = radius - 0.43,
      z = positions$z[[i]] + jitter[[2]],
      radius = radius,
      material = material
    )
  }

  scene = rayrender::add_object(scene, do.call(rbind, objects))
  scene = rayrender::add_object(
    scene,
    rayrender::sphere(
      x = -3.5,
      y = 6,
      z = -4,
      radius = 1.2,
      material = rayrender::light(intensity = 35, importance_sample = TRUE)
    )
  )
  rayrender::add_object(
    scene,
    rayrender::sphere(
      x = 4,
      y = 5,
      z = 2,
      radius = 0.8,
      material = rayrender::light(intensity = 20, importance_sample = TRUE)
    )
  )
}

render_many_spheres_scene = function(scene, settings) {
  render_args = list(
    scene = scene,
    width = settings$width,
    height = settings$height,
    samples = settings$samples,
    preview = FALSE,
    plot_scene = FALSE,
    progress = FALSE,
    denoise = FALSE,
    parallel = settings$threads > 1,
    ambient_light = FALSE,
    lookfrom = c(7, 3.2, -10),
    lookat = c(0, 0, 0),
    aperture = 0,
    fov = 34,
    clamp_value = 8,
    sample_method = "random",
    integrator_type = "nee",
    new_page = FALSE
  )
  call_with_supported_args(rayrender::render_scene, render_args)
}

run_benchmark = function(settings) {
  if (isTRUE(settings$time_build)) {
    scene_build_seconds = system.time({
      scene = build_many_spheres_scene(settings)
    })[["elapsed"]]
    render_seconds = system.time({
      image = render_many_spheres_scene(scene, settings)
    })[["elapsed"]]
    total_seconds = scene_build_seconds + render_seconds
  } else {
    scene = build_many_spheres_scene(settings)
    render_seconds = system.time({
      image = render_many_spheres_scene(scene, settings)
    })[["elapsed"]]
    scene_build_seconds = NA_real_
    total_seconds = render_seconds
  }
  attr(image, "render_seconds") = as.numeric(render_seconds)
  attr(image, "scene_build_seconds") = as.numeric(scene_build_seconds)
  attr(image, "bvh_build_seconds") = NA_real_
  attr(image, "total_seconds") = as.numeric(total_seconds)
  image
}

`%||%` = function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) {
    y
  } else {
    x
  }
}

call_with_supported_args = function(fun, args) {
  supported = names(formals(fun))
  do.call(fun, args[names(args) %in% supported])
}
