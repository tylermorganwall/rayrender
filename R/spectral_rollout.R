#' Spectral Renderer Capabilities and Rollout Policy
#'
#' @return A list describing the current spectral renderer rollout contract.
#' @export
spectral_render_capabilities = function() {
  list(
    status = "compiler_input_preview",
    default_render_mode = "rgb_legacy",
    default_policy = list(
      decision = "keep_rgb_legacy_default",
      earliest_spectral_default = "after one stable opt-in release and accepted PR24 gates",
      legacy_escape_hatch = 'render_mode = "rgb_legacy"',
      remove_legacy_renderer_in_default_switch = FALSE,
      future_deprecation_plan = "separate post-stable-release plan"
    ),
    compiler_input = list(
      still = TRUE,
      animation_frame_zero = TRUE,
      direct_runtime_bridge = FALSE
    ),
    runtime_rendering = list(
      still = FALSE,
      animation = FALSE,
      preview = FALSE,
      progress = FALSE,
      cancellation = FALSE,
      adaptive_sampling = FALSE,
      denoising = FALSE,
      file_output = FALSE,
      debug_channels = FALSE
    ),
    aovs = spectral_aov_contract(),
    denoiser = spectral_denoiser_contract(),
    adaptive_sampling = spectral_adaptive_contract(),
    alpha = spectral_alpha_contract(),
    still_animation_shared = spectral_still_animation_contract(),
    warning_schedule = spectral_warning_schedule()
  )
}

spectral_aov_contract = function() {
  list(
    space = "Film/PixelSensor output-linear RGB unless noted",
    available_in_runtime = FALSE,
    required = c(
      "beauty",
      "alpha",
      "albedo",
      "normal",
      "depth",
      "emission",
      "variance"
    ),
    semantics = list(
      beauty = "sensor-integrated transported radiance",
      alpha = "coverage/transmittance, not wavelength-packet opacity",
      albedo = "visible-surface spectral albedo integrated through the selected sensor/output transform",
      normal = "geometric or shading normal in a documented coordinate system",
      depth = "camera-space or ray distance with documented units",
      emission = "sensor-integrated emitted contribution where supported",
      variance = "Film estimator metadata"
    )
  )
}

spectral_denoiser_contract = function() {
  list(
    available_in_runtime = FALSE,
    input_space = "output-linear RGB",
    beauty = "Film/PixelSensor output-linear RGB before display encoding",
    albedo = "same output-linear RGB basis as beauty",
    normal = "geometric vectors, not spectral samples",
    prohibited = c(
      "packet components",
      "wavelength PDF reweighting after Film accumulation",
      "sRGB-encoded denoiser inputs",
      "unrecorded clipping of negative sensor values"
    )
  )
}

spectral_adaptive_contract = function() {
  list(
    available_in_runtime = FALSE,
    decision_space = "linear sensor or output-linear RGB Film channels",
    uses_packet_components = FALSE,
    minimum_samples_required = TRUE,
    validation_required = "adaptive estimates must be statistically consistent with fixed sampling"
  )
}

spectral_alpha_contract = function() {
  list(
    available_in_runtime = FALSE,
    policy = "stochastic alpha cutout before Material closure and region transition",
    alpha_rejected_hit = list(
      changes_region = FALSE,
      changes_medium = FALSE,
      increments_depth = FALSE,
      preserves_mis_state = TRUE
    ),
    transmissive_dielectric_is_alpha = FALSE
  )
}

spectral_still_animation_contract = function() {
  list(
    shared = c(
      "scene schema/compiler",
      "spectrum and image caches",
      "RGB-to-spectrum tables",
      "Film/PixelSensor implementation",
      "Material/Light/Medium registries",
      "integrator code",
      "output color conversion",
      "region initialization logic"
    ),
    frame_varying = c(
      "camera transforms",
      "animated object transforms",
      "time-dependent lights/materials when supported",
      "per-frame output"
    )
  )
}

spectral_warning_schedule = function() {
  list(
    current_release = "spectral mode is opt-in and returns compiler input only",
    next_stable_release = "continue rgb_legacy default while collecting acceptance data",
    default_switch_candidate = "requires accepted PR24 gates and release notes",
    deprecation = "open a separate future plan after at least one stable spectral release"
  )
}

spectral_runtime_request = function(
  context,
  preview,
  progress,
  denoise,
  adaptive,
  debug_channel,
  filename,
  transparent_background
) {
  requested = list(
    context = context,
    preview = isTRUE(preview),
    progress = isTRUE(progress),
    denoise = isTRUE(denoise),
    adaptive_sampling = isTRUE(adaptive),
    debug_channel = debug_channel,
    file_output = !is.null(filename) && !is.na(filename),
    transparent_background = isTRUE(transparent_background)
  )
  unsupported = character()
  for (name in setdiff(names(requested), c("context", "debug_channel"))) {
    if (isTRUE(requested[[name]])) {
      unsupported = c(unsupported, name)
    }
  }
  if (!identical(debug_channel, "none")) {
    unsupported = c(unsupported, "debug_channel")
  }
  list(
    request = requested,
    unsupported_runtime_features = unsupported,
    runtime_bridge = "not_connected"
  )
}

spectral_runtime_unavailable = function(context) {
  stop(
    sprintf(
      "%s(render_mode = \"spectral\") builds schema-v2 compiler input only; the spectral runtime renderer is not connected yet. Set return_result = TRUE to inspect the compiler input.",
      context
    ),
    call. = FALSE
  )
}

spectral_default_integrator = function(integrator_type, max_depth) {
  spectral_max_depth = if (is.na(max_depth)) 50 else max_depth
  if (integrator_type %in% c("randomwalk", "random_walk")) {
    random_walk_integrator(max_depth = spectral_max_depth)
  } else {
    path_integrator(max_depth = spectral_max_depth)
  }
}

spectral_default_camera = function(
  fov,
  lookfrom,
  lookat,
  camera_up,
  aperture,
  focal_distance,
  shutteropen,
  shutterclose,
  ortho_dimensions
) {
  if (identical(fov, 0) || isTRUE(fov == 0)) {
    orthographic_camera(
      lookfrom = lookfrom,
      lookat = lookat,
      up = camera_up,
      ortho_dimensions = ortho_dimensions,
      shutteropen = shutteropen,
      shutterclose = shutterclose,
      initial_regions = character()
    )
  } else {
    perspective_camera(
      lookfrom = lookfrom,
      lookat = lookat,
      up = camera_up,
      fov = fov,
      aperture = aperture,
      focal_distance = focal_distance,
      shutteropen = shutteropen,
      shutterclose = shutterclose,
      initial_regions = character()
    )
  }
}

spectral_camera_motion_frame = function(camera_motion, frame) {
  required_columns = c(
    "x",
    "y",
    "z",
    "dx",
    "dy",
    "dz",
    "aperture",
    "fov",
    "focal",
    "orthox",
    "orthoy",
    "upx",
    "upy",
    "upz"
  )
  missing_columns = setdiff(required_columns, names(camera_motion))
  if (length(missing_columns) > 0) {
    schema_stop(
      "render_animation(camera_motion)",
      paste(
        "must contain columns",
        paste(sprintf("'%s'", required_columns), collapse = ", ")
      )
    )
  }
  if (!is.numeric(frame) || length(frame) != 1 || is.na(frame)) {
    schema_stop("render_animation(start_frame)", "must be a frame number")
  }
  frame = as.integer(frame)
  if (frame < 1 || frame > nrow(camera_motion)) {
    schema_stop("render_animation(start_frame)", "is outside camera_motion")
  }
  row = camera_motion[frame, , drop = FALSE]
  list(
    lookfrom = c(row$x, row$y, row$z),
    lookat = c(row$dx, row$dy, row$dz),
    camera_up = c(row$upx, row$upy, row$upz),
    aperture = row$aperture,
    fov = row$fov,
    focal_distance = row$focal,
    ortho_dimensions = c(row$orthox, row$orthoy)
  )
}

spectral_validate_render_descriptors = function(
  integrator,
  sampler,
  camera,
  film,
  spectral,
  validation,
  environment,
  context
) {
  if (!inherits(integrator, "ray_integrator")) {
    schema_stop(
      paste0(context, "(integrator)"),
      "must be a ray_integrator descriptor"
    )
  }
  if (!inherits(sampler, "ray_sampler")) {
    schema_stop(
      paste0(context, "(sampler)"),
      "must be a ray_sampler descriptor"
    )
  }
  if (!inherits(camera, "ray_camera")) {
    schema_stop(paste0(context, "(camera)"), "must be a ray_camera descriptor")
  }
  if (!inherits(film, "ray_film")) {
    schema_stop(paste0(context, "(film)"), "must be a ray_film descriptor")
  }
  if (!inherits(spectral, "ray_spectral_options")) {
    schema_stop(
      paste0(context, "(spectral)"),
      "must be a ray_spectral_options descriptor"
    )
  }
  if (!inherits(validation, "ray_scene_validation")) {
    schema_stop(
      paste0(context, "(validation)"),
      "must be a ray_scene_validation descriptor"
    )
  }
  if (
    !is.null(environment) &&
      (!inherits(environment, "ray_light") ||
        !isTRUE(attr(environment, "infinite")))
  ) {
    schema_stop(
      paste0(context, "(environment)"),
      "must be an infinite ray_light descriptor or NULL"
    )
  }
}

spectral_compile_render_input = function(
  scene,
  width,
  height,
  samples,
  fov,
  lookfrom,
  lookat,
  camera_up,
  aperture,
  focal_distance,
  shutteropen,
  shutterclose,
  ortho_dimensions,
  max_depth,
  integrator_type,
  filename,
  environment_light,
  rotate_env,
  intensity_env,
  preview,
  progress,
  denoise,
  min_variance,
  debug_channel,
  transparent_background,
  integrator = NULL,
  sampler = NULL,
  camera = NULL,
  film = NULL,
  spectral = spectral_options(),
  environment = NULL,
  validation = scene_validation(),
  context = "render_scene",
  frame_index = NULL
) {
  if (is.null(integrator)) {
    integrator = spectral_default_integrator(integrator_type, max_depth)
  }
  if (is.null(sampler)) {
    sampler = sobol_sampler(pixel_samples = samples)
  }
  if (is.null(camera)) {
    camera = spectral_default_camera(
      fov = fov,
      lookfrom = lookfrom,
      lookat = lookat,
      camera_up = camera_up,
      aperture = aperture,
      focal_distance = focal_distance,
      shutteropen = shutteropen,
      shutterclose = shutterclose,
      ortho_dimensions = ortho_dimensions
    )
  }
  if (is.null(film)) {
    film = rgb_film(
      width = width,
      height = height,
      filename = if (is.null(filename) || is.na(filename)) NULL else filename
    )
  }
  if (is.null(environment) && !is.null(environment_light)) {
    environment = image_infinite_light(
      filename = environment_light,
      scale = intensity_env,
      rotation = rotate_env
    )
  }
  spectral_validate_render_descriptors(
    integrator,
    sampler,
    camera,
    film,
    spectral,
    validation,
    environment,
    context
  )
  scene = legacy_scene_to_schema_v2(scene, validation = validation)
  scene_attrs = ray_scene_attrs(scene)
  scene_attrs$render_defaults = list(
    integrator = integrator,
    sampler = sampler,
    camera = camera,
    film = film,
    spectral = spectral,
    spectral_policy = spectral_render_capabilities(),
    runtime = spectral_runtime_request(
      context = context,
      preview = preview,
      progress = progress,
      denoise = denoise,
      adaptive = min_variance > 0,
      debug_channel = debug_channel,
      filename = filename,
      transparent_background = transparent_background
    ),
    frame_index = frame_index
  )
  if (!is.null(environment)) {
    scene_attrs$environment = environment
  }
  scene = restore_ray_scene_attrs(scene, scene_attrs)
  as_scene_compiler_input(scene, validation = validation)
}
