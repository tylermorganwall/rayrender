#' @keywords internal
sky_light_celestial_settings = function(light) {
  defaults = list(
    sun = TRUE,
    moon = TRUE,
    sun_resolution = 256,
    moon_resolution = 256,
    earthshine = TRUE,
    earthshine_albedo = 0.19,
    solar_irradiance_w_m2 = 1300,
    stars = FALSE,
    star_width = 1,
    stars_exposure = 0,
    planets = FALSE,
    celestial_resolution = 2048,
    number_cores = 1
  )
  values = utils::modifyList(
    defaults,
    light[intersect(names(defaults), names(light))],
    keep.null = TRUE
  )
  for (field in c("sun", "moon", "stars", "planets")) {
    value = values[[field]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(field, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!values$sun && identical(light$sky_args$render_mode, "sun")) {
    stop('render_mode = "sun" requires sun = TRUE.', call. = FALSE)
  }
  for (field in c(
    "sun_resolution",
    "moon_resolution",
    "celestial_resolution",
    "number_cores"
  )) {
    value = values[[field]]
    minimum = if (field == "number_cores") 1 else 16
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        value < minimum ||
        value != floor(value) ||
        value > .Machine$integer.max / 2
    ) {
      stop(
        field,
        " must be an integer of at least ",
        minimum,
        ".",
        call. = FALSE
      )
    }
  }
  for (field in c("star_width", "stars_exposure")) {
    value = values[[field]]
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        (field == "star_width" && value <= 0) ||
        (field == "stars_exposure" && !is.finite(2^value))
    ) {
      stop("Invalid ", field, ".", call. = FALSE)
    }
  }
  celestial_moon_args(values[c(
    "earthshine",
    "earthshine_albedo",
    "solar_irradiance_w_m2"
  )])
  values
}

#' @keywords internal
sky_light_celestial_lights = function(light) {
  controls = sky_light_celestial_settings(light)
  settings = prague_sky_settings(light$sky_args)
  sky_args = settings[c(
    "altitude",
    "visibility",
    "albedo",
    "prague_rgb_correction",
    "prague_rgb_correction_strength",
    "prague_rgb_correction_gain"
  )]
  sky_args$number_cores = controls$number_cores
  common = list(
    lat = light$lat,
    long = light$long,
    datetime = light$datetime,
    sky_args = sky_args,
    intensity = light$intensity,
    rotation = light$rotation
  )
  lights = list()
  if (controls$sun && settings$render_mode != "atmosphere") {
    lights$sun = do.call(
      sun_light,
      c(
        common,
        list(
          resolution = controls$sun_resolution,
          name = paste0(light$name, "::sun")
        )
      )
    )
  }
  if (controls$moon) {
    lights$moon = do.call(
      moon_light,
      c(
        common,
        list(
          resolution = controls$moon_resolution,
          name = paste0(light$name, "::moon"),
          moon_args = controls[c(
            "earthshine",
            "earthshine_albedo",
            "solar_irradiance_w_m2"
          )]
        )
      )
    )
  }
  lights
}

#' @keywords internal
prepare_sky_celestial_background = function(light) {
  # Stars and planets are distant emitters. Generate their full sphere without
  # extinction or a horizon mask; the native atmosphere applies both locally.
  controls = sky_light_celestial_settings(light)
  args = list(
    lat = light$lat,
    lon = light$long,
    datetime = light$datetime,
    altitude = prague_sky_settings(light$sky_args)$altitude,
    resolution = controls$celestial_resolution,
    number_cores = controls$number_cores,
    atmosphere_effects = FALSE,
    upper_hemisphere_only = FALSE
  )
  cache = file.path(tempdir(), "rayrender-sky-celestial-backgrounds")
  dir.create(cache, showWarnings = FALSE)
  key_file = tempfile(tmpdir = cache)
  on.exit(unlink(key_file), add = TRUE)
  saveRDS(
    list(
      args = args,
      controls = controls[c(
        "stars",
        "planets",
        "star_width",
        "stars_exposure"
      )],
      skymodelr = as.character(utils::packageVersion("skymodelr"))
    ),
    key_file,
    version = 2
  )
  filename = file.path(cache, paste0(unname(tools::md5sum(key_file)), ".exr"))
  if (!file.exists(filename)) {
    pixels = NULL
    if (controls$stars) {
      pixels = do.call(
        skymodelr::generate_stars,
        c(args, list(star_width = controls$star_width))
      )[,, 1:3] *
        2^controls$stars_exposure
    }
    if (controls$planets) {
      planets = do.call(
        skymodelr::generate_planets,
        c(args, list(planet_width = controls$star_width))
      )[,, 1:3]
      pixels = if (is.null(pixels)) planets else pixels + planets
    }
    pixels = rayimage::ray_read_image(
      pixels,
      normalize = FALSE,
      assume_colorspace = rayimage::CS_SRGB,
      source_linear = TRUE
    )
    success = FALSE
    on.exit(if (!success) unlink(filename), add = TRUE)
    rayimage::ray_write_image(pixels, filename, clamp = FALSE)
    success = TRUE
  }
  infinite_light(
    filename,
    intensity = light$intensity,
    rotation = light$rotation,
    name = paste0(light$name, "::background")
  )
}

#' @keywords internal
split_hosek_sky_lights = function(
  lights,
  index = NULL,
  elevation = NULL,
  azimuth = NULL
) {
  # Keep the editable sky in its original slot and append automatic emitters.
  # Explicit Sun lights replace the automatic disk, just as for native Prague.
  explicit_sun = any(vapply(
    lights,
    function(light) light$type == "sun",
    logical(1)
  ))
  indices = if (is.null(index)) seq_along(lights) else index
  for (i in indices) {
    sky = lights[[i]]
    if (
      !identical(sky$type, "sky_image") || identical(sky$sky_args$hosek, FALSE)
    ) {
      next
    }
    mode = sky$sky_args$render_mode
    if (is.null(mode)) {
      mode = "all"
    }
    include_sun = !identical(sky$sun, FALSE) && mode != "atmosphere"
    if (!include_sun) {
      next
    }
    if (!explicit_sun) {
      args = sky$sky_args[intersect(
        names(sky$sky_args),
        c("altitude", "albedo", "turbidity", "number_cores")
      )]
      args$hosek = TRUE
      if (!is.null(elevation)) {
        args$elevation = clamp_sky_sun_elevation(elevation)
        args$azimuth = azimuth
      }
      sun = sun_light(
        sky$lat,
        sky$long,
        sky$datetime,
        sky_args = args,
        resolution = 64,
        intensity = sky$intensity,
        rotation = sky$rotation,
        name = paste0(sky$name, "::sun")
      )
      # Apply the same optional white adaptation to the sky and its solar disk.
      sun$environment_light_bake_white = sky$environment_light_bake_white
      sun$environment_light_bake_white_target = sky$environment_light_bake_white_target
      sun$exr_adopted_white = sky$sky_args$exr_adopted_white
      lights = c(lights, list(sun))
    }
    sky$sun = FALSE
    sky$separate_sun = TRUE
    sky$omit_solar_atmosphere = identical(mode, "sun")
    lights[[i]] = sky
  }
  lights
}
