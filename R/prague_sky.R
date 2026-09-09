#' @importFrom skymodelr get_prague_sky_metadata
#' @keywords internal
prague_sky_settings = function(args) {
  defaults = list(
    altitude = 0,
    visibility = 50,
    albedo = 0.5,
    resolution = 64,
    render_mode = "all",
    hosek = FALSE,
    wide_spectrum = FALSE,
    number_cores = 1,
    prague_rgb_correction = TRUE,
    prague_rgb_correction_strength = 1,
    prague_rgb_correction_gain = "auto",
    below_horizon = TRUE,
    stars = FALSE,
    moon = FALSE,
    planets = FALSE,
    verbose = FALSE
  )
  unknown = setdiff(names(args), names(defaults))
  if (length(unknown)) {
    stop(
      "Unsupported atmospheric sky_args: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  result = utils::modifyList(defaults, args, keep.null = TRUE)
  for (flag in c("hosek", "wide_spectrum", "stars", "moon", "planets")) {
    if (!identical(result[[flag]], FALSE)) {
      stop(
        "Atmospheric sky lights require sky_args$",
        flag,
        " = FALSE.",
        call. = FALSE
      )
    }
  }
  if (!identical(result$below_horizon, TRUE)) {
    stop(
      "Atmospheric skies require below_horizon = TRUE for downward queries.",
      call. = FALSE
    )
  }
  for (field in c(
    "altitude",
    "visibility",
    "albedo",
    "resolution",
    "number_cores"
  )) {
    bounds = switch(
      field,
      altitude = c(0, 15000),
      visibility = c(20, 131.8),
      albedo = c(0, 1),
      resolution = c(16, 2048),
      number_cores = c(1, .Machine$integer.max)
    )
    value = result[[field]]
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.finite(value) ||
        value < bounds[1] ||
        value > bounds[2] ||
        (field %in% c("resolution", "number_cores") && value != floor(value))
    ) {
      stop("Invalid atmospheric sky_args$", field, ".", call. = FALSE)
    }
  }
  if (
    !is.character(result$render_mode) ||
      length(result$render_mode) != 1L ||
      !result$render_mode %in% c("all", "atmosphere", "sun")
  ) {
    stop(
      'sky_args$render_mode must be "all", "atmosphere", or "sun".',
      call. = FALSE
    )
  }
  result
}

#' @keywords internal
validate_prague_sky_light = function(light) {
  if (!identical(light$type, "sky")) {
    stop("The atmosphere option belongs to sky_light().", call. = FALSE)
  }
  scale = light$meters_per_unit
  if (
    !is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0
  ) {
    stop("meters_per_unit must be a finite positive number.", call. = FALSE)
  }
  origin = light$atmosphere_origin
  if (!is.numeric(origin) || length(origin) != 3L || any(!is.finite(origin))) {
    stop(
      "atmosphere_origin must be a finite world-space c(x, y, z) point.",
      call. = FALSE
    )
  }
  prague_sky_settings(light$sky_args)
  invisible(TRUE)
}

#' @keywords internal
prepare_prague_sky_light = function(light) {
  if (!requireNamespace("skymodelr", quietly = TRUE)) {
    stop("Atmospheric sky lights require the skymodelr package.", call. = FALSE)
  }
  if (!"get_prague_sky_metadata" %in% getNamespaceExports("skymodelr")) {
    stop(
      "Update skymodelr to a version exporting get_prague_sky_metadata().",
      call. = FALSE
    )
  }
  settings = prague_sky_settings(light$sky_args)
  metadata = do.call(
    skymodelr::get_prague_sky_metadata,
    c(
      list(
        datetime = light$datetime,
        lat = light$lat,
        lon = light$long
      ),
      settings[c(
        "altitude",
        "visibility",
        "albedo",
        "prague_rgb_correction",
        "prague_rgb_correction_strength",
        "prague_rgb_correction_gain"
      )]
    )
  )
  list(
    type = "prague",
    filename = metadata$filename,
    name = light$name,
    intensity = light$intensity,
    rotation = light$rotation,
    origin = unname(light$atmosphere_origin),
    meters_per_unit = light$meters_per_unit,
    altitude = settings$altitude,
    visibility = settings$visibility,
    albedo = settings$albedo,
    elevation = metadata$elevation_deg,
    azimuth = metadata$azimuth_deg,
    angular_diameter = metadata$angular_diameter_deg,
    rgb_gain = unname(metadata$rgb_gain),
    resolution = as.integer(settings$resolution),
    render_mode = settings$render_mode
  )
}

#' @keywords internal
prepare_scene_infinite_lights = function(lights) {
  atmospheric = which(vapply(
    lights,
    function(light) isTRUE(light$atmosphere),
    logical(1)
  ))
  if (!length(atmospheric)) {
    return(unname(lapply(lights, prepare_infinite_light)))
  }
  settings = prague_sky_settings(lights[[atmospheric]]$sky_args)
  explicit_sun = any(vapply(
    lights,
    function(light) light$type == "sun",
    logical(1)
  ))
  unname(lapply(lights, function(light) {
    if (light$type %in% c("sun", "moon")) {
      if (is.null(light$sky_args$altitude)) {
        light$sky_args$altitude = settings$altitude
      }
      return(prepare_celestial_light(light, atmospheric_attenuation = FALSE))
    }
    result = prepare_infinite_light(light)
    if (identical(result$type, "prague") && explicit_sun) {
      result$include_sun = FALSE
    }
    result
  }))
}
