#' Sun and Moon Infinite Lights
#' @md
#'
#' @description
#' Create detailed celestial disks without embedding them in a latitude-longitude
#' environment map. skymodelr supplies the position and apparent angular size for
#' the observer, date, and time. The Sun uses Prague solar radiance queries; the
#' Moon uses skymodelr's surface texture, phase, earthshine, and radiometry routines.
#' Add these descriptions to a scene with [add_infinite_light()].
#'
#' @param lat Latitude in degrees, between -90 and 90.
#' @param long Longitude in degrees, between -180 and 180.
#' @param datetime A single `POSIXct` date and time, with an explicit time zone.
#' @param sky_args Default `list()`. Named atmospheric settings for skymodelr's
#' Prague model: `altitude`, `visibility`, `albedo`, `wide_spectrum`,
#' `number_cores`, and the `prague_rgb_correction` options accepted by
#' `skymodelr::calculate_sky_values()`. `hosek = FALSE` is accepted;
#' `hosek = TRUE` is unsupported because these lights use per-direction queries.
#' Defaults match that function. Without a native atmospheric sky, altitude
#' describes one observer for the whole light. With [sky_light()]'s
#' `atmosphere = TRUE`, a missing altitude uses that sky's reference altitude
#' for ephemeris placement; atmospheric filtering uses each interaction's position.
#' @param resolution Default `256` for the Sun and `1024` for the Moon. Target
#' disk image width and height in pixels, at least 16. The Moon's padded image is
#' cropped without downsampling; edge coverage can add a few pixels. Texture
#' detail is evaluated directly at this resolution, independently of the sky map.
#' @param intensity Default `1`. Nonnegative multiplier on physical radiance.
#' @param rotation Default `0`. Rotation in degrees around world Y, with the
#' same convention as [infinite_light()].
#' @param name Default `"sun"` or `"moon"`. Unique light name within the scene.
#' @param moon_args Default `list()`. Named options for skymodelr's Moon routines:
#' `earthshine = TRUE`, `earthshine_albedo = 0.19`,
#' `solar_irradiance_w_m2 = 1300`, and `moon_extinction_kV = 0.172`.
#' Native atmospheric scenes ignore `moon_extinction_kV`, using Prague
#' transmission instead.
#'
#' @details Requires skymodelr and its Prague data. Install datasets with
#' `skymodelr::download_sky_data()` before rendering. Requires the public
#' `skymodelr::generate_sun_disk()` and `skymodelr::generate_moon_disk()` exports.
#'
#' Pair a Sun disk with `sky_light(..., sky_args = list(hosek = FALSE,
#' render_mode = "atmosphere"))` to exclude the sky map's rasterized Sun. Leave
#' `moon = FALSE` in that sky when adding a Moon disk. Lights add radiance; they
#' do not eclipse or occlude each other, and adding a second copy doubles its light.
#' With `sky_light(atmosphere = TRUE)`, an explicit Sun automatically replaces
#' the built-in Sun; there is no need to change `render_mode`. The clear-air sky
#' and haze retain that sky light's own location, time, rotation, and intensity.
#'
#' Both lights use north at world +Z, east at world -X, and up at world +Y,
#' matching [sky_light()]. They have no parallax and add no scene geometry.
#' The geometric horizon clips their emission. Preparation generates and caches
#' linear EXRs before worker threads start; phase, time, and position stay fixed
#' during animation. Direct illumination samples each disk's solid angle, so a
#' small apparent diameter does not depend on a high-resolution sky sampler.
#'
#' The Sun texture uses Prague's solar disk profile, mapped to the ephemeris
#' angular diameter. Moon radiance preserves skymodelr's phase-dependent total
#' irradiance and atmospheric extinction, with its spectral RGB and atmospheric
#' tint. Earthshine is part of the generated phase texture. These are radiance
#' images; exposure and tone mapping are applied only by the final render.
#'
#' With a native atmospheric sky, preparation automatically requests
#' `atmospheric_attenuation = FALSE` through skymodelr's public disk generators.
#' Sun radiance and lunar phase, surface detail, and earthshine are generated
#' before atmospheric filtering. The renderer applies the Sun's spectrum or the
#' Moon's reference spectrum to Prague transmission along each light direction.
#' Spherical Earth visibility clips each direction, so a partially set disk can
#' illuminate high clouds while being hidden from the ground. Texture detail
#' remains independent of the sky sampling resolution. Skymodelr must support
#' the new argument. Without a native atmosphere, the existing fixed-observer
#' attenuation and geometric horizon clipping remain in effect.
#' The precomputed clear-air haze is Sun-driven; lunar atmospheric in-scattering
#' and halos are not modeled. See [sky_light()] for the model's supported domain.
#'
#' @return A `ray_infinite_light` description.
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' time = as.POSIXct("2026-01-28 21:00:00", tz = "Pacific/Auckland")
#' scene = sphere(material = diffuse()) |>
#'   add_infinite_light(moon_light(-36.87593, 174.7647, time)) |>
#'   add_camera(camera(aperture = 0))
#' render_scene(scene, integrator_type = "nee", auto_exposure = TRUE)
sun_light = function(
  lat,
  long,
  datetime,
  sky_args = list(),
  resolution = 256,
  intensity = 1,
  rotation = 0,
  name = "sun"
) {
  light = sky_light(lat, long, datetime, sky_args, intensity, rotation, name)
  light$type = "sun"
  light$resolution = resolution
  validate_infinite_light(light)
  light
}

#' @rdname sun_light
#' @export
moon_light = function(
  lat,
  long,
  datetime,
  sky_args = list(),
  moon_args = list(),
  resolution = 1024,
  intensity = 1,
  rotation = 0,
  name = "moon"
) {
  light = sky_light(lat, long, datetime, sky_args, intensity, rotation, name)
  light$type = "moon"
  light$resolution = resolution
  light$moon_args = moon_args
  validate_infinite_light(light)
  light
}

#' @keywords internal
celestial_sky_args = function(args) {
  if (!is.null(args$hosek) && !identical(args$hosek, FALSE)) {
    stop(
      "Sun and Moon disk lights require Prague queries: use sky_args = list(hosek = FALSE).",
      call. = FALSE
    )
  }
  args$hosek = NULL
  defaults = list(
    altitude = 0,
    visibility = 50,
    albedo = 0.5,
    wide_spectrum = FALSE,
    number_cores = 1,
    prague_rgb_correction = TRUE,
    prague_rgb_correction_strength = 1,
    prague_rgb_correction_gain = "auto"
  )
  unknown = setdiff(names(args), names(defaults))
  if (length(unknown)) {
    stop(
      "Unsupported celestial sky_args: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  values = utils::modifyList(defaults, args, keep.null = TRUE)
  for (field in c("altitude", "visibility", "albedo", "number_cores")) {
    v = values[[field]]
    bounds = switch(
      field,
      altitude = c(0, 15000),
      visibility = c(20, 131.8),
      albedo = c(0, 1),
      number_cores = c(1, .Machine$integer.max)
    )
    if (
      !is.numeric(v) ||
        length(v) != 1 ||
        !is.finite(v) ||
        v < bounds[1] ||
        v > bounds[2] ||
        (field == "number_cores" && v != floor(v))
    ) {
      stop("Invalid celestial sky_args$", field, ".", call. = FALSE)
    }
  }
  if (
    !is.logical(values$wide_spectrum) ||
      length(values$wide_spectrum) != 1 ||
      is.na(values$wide_spectrum)
  ) {
    stop("sky_args$wide_spectrum must be TRUE or FALSE.", call. = FALSE)
  }
  values
}

#' @keywords internal
celestial_moon_args = function(args) {
  defaults = list(
    earthshine = TRUE,
    earthshine_albedo = 0.19,
    solar_irradiance_w_m2 = 1300,
    moon_extinction_kV = 0.172
  )
  if (
    !is.list(args) ||
      (length(args) &&
        (is.null(names(args)) ||
          anyNA(names(args)) ||
          any(!nzchar(names(args))) ||
          anyDuplicated(names(args))))
  ) {
    stop("moon_args must be a uniquely named list.", call. = FALSE)
  }
  unknown = setdiff(names(args), names(defaults))
  if (length(unknown)) {
    stop(
      "Unsupported moon_args: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  values = utils::modifyList(defaults, args, keep.null = TRUE)
  if (
    !is.logical(values$earthshine) ||
      length(values$earthshine) != 1 ||
      is.na(values$earthshine)
  ) {
    stop("moon_args$earthshine must be TRUE or FALSE.", call. = FALSE)
  }
  for (field in setdiff(names(defaults), "earthshine")) {
    v = values[[field]]
    if (
      !is.numeric(v) ||
        length(v) != 1 ||
        !is.finite(v) ||
        v < 0 ||
        (field == "solar_irradiance_w_m2" && v == 0)
    ) {
      stop("Invalid moon_args$", field, ".", call. = FALSE)
    }
  }
  values
}

#' @keywords internal
validate_celestial_light = function(light) {
  n = light$resolution
  if (
    !is.numeric(n) ||
      length(n) != 1 ||
      !is.finite(n) ||
      n < 16 ||
      n != floor(n) ||
      n > .Machine$integer.max / 2
  ) {
    stop(
      "Celestial resolution must be an integer of at least 16.",
      call. = FALSE
    )
  }
  celestial_sky_args(light$sky_args)
  if (light$type == "moon") {
    celestial_moon_args(light$moon_args)
  }
  invisible(TRUE)
}

#' @keywords internal
validate_celestial_disk = function(light) {
  d = light$direction
  if (
    !is.numeric(d) || length(d) != 3 || any(!is.finite(d)) || max(abs(d)) == 0
  ) {
    stop(
      "Celestial disk direction must be a finite nonzero three-vector.",
      call. = FALSE
    )
  }
  a = light$angular_diameter
  if (!is.numeric(a) || length(a) != 1 || !is.finite(a) || a <= 0 || a >= 180) {
    stop(
      "Celestial disk angular diameter must be between 0 and 180 degrees.",
      call. = FALSE
    )
  }
  if (
    !is.logical(light$clip_horizon) ||
      length(light$clip_horizon) != 1 ||
      is.na(light$clip_horizon)
  ) {
    stop("Celestial disk clip_horizon must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
celestial_disk = function(
  filename,
  direction,
  angular_diameter,
  intensity = 1,
  rotation = 0,
  name = "disk",
  clip_horizon = TRUE
) {
  light = infinite_light(filename, intensity, rotation, name)
  light$type = "disk"
  light$direction = direction / max(abs(direction))
  light$direction = light$direction / sqrt(sum(light$direction^2))
  light$angular_diameter = angular_diameter
  light$clip_horizon = clip_horizon
  validate_infinite_light(light)
  light
}

#' @keywords internal
prepare_celestial_light = function(light, atmospheric_attenuation = TRUE) {
  if (!requireNamespace("skymodelr", quietly = TRUE)) {
    stop(
      "Sun and Moon lights require skymodelr. Install it with install.packages('skymodelr').",
      call. = FALSE
    )
  }
  generator = paste0("generate_", light$type, "_disk")
  if (!generator %in% getNamespaceExports("skymodelr")) {
    stop(
      "Update skymodelr to a version exporting ",
      generator,
      "() to use this light.",
      call. = FALSE
    )
  }
  supported = names(formals(getExportedValue("skymodelr", generator)))
  if (
    !atmospheric_attenuation &&
      !any(c("atmospheric_attenuation", "...") %in% supported)
  ) {
    stop(
      "Update skymodelr to a version supporting atmospheric_attenuation = FALSE in ",
      generator,
      "().",
      call. = FALSE
    )
  }
  settings = celestial_sky_args(light$sky_args)
  if (!atmospheric_attenuation) {
    settings$atmospheric_attenuation = FALSE
  }
  cache = file.path(tempdir(), "rayrender-celestial-disk-api-v2")
  dir.create(cache, showWarnings = FALSE)
  key = light
  key$atmospheric_attenuation = atmospheric_attenuation
  key[c("intensity", "rotation", "name")] = NULL
  key_file = tempfile(tmpdir = cache)
  on.exit(unlink(key_file), add = TRUE)
  saveRDS(
    list(
      light = key,
      skymodelr = as.character(utils::packageVersion("skymodelr")),
      rayvertex = as.character(utils::packageVersion("rayvertex"))
    ),
    key_file,
    version = 2
  )
  stem = file.path(cache, unname(tools::md5sum(key_file)))
  filename = paste0(stem, ".exr")
  metadata = paste0(stem, ".rds")
  if (!file.exists(filename) || !file.exists(metadata)) {
    success = FALSE
    on.exit(if (!success) unlink(c(filename, metadata)), add = TRUE)
    args = c(
      list(
        datetime = light$datetime,
        lat = light$lat,
        lon = light$long,
        resolution = light$resolution
      ),
      settings
    )
    disk = if (light$type == "sun") {
      do.call(skymodelr::generate_sun_disk, args)
    } else {
      do.call(
        skymodelr::generate_moon_disk,
        c(args, celestial_moon_args(light$moon_args))
      )
    }
    validate_skymodelr_disk(disk)
    if (
      !atmospheric_attenuation &&
        !identical(disk$atmospheric_attenuation, FALSE)
    ) {
      stop(
        "skymodelr did not return an unattenuated celestial disk. Update skymodelr.",
        call. = FALSE
      )
    }
    azimuth = disk$azimuth_deg * pi / 180
    elevation = disk$elevation_deg * pi / 180
    direction = c(
      -sin(azimuth) * cos(elevation),
      sin(elevation),
      cos(azimuth) * cos(elevation)
    )
    diameter = disk$angular_diameter_deg
    pixels = disk$image
    # Explicit linear sRGB tagging avoids a color conversion or exposure change.
    pixels = rayimage::ray_read_image(
      pixels,
      normalize = FALSE,
      assume_colorspace = rayimage::CS_SRGB,
      source_linear = TRUE
    )
    rayimage::ray_write_image(pixels, filename, clamp = FALSE)
    saveRDS(list(direction = direction, diameter = diameter), metadata)
    success = TRUE
  }
  info = readRDS(metadata)
  result = celestial_disk(
    filename,
    info$direction,
    info$diameter,
    light$intensity,
    light$rotation,
    light$name,
    clip_horizon = atmospheric_attenuation
  )
  result$radiance_spectrum = if (atmospheric_attenuation) "rgb" else light$type
  result$atmospheric_attenuation = atmospheric_attenuation
  result
}

#' @keywords internal
validate_skymodelr_disk = function(disk) {
  fail = function() {
    stop(
      "skymodelr returned invalid disk image or geometry metadata.",
      call. = FALSE
    )
  }
  if (!is.list(disk) || !identical(disk$projection, "rectilinear")) {
    fail()
  }
  pixels = disk$image
  dims = dim(pixels)
  if (
    !is.numeric(pixels) ||
      length(dims) != 3L ||
      dims[3] != 3L ||
      any(dims[1:2] < 1) ||
      any(!is.finite(pixels))
  ) {
    fail()
  }
  for (field in c("azimuth_deg", "elevation_deg", "angular_diameter_deg")) {
    value = disk[[field]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) fail()
  }
  if (
    disk$azimuth_deg < 0 ||
      disk$azimuth_deg >= 360 ||
      abs(disk$elevation_deg) > 90 ||
      disk$angular_diameter_deg <= 0 ||
      disk$angular_diameter_deg >= 180
  ) {
    fail()
  }
}
