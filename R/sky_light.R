#' Location and Time Sky Light
#' @md
#'
#' @description
#' Creates a location-and-time sky light. Add it with [add_infinite_light()].
#' By default, `skymodelr::generate_sky_latlong()` generates a cached EXR.
#' Set `atmosphere = TRUE` to evaluate the Prague model throughout the scene,
#' including finite-distance atmospheric attenuation and in-scattering.
#'
#' @param lat Latitude in degrees, between -90 and 90.
#' @param long Longitude in degrees, between -180 and 180.
#' @param datetime A single `POSIXct` date and time. Specify its time zone when
#' constructing it with `as.POSIXct()`.
#' @param sky_args Default `list()`. Named arguments passed to
#' `skymodelr::generate_sky_latlong()`, such as `list(hosek = FALSE)` for the
#' Prague model, `resolution`, `altitude`, `moon`, or `stars`. Location, datetime,
#' and the generated filename are supplied by this light and cannot be overridden.
#' Atmospheric mode accepts the Prague settings described below.
#' @param intensity Default `1`. Nonnegative multiplier for this light's radiance.
#' @param rotation Default `0`. Additional rotation in degrees around the world
#' Y axis, using the same convention as [infinite_light()].
#' @param name Default `"sky"`. Unique light name within the scene.
#' @param atmosphere Default `FALSE`. Use native Prague queries at each path
#' location and include atmospheric attenuation and in-scattering over finite
#' distances. Requires `integrator_type = "nee"` and the full-altitude dataset.
#' @param meters_per_unit Default `1`. Physical meters per world-space unit when
#' `atmosphere = TRUE`. Applies to distances along every axis.
#' @param atmosphere_origin Default `c(0, 0, 0)`. World-space location of the
#' geographic reference point, at `sky_args$altitude` meters above sea level.
#' World +Y is up; horizontal offsets follow the model's spherical Earth.
#'
#' @details skymodelr is a required package dependency. Install any required
#' Prague datasets with `skymodelr::download_sky_data()` before rendering.
#' Atmospheric queries run through skymodelr's registered native API; its
#' implementation and coefficients are shared without calling R from workers.
#'
#' By default this is a static image-based sky: `sky_args$altitude` selects one
#' observer altitude for the entire environment. With `atmosphere = TRUE`, the
#' sky and Sun change with the location of each surface or cloud interaction.
#' Date and time stay fixed during an animation in both modes.
#' The image uses skymodelr's orientation: north at the image seam, east one
#' quarter across. With zero rotation, north is world +Z and east is world -X.
#' Other infinite lights add to the sky. Use `auto_exposure = TRUE` in
#' [render_scene()] when appropriate for the sky's physical radiance scale.
#'
#' Atmospheric mode uses the Prague model regardless of the static sky default.
#' Install its data with `skymodelr::download_sky_data(sea_level = FALSE)`.
#' A scene can have one atmospheric sky. Do not combine its haze with a medium
#' modeling the same atmospheric scattering or absorption: that would count the
#' atmosphere twice. Separate cloud volumes can be added normally.
#' The Sun is sampled as a disk independently of the sky sampling resolution.
#' Its visibility accounts for Earth's curvature, allowing high clouds to remain
#' sunlit after ground-level sunset. Atmospheric refraction is not modeled.
#'
#' Supported atmospheric `sky_args` are `altitude` (default 0 m, range
#' 0--15000), `visibility` (default 50 km, range 20--131.8), `albedo` (default
#' 0.5, range 0--1), `render_mode` (`"all"`, `"atmosphere"`, or `"sun"`),
#' and skymodelr's `prague_rgb_correction`, `prague_rgb_correction_strength`,
#' and `prague_rgb_correction_gain`. `resolution` (default 64) controls the
#' height of directional importance-sampling tables, not the rendered sky's
#' detail. `hosek`, `wide_spectrum`, `moon`, `stars`, and `planets` must be
#' `FALSE` when supplied; `below_horizon` must be `TRUE`. Solar elevations
#' outside -4.2--90 degrees are rejected. Queries outside the altitude range
#' use the nearest modeled altitude; keep scene interactions within that range.
#'
#' The model precomputes clear-air multiple scattering over a spherical Earth
#' with uniform ground albedo. Local geometry and clouds block direct Sun and
#' sky lighting, but do not cast shadows into this precomputed atmospheric
#' in-scattering field. Haze is disabled inside dielectric solids.
#' Radiance is integrated spectrally and converted to renderer RGB; attenuation
#' of RGB materials uses a broadband approximation. Finite-distance fitted
#' transmittance is normalized at zero distance and interpolated in optical
#' depth over the first 100 m. Ray-anchored cumulative transport prevents those
#' fitting errors from accumulating at each cloud null event or boundary.
#'
#' Additional image lights represent radiance outside the atmosphere and receive
#' atmospheric attenuation. Do not supply an image that already includes the
#' same haze. In atmospheric scenes, [sun_light()] and [moon_light()] automatically
#' request unattenuated textures from skymodelr. The renderer applies spectral
#' atmospheric filtering and Earth occlusion at each interaction, including
#' partial disks and the depressed horizon at altitude. An explicit Sun light
#' replaces this sky's built-in solar disk, preserving the sky and haze.
#' Without an explicit disk altitude, celestial ephemerides use this sky's
#' reference altitude. Light rotations and intensities remain independent;
#' match them when the disk should correspond to this sky's solar illumination.
#' The precomputed haze remains Sun-driven: a Moon disk illuminates surfaces and
#' clouds but does not add moonlit atmospheric in-scattering or a lunar halo.
#' With a transparent background,
#' atmospheric in-scattering remains foreground radiance, with scalar opacity
#' derived from primary-ray atmospheric transmittance. RGB transmission into
#' an arbitrary compositing background remains an approximation.
#'
#' @return A `ray_infinite_light` object containing the sky description.
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' if (requireNamespace("skymodelr", quietly = TRUE)) {
#'   scene = sphere(material = diffuse(color = "white")) |>
#'     add_infinite_light(sky_light(
#'       lat = 40.7, long = -74,
#'       datetime = as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York"),
#'       sky_args = list(resolution = 256)
#'     )) |>
#'     add_camera(camera(aperture = 0))
#'   render_scene(scene, integrator_type = "nee", auto_exposure = TRUE)
#' }
sky_light = function(
  lat,
  long,
  datetime,
  sky_args = list(),
  intensity = 1,
  rotation = 0,
  name = "sky",
  atmosphere = FALSE,
  meters_per_unit = 1,
  atmosphere_origin = c(0, 0, 0)
) {
  result = structure(
    list(
      type = "sky",
      lat = lat,
      long = long,
      datetime = datetime,
      sky_args = sky_args,
      intensity = intensity,
      rotation = rotation,
      name = name,
      atmosphere = atmosphere,
      meters_per_unit = meters_per_unit,
      atmosphere_origin = atmosphere_origin
    ),
    class = "ray_infinite_light"
  )
  validate_infinite_light(result)
  result
}

#' @keywords internal
validate_sky_light = function(light) {
  # Descriptions serialized before the atmosphere option retain the image path.
  if ("atmosphere" %in% names(light)) {
    if (
      !is.logical(light$atmosphere) ||
        length(light$atmosphere) != 1L ||
        is.na(light$atmosphere)
    ) {
      stop("atmosphere must be TRUE or FALSE.", call. = FALSE)
    }
  }
  for (field in c("lat", "long")) {
    value = light[[field]]
    limit = if (field == "lat") 90 else 180
    if (
      !is.numeric(value) ||
        length(value) != 1 ||
        !is.finite(value) ||
        abs(value) > limit
    ) {
      stop(
        "Sky light ",
        field,
        " must be a finite number between -",
        limit,
        " and ",
        limit,
        ".",
        call. = FALSE
      )
    }
  }
  if (
    !inherits(light$datetime, "POSIXct") ||
      length(light$datetime) != 1 ||
      !is.finite(as.numeric(light$datetime))
  ) {
    stop(
      "Sky light datetime must be one finite POSIXct date and time.",
      call. = FALSE
    )
  }
  args = light$sky_args
  if (
    !is.list(args) ||
      (length(args) &&
        (is.null(names(args)) ||
          anyNA(names(args)) ||
          any(!nzchar(names(args))) ||
          anyDuplicated(names(args))))
  ) {
    stop("sky_args must be a uniquely named list.", call. = FALSE)
  }
  reserved = intersect(
    names(args),
    c("lat", "lon", "long", "datetime", "filename", "allow_download")
  )
  if (length(reserved)) {
    stop(
      "sky_args cannot override: ",
      paste(reserved, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (isTRUE(light$atmosphere)) {
    validate_prague_sky_light(light)
  }
  invisible(TRUE)
}

#' @keywords internal
prepare_infinite_light = function(light) {
  validate_infinite_light(light)
  if (isTRUE(light$atmosphere)) {
    return(prepare_prague_sky_light(light))
  }
  if (light$type %in% c("image", "disk")) {
    return(light)
  }
  if (light$type %in% c("sun", "moon")) {
    return(prepare_celestial_light(light))
  }
  if (!requireNamespace("skymodelr", quietly = TRUE)) {
    stop(
      "sky_light() requires skymodelr. Install it with install.packages('skymodelr').",
      call. = FALSE
    )
  }
  args = c(
    list(lat = light$lat, lon = light$long, datetime = light$datetime),
    light$sky_args
  )
  supported = names(formals(skymodelr::generate_sky_latlong))
  unknown = setdiff(names(args), supported)
  if (length(unknown) && !"..." %in% supported) {
    stop(
      "Unsupported sky_args: ",
      paste(unknown, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if ("allow_download" %in% supported) {
    args$allow_download = FALSE
  }
  cache = file.path(tempdir(), "rayrender-skies")
  dir.create(cache, showWarnings = FALSE)
  key_file = tempfile(tmpdir = cache)
  on.exit(unlink(key_file), add = TRUE)
  saveRDS(
    list(
      args = args,
      version = as.character(utils::packageVersion("skymodelr"))
    ),
    key_file,
    version = 2
  )
  filename = file.path(cache, paste0(unname(tools::md5sum(key_file)), ".exr"))
  if (!file.exists(filename)) {
    success = FALSE
    on.exit(if (!success) unlink(filename), add = TRUE)
    args$filename = filename
    do.call(skymodelr::generate_sky_latlong, args)
    success = TRUE
  }
  infinite_light(
    filename,
    intensity = light$intensity,
    rotation = light$rotation,
    name = light$name
  )
}
