#' Image-Based Location and Time Sky Light
#' @md
#'
#' @description
#' Generate a cached sky EXR with `skymodelr::generate_sky_latlong()` and use it
#' as an infinite light. Add it to a scene with [add_infinite_light()]. The image
#' represents one observer and illuminates every scene position with the same
#' sky. Use [sky_light()] for altitude-dependent lighting and finite-distance haze.
#' Both constructors include Sun and Moon by default, controlled with `sun` and
#' `moon`, and support optional stars and planets.
#'
#' @inheritParams sky_light
#' @param altitude Default `0`. Observer altitude in meters above sea level for
#'   the entire sky image. Prague supports 0--15000 m with its full-altitude data.
#' @param resolution Default `2048`. Height of the cached image in pixels;
#'   its width is twice the height.
#' @param hosek Default `TRUE`. Generate a Hosek sky. Set `FALSE` to use Prague.
#' @param turbidity Default `3`. Hosek turbidity, from 1.7 to 10.
#' @param wide_spectrum Default `FALSE`. Use Prague's 55-channel sea-level data.
#' @param below_horizon Default `TRUE`. Include atmospheric radiance below the horizon.
#' @param stars Default `FALSE`. Composite stars into the sky image.
#' @param star_width Default `1`. Stellar point-spread size, passed to
#'   `skymodelr::generate_stars()`.
#' @param stars_exposure Default `0`. Artistic exposure adjustment for stars, in stops.
#' @param planets Default `FALSE`. Composite bright planets into the sky image.
#' @param sun Default `TRUE`. Include the solar disk when selected by `render_mode`.
#'   Set `FALSE` to omit it from the image.
#' @param moon Default `TRUE`. Composite a Moon image into the sky.
#'   Set `FALSE` when adding a separate [moon_light()].
#' @param moon_atmosphere Default `FALSE`. Include atmospheric scattering of moonlight.
#' @param moon_hosek Default `TRUE`. Use Hosek for moonlight scattering.
#'   Set `FALSE` to use Prague.
#' @param exr_adopted_white Default `"D60"`. Adopted white for EXR metadata:
#'   `"D60"`, `"D65"`, or numeric XYZ with Y = 1. Does not change image pixels.
#' @param exr_metadata Default `TRUE`. Attach skymodelr color metadata to the EXR.
#' @param number_cores Default `1`. CPU threads used to generate the cached image.
#' @param verbose Default `FALSE`. Print sky-generation progress information.
#' @param ... Additional named arguments forwarded by
#'   `skymodelr::generate_sky_latlong()` to its star, planet, and Moon generators.
#'   Pass settings directly; location, datetime, and the cached filename are
#'   managed by this light. Native atmospheric controls belong to [sky_light()].
#'
#' @details Image generation happens before rendering and is cached for the R
#' session. Changing only `intensity`, `rotation`, or `name` reuses the image.
#' Install any required Prague data with `skymodelr::download_sky_data()` first;
#' rendering does not download datasets.
#'
#' With zero rotation, north is world +Z and east is world -X. Date, time, and
#' observer altitude remain fixed throughout the image. This light supports the
#' same integrators as [infinite_light()] and adds no finite atmospheric haze.
#' Use [render_scene()]'s `iso` to adjust exposure.
#'
#' For a separately sampled Sun, set `render_mode = "atmosphere"` and add
#' [sun_light()] with matching location, time, and atmospheric settings. This
#' avoids relying on the environment image to resolve the small solar disk.
#'
#' @return A `ray_infinite_light` containing a cached-image sky description.
#' @seealso [sky_light()], [infinite_light()], [sun_light()], [moon_light()]
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' time = as.POSIXct("2026-06-21 18:00:00", tz = "America/New_York")
#' scene = sphere(material = diffuse("white")) |>
#'   add_object(generate_ground())
#'
#' scene |>
#'   add_infinite_light(sky_light_image(40.7, -74, time)) |>
#'   render_scene(
#'     lookfrom = c(0, 1, 10),
#'     lookat = c(0, 0, 0),
#'     width = 600,
#'     height = 300,
#'     samples = 32,
#'     iso = 3
#'   )
#'
#' # A Prague image sky with an independently sampled solar disk.
#' # Install Prague data with skymodelr::download_sky_data() before rendering.
#' scene |>
#'   add_infinite_light(sky_light_image(
#'     40.7,
#'     -74,
#'     time,
#'     hosek = FALSE,
#'     render_mode = "atmosphere",
#'     altitude = 0,
#'     visibility = 50,
#'     albedo = 0.3
#'   )) |>
#'   add_infinite_light(sun_light(
#'     40.7,
#'     -74,
#'     time,
#'     sky_args = list(altitude = 0, visibility = 50, albedo = 0.3)
#'   )) |>
#' render_scene(
#'   lookfrom = c(0, 1, 10),
#'   lookat = c(0, 0, 0),
#'   aperture = 0,
#'   width = 600,
#'   height = 300,
#'   samples = 32,
#'   iso = 3
#' )
sky_light_image = function(
  lat,
  long,
  datetime,
  intensity = 1,
  rotation = 0,
  name = "sky",
  altitude = 0,
  visibility = 131.8,
  albedo = 0.5,
  resolution = 2048,
  hosek = TRUE,
  render_mode = "all",
  turbidity = 3,
  wide_spectrum = FALSE,
  below_horizon = TRUE,
  prague_rgb_correction = TRUE,
  prague_rgb_correction_strength = 1,
  prague_rgb_correction_gain = "auto",
  stars = FALSE,
  star_width = 1,
  stars_exposure = 0,
  planets = FALSE,
  moon = TRUE,
  moon_atmosphere = FALSE,
  moon_hosek = TRUE,
  exr_adopted_white = "D60",
  exr_metadata = TRUE,
  number_cores = 1,
  verbose = FALSE,
  sun = TRUE,
  ...
) {
  for (field in c("sun", "moon")) {
    value = get(field)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(field, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (!sun && identical(render_mode, "sun")) {
    stop('render_mode = "sun" requires sun = TRUE.', call. = FALSE)
  }
  extra = list(...)
  if ("sky_args" %in% names(extra)) {
    stop(
      "Pass sky settings directly to sky_light_image(), not in sky_args.",
      call. = FALSE
    )
  }
  native = intersect(
    names(extra),
    c(
      "atmosphere",
      "meters_per_unit",
      "atmosphere_origin",
      "attenuation",
      "query_altitude",
      "haze_in_volumes",
      "deferred_haze",
      "haze_correction_probability",
      "cache_spectra",
      "transmission_table",
      "transmission_table_max_mb"
    )
  )
  if (length(native)) {
    stop(
      "Use sky_light() for native atmospheric controls: ",
      paste(native, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  sky_args = c(
    list(
      altitude = altitude,
      visibility = visibility,
      albedo = albedo,
      resolution = resolution,
      hosek = hosek,
      render_mode = render_mode,
      turbidity = turbidity,
      wide_spectrum = wide_spectrum,
      below_horizon = below_horizon,
      prague_rgb_correction = prague_rgb_correction,
      prague_rgb_correction_strength = prague_rgb_correction_strength,
      prague_rgb_correction_gain = prague_rgb_correction_gain,
      stars = stars,
      star_width = star_width,
      stars_exposure = stars_exposure,
      planets = planets,
      moon = moon,
      moon_atmosphere = moon_atmosphere,
      moon_hosek = moon_hosek,
      exr_adopted_white = exr_adopted_white,
      exr_metadata = exr_metadata,
      number_cores = number_cores,
      verbose = verbose
    ),
    extra
  )
  result = new_sky_light(
    "sky_image",
    lat,
    long,
    datetime,
    sky_args,
    intensity,
    rotation,
    name
  )
  result["sun"] = list(sun)
  validate_infinite_light(result)
  result
}

#' @keywords internal
prepare_sky_light_image = function(light) {
  if (!requireNamespace("skymodelr", quietly = TRUE)) {
    stop(
      "sky_light_image() requires skymodelr. Install it with install.packages('skymodelr').",
      call. = FALSE
    )
  }
  args = c(
    list(lat = light$lat, lon = light$long, datetime = light$datetime),
    light$sky_args
  )
  if (identical(light$sun, FALSE)) {
    args$render_mode = "atmosphere"
  }
  supported = names(formals(skymodelr::generate_sky_latlong))
  unknown = setdiff(names(args), supported)
  if (length(unknown) && !"..." %in% supported) {
    stop(
      "Unsupported sky arguments: ",
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
