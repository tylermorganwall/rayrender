#' Image-Based Sky Light
#' @md
#'
#' @description
#' Generate a cached sky EXR with `skymodelr::generate_sky_latlong()` (location
#' and time) or `skymodelr::generate_sky()` (direct Sun angles) and use it
#' as an infinite light. Add it to a scene with [add_infinite_light()]. The image
#' represents one observer and illuminates every scene position with the same
#' sky. Use [sky_light()] for altitude-dependent lighting and finite-distance haze.
#' Location/time skies include Sun and Moon by default, controlled with `sun` and
#' `moon`, and support optional stars and planets. Direct-angle skies include
#' the Sun and have no ephemeris-dependent celestial components.
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
#' @param moon Default `TRUE` for location/time skies, `FALSE` for direct skies. Composite a Moon image into the sky.
#'   Set `FALSE` when adding a separate [moon_light()].
#' @param moon_atmosphere Default `FALSE`. Include atmospheric scattering of moonlight.
#' @param moon_hosek Default `TRUE`. Use Hosek for moonlight scattering.
#'   Set `FALSE` to use Prague.
#' @param exr_adopted_white Default `"D60"`. Adopted white for EXR metadata:
#'   `"D60"`, `"D65"`, or numeric XYZ with Y = 1. Does not change image pixels.
#' @param exr_metadata Default `TRUE`. Attach skymodelr color metadata to the EXR.
#' @param environment_light_bake_white Default `FALSE`. Bake chromatic adaptation
#'   from the generated EXR's `white_current` metadata into its RGB pixels before
#'   using it for lighting. Requires `exr_metadata = TRUE`. The adapted image is
#'   cached and reused for still images and animations.
#' @param environment_light_bake_white_target Default `"D65"`. Target white point
#'   when baking: `"D50"`, `"D55"`, `"D60"`, `"D65"`, `"D75"`, `"E"`, or a finite
#'   numeric XYZ vector with positive Y (normalized to Y = 1). Unlike
#'   `exr_adopted_white`, this changes the image pixels when baking is enabled.
#' @param number_cores Default `1`. CPU threads used to generate the cached image.
#' @param verbose Default `FALSE`. Print sky-generation progress information.
#' @param ... Additional named arguments forwarded by
#'   `skymodelr::generate_sky_latlong()` to its star, planet, and Moon generators.
#'   These extra celestial settings apply only to location/time skies.
#'   Pass settings directly; location, datetime, and the cached filename are
#'   managed by this light. Native atmospheric controls belong to [sky_light()].
#'
#' @details Image generation happens before rendering and is cached for the R
#' session. Changing only `intensity`, `rotation`, or `name` reuses the image.
#' Changing the baked target white point reuses the generated sky and caches a
#' separate adapted image. White-balance controls belong to this light; pass the
#' light to [add_infinite_light()] before calling [render_scene()] or
#' [render_animation()].
#' Install any required Prague data with `skymodelr::download_sky_data()` first;
#' ordinary rendering does not download datasets. Pkgdown builds download
#' missing datasets automatically.
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
#' # Direct Sun angles also work with the other atmospheric controls.
#' generate_ground() |>
#'   add_object(sphere(y = 1)) |>
#'   add_infinite_light(sky_light_image(
#'     elevation = 25, azimuth = 135, turbidity = 4, resolution = 512
#'   )) |>
#'   render_scene(lookfrom = c(0, 2, -8), lookat = c(0, 1, 0),
#'                aperture = 0, samples = 32, iso = 5)
#'
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
  lat = NULL,
  long = NULL,
  datetime = NULL,
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
  environment_light_bake_white = FALSE,
  environment_light_bake_white_target = "D65",
  elevation = NULL,
  azimuth = NULL,
  ...
) {
  if (missing(moon) && (!is.null(elevation) || !is.null(azimuth))) {
    moon = FALSE
  }
  for (field in c("sun", "moon", "environment_light_bake_white")) {
    value = get(field)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(field, " must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (environment_light_bake_white) {
    if (!isTRUE(exr_metadata)) {
      stop(
        "environment_light_bake_white requires exr_metadata = TRUE.",
        call. = FALSE
      )
    }
    environment_light_white_xyz(environment_light_bake_white_target)
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
      "haze",
      "query_altitude",
      "haze_in_volumes",
      "deferred_haze",
      "haze_filter",
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
  result[c("elevation", "azimuth")] = list(elevation, azimuth)
  result["sun"] = list(sun)
  result$environment_light_bake_white = environment_light_bake_white
  result$environment_light_bake_white_target = environment_light_bake_white_target
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
  direct = !is.null(light$elevation)
  if (
    identical(light$sky_args$hosek, FALSE) ||
      (isTRUE(light$sky_args$moon_atmosphere) &&
        identical(light$sky_args$moon_hosek, FALSE))
  ) {
    prepare_pkgdown_sky_data(
      light$sky_args$altitude,
      light$sky_args$wide_spectrum
    )
  }
  generator = if (direct) {
    skymodelr::generate_sky
  } else {
    skymodelr::generate_sky_latlong
  }
  if (direct) {
    settings = light$sky_args
    settings[c(
      "moon",
      "stars",
      "planets",
      "moon_atmosphere",
      "moon_hosek",
      "star_width",
      "stars_exposure"
    )] = NULL
    args = c(
      list(elevation = light$elevation, azimuth = light$azimuth),
      settings
    )
  } else {
    args = c(
      list(lat = light$lat, lon = light$long, datetime = light$datetime),
      light$sky_args
    )
  }
  if (identical(light$sun, FALSE)) {
    args$render_mode = "atmosphere"
  }
  supported = names(formals(generator))
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
    do.call(generator, args)
    success = TRUE
  }

  # Preserve the generated sky so each white point can reuse it. Adapted pixels
  # have their own session cache and stay available across frames and renders.
  if (isTRUE(light$environment_light_bake_white)) {
    target_white = environment_light_white_xyz(
      light$environment_light_bake_white_target
    )
    saveRDS(
      list(
        image = basename(filename),
        target_white = target_white,
        rayimage_version = as.character(utils::packageVersion("rayimage"))
      ),
      key_file,
      version = 2
    )
    baked_filename = file.path(
      cache,
      paste0(unname(tools::md5sum(key_file)), "-white.exr")
    )
    if (!file.exists(baked_filename)) {
      white_balance = prepare_environment_light_white_balance(
        filename,
        environment_light_bake_white = TRUE,
        environment_light_bake_white_target = target_white
      )
      on.exit(unlink(white_balance$cleanup), add = TRUE)
      if (!file.copy(white_balance$environment_light, baked_filename)) {
        unlink(baked_filename)
        stop("Could not cache the white-balanced sky image.", call. = FALSE)
      }
    }
    filename = baked_filename
  }
  infinite_light(
    filename,
    intensity = light$intensity,
    rotation = light$rotation,
    name = light$name
  )
}
