#' @keywords internal
native_sky_controls = function(lights) {
  index = which(vapply(
    lights,
    function(light) {
      light$type %in% c("sky", "sky_image")
    },
    logical(1)
  ))
  if (length(index) != 1L) {
    return(NULL)
  }
  sky = lights[[index]]
  model = if (isTRUE(sky$atmosphere) || identical(sky$sky_args$hosek, FALSE)) {
    1L
  } else {
    0L
  }
  initial_model = model
  position = native_sky_position(sky)
  list(
    index = as.integer(index - 1L),
    model = model,
    latitude = sky$lat,
    longitude = sky$long,
    datetime = format(sky$datetime, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    elevation = position$elevation,
    azimuth = position$azimuth,
    base_altitude = if (is.null(sky$sky_args$altitude)) {
      0
    } else {
      sky$sky_args$altitude
    },
    meters_per_unit = if (is.null(sky$meters_per_unit)) {
      1
    } else {
      sky$meters_per_unit
    },
    update = function(
      latitude,
      longitude,
      datetime,
      model = initial_model,
      elevation = NULL,
      azimuth = NULL,
      base_altitude = NULL,
      meters_per_unit = NULL,
      fast = FALSE
    ) {
      if (
        !is.finite(latitude) ||
          abs(latitude) > 90 ||
          !is.finite(longitude) ||
          abs(longitude) > 180
      ) {
        return(list(
          error = "Latitude must be -90 to 90; longitude -180 to 180."
        ))
      }
      time = as.POSIXct(strptime(datetime, "%Y-%m-%d %H:%M:%S", tz = "UTC"))
      if (
        is.na(time) ||
          !identical(format(time, "%Y-%m-%d %H:%M:%S", tz = "UTC"), datetime)
      ) {
        return(list(error = "Use a valid UTC date/time: YYYY-MM-DD HH:MM:SS."))
      }
      if (!model %in% 0:1) {
        return(list(error = "Choose Hosek or Prague."))
      }
      if (
        !is.null(base_altitude) &&
          (!is.numeric(base_altitude) ||
            length(base_altitude) != 1L ||
            !is.finite(base_altitude) ||
            base_altitude < 0 ||
            base_altitude > 15000)
      ) {
        return(list(
          error = "Base altitude must be between 0 and 15000 meters."
        ))
      }
      if (
        !is.null(meters_per_unit) &&
          (!is.numeric(meters_per_unit) ||
            length(meters_per_unit) != 1L ||
            !is.finite(meters_per_unit) ||
            meters_per_unit <= 0 ||
            meters_per_unit > 1e12)
      ) {
        return(list(
          error = "Atmosphere scale must be positive and at most 1e12 meters per unit."
        ))
      }
      manual = !is.null(elevation) || !is.null(azimuth)
      if (
        manual &&
          (length(elevation) != 1L ||
            length(azimuth) != 1L ||
            !is.finite(elevation) ||
            abs(elevation) > 90 ||
            !is.finite(azimuth) ||
            azimuth < 0 ||
            azimuth > 360)
      ) {
        return(list(error = "Invalid Sun elevation or azimuth."))
      }
      # Prepare replacements without changing the original scene or committed GUI
      # settings. File/data errors are shown in the inspector, leaving rendering live.
      tryCatch(
        {
          updated = lights
          next_sky = native_sky_model(sky, model)
          next_sky$lat = latitude
          next_sky$long = longitude
          next_sky$datetime = time
          # Reduce only the temporary sky map/proposal during a sun gesture.
          # The original description is unchanged, so release, undo and export
          # always rebuild with the scene's requested resolution.
          if (isTRUE(fast)) {
            resolution = next_sky$sky_args$resolution
            if (is.null(resolution)) {
              resolution = 512L
            }
            limit = if (isTRUE(next_sky$atmosphere)) 16L else 128L
            next_sky$sky_args$resolution = min(resolution, limit)
          }
          # Altitude affects ephemerides and Prague's RGB calibration as well as
          # transport. Set it before preparing any component of the replacement.
          if (!is.null(base_altitude)) {
            next_sky$sky_args$altitude = base_altitude
          }
          if (isTRUE(next_sky$atmosphere) && !is.null(meters_per_unit)) {
            next_sky$meters_per_unit = meters_per_unit
          }
          position = if (manual) {
            list(
              elevation = clamp_sky_sun_elevation(elevation),
              azimuth = azimuth
            )
          } else {
            native_sky_position(next_sky)
          }
          updated[[index]] = next_sky
          if (manual && !isTRUE(next_sky$atmosphere)) {
            updated = split_hosek_sky_lights(
              updated,
              index,
              position$elevation,
              position$azimuth
            )
            updated[[index]] = native_sky_manual_image(
              updated[[index]],
              position$elevation,
              position$azimuth
            )
          }
          prepared = prepare_scene_infinite_lights(updated)
          if (manual) {
            prepared = native_sky_apply_direction(
              prepared,
              index,
              position$elevation,
              position$azimuth
            )
          }
          list(
            error = "",
            lights = prepared,
            elevation = position$elevation,
            azimuth = position$azimuth
          )
        },
        error = function(error) list(error = conditionMessage(error))
      )
    }
  )
}

#' @keywords internal
native_sky_position = function(sky) {
  # skymodelr's public Prague metadata helper requires its coefficient dataset.
  # Use the same ephemeris routine as generate_sky_latlong so Hosek needs no data.
  ephemeris = utils::getFromNamespace("swe_dirs_topo_moon_sun", "skymodelr")
  altitude = if (is.null(sky$sky_args$altitude)) 0 else sky$sky_args$altitude
  direction = ephemeris(
    sky$datetime,
    sky$lat,
    sky$long,
    elev_m = altitude
  )$sun_dir_topo
  list(
    elevation = clamp_sky_sun_elevation(
      asin(pmax(-1, pmin(1, direction[3]))) * 180 / pi
    ),
    azimuth = (180 + atan2(direction[1], direction[2]) * 180 / pi) %% 360
  )
}

#' @keywords internal
native_sky_model = function(sky, model) {
  if (model == 1L) {
    if (isTRUE(sky$atmosphere)) {
      return(sky)
    }
    # The editor's Prague choice needs native transport for finite haze and
    # position-dependent lighting. Image resolution and Hosek-only parameters
    # do not apply to its importance-sampling tables or atmospheric queries.
    shared = intersect(names(sky$sky_args), names(formals(sky_light)))
    args = c(
      list(
        lat = sky$lat,
        long = sky$long,
        datetime = sky$datetime,
        intensity = sky$intensity,
        rotation = sky$rotation,
        name = sky$name
      ),
      sky$sky_args[shared]
    )
    if (!is.null(sky$sun)) {
      args$sun = sky$sun
    }
    return(do.call(sky_light, args))
  }
  if (!isTRUE(sky$atmosphere)) {
    sky$sky_args$hosek = model == 0L
    return(sky)
  }
  # A native Prague sky can temporarily use Hosek's image model. Switching back
  # restores its original atmosphere configuration, including volume boundaries.
  args = c(
    list(
      lat = sky$lat,
      long = sky$long,
      datetime = sky$datetime,
      intensity = sky$intensity,
      rotation = sky$rotation,
      name = sky$name
    ),
    sky$sky_args
  )
  args$resolution = 512L
  for (field in c(
    "sun",
    "moon",
    "stars",
    "star_width",
    "stars_exposure",
    "planets",
    "number_cores"
  )) {
    if (!is.null(sky[[field]])) args[[field]] = sky[[field]]
  }
  args$hosek = TRUE
  do.call(sky_light_image, args)
}

#' @keywords internal
clamp_sky_sun_elevation = function(elevation) {
  # Match MaxSunElevationDegrees in lights/sun_direction.h before generating
  # image skies, so neither sky model receives the singular zenith endpoint.
  pmin(elevation, 89.9)
}

#' @keywords internal
native_sky_manual_image = function(sky, elevation, azimuth) {
  elevation = clamp_sky_sun_elevation(elevation)
  args = sky$sky_args
  args$render_mode = if (identical(sky$sun, FALSE)) {
    "atmosphere"
  } else {
    args$render_mode
  }
  solar_args = args[intersect(
    names(args),
    names(formals(skymodelr::generate_sky))
  )]
  solar_args$elevation = elevation
  solar_args$azimuth = azimuth
  cache = file.path(tempdir(), "rayrender-editor-skies")
  dir.create(cache, showWarnings = FALSE)
  key = tempfile(tmpdir = cache)
  on.exit(unlink(key), add = TRUE)
  saveRDS(
    list(
      sky = sky,
      elevation = elevation,
      azimuth = azimuth,
      version = as.character(utils::packageVersion("skymodelr"))
    ),
    key,
    version = 2
  )
  filename = file.path(cache, paste0(unname(tools::md5sum(key)), ".exr"))
  if (!file.exists(filename)) {
    success = FALSE
    on.exit(if (!success) unlink(filename), add = TRUE)
    pixels = do.call(skymodelr::generate_sky, solar_args)
    if (isTRUE(sky$omit_solar_atmosphere)) {
      # A sun-only sky still keeps its own slot and optional celestial background.
      pixels[,, 1:3] = 0
    }
    if (any(vapply(args[c("moon", "stars", "planets")], isTRUE, logical(1)))) {
      # Move only the Sun and its scattering. Retain the Moon, stars and planets
      # from the scene's location/time instead of dropping them during a drag.
      original = native_sky_position(sky)
      original_args = solar_args
      original_args$elevation = original$elevation
      original_args$azimuth = original$azimuth
      old_sky = do.call(skymodelr::generate_sky, original_args)
      full_sky = do.call(
        skymodelr::generate_sky_latlong,
        c(list(datetime = sky$datetime, lat = sky$lat, lon = sky$long), args)
      )
      pixels[,, 1:3] = pixels[,, 1:3] + full_sky[,, 1:3] - old_sky[,, 1:3]
    }
    # Keep color metadata while omitting skymodelr-specific diagnostic tags
    # that libopenexr does not accept as file attributes.
    metadata = attr(pixels, "exr", exact = TRUE)
    attr(pixels, "exr") = metadata[intersect(
      names(metadata),
      c("chromaticities", "adoptedNeutral", "whiteLuminance", "envmap")
    )]
    rayimage::ray_write_image(pixels, filename, clamp = FALSE)
    if (isTRUE(sky$environment_light_bake_white)) {
      balance = prepare_environment_light_white_balance(
        filename,
        environment_light_bake_white = TRUE,
        environment_light_bake_white_target = environment_light_white_xyz(
          sky$environment_light_bake_white_target
        )
      )
      on.exit(unlink(balance$cleanup), add = TRUE)
      if (
        !identical(balance$environment_light, filename) &&
          !file.copy(balance$environment_light, filename, overwrite = TRUE)
      ) {
        stop("Could not cache the adjusted sky image.", call. = FALSE)
      }
    }
    success = TRUE
  }
  infinite_light(
    filename,
    intensity = sky$intensity,
    rotation = sky$rotation,
    name = sky$name
  )
}

#' @keywords internal
native_sky_apply_direction = function(lights, index, elevation, azimuth) {
  if (!identical(lights[[index]]$type, "prague")) {
    return(lights)
  }
  elevation = clamp_sky_sun_elevation(elevation)
  lights[[index]]$elevation = elevation
  lights[[index]]$azimuth = azimuth
  radians = pi / 180
  direction = c(
    -sin(azimuth * radians) * cos(elevation * radians),
    sin(elevation * radians),
    cos(azimuth * radians) * cos(elevation * radians)
  )
  for (i in seq_along(lights)) {
    if (
      identical(lights[[i]]$type, "disk") &&
        identical(lights[[i]]$radiance_spectrum, "sun")
    ) {
      lights[[i]]$direction = direction
      lights[[i]]$rotation = lights[[index]]$rotation
    }
  }
  lights
}

#' @keywords internal
native_sky_restore = function(lights, state) {
  controls = native_sky_controls(lights)
  if (is.null(controls)) {
    stop("Saved sky edits need exactly one editable sky.", call. = FALSE)
  }
  required = c(
    "model",
    "latitude",
    "longitude",
    "datetime",
    "manual",
    "elevation",
    "azimuth",
    "haze",
    "altitude"
  )
  if (!all(required %in% names(state))) {
    stop("Incomplete saved sky settings.", call. = FALSE)
  }
  result = controls$update(
    state$latitude,
    state$longitude,
    state$datetime,
    state$model,
    if (isTRUE(state$manual)) state$elevation else NULL,
    if (isTRUE(state$manual)) state$azimuth else NULL,
    base_altitude = state$base_altitude,
    meters_per_unit = state$meters_per_unit
  )
  if (nzchar(result$error)) {
    stop(result$error, call. = FALSE)
  }
  index = controls$index + 1L
  if (identical(result$lights[[index]]$type, "prague")) {
    result$lights[[index]]$haze = isTRUE(state$haze)
    result$lights[[index]]$query_altitude = isTRUE(state$altitude)
  }
  # Older exports have no numeric atmosphere controls; preserve their defaults.
  for (name in c(
    required,
    intersect(c("base_altitude", "meters_per_unit"), names(state))
  )) {
    controls[[name]] = state[[name]]
  }
  controls$elevation = result$elevation
  controls$azimuth = result$azimuth
  list(lights = result$lights, controls = controls)
}
