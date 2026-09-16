#' @keywords internal
native_sky_controls = function(lights) {
  index = which(vapply(
    lights,
    function(light) isTRUE(light$atmosphere),
    logical(1)
  ))
  if (length(index) != 1L) {
    return(NULL)
  }
  sky = lights[[index]]
  list(
    latitude = sky$lat,
    longitude = sky$long,
    datetime = format(sky$datetime, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    update = function(latitude, longitude, datetime) {
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
      updated = lights
      updated[[index]]$lat = latitude
      updated[[index]]$long = longitude
      updated[[index]]$datetime = time
      list(error = "", lights = prepare_scene_infinite_lights(updated))
    }
  )
}
