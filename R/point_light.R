#' Point and Spot Lights
#'
#' Geometry-free lights for next-event estimation. Point lights emit in every
#' direction; spot lights use PBRT v4's cubic smoothstep falloff in cosine space.
#' Illumination follows the inverse-square law and produces hard shadows.
#'
#' @param position Default `c(0, 1, 0)`. World-space light position.
#' @param color Default `"white"`. RGB light color.
#' @param intensity Default `1`. Radiant intensity multiplier, per steradian.
#'   A point light emits total RGB power `4 * pi * color * intensity`.
#' @param name Default `"point"` for point lights or `"spot"` for spot lights.
#'   Unique scene light name.
#' @param direction Default `c(0, -1, 0)`. World-space spotlight axis, pointing
#'   from the light toward the illuminated region.
#' @param cone_angle Default `30`. Outer cone half-angle in degrees, in `(0, 180]`.
#' @param falloff_angle Default `5`. Angular width of the falloff region inside
#'   `cone_angle`. Zero gives a hard cone edge.
#' @return A `ray_light` description for [add_light()]. Lights have no visible
#'   surface and automatically select the NEE integrator. Surface and volume
#'   scattering receive their illumination; specular caustics require a transport
#'   method that can connect through specular paths.
#' @details Attached participating media support point/spot illumination.
#'   The analytic sky atmosphere attenuates their shadow connections; its
#'   atmospheric in-scattering remains driven by its infinite lights.
#' @export
#' @md
#' @examples
#' scene = sphere() |>
#'   add_light(point_light(position = c(-3, 4, -2), intensity = 40))
#' scene = scene |>
#'   add_light(spot_light(position = c(3, 4, -2), direction = c(-3, -4, 2),
#'                       color = "lightblue", intensity = 80))
point_light = function(
  position = c(0, 1, 0),
  color = "white",
  intensity = 1,
  name = "point"
) {
  result = structure(
    list(
      type = "point",
      position = position,
      color = convert_color(color),
      intensity = intensity,
      name = name
    ),
    class = "ray_light"
  )
  validate_point_light(result)
  result
}

#' @rdname point_light
#' @export
spot_light = function(
  position = c(0, 1, 0),
  direction = c(0, -1, 0),
  cone_angle = 30,
  falloff_angle = 5,
  color = "white",
  intensity = 1,
  name = "spot"
) {
  result = structure(
    list(
      type = "spot",
      position = position,
      direction = direction,
      cone_angle = cone_angle,
      falloff_angle = falloff_angle,
      color = convert_color(color),
      intensity = intensity,
      name = name
    ),
    class = "ray_light"
  )
  validate_point_light(result)
  direction = direction / max(abs(direction))
  result$direction = direction / sqrt(sum(direction^2))
  result
}

#' Attach a Point or Spot Light
#'
#' Attach a named world-space light without adding emissive geometry. Use
#' [add_infinite_light()] for environment and sky lights. Adding objects preserves
#' attached lights; geometry transforms do not move these world-space lights.
#'
#' @param scene Scene to modify.
#' @param light Description from [point_light()] or `spot_light()`.
#' @param name Default `NULL`. Optional override for the light name.
#' @param replace Default `FALSE`. Replace an existing light with this name.
#' @return The modified scene.
#' @export
#' @md
add_light = function(scene, light, name = NULL, replace = FALSE) {
  if (!inherits(scene, "ray_scene")) {
    stop("`scene` must be a ray_scene.", call. = FALSE)
  }
  if (!is.null(name)) {
    light$name = name
  }
  validate_point_light(light)
  if (!is.logical(replace) || length(replace) != 1L || is.na(replace)) {
    stop("`replace` must be TRUE or FALSE.", call. = FALSE)
  }
  lights = ray_scene_point_lights(scene)
  if (!replace && light$name %in% names(lights)) {
    stop(
      "Light name already exists: ",
      light$name,
      ". Use replace = TRUE.",
      call. = FALSE
    )
  }
  lights[[light$name]] = light
  attr(scene, "ray_lights") = lights
  scene
}

#' Inspect or Remove Point and Spot Lights
#' @param scene Scene containing attached lights.
#' @param name Light name.
#' @return `get_light()` returns one light, `list_lights()` returns the named light
#'   list, and `remove_light()` returns the modified scene.
#' @export
#' @md
get_light = function(scene, name) {
  lights = ray_scene_point_lights(scene)
  if (
    !is.character(name) ||
      length(name) != 1L ||
      is.na(name) ||
      !name %in% names(lights)
  ) {
    stop("Unknown light name.", call. = FALSE)
  }
  lights[[name]]
}

#' @rdname get_light
#' @export
list_lights = function(scene) {
  ray_scene_point_lights(scene)
}

#' @rdname get_light
#' @export
remove_light = function(scene, name) {
  get_light(scene, name)
  lights = ray_scene_point_lights(scene)
  lights[[name]] = NULL
  attr(scene, "ray_lights") = if (length(lights)) lights else NULL
  scene
}

#' @param light Point/spot description.
#' @return NULL after validation.
#' @keywords internal
#' @noRd
validate_point_light = function(light) {
  if (!inherits(light, "ray_light") || !light$type %in% c("point", "spot")) {
    stop("Expected a point_light() or spot_light().", call. = FALSE)
  }
  for (field in c("position", "color", if (light$type == "spot") "direction")) {
    value = light[[field]]
    if (
      !is.numeric(value) ||
        length(value) != 3L ||
        any(!is.finite(value)) ||
        any(abs(value) > 1e30)
    ) {
      stop("Light ", field, " must contain three finite values.", call. = FALSE)
    }
  }
  if (any(light$color < 0 | light$color > 1)) {
    stop("Light color must lie in [0,1].", call. = FALSE)
  }
  if (
    !is.numeric(light$intensity) ||
      length(light$intensity) != 1L ||
      !is.finite(light$intensity) ||
      light$intensity < 0 ||
      light$intensity > 1e30
  ) {
    stop(
      "Light intensity must be a finite nonnegative scalar at most 1e30.",
      call. = FALSE
    )
  }
  if (
    !is.character(light$name) ||
      length(light$name) != 1L ||
      is.na(light$name) ||
      !nzchar(trimws(light$name))
  ) {
    stop("Light name must be a nonempty string.", call. = FALSE)
  }
  if (light$type == "spot") {
    if (max(abs(light$direction)) == 0) {
      stop("Spot direction cannot be zero.", call. = FALSE)
    }
    for (field in c("cone_angle", "falloff_angle")) {
      value = light[[field]]
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop("Spot angles must be finite scalars.", call. = FALSE)
      }
    }
    if (
      light$cone_angle <= 0 ||
        light$cone_angle > 180 ||
        light$falloff_angle < 0 ||
        light$falloff_angle > light$cone_angle
    ) {
      stop(
        "Require 0 < cone_angle <= 180 and 0 <= falloff_angle <= cone_angle.",
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

#' @param scene Scene, or NULL.
#' @return Validated named point/spot light list.
#' @keywords internal
#' @noRd
ray_scene_point_lights = function(scene) {
  lights = attr(scene, "ray_lights", exact = TRUE)
  if (is.null(lights)) {
    return(list())
  }
  if (
    !is.list(lights) ||
      is.null(names(lights)) ||
      anyNA(names(lights)) ||
      any(!nzchar(names(lights))) ||
      anyDuplicated(names(lights))
  ) {
    stop("ray_lights must be a uniquely named list.", call. = FALSE)
  }
  for (i in seq_along(lights)) {
    validate_point_light(lights[[i]])
    if (!identical(names(lights)[i], lights[[i]]$name)) {
      stop("Light names must match their scene entries.", call. = FALSE)
    }
  }
  lights
}
