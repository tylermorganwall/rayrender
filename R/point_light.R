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
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # point_light(): a miniature sundial with a starburst of hard shadows.
#' # One geometry-free source illuminates every direction; nearby pegs receive more light.
#' scene = generate_ground(depth = -.12, material = diffuse("#20354b")) |>
#'   add_object(cylinder(
#'     y = -.05,
#'     radius = 1.9,
#'     length = .1,
#'     material = diffuse("#e8d6af")
#'   ))
#' for (i in 0:11) {
#'   phi = i * pi / 6
#'   scene = add_object(
#'     scene,
#'     cylinder(
#'       x = 1.05 * cos(phi),
#'       z = 1.05 * sin(phi),
#'       y = .25,
#'       radius = .055,
#'       length = .5,
#'       material = diffuse("#b56d4c")
#'     )
#'   )
#' }
#' scene = add_light(
#'   scene,
#'   point_light(position = c(.2, 1.3, .1), intensity = 9, color = "#ffe0ac")
#' )
#' render_scene(
#'   scene,
#'   lookfrom = c(3, 5, 5),
#'   lookat = c(0, 0, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.6, 3.8),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
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
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # spot_light(): a hanging mobile makes a quacking duck only in shadow.
#' # A blocker at fraction t of the light-to-wall distance needs t times the
#' # desired shadow dimensions. Different t values scatter the toys in depth.
#' # Keep them near the projector, and view from the side to separate toys and shadow.
#' projector = c(-.5, 2.4, 4.5)
#' wall_z = -2
#' thread = diffuse("#c4c7cf")
#' hanging_height = 3.9
#' parts = data.frame(
#'   x = c(-.15, .48, .28, -.79, -.76),
#'   y = c(1.15, 1.9, 1.55, 1.46, 1.49),
#'   a = c(.75, .28, .18, .34, .30),
#'   b = c(.38, .28, .4, .09, .075),
#'   tilt = c(0, 0, 0, -25, -48),
#'   fraction = c(.28, .44, .36, .10, .18)
#' )
#' scene = generate_ground(depth = 0, material = diffuse("#39475c")) |>
#'   add_object(xy_rect(
#'     x = -.9,
#'     y = 1.8,
#'     z = wall_z,
#'     xwidth = 6.8,
#'     ywidth = 3.6,
#'     material = diffuse("#f5e7ce")
#'   ))
#' for (i in seq_len(nrow(parts))) {
#'   t = parts$fraction[i]
#'   p = projector + t * (c(parts$x[i], parts$y[i], wall_z) - projector)
#'   scene = add_object(
#'     scene,
#'     ellipsoid(
#'       x = p[1],
#'       y = p[2],
#'       z = p[3],
#'       a = t * parts$a[i],
#'       b = t * parts$b[i],
#'       c = .06,
#'       angle = c(0, 0, parts$tilt[i]),
#'       material = diffuse(c(
#'         "#e5ad48",
#'         "#e27160",
#'         "#639d99",
#'         "#617899",
#'         "#789cb5"
#'       )[i])
#'     )
#'   )
#'   scene = add_object(
#'     scene,
#'     segment(
#'       start = p,
#'       end = c(p[1], hanging_height, p[3]),
#'       radius = .0018,
#'       material = thread
#'     )
#'   )
#' }
#' # Two tilted blocks project the upper and lower halves of an open bill.
#' for (i in 1:2) {
#'   t = c(.37, .50)[i]
#'   p = projector + t * (c(.92, c(2.02, 1.76)[i], wall_z) - projector)
#'   scene = add_object(
#'     scene,
#'     cube(
#'       x = p[1],
#'       y = p[2],
#'       z = p[3],
#'       xwidth = .64 * t,
#'       ywidth = .10 * t,
#'       zwidth = .03,
#'       angle = c(0, 0, c(20, -20)[i]),
#'       material = diffuse("#efb947")
#'     )
#'   )
#'   scene = add_object(
#'     scene,
#'     segment(
#'       start = p,
#'       end = c(p[1], hanging_height, p[3]),
#'       radius = .0018,
#'       material = thread
#'     )
#'   )
#' }
#' # Hang the glass first, nearest the projector. Scale its dimensions with
#' # distance to preserve the water projection, and preserve its optical depth.
#' water_fraction = .04
#' water_scale = water_fraction / .30
#' water_width = 1.4 * water_scale
#' water_height = .22 * water_scale
#' water_target = c(.12, .83 - water_height / (2 * water_fraction), wall_z)
#' water = projector + water_fraction * (water_target - projector)
#' # Give the glass a closed, thin volume so absorption colors transmitted light.
#' # Match air's IOR to project a blue filter without bending the spotlight rays.
#' scene = add_object(
#'   scene,
#'   cube(
#'     x = water[1],
#'     y = water[2],
#'     z = water[3],
#'     xwidth = water_width,
#'     ywidth = water_height,
#'     zwidth = .035 * water_scale,
#'     material = dielectric(
#'       refraction = 1,
#'       attenuation = c(30, 10, 1) / water_scale
#'     )
#'   )
#' )
#' for (x in water[1] + c(-1, 1) * water_width / 2) {
#'   scene = add_object(
#'     scene,
#'     segment(
#'       start = c(x, water[2] + water_height / 2, water[3]),
#'       end = c(x, hanging_height, water[3]),
#'       radius = .0025 * water_scale,
#'       material = thread
#'     )
#'   )
#' }
#' scene = add_light(
#'   scene,
#'   spot_light(
#'     position = projector,
#'     direction = c(0, 1.5, wall_z) - projector,
#'     cone_angle = 36,
#'     falloff_angle = 5,
#'     intensity = 65,
#'     color = "#fff1d4"
#'   )
#' )
#' render_scene(
#'   scene,
#'   lookfrom = c(4.5, 4.1, 8),
#'   lookat = c(-.2, 1.5, -.4),
#'   fov = 0,
#'   ortho_dimensions = c(6.3, 4.7),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
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
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # add_light(): three colored lamps give one rocket three colored shadows.
#' # Each named light adds to the scene; it does not replace the previous one.
#' scene = generate_ground(depth = 0, material = diffuse("white")) |>
#'   add_object(cylinder(
#'     y = .8,
#'     radius = .27,
#'     length = 1.1,
#'     material = diffuse("ivory")
#'   )) |>
#'   add_object(cone(
#'     start = c(0, 1.35, 0),
#'     end = c(0, 1.9, 0),
#'     radius = .27,
#'     material = diffuse("ivory")
#'   ))
#' for (i in 1:3) {
#'   scene = add_light(
#'     scene,
#'     point_light(
#'       position = c(c(-2, 0, 2)[i], 3, 2),
#'       color = c("#ff574a", "#64ff91", "#7298ff")[i],
#'       intensity = 18
#'     ),
#'     name = c("red", "green", "blue")[i]
#'   )
#' }
#' render_scene(
#'   scene,
#'   lookfrom = c(0, 7, 4),
#'   lookat = c(0, 0, -1.7),
#'   fov = 0,
#'   ortho_dimensions = c(8, 6.5),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
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
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # get_light(): hand the followspot from one paper star to the other.
#' # Retrieve the descriptor, edit its direction, and replace the same named light.
#' a = pi / 2 + (0:9) * pi / 5
#' star = cbind(cos(a), sin(a)) * rep(c(.58, .26), 5)
#' scene = generate_ground(depth = 0, material = diffuse("#272d3a")) |>
#'   add_object(xy_rect(
#'     y = 1.2,
#'     z = -.2,
#'     xwidth = 4.5,
#'     ywidth = 2.4,
#'     material = diffuse("#687888")
#'   ))
#' for (i in 1:2) {
#'   scene = add_object(
#'     scene,
#'     extruded_polygon(
#'       star,
#'       x = c(-1.1, 1.1)[i],
#'       y = 1.1,
#'       plane = "xy",
#'       bottom = 0,
#'       top = .06,
#'       material = diffuse(c("#efa24b", "#479ec4")[i])
#'     )
#'   )
#' }
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1, 1),
#'     angular_diameter = 40,
#'     intensity = .2,
#'     name = "fill"
#'   )) |>
#'   add_light(
#'     spot_light(
#'       position = c(0, 2.6, 4),
#'       direction = c(-1.1, -1.5, -4),
#'       cone_angle = 12,
#'       falloff_angle = 3,
#'       intensity = 18
#'     ),
#'     name = "followspot"
#'   )
#' # Before: the amber star has the spotlight.
#'
#' render_scene(
#'   scene,
#'   lookfrom = c(0, 2, 7),
#'   lookat = c(0, 1, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.9, 3.4),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE
#' )
#'
#' beam = get_light(scene, "followspot")
#' beam$direction = c(1.1, 1.1, 0) - beam$position
#' scene = add_light(scene, beam, name = "followspot", replace = TRUE)
#' # After: the blue star has the very same light, aimed in a new direction.
#'
#' render_scene(
#'   scene,
#'   lookfrom = c(0, 2, 7),
#'   lookat = c(0, 1, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.9, 3.4),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE
#' )
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
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # list_lights(): turn the lighting inventory into a string of visible lanterns.
#' # The returned named list contains positions and colors, not just light names.
#' scene = generate_ground(depth = 0, material = diffuse("#374054"))
#' for (i in 1:5) {
#'   scene = add_light(
#'     scene,
#'     point_light(
#'       position = c((i - 3) * .65, 1.2 + .12 * (i - 3)^2, 0),
#'       color = c("#ed8b68", "#f4d16f", "#a9d899", "#7abcc4", "#c5a5df")[i],
#'       intensity = 3
#'     ),
#'     name = paste0("lantern", i)
#'   )
#' }
#' lamps = list_lights(scene)
#' names(lamps)
#' for (lamp in lamps) {
#'   p = lamp$position
#'   # Put the visible bulb above its point source so it does not block the light.
#'   p[2] = p[2] + .18
#'   scene = scene |>
#'     add_object(sphere(
#'       x = p[1],
#'       y = p[2],
#'       z = p[3],
#'       radius = .13,
#'       material = light(
#'         color = lamp$color,
#'         intensity = 1,
#'         importance_sample = FALSE
#'       )
#'     )) |>
#'     add_object(cylinder(
#'       x = p[1],
#'       y = (p[2] + 2.1) / 2,
#'       radius = .008,
#'       length = 2.1 - p[2],
#'       material = diffuse("#172b39")
#'     ))
#' }
#' render_scene(
#'   scene,
#'   lookfrom = c(1, 2, 7),
#'   lookat = c(0, 1, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.2, 3.4),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
list_lights = function(scene) {
  ray_scene_point_lights(scene)
}

#' @rdname get_light
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # remove_light(): the campfire has gone out, but the moon still lights the logs.
#' # Removing a named light keeps the geometry and every other light intact.
#' scene = generate_ground(depth = 0, material = diffuse("#55554c"))
#' for (i in 1:3) {
#'   scene = add_object(
#'     scene,
#'     cylinder(
#'       y = .12 + .12 * (i - 1),
#'       radius = .12,
#'       length = 1.8,
#'       angle = c(90, (i - 1) * 60, 0),
#'       material = diffuse("#885538")
#'     )
#'   )
#' }
#' scene = scene |>
#'   add_light(
#'     point_light(position = c(0, 1, 0), color = "#ff9e43", intensity = 12),
#'     name = "fire"
#'   ) |>
#'   add_light(
#'     point_light(position = c(-3, 4, 2), color = "#8caaff", intensity = 25),
#'     name = "moon"
#'   )
#'
#' # Before: warm firelight and cool moonlight together.
#' render_scene(
#'   scene,
#'   lookfrom = c(3, 3, 5),
#'   lookat = c(0, 0.2, 0),
#'   fov = 0,
#'   ortho_dimensions = c(3.8, 3),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
#' # After: remove only the firelight.
#' scene = remove_light(scene, "fire")
#' names(list_lights(scene)) # Only "moon" remains.
#' render_scene(
#'   scene,
#'   lookfrom = c(3, 3, 5),
#'   lookat = c(0, 0.2, 0),
#'   fov = 0,
#'   ortho_dimensions = c(3.8, 3),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
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
