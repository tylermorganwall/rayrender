#' Image-Based Infinite Light
#' @md
#'
#' @description
#' Creates an infinite light from an equirectangular environment image. Attach
#' lights with [add_infinite_light()]. Their radiance adds together in the
#' background, reflections, and illumination of surfaces and volumes.
#'
#' @param filename Environment image filename. Supports EXR, HDR, PNG, and JPEG,
#' using the same image loading and color conventions as `environment_light`
#' in [render_scene()]. HDR and EXR preserve linear high dynamic range values.
#' @param intensity Default `1`. Nonnegative multiplier for this light's radiance.
#' @param rotation Default `0`. Rotation in degrees around the world Y axis,
#' with the same direction as `rotate_env` in [render_scene()].
#' @param name Default `"environment"`. Unique name within the scene.
#'
#' @details Infinite lights are scene metadata, like cameras, and do not add
#' geometry or change scene bounds. Grouping geometry does not transform them.
#' `add_object()` preserves lights when combining scenes and rejects duplicate
#' names. Scene lights suppress automatic fallback illumination, including when
#' all their intensities are zero. An explicit `environment_light` render argument
#' adds another image light; `intensity_env` controls that legacy light only.
#' `rotate_env` and interactive environment rotation rotate all infinite lights
#' together, in addition to each light's own rotation.
#'
#' Additive lighting lets you keep an existing HDRI and attach independently
#' adjustable fill lights. At each direction, the background radiance is the
#' sum of the lights after applying their intensities and rotations. For example,
#' a cool fill can reveal detail in blue materials under a warm environment.
#' Keep exposure fixed when comparing the result to see the added illumination.
#'
#' @return A `ray_infinite_light` object.
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # Write explicitly linear RGB radiance to EXR, including values greater than 1.
#' write_environment = function(image) {
#'   file = tempfile(fileext = ".exr")
#'   image = rayimage::ray_read_image(
#'     image,
#'     normalize = FALSE,
#'     source_linear = TRUE,
#'     assume_colorspace = rayimage::CS_SRGB
#'   )
#'   rayimage::ray_write_image(image, file, clamp = FALSE)
#'   file
#' }
#'
#' # 1. Add an independently adjustable fill to an existing environment.
#' # These constant-color maps make the color and brightness changes easy to see.
#' warm = cool = array(0, c(32, 64, 3))
#' warm_rgb = c(0.8, 0.3, 0.1)
#' cool_rgb = c(0.1, 0.3, 0.8)
#' for (channel in 1:3) {
#'   warm[,, channel] = warm_rgb[channel]
#'   cool[,, channel] = cool_rgb[channel]
#' }
#' warm_file = write_environment(warm)
#' cool_file = write_environment(cool)
#'
#' # The included R logo retains its blue lettering and gray surround.
#' # The ground and oblique camera reveal contact shadows and the beveled edges.
#' base_scene = obj_model(r_obj(), scale = 2.5) |>
#'   add_object(generate_ground(
#'     depth = -0.96,
#'     material = diffuse(color = "grey20")
#'   )) |>
#'   add_camera(camera(
#'     lookfrom = c(2.5, 1.2, -6),
#'     lookat = c(0, -0.05, 0),
#'     fov = 28,
#'     aperture = 0
#'   ))
#' warm_scene = base_scene |>
#'   add_infinite_light(infinite_light(warm_file, name = "warm"))
#' set.seed(724)
#' render_scene(
#'   warm_scene,
#'   width = 600,
#'   height = 500,
#'   samples = 128,
#'   integrator_type = "nee",
#'   bloom = FALSE
#' )
#'
#' # The two lights contribute warm + 0.5 * cool. The extra cool radiance
#' # lifts the blue lettering and shifts the background color.
#' # Keep exposure fixed between renders to see both the color and brightness change.
#' # The same approach adds fill to a photographic HDRI without editing that image.
#' additive_scene = warm_scene |>
#'   add_infinite_light(infinite_light(cool_file, intensity = 0.5, name = "cool"))
#' set.seed(724)
#' render_scene(
#'   additive_scene,
#'   width = 600,
#'   height = 500,
#'   samples = 128,
#'   integrator_type = "nee",
#'   bloom = FALSE
#' )
#'
#' # 2. Warm and cool fills on opposite sides, like two broad studio softboxes.
#' # Each EXR is dark except for a bright patch above the horizon. The patch is
#' # one-quarter across the map. Rotations of 30 and 210 degrees put the patches
#' # on opposite sides. A dim neutral environment keeps the front readable.
#' u = (seq_len(256) - 0.5) / 256
#' v = (seq_len(128) - 0.5) / 128
#' # Longitude wraps around: keep the softbox continuous at the EXR's edges.
#' du = ((u - 0.25 + 0.5) %% 1) - 0.5
#' softbox = outer(
#'   exp(-((v - 0.32) / 0.18)^2),
#'   exp(-(du / 0.12)^2)
#' )
#' warm_side = cool_side = array(0, c(128, 256, 3))
#' for (channel in 1:3) {
#'   warm_side[,, channel] = 8 * warm_rgb[channel] * softbox
#'   cool_side[,, channel] = 8 * cool_rgb[channel] * softbox
#' }
#' neutral_file = write_environment(array(0.08, c(32, 64, 3)))
#' warm_side_file = write_environment(warm_side)
#' cool_side_file = write_environment(cool_side)
#' fill_scene = base_scene |>
#'   add_infinite_light(infinite_light(neutral_file, name = "neutral")) |>
#'   add_infinite_light(infinite_light(
#'     warm_side_file,
#'     rotation = 30,
#'     name = "warm_fill"
#'   )) |>
#'   add_infinite_light(infinite_light(
#'     cool_side_file,
#'     rotation = 210,
#'     intensity = 0.75,
#'     name = "cool_fill"
#'   ))
#' set.seed(724)
#' render_scene(
#'   fill_scene,
#'   width = 600,
#'   height = 500,
#'   samples = 128,
#'   integrator_type = "nee",
#'   bloom = FALSE
#' )
#'
#' # Change either fill's intensity or rotation independently to shape the lighting.
#' # This gives independent control over both sides while retaining the base light.
infinite_light = function(
  filename,
  intensity = 1,
  rotation = 0,
  name = "environment"
) {
  result = structure(
    list(
      type = "image",
      filename = filename,
      intensity = intensity,
      rotation = rotation,
      name = name
    ),
    class = "ray_infinite_light"
  )
  validate_infinite_light(result)
  result$filename = path.expand(filename)
  result
}

#' Add an Infinite Light
#' @md
#'
#' @param scene Scene to modify.
#' @param light Light created with [infinite_light()], [sky_light()],
#' [sun_light()], or [moon_light()].
#' @param name Default `NULL`. Optional name overriding the light's name.
#' @param replace Default `FALSE`. Replace an existing light with the same name.
#' @return A modified scene.
#' @export
add_infinite_light = function(scene, light, name = NULL, replace = FALSE) {
  if (!inherits(light, "ray_infinite_light")) {
    stop("light must inherit from class 'ray_infinite_light'.", call. = FALSE)
  }
  if (!is.null(name)) {
    light$name = name
  }
  validate_infinite_light(light)
  if (!is.logical(replace) || length(replace) != 1 || is.na(replace)) {
    stop("replace must be TRUE or FALSE.", call. = FALSE)
  }
  lights = ray_scene_infinite_lights(scene)
  if (!replace && light$name %in% names(lights)) {
    stop(
      "Infinite light name '",
      light$name,
      "' already exists. Use replace = TRUE to replace it.",
      call. = FALSE
    )
  }
  lights[[light$name]] = light
  attr(scene, "ray_infinite_lights") = lights
  scene
}

#' Get an Infinite Light
#'
#' @param scene Scene containing infinite lights.
#' @param name Name of the infinite light.
#' @return A `ray_infinite_light` object.
#' @export
get_infinite_light = function(scene, name) {
  lights = ray_scene_infinite_lights(scene)
  if (
    !is.character(name) ||
      length(name) != 1 ||
      is.na(name) ||
      !name %in% names(lights)
  ) {
    stop("name must identify an infinite light in the scene.", call. = FALSE)
  }
  lights[[name]]
}

#' List Infinite Lights
#'
#' @param scene Scene containing infinite lights.
#' @return A named list of `ray_infinite_light` objects.
#' @export
list_infinite_lights = function(scene) {
  ray_scene_infinite_lights(scene)
}

#' Remove an Infinite Light
#'
#' @param scene Scene to modify.
#' @param name Name of the infinite light to remove.
#' @return A modified scene.
#' @export
remove_infinite_light = function(scene, name) {
  get_infinite_light(scene, name)
  lights = ray_scene_infinite_lights(scene)
  lights[[name]] = NULL
  attr(scene, "ray_infinite_lights") = if (length(lights)) lights else NULL
  scene
}

#' @export
print.ray_infinite_light = function(x, ...) {
  if (x$type %in% c("sky", "sun", "moon")) {
    cat(sprintf(
      "Infinite light '%s' (%s)\n  location: %g, %g\n  datetime: %s\n  intensity: %g\n  rotation: %g degrees\n",
      x$name,
      x$type,
      x$lat,
      x$long,
      format(x$datetime, usetz = TRUE),
      x$intensity,
      x$rotation
    ))
    return(invisible(x))
  }
  cat(sprintf(
    "Infinite light '%s' (image)\n  file: %s\n  intensity: %g\n  rotation: %g degrees\n",
    x$name,
    x$filename,
    x$intensity,
    x$rotation
  ))
  invisible(x)
}

#' @keywords internal
validate_infinite_light = function(light) {
  if (
    !inherits(light, "ray_infinite_light") ||
      !is.character(light$type) ||
      length(light$type) != 1 ||
      !light$type %in% c("image", "sky", "sun", "moon", "disk")
  ) {
    stop(
      "Expected an image, sky, sun, or moon ray_infinite_light.",
      call. = FALSE
    )
  }
  if (
    !is.character(light$name) ||
      length(light$name) != 1 ||
      is.na(light$name) ||
      !nzchar(trimws(light$name))
  ) {
    stop("Infinite light name must be a nonempty string.", call. = FALSE)
  }
  if (light$type %in% c("sky", "sun", "moon")) {
    validate_sky_light(light)
    if (light$type != "sky") validate_celestial_light(light)
  } else {
    if (
      !is.character(light$filename) ||
        length(light$filename) != 1 ||
        is.na(light$filename) ||
        !nzchar(light$filename)
    ) {
      stop("Infinite light filename must be a nonempty string.", call. = FALSE)
    }
    filename = path.expand(light$filename)
    if (!file.exists(filename) || dir.exists(filename)) {
      stop(
        "Infinite light file does not exist or is a directory: ",
        filename,
        call. = FALSE
      )
    }
  }
  if (light$type == "disk") {
    validate_celestial_disk(light)
  }
  for (field in c("intensity", "rotation")) {
    value = light[[field]]
    if (
      !is.numeric(value) ||
        length(value) != 1 ||
        !is.finite(value) ||
        (field == "intensity" && value < 0)
    ) {
      stop(
        "Infinite light ",
        field,
        " must be a finite numeric scalar",
        if (field == "intensity") " greater than or equal to zero" else "",
        ".",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

#' @keywords internal
ray_scene_infinite_lights = function(scene) {
  lights = attr(scene, "ray_infinite_lights", exact = TRUE)
  if (is.null(lights)) {
    return(list())
  }
  if (
    !is.list(lights) ||
      (length(lights) &&
        (is.null(names(lights)) ||
          anyNA(names(lights)) ||
          any(!nzchar(names(lights))) ||
          anyDuplicated(names(lights))))
  ) {
    stop("ray_infinite_lights must be a uniquely named list.", call. = FALSE)
  }
  for (i in seq_along(lights)) {
    validate_infinite_light(lights[[i]])
    if (!identical(names(lights)[i], lights[[i]]$name)) {
      stop(
        "Infinite light names must match their scene entries.",
        call. = FALSE
      )
    }
  }
  lights
}
