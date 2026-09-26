#' Uniform Infinite Disk Light
#' @md
#'
#' @description
#' Create a uniformly colored disk at infinite distance. Add it to a scene with
#' [add_infinite_light()] to illuminate surfaces and volumes and show the disk
#' in the background and reflections. No image files or sky datasets are needed.
#'
#' @param color Default `"white"`. Light color, as an R color name, hexadecimal
#' string, or numeric RGB vector of three values between 0 and 1, as in [light()].
#' @param intensity Default `1`. Nonnegative radiance multiplier for the color.
#' @param angular_diameter Default `0.53`. Full angular diameter in degrees,
#' greater than 0 and less than 180. Larger disks produce softer shadows and
#' more total illumination at the same intensity.
#' @param direction Default `c(0, 1, 0)`. Nonzero vector pointing from the scene
#' toward the disk center. Normalized automatically. World +Y is up, +Z is
#' north, and -X is east, matching [sun_light()] and [moon_light()].
#' @param rotation Default `0`. Rotation in degrees around world Y, with the
#' same convention as [infinite_light()].
#' @param name Default `"disk"`. Unique light name within the scene.
#'
#' @details The light has no distance falloff or parallax and adds no geometry.
#' Its color and radiance are constant across its angular extent. Intensity is
#' not normalized by disk area. The disk can be placed below the horizon;
#' it has no automatic horizon clipping or celestial atmospheric filtering.
#' It adds to other infinite lights without replacing automatic Sun or Moon disks.
#' Global `rotate_env` rotation applies in addition to its own rotation.
#'
#' @return A `ray_infinite_light` description.
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' scene = generate_ground(material = diffuse("grey60")) |>
#'   add_object(sphere(y = 1, material = diffuse("coral"))) |>
#'   add_infinite_light(disk_light(
#'     color = "#fff0dd", intensity = 100, angular_diameter = 10,
#'     direction = c(-1, 2, 1)
#'   ))
#' render_scene(scene, lookfrom = c(4, 3, -8), lookat = c(0, 1, 0),
#'              samples = 64, aperture = 0, ambient_light = FALSE)
disk_light = function(
  color = "white",
  intensity = 1,
  angular_diameter = 0.53,
  direction = c(0, 1, 0),
  rotation = 0,
  name = "disk"
) {
  result = structure(
    list(
      type = "uniform_disk",
      color = convert_color(color),
      intensity = intensity,
      angular_diameter = angular_diameter,
      direction = direction,
      rotation = rotation,
      name = name,
      clip_horizon = FALSE
    ),
    class = "ray_infinite_light"
  )
  validate_infinite_light(result)
  # Scale first so normalizing very large or very small vectors stays finite.
  result$direction = direction / max(abs(direction))
  result$direction = result$direction / sqrt(sum(result$direction^2))
  result
}
