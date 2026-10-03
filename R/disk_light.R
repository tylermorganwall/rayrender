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
#' # disk_light(): a toy skyline in the long shadows of a low sun.
#' # A broad angular diameter softens shadows even though the light is infinitely far away.
#' scene = generate_ground(depth = 0, material = diffuse("#edd9ad"))
#' for (i in 1:5) {
#'   h = c(.6, 1.1, 1.6, .9, 1.3)[i]
#'   scene = add_object(
#'     scene,
#'     cube(
#'       x = (i - 3) * .6,
#'       y = h / 2,
#'       xwidth = .4,
#'       ywidth = h,
#'       zwidth = .45,
#'       material = diffuse(c(
#'         "#b96853",
#'         "#cf9466",
#'         "#708c8e",
#'         "#b96853",
#'         "#cf9466"
#'       )[i])
#'     )
#'   )
#' }
#' scene = add_infinite_light(
#'   scene,
#'   disk_light(
#'     direction = c(-1.82, .55, -1.26),
#'     angular_diameter = 12,
#'     intensity = 110,
#'     color = "#ffce93"
#'   )
#' )
#' render_scene(
#'   scene,
#'   lookfrom = c(2.75, 2.98, 4.81),
#'   lookat = c(.12, .84, .22),
#'   rotate_env = -72,
#'   fov = 0,
#'   ortho_dimensions = c(5.5, 4.4),
#'   width = 420,
#'   height = 340,
#'   samples = 16,
#'   denoise = TRUE
#' )
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
