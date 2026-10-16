#' Project an Image into a Scene
#'
#' Create a slide-projector light that paints an image onto objects and through
#' scattering media. Use it for stained-glass lighting, stage patterns, or a
#' miniature cinema. Attach it to a scene with [add_light()].
#'
#' @param image Image filename or numeric RGB array with dimensions height,
#' width, channels. The first row is the top of the projected image. Alpha is ignored.
#' @param position Default `c(0, 1, 0)`. Projector position in world coordinates.
#' @param direction Default `c(0, 0, -1)`. Nonzero direction pointing into the scene.
#' @param up Default `c(0, 1, 0)`. Image-up direction; cannot be parallel to `direction`.
#' @param fov Default `45`. Full field of view in degrees along the shorter image
#' dimension, greater than zero and less than 180. Aspect ratio follows the image.
#' @param intensity Default `1`. Nonnegative radiant-intensity multiplier.
#' Doubling the distance reduces illumination by a factor of four.
#' @param color Default `"white"`. RGB tint, as an R color or three values in [0,1].
#' @param encoding Default `"auto"`. `"linear"` or `"srgb"`, or `"auto"` to
#' treat arrays and HDR/EXR files as linear and other image files as sRGB.
#' @param name Default `"projection"`. Unique light name within the scene.
#' @return A `ray_light` description containing the decoded image.
#' @export
#' @details The image is read when the light is created, so subsequent rendering
#' does not depend on a temporary file. Values are linear RGB after decoding;
#' floating-point images may exceed one. Images use bilinear filtering, a hard
#' rectangular boundary, and no repetition. This is an ideal point projector:
#' shadows are hard and there is no lens defocus. CPU and Metal use the same
#' projection. Looking from the projector toward the scene, image right is
#' `direction` crossed with `up`. Geometry transforms do not move the light.
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # A tiny moon cinema: the crescent bends around a white ball and lands on the wall.
#' axis = seq(-1, 1, length.out = 128)
#' x = matrix(rep(axis, each = 128), 128)
#' y = -t(x)
#' crescent = x^2 + y^2 < .65^2 & (x - .28)^2 + (y - .1)^2 > .57^2
#' slide = array(0, c(128, 128, 3))
#' for (channel in 1:3) {
#'   slide[,, channel] = crescent * c(1, .7, .3)[channel]
#' }
#' scene = xy_rect(y = 1, z = -1, xwidth = 6, ywidth = 4,
#'                 material = diffuse("white")) |>
#'   add_object(sphere(x = -.35, y = 1, radius = .55,
#'                     material = diffuse("white"))) |>
#'   add_light(projection_light(slide, position = c(0, 1, 4), intensity = 40)) |>
#'   add_infinite_light(disk_light(direction = c(-1, 2, 2),
#'                                 angular_diameter = 20, intensity = .4))
#' render_scene(scene, lookfrom = c(3, 2.1, 7), lookat = c(0, 1, 0),
#'              fov = 35, width = 500, height = 400, samples = 16, denoise = TRUE)
projection_light = function(
  image,
  position = c(0, 1, 0),
  direction = c(0, 0, -1),
  up = c(0, 1, 0),
  fov = 45,
  intensity = 1,
  color = "white",
  encoding = "auto",
  name = "projection"
) {
  encoding = match.arg(encoding, c("auto", "linear", "srgb"))
  if (is.character(image)) {
    if (
      length(image) != 1L || is.na(image) || !file.exists(path.expand(image))
    ) {
      stop("`image` must name an existing image file.", call. = FALSE)
    }
    image = path.expand(image)
  } else if (
    !is.numeric(image) ||
      length(dim(image)) != 3L ||
      !dim(image)[3] %in% c(3L, 4L) ||
      any(!is.finite(image)) ||
      any(image < 0)
  ) {
    stop(
      "`image` must be a finite nonnegative RGB or RGBA array.",
      call. = FALSE
    )
  }
  pixels = rayimage::ray_read_image(
    image,
    convert_to_array = TRUE,
    normalize = FALSE,
    source_linear = switch(encoding, auto = NA, linear = TRUE, srgb = FALSE),
    assume_colorspace = rayimage::CS_SRGB
  )
  pixels = array(
    as.numeric(pixels[,, 1:3, drop = FALSE]),
    c(dim(pixels)[1:2], 3L)
  )
  result = structure(
    list(
      type = "projection",
      position = position,
      direction = direction,
      up = up,
      fov = fov,
      image = pixels,
      color = convert_color(color),
      intensity = intensity,
      name = name
    ),
    class = "ray_light"
  )
  validate_point_light(result)
  result$direction = direction / max(abs(direction))
  result$direction = result$direction / sqrt(sum(result$direction^2))
  result$up = up / max(abs(up))
  result$up = result$up / sqrt(sum(result$up^2))
  result
}
