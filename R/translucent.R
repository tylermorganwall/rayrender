#' Translucent Sheet Material
#'
#' Lambertian reflection and diffuse transmission for thin surfaces such as
#' paper, lampshades, or leaves. Light scatters at the surface without refraction or a
#' volumetric walk. Uses geometric normals; bump mapping is not supported.
#' Reflected plus transmitted energy is limited to one per color channel;
#' overbright texture values are proportionally reduced at each lookup.
#' @md
#'
#' @param reflectance Default `c(0.25, 0.25, 0.25)`. Reflected RGB fraction or R color.
#' @param transmittance Default `c(0.25, 0.25, 0.25)`. Transmitted RGB fraction or R color.
#' @param image_texture Default `""`. RGB array or image filename replacing reflectance.
#' @param transmission_texture Default `""`. RGB array or image filename replacing transmittance.
#' @param image_repeat Default `1`. Scalar or length-two UV repeat for both images.
#' @param image_offset Default `c(0, 0)`. Length-two UV offset, applied after repeat.
#' @param alpha_texture Default `""`. Surface coverage image, separate from diffuse transmission.
#' @return A rayrender material.
#' @export
#' @examples
#' translucent(reflectance = c(.25, .25, .25), transmittance = c(.6, .5, .3))
translucent = function(
  reflectance = c(.25, .25, .25),
  transmittance = c(.25, .25, .25),
  image_texture = "",
  transmission_texture = "",
  image_repeat = 1,
  image_offset = c(0, 0),
  alpha_texture = ""
) {
  reflectance = convert_color(reflectance)
  transmittance = convert_color(transmittance)
  if (
    !length(reflectance) %in% c(1L, 3L) ||
      !length(transmittance) %in% c(1L, 3L) ||
      any(!is.finite(c(reflectance, transmittance)))
  ) {
    stop(
      "Reflectance and transmittance must be finite scalar or RGB colors.",
      call. = FALSE
    )
  }
  reflectance = rep(reflectance, length.out = 3L)
  transmittance = rep(transmittance, length.out = 3L)
  if (any(reflectance + transmittance > 1 + 1e-8)) {
    stop(
      "`reflectance + transmittance` must not exceed one per channel.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(image_repeat) ||
      !length(image_repeat) %in% c(1L, 2L) ||
      any(!is.finite(image_repeat))
  ) {
    stop(
      "`image_repeat` must contain one or two finite numbers.",
      call. = FALSE
    )
  }
  out = diffuse(
    color = reflectance,
    image_texture = image_texture,
    image_repeat = image_repeat,
    image_offset = image_offset,
    alpha_texture = alpha_texture
  )
  out[[1]]$type = get_material_enum("translucent")
  out[[1]]$transmittance = transmittance
  out[[1]]$transmission_texture = check_image_texture(transmission_texture)
  out[[1]]$transmission_repeat = rep(image_repeat, length.out = 2L)
  out[[1]]$transmission_offset = check_image_offset(image_offset)
  out
}
