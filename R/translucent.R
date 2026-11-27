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
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # translucent(): a folded-paper lantern lit from within.
#' # Reflectance controls light returning from the paper; transmittance controls
#' # the warm glow passing through it. Their sum stays below one in each channel.
#' paper = translucent(
#'   reflectance = c(.25, .18, .1),
#'   transmittance = c(.65, .45, .18)
#' )
#' scene = generate_ground(depth = 0, material = diffuse("#253541")) |>
#'   add_object(xy_rect(
#'     y = .8,
#'     z = -.5,
#'     xwidth = 1,
#'     ywidth = 1.5,
#'     material = paper
#'   )) |>
#'   add_object(xy_rect(
#'     y = .8,
#'     z = .5,
#'     xwidth = 1,
#'     ywidth = 1.5,
#'     material = paper
#'   )) |>
#'   add_object(yz_rect(
#'     x = -.5,
#'     y = .8,
#'     ywidth = 1.5,
#'     zwidth = 1,
#'     material = paper
#'   )) |>
#'   add_object(yz_rect(
#'     x = .5,
#'     y = .8,
#'     ywidth = 1.5,
#'     zwidth = 1,
#'     material = paper
#'   )) |>
#'   add_light(point_light(
#'     position = c(0, .8, 0),
#'     color = "#ffdb98",
#'     intensity = 2
#'   ))
#' for (x in c(-.52, .52)) {
#'   for (z in c(-.52, .52)) {
#'     scene = add_object(
#'       scene,
#'       cylinder(
#'         x = x,
#'         y = .8,
#'         z = z,
#'         radius = .025,
#'         length = 1.65,
#'         material = diffuse("#4b2633")
#'       )
#'     )
#'   }
#' }
#' render_scene(
#'   scene,
#'   lookfrom = c(3, 2, 5),
#'   lookat = c(0, .7, 0),
#'   fov = 0,
#'   ortho_dimensions = c(2.7, 2.7),
#'   width = 360,
#'   height = 360,
#'   samples = 16,
#'   denoise = TRUE,
#'   max_depth = 12
#' )
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
