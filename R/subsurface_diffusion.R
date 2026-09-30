#' Normalized Diffusion Subsurface Scattering
#'
#' A fast Christensen--Burley diffusion material for optically thick solids,
#' with a neutral dielectric surface. It samples a surface-to-surface scattering
#' event instead of following an internal random walk. Closed geometry and
#' dielectric priorities work as for [subsurface()].
#'
#' @param color Default `"white"`. Apparent body reflectance, as a color name or
#'   RGB triple in `[0, 1]`, using the same color convention as [diffuse()].
#' @param radius Default `1`. Positive scalar or RGB diffusion-profile scale
#'   in world units. This is not the extinction mean free path in [subsurface()].
#' @param scale Default `1`. Positive multiplier applied to `radius`.
#' @param refraction Default `1.3`. Positive finite interior index of refraction.
#' @param roughness Default `0`. Entry boundary roughness, zero or from `0.0001`
#'   to one. The diffuse exit lobe uses normalized smooth Fresnel transmission.
#' @param priority Default `0`. Nonnegative integer dielectric priority; lower
#'   values win. Give glass a lower value than the liquid and overlap their walls
#'   slightly to avoid an air gap. Use distinct values for non-nested overlaps.
#'
#' @details
#' For reflectance A and scaled radius d, the radial area profile is
#' `A * (exp(-r/d) + exp(-r/(3*d))) / (8*pi*d*r)`. Its integral over a plane is A.
#' `radius` sets the spread, not the object's opacity: this model has no ballistic
#' transmission. Black still retains the neutral surface reflection.
#'
#' This is an empirical, separable diffusion approximation.
#' The profile is normalized on a plane;
#' finite-geometry energy conservation is not enforced. Thin, sharply curved,
#' concave, or disconnected geometry can gain energy or leak light. Prefer
#' [subsurface()] for finite-geometry energy conservation, accurate thin transmission,
#' anisotropy, cameras inside the material, or opaque objects embedded within it.
#' Rays starting inside an active diffusion region are terminated and recorded
#' in `path_warnings`. A console warning is printed once per render; interactive
#' rendering continues and recovers when the camera leaves the region. Surface textures and
#' internal occlusion by ordinary opaque objects are not modeled.
#'
#' Surface samples belong to the same object placement, including instances.
#' Priority-winning dielectrics clip the diffusion region: glass may supply the
#' effective liquid boundary. Hidden interfaces do not add refraction or scattering.
#' Exit directions use the adjacent winning dielectric's index of refraction.
#' The spatial kernel remains an approximation for domains cut by other solids.
#' Straight light connections cannot cross a refractive glass wall; use finite
#' emitters so continued paths can sample illumination through glass.
#'
#' Each entry interface and diffusion event consumes an ordinary depth event.
#' There are no internal random-walk collisions. Atmospheric haze is excluded.
#' Explicit [set_medium()] attachments on the same object are an error.
#'
#' @return A `ray_material` descriptor.
#' @references Christensen and Burley (2015), Approximate Reflectance Profiles
#'   for Efficient Subsurface Scattering, Pixar Technical Memo 15-04.
#'   <https://www.seanet.com/~myandper/abstract/memo1504.htm>.
#'   Pharr, Jakob, and Humphreys, Physically Based Rendering, third edition,
#'   section 15.4, Sampling Subsurface Reflection Functions.
#' @seealso [subsurface()], [dielectric()]
#' @md
#' @export
#' @examples
#' milk = subsurface_diffusion(color = "white", radius = 0.02,
#'                             refraction = 1.333, priority = 1)
#' sphere(material = milk)
subsurface_diffusion = function(
  color = "white",
  radius = 1,
  scale = 1,
  refraction = 1.3,
  roughness = 0,
  priority = 0
) {
  # Reuse the existing authoring validation and automatic interior ownership.
  # Replace its transport descriptor: diffusion has no physical extinction.
  out = subsurface(
    color = color,
    radius = radius,
    scale = scale,
    refraction = refraction,
    roughness = roughness,
    priority = priority
  )
  medium = out[[1]]$subsurface
  medium$sigma_a = rep(0, 3)
  medium$sigma_s = rep(0, 3)
  medium$subsurface = list(
    method = "diffusion",
    color = convert_color(color),
    radius = rep(radius, length.out = 3) * scale,
    refraction = refraction,
    roughness = roughness
  )
  out[[1]]$subsurface = medium
  out
}
