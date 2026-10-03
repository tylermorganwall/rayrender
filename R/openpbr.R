#' OpenPBR Surface Material
#'
#' A separate layered uber shader implementing OpenPBR Surface 1.1.1 using
#' Adobe's reference BSDF. Includes metal and rough diffuse bases, refraction,
#' subsurface scattering, coat, fuzz, thin film, and emission.
#'
#' @param base_weight Default `1.0`. Multiplier on the intensity of the reflection from the diffuse and metallic base.
#' @param base_color Default `c(0.8, 0.8, 0.8)`. Color of the reflection from the diffuse and metallic base.
#' @param base_diffuse_roughness Default `0.0`. Roughness of the diffuse reflection. Higher values cause the surface to appear flatter.
#' @param base_metalness Default `0.0`. Specifies how metallic the base material appears (dials the base from pure dielectric to pure metal).
#' @param specular_weight Default `1.0`. Multiplies the specular reflectivity.
#' @param specular_color Default `c(1, 1, 1)`. Color of the specular reflection (controls the physical edge-tint for metals, and a non-physical overall tint for dielectrics).
#' @param specular_roughness Default `0.3`. The roughness of the specular reflection. Lower numbers produce sharper reflections, higher numbers produce blurrier reflections.
#' @param specular_ior Default `1.5`. Index of refraction of the dielectric base.
#' @param specular_roughness_anisotropy Default `0.0`. The directional bias of the roughness of the metal/dielectric base, resulting in increasingly stretched highlights along the tangent direction.
#' @param transmission_weight Default `0.0`. Mixture weight between the transparent and opaque dielectric base. The greater the value the more transparent the material.
#' @param transmission_color Default `c(1, 1, 1)`. Controls color of the transparent base due to Beer's law volumetric absorption under the surface (reverts to a non-physical tint when transmission_depth is zero).
#' @param transmission_depth Default `0.0`. Specifies the distance light travels inside the transparent base before it becomes exactly the transmission_color according to Beer's law.
#' @param transmission_scatter Default `c(0, 0, 0)`. Controls the color of light volumetrically scattered inside the transparent base. Suitable for materials with visually significant scattering such as honey, fruit juice, murky water, opalescent glass, or milky glass.
#' @param transmission_scatter_anisotropy Default `0.0`. The amount of directional bias, or anisotropy, of the volumetric scattering in the transparent base.
#' @param transmission_dispersion_scale Default `0.0`. Linearly scales the amount of dispersion.
#' @param transmission_dispersion_abbe_number Default `20.0`. Physical Abbe number of the dielectric medium, describing how much the dielectric index of refraction varies across wavelengths.
#' @param subsurface_weight Default `0`. Mixture weight which dials the opaque dielectric base between diffuse reflection and subsurface scattering. A value of 1.0 indicates full subsurface scattering and a value 0 for diffuse reflection only.
#' @param subsurface_color Default `c(0.8, 0.8, 0.8)`. The observed reflection color of the subsurface scattering medium.
#' @param subsurface_radius Default `1.0`. Length scale of the subsurface scattering mean free path.
#' @param subsurface_radius_scale Default `c(1.0, 0.5, 0.25)`. RGB multiplier to subsurface_radius, giving the per-channel scattering mean-free-paths.
#' @param subsurface_scatter_anisotropy Default `0.0`. Controls the phase-function of subsurface scattering, where zero scatters light evenly, positive values scatter forwards, and negative values scatter backwards.
#' @param fuzz_weight Default `0.0`. The presence weight of a fuzz layer that can be used to approximate microfibers, for fabrics such as velvet and satin as well as dust grains.
#' @param fuzz_color Default `c(1, 1, 1)`. The color of the fuzz layer.
#' @param fuzz_roughness Default `0.5`. The roughness of the fuzz layer.
#' @param coat_weight Default `0.0`. The presence weight of a reflective clear-coat layer on top of the material. Use for materials such as car paint or an oily layer.
#' @param coat_color Default `c(1, 1, 1)`. The color of the clear-coat layer's transparency, due to absorption in the coat.
#' @param coat_roughness Default `0.0`. The roughness of the clear-coat reflections. The lower the value, the sharper the reflection.
#' @param coat_roughness_anisotropy Default `0.0`. The directional bias of the roughness of the clear-coat layer, resulting in increasingly stretched highlights along the coat tangent direction.
#' @param coat_ior Default `1.6`. The index of refraction of the clear-coat layer.
#' @param coat_darkening Default `1.0`. Modulates the physical coat darkening effect.
#' @param thin_film_weight Default `0`. Coverage weight of the thin-film. Use for materials such as multi-tone car paint or soap bubbles.
#' @param thin_film_thickness Default `0.5`. The thickness of the thin-film layer on the base (in micrometers).
#' @param thin_film_ior Default `1.4`. The index of refraction of the thin-film.
#' @param emission_luminance Default `0.0`. The amount of emitted light, as a luminance in nits.
#' @param emission_color Default `c(1, 1, 1)`. The color of the emitted light.
#' @param geometry_opacity Default `1`. The opacity of the entire material.
#' @param geometry_thin_walled Default `FALSE`. If true the surface is double-sided and represents an infinitesimally thin shell. Suitable for extremely geometrically thin objects such as leaves or paper.
#' @param geometry_normal Default `NULL`. Input geometric normal
#' @param geometry_coat_normal Default `NULL`. Input normal for coat layer
#' @param geometry_tangent Default `NULL`. Input geometric tangent
#' @param geometry_coat_tangent Default `NULL`. Input geometric tangent for coat layer
#' @param priority Default `0`. Nonnegative dielectric priority; lower values win.
#' @param image_texture Default `""`. Base-color texture filename or RGB array; replaces base_color.
#' @param image_repeat Default `1`. One or two positive UV repeat factors.
#' @param image_offset Default `c(0, 0)`. Finite length-two UV translation applied after repeat, before wrapping, to base-color, bump, and roughness textures.
#' @param bump_texture Default `""`. Height-map filename, matrix, or array.
#' @param bump_intensity Default `1`. Height-map scale. Slopes are measured per UV unit, including texture repeats, independently of image resolution. High values may lead to unphysical results.
#' @param roughness_texture Default `""`. Scalar image replacing specular_roughness, without gamma decoding.
#' @param importance_sample Default `TRUE`. Include emissive surfaces in direct-light sampling.
#'
#' @details
#' Colors use rayrender's color convention (numeric RGB triples or R color names).
#' Numeric colors are passed directly to the linear BSDF, as for diffuse().
#' `base_color` and `specular_roughness` also accept composable [textures].
#' Each graph has its own coordinate mapping; legacy image arguments cannot
#' be combined with a graph for the same input.
#' subsurface_radius_scale is a numeric RGB distance multiplier, not a color.
#' Geometry vectors default to the surface normal and UV tangent. Explicit vectors
#' are in world coordinates. The coat uses the base frame unless overridden.
#'
#' Solid transmission and subsurface materials automatically own a homogeneous
#' interior and require closed spheres, cubes, ellipsoids, or triangle meshes.
#' They cannot be combined with set_medium() or opacity cutouts. Thin-walled
#' materials instead scatter at the surface and can be used on open geometry.
#' Distances are in world units, except thin_film_thickness, which is in micrometers.
#' The interior uses physical random walks with dielectric-priority handling.
#' OpenPBR scenes select the NEE integrator, including inside instances.
#'
#' The reference implementation clamps microfacet roughness to at least 0.001
#' and regularizes zero volume distances and perfectly directional scattering.
#' Dielectric priority contacts use the nominal specular_ior; keep specular_weight
#' at one for physical IOR matching across different material types.
#' It uses representative RGB wavelengths for dispersion and thin film; this
#' is not a spectral renderer. Thin-film iridescence does not affect thin-walled
#' transmission in the reference implementation. Coat/fuzz use its default
#' energy-conserving camera-to-light layering, which is not fully reciprocal.
#' Shading-normal hemisphere clipping may lose energy at extreme bump slopes.
#'
#' @return A ray_material descriptor.
#' @references OpenPBR Surface 1.1.1, <https://academysoftwarefoundation.github.io/OpenPBR/>.
#' Adobe OpenPBR BSDF, <https://github.com/adobe/openpbr-bsdf>.
#' @seealso [subsurface()], [dielectric()], [microfacet()]
#' @export
#' @md
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # openpbr(): six eggs from a very unusual Easter basket.
#' # The same geometry shows fuzz, metal, coat, thin film, transmission, and SSS.
#' materials = list(
#'   Velvet = openpbr(
#'     base_color = "#531c46",
#'     specular_roughness = .6,
#'     fuzz_weight = 1,
#'     fuzz_color = "#f5a5de",
#'     fuzz_roughness = .65
#'   ),
#'   Copper = openpbr(
#'     base_color = "#df783c",
#'     base_metalness = 1,
#'     specular_roughness = .23,
#'     specular_roughness_anisotropy = .8
#'   ),
#'   Lacquer = openpbr(
#'     base_color = "#961f29",
#'     coat_weight = 1,
#'     coat_roughness = .08,
#'     specular_roughness = .45
#'   ),
#'   Iridescent = openpbr(
#'     base_color = "#bccac9",
#'     base_metalness = 1,
#'     specular_roughness = .17,
#'     thin_film_weight = 1,
#'     thin_film_thickness = .65
#'   ),
#'   Glass = openpbr(
#'     transmission_weight = 1,
#'     transmission_color = "#79cedd",
#'     transmission_depth = .5,
#'     specular_roughness = .035
#'   ),
#'   Wax = openpbr(
#'     subsurface_weight = 1,
#'     subsurface_color = "#ffe2b0",
#'     subsurface_radius = .06,
#'     specular_roughness = .3
#'   )
#' )
#' scene = generate_ground(depth = 0, material = diffuse("#242732"))
#' x = rep(c(-1.3, 0, 1.3), 2)
#' z = rep(c(-.85, .85), each = 3)
#' for (i in seq_along(materials)) {
#'   scene = scene |>
#'     add_object(cylinder(
#'       x = x[i],
#'       y = .06,
#'       z = z[i],
#'       radius = .52,
#'       length = .12,
#'       material = diffuse("#47505c")
#'     )) |>
#'     add_object(ellipsoid(
#'       x = x[i],
#'       y = .74,
#'       z = z[i],
#'       a = .44,
#'       b = .62,
#'       c = .44,
#'       angle = c(0, 0, -10),
#'       material = materials[[i]]
#'     ))
#' }
#' scene = scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1, 1),
#'     angular_diameter = 38,
#'     intensity = 12,
#'     name = "softbox"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .6, -1),
#'     angular_diameter = 25,
#'     intensity = 8,
#'     color = "#b7d9ff",
#'     name = "edge"
#'   ))
#' labels = screen_text(
#'   names(materials),
#'   x = x,
#'   y = rep(c(1.65, .06), each = 3),
#'   z = z + rep(c(0, .6), each = 3),
#'   size = 15,
#'   hjust = .5,
#'   color = "white"
#' )
#'
#' render_scene(
#'   scene,
#'   lookfrom = c(0, 4.5, 8),
#'   lookat = c(0, 0.6, 0),
#'   fov = 0,
#'   ortho_dimensions = c(4.9, 3.9),
#'   width = 500,
#'   height = 400,
#'   samples = 16,
#'   denoise = TRUE,
#'   screen_text = labels
#' )
openpbr = function(
  base_weight = 1.0,
  base_color = c(0.8, 0.8, 0.8),
  base_diffuse_roughness = 0.0,
  base_metalness = 0.0,
  specular_weight = 1.0,
  specular_color = c(1, 1, 1),
  specular_roughness = 0.3,
  specular_ior = 1.5,
  specular_roughness_anisotropy = 0.0,
  transmission_weight = 0.0,
  transmission_color = c(1, 1, 1),
  transmission_depth = 0.0,
  transmission_scatter = c(0, 0, 0),
  transmission_scatter_anisotropy = 0.0,
  transmission_dispersion_scale = 0.0,
  transmission_dispersion_abbe_number = 20.0,
  subsurface_weight = 0,
  subsurface_color = c(0.8, 0.8, 0.8),
  subsurface_radius = 1.0,
  subsurface_radius_scale = c(1.0, 0.5, 0.25),
  subsurface_scatter_anisotropy = 0.0,
  fuzz_weight = 0.0,
  fuzz_color = c(1, 1, 1),
  fuzz_roughness = 0.5,
  coat_weight = 0.0,
  coat_color = c(1, 1, 1),
  coat_roughness = 0.0,
  coat_roughness_anisotropy = 0.0,
  coat_ior = 1.6,
  coat_darkening = 1.0,
  thin_film_weight = 0,
  thin_film_thickness = 0.5,
  thin_film_ior = 1.4,
  emission_luminance = 0.0,
  emission_color = c(1, 1, 1),
  geometry_opacity = 1,
  geometry_thin_walled = FALSE,
  geometry_normal = NULL,
  geometry_coat_normal = NULL,
  geometry_tangent = NULL,
  geometry_coat_tangent = NULL,
  priority = 0,
  image_texture = "",
  image_repeat = 1,
  image_offset = c(0, 0),
  bump_texture = "",
  bump_intensity = 1,
  roughness_texture = "",
  importance_sample = TRUE
) {
  if (inherits(base_color, "ray_texture") && base_color$op == "constant") {
    base_color = rep(base_color$value, length.out = 3L)
  }
  if (
    inherits(specular_roughness, "ray_texture") &&
      specular_roughness$op == "constant"
  ) {
    specular_roughness = texture_scalar(specular_roughness)$value
  }
  color_graph = if (inherits(base_color, "ray_texture")) base_color else NULL
  roughness_graph = if (inherits(specular_roughness, "ray_texture")) {
    texture_scalar(specular_roughness)
  } else {
    NULL
  }
  if (!is.null(color_graph)) {
    if (!identical(image_texture, "")) {
      stop(
        "Choose a base-color graph or image_texture, not both.",
        call. = FALSE
      )
    }
    base_color = c(1, 1, 1)
  }
  if (!is.null(roughness_graph)) {
    if (!identical(roughness_texture, "")) {
      stop(
        "Choose a roughness graph or roughness_texture, not both.",
        call. = FALSE
      )
    }
    specular_roughness = 0.3
  }
  parameters = as.list(environment())[c(
    "base_weight",
    "base_color",
    "base_diffuse_roughness",
    "base_metalness",
    "specular_weight",
    "specular_color",
    "specular_roughness",
    "specular_ior",
    "specular_roughness_anisotropy",
    "transmission_weight",
    "transmission_color",
    "transmission_depth",
    "transmission_scatter",
    "transmission_scatter_anisotropy",
    "transmission_dispersion_scale",
    "transmission_dispersion_abbe_number",
    "subsurface_weight",
    "subsurface_color",
    "subsurface_radius",
    "subsurface_radius_scale",
    "subsurface_scatter_anisotropy",
    "fuzz_weight",
    "fuzz_color",
    "fuzz_roughness",
    "coat_weight",
    "coat_color",
    "coat_roughness",
    "coat_roughness_anisotropy",
    "coat_ior",
    "coat_darkening",
    "thin_film_weight",
    "thin_film_thickness",
    "thin_film_ior",
    "emission_luminance",
    "emission_color",
    "geometry_opacity",
    "geometry_thin_walled",
    "geometry_normal",
    "geometry_coat_normal",
    "geometry_tangent",
    "geometry_coat_tangent",
    "priority"
  )]
  # Validate the complete descriptor before allocating textures or attaching an interior.
  for (name in names(parameters)) {
    value = parameters[[name]]
    if (inherits(value, "ray_texture")) {
      stop(
        sprintf("Texture graphs are not supported for `%s` yet.", name),
        call. = FALSE
      )
    }
    if (name == "geometry_thin_walled") {
      if (
        !is.logical(value) ||
          length(value) != 1L ||
          is.na(value) ||
          !is.null(dim(value))
      ) {
        stop("`geometry_thin_walled` must be TRUE or FALSE.", call. = FALSE)
      }
      next
    }
    if (
      name %in%
        c(
          "geometry_normal",
          "geometry_tangent",
          "geometry_coat_normal",
          "geometry_coat_tangent"
        )
    ) {
      if (is.null(value)) {
        next
      }
      if (
        !is.numeric(value) ||
          length(value) != 3L ||
          !is.null(dim(value)) ||
          any(!is.finite(value)) ||
          !any(value != 0)
      ) {
        stop(
          sprintf("`%s` must be NULL or a finite nonzero three-vector.", name),
          call. = FALSE
        )
      }
      parameters[[name]] = value /
        sqrt(sum((value / max(abs(value)))^2)) /
        max(abs(value))
      next
    }
    if (
      grepl("_color$", name) ||
        name %in% c("transmission_scatter", "subsurface_radius_scale")
    ) {
      if (
        name != "subsurface_radius_scale" &&
          is.character(value) &&
          length(value) == 1L &&
          !is.na(value)
      ) {
        value = convert_color(value)
      }
      if (
        !is.numeric(value) ||
          length(value) != 3L ||
          !is.null(dim(value)) ||
          any(!is.finite(value)) ||
          any(value < 0 | value > 1)
      ) {
        stop(
          sprintf("`%s` must be an RGB triple in [0, 1].", name),
          call. = FALSE
        )
      }
      parameters[[name]] = value
      next
    }
    lower = if (grepl("scatter_anisotropy$", name)) -1 else 0
    upper = if (
      name %in%
        c(
          "specular_weight",
          "specular_ior",
          "transmission_depth",
          "transmission_dispersion_abbe_number",
          "subsurface_radius",
          "coat_ior",
          "thin_film_thickness",
          "thin_film_ior",
          "emission_luminance"
        )
    ) {
      1e30
    } else {
      1
    }
    if (name == "priority") {
      upper = .Machine$integer.max
    }
    positive = name %in%
      c(
        "specular_ior",
        "coat_ior",
        "thin_film_ior",
        "transmission_dispersion_abbe_number"
      )
    if (
      !is.numeric(value) ||
        length(value) != 1L ||
        !is.null(dim(value)) ||
        !is.finite(value) ||
        value < lower ||
        value > upper ||
        (positive && value <= 0) ||
        (name == "priority" && value != floor(value))
    ) {
      stop(
        sprintf(
          "Invalid `%s`: expected a finite scalar in [%s, %s]%s.",
          name,
          lower,
          upper,
          if (positive) " greater than zero" else ""
        ),
        call. = FALSE
      )
    }
  }
  owns_interior = !geometry_thin_walled &&
    base_metalness < 1 &&
    (transmission_weight > 0 || subsurface_weight > 0)
  if (owns_interior && geometry_opacity != 1) {
    stop(
      "Solid OpenPBR interiors require `geometry_opacity = 1`.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(image_repeat) ||
      !length(image_repeat) %in% 1:2 ||
      any(!is.finite(image_repeat)) ||
      any(image_repeat <= 0)
  ) {
    stop(
      "`image_repeat` must contain one or two positive finite values.",
      call. = FALSE
    )
  }
  if (
    !is.numeric(bump_intensity) ||
      length(bump_intensity) != 1L ||
      !is.finite(bump_intensity)
  ) {
    stop("`bump_intensity` must be a finite scalar.", call. = FALSE)
  }
  if (
    !is.logical(importance_sample) ||
      length(importance_sample) != 1L ||
      is.na(importance_sample)
  ) {
    stop("`importance_sample` must be TRUE or FALSE.", call. = FALSE)
  }
  out = diffuse(
    color = parameters$base_color,
    image_texture = image_texture,
    image_repeat = image_repeat,
    image_offset = image_offset,
    bump_texture = bump_texture,
    bump_intensity = bump_intensity,
    importance_sample = importance_sample && emission_luminance > 0
  )
  out[[1]]$type = get_material_enum("openpbr")
  out[[1]]$texture_graphs = list(
    color = color_graph,
    roughness = roughness_graph
  )
  out[[1]]$openpbr = parameters
  out[[1]]$roughness_texture = check_image_texture(roughness_texture)
  if (owns_interior) {
    interior = homogeneous_medium(sigma_a = 0, sigma_s = 0, haze = FALSE)
    interior$openpbr = parameters
    # The existing automatic-interior preparation validates closed shapes,
    # imported face overrides, instances, and conflicts with set_medium().
    out[[1]]$subsurface = interior
  }
  out
}
