#' Homogeneous Subsurface Scattering
#'
#' A neutral dielectric boundary containing a homogeneous RGB scattering medium.
#' Pass this material directly to a closed sphere, cube, ellipsoid, or watertight
#' triangle mesh. An interior medium is attached automatically during preparation.
#' Explicit [set_medium()] attachments on the same object are an error.
#'
#' @param color Default `"white"`. Target body reflectance in artist mode, using
#'   the same color convention as [diffuse()]. A color name or RGB triple.
#' @param radius Default `1`. Extinction mean free path in world units before
#'   scaling; a positive scalar or RGB vector.
#' @param scale Default `1`. Positive multiplier applied to `radius` in artist mode.
#' @param g Default `0`. Physical HG anisotropy, strictly between -1 and 1.
#'   Positive values favor forward propagation.
#' @param refraction Default `1.3`. Positive finite interior index of refraction.
#' @param roughness Default `0`. Boundary roughness: zero for smooth, or from
#'   `0.0001` to one for a single-scattering GGX dielectric (`alpha = roughness^2`).
#' @param sigma_a Default `NULL`. Absorption coefficients in inverse world units
#'   for physical mode; supply with `sigma_s`.
#' @param sigma_s Default `NULL`. Scattering coefficients in inverse world units
#'   for physical mode; supply with `sigma_a`. Scalars broadcast to RGB.
#' @param method Default `c("random_walk", "guided")`. Sampling method; both
#'   methods use the same transport model. The reference method remains the
#'   default pending broader production benchmarks.
#' @param priority Default `0`. Nonnegative integer dielectric priority. Lower
#'   values win in overlaps with dielectrics or other subsurface bodies. Equal
#'   values support innermost selection for properly nested bodies. Use distinct
#'   priorities for non-nested overlaps. For liquid in glass, give the glass a
#'   lower value and extend the liquid slightly into the glass wall.
#'
#' @details
#' In physical mode coefficients are used directly without color conversion.
#' Explicit `color`, `radius`, or `scale` arguments are rejected in that mode.
#' Vacuum, pure absorption, conservative scattering, and zero-extinction channels
#' are supported without adding artificial absorption.
#'
#' In artist mode, `sigma_t = 1 / (scale * radius)`. The Chiang-Kutz-Burley
#' reflectance fit splits extinction into scattering and absorption. Its
#' calibration used isotropic scattering, a semi-infinite slab, a diffuse
#' interface, and a white environment. It is a convenient authoring approximation,
#' not exact color matching for dielectric boundaries, arbitrary geometry, or
#' nonzero `g`. Black maps to zero scattering; white maps to zero absorption.
#' This distance is not a reduced mean free path or a diffusion-profile radius.
#'
#' Closed subsurface and dielectric boundaries may overlap; only the winning
#'   material scatters or absorbs in their overlap, and hidden interfaces do not
#'   reflect or refract. Ordinary explicit media still require proper nesting.
#'   Interior objects and
#' nested media use ordinary scene intersections. Body textures are not supported.
#' Atmospheric haze is excluded from the solid. Refractive interfaces block
#' straight next-event light connections; ideal delta lights behind specular
#' chains need a dedicated connection method. Use finite emitters for comparisons.
#'
#' Internal collisions and boundary reflections do not consume `max_depth`.
#' Entry consumes one ordinary event; exit is part of that event. Internal walks
#' use compensated roulette after five events, with no finite internal cap.
#' Embedded surfaces and nested entries still consume ordinary depth. The outer
#' depth limit remains a truncation. Camera-inside paths use unguided sampling.
#' The current pole guide is restricted to `g = 0`; nonzero anisotropy uses the
#' same tracked unguided fallback without changing the physical phase function.
#'
#' @references Chiang, Kutz, and Burley (2016), Practical and Controllable
#'   Subsurface Scattering for Production Path Tracing. DOI: 10.1145/2897839.2927433.
#'   d'Eon and Krivanek (2020), Zero-Variance Theory for Efficient Subsurface
#'   Scattering, section 6.6, <https://eugenedeon.com/pdfs/zv2020.pdf>.
#' @return A `ray_material` descriptor.
#' @md
#' @export
#' @examples
#' wax = subsurface(color = "ivory", radius = c(0.4, 0.2, 0.1))
#' physical = subsurface(sigma_a = c(0.1, 0.3, 0.8), sigma_s = 4)
#' thin = cube(zwidth = 0.1, material = physical)
#' reference = sphere(material = subsurface(sigma_a = 0.1, sigma_s = 4,
#'                                          method = "random_walk"))
#' guided = sphere(material = subsurface(sigma_a = 0.1, sigma_s = 4,
#'                                       method = "guided"))
subsurface = function(
  color = "white",
  radius = 1,
  scale = 1,
  g = 0,
  refraction = 1.3,
  roughness = 0,
  sigma_a = NULL,
  sigma_s = NULL,
  method = c("random_walk", "guided"),
  priority = 0
) {
  method = match.arg(method)
  scalar = function(x, name, lower, upper, strict = FALSE) {
    if (
      !is.numeric(x) ||
        length(x) != 1L ||
        !is.null(dim(x)) ||
        !is.finite(x) ||
        if (strict) x <= lower || x >= upper else x < lower || x > upper
    ) {
      stop(
        sprintf(
          "`%s` must be a finite scalar %s (%s, %s).",
          name,
          if (strict) "strictly within" else "within the closed bounds",
          lower,
          upper
        ),
        call. = FALSE
      )
    }
  }
  scalar(g, "g", -1, 1, TRUE)
  scalar(refraction, "refraction", 0, Inf, TRUE)
  scalar(roughness, "roughness", 0, 1)
  scalar(priority, "priority", 0, .Machine$integer.max)
  if (priority != floor(priority)) {
    stop("`priority` must be an integer.", call. = FALSE)
  }
  if (roughness > 0 && roughness < 0.0001) {
    stop(
      "`roughness` must be zero or at least 0.0001 for a resolvable rough boundary.",
      call. = FALSE
    )
  }
  physical = !is.null(sigma_a) || !is.null(sigma_s)
  if (physical) {
    if (is.null(sigma_a) || is.null(sigma_s)) {
      stop("Supply both `sigma_a` and `sigma_s` together.", call. = FALSE)
    }
    if (!missing(color) || !missing(radius) || !missing(scale)) {
      stop(
        "Physical coefficients cannot be combined with explicit `color`, `radius`, or `scale`.",
        call. = FALSE
      )
    }
    sigma_a = medium_rgb(sigma_a, "sigma_a")
    sigma_s = medium_rgb(sigma_s, "sigma_s")
  } else {
    radius = medium_rgb(radius, "radius")
    if (any(radius <= 0)) {
      stop("`radius` must be positive.", call. = FALSE)
    }
    scalar(scale, "scale", 0, Inf, TRUE)
    if (is.character(color) && length(color) == 1L && !is.na(color)) {
      color = convert_color(color)
    }
    if (
      !is.numeric(color) ||
        length(color) != 3L ||
        !is.null(dim(color)) ||
        any(!is.finite(color)) ||
        any(color < 0 | color > 1)
    ) {
      stop(
        "`color` must be one color name or a finite RGB triple in [0, 1].",
        call. = FALSE
      )
    }
    color = convert_color(color)
    distance = scale * radius
    if (
      any(!is.finite(distance)) ||
        any(distance <= 0) ||
        any(!is.finite(1 / distance))
    ) {
      stop(
        "`scale * radius` and its reciprocal must be finite and positive.",
        call. = FALSE
      )
    }
    albedo = -expm1(-5.09406 * color + 2.61188 * color^2 - 4.31805 * color^3)
    albedo[color == 0] = 0
    albedo[color == 1] = 1
    sigma_s = albedo / distance
    sigma_a = (1 - albedo) / distance
  }
  medium = homogeneous_medium(
    sigma_a = sigma_a,
    sigma_s = sigma_s,
    g = g,
    haze = FALSE
  )
  medium$subsurface = list(
    method = method,
    refraction = refraction,
    roughness = roughness
  )
  out = dielectric(refraction = refraction, priority = priority)
  out[[1]]$subsurface = medium
  out
}

#' Expand auto-owned subsurface interiors before medium feature detection
#' @param scene A scene, including processed children stored in instances.
#' @keywords internal
prepare_subsurface = function(scene) {
  for (i in seq_len(nrow(scene))) {
    info = scene$shape_info[[i]]
    body = scene$material[[i]]$subsurface
    if (identical(info$medium_owner, "subsurface")) {
      info$medium = NULL
      info$medium_keep_surface = NULL
      info$medium_owner = NULL
    }
    if (!is.null(body)) {
      if (!is.null(info$medium)) {
        stop(
          "A subsurface material conflicts with an explicit medium. Remove it with set_medium(scene, NULL) or replace the material.",
          call. = FALSE
        )
      }
      if (
        !scene$shape[i] %in%
          c(
            "sphere",
            "box",
            "ellipsoid",
            "obj",
            "ply",
            "mesh3d",
            "raymesh",
            1,
            5,
            6,
            9,
            12,
            13,
            14
          )
      ) {
        stop(
          "Subsurface boundaries require spheres, cubes, ellipsoids, or closed triangle meshes. Apply the material before creating instances.",
          call. = FALSE
        )
      }
      info$medium = body
      info$medium_keep_surface = TRUE
      info$medium_owner = "subsurface"
      # A body's boundary is neutral even when an imported mesh has materials.
      if (!is.null(info$shape_properties$load_material)) {
        info$shape_properties$load_material = FALSE
      }
      if (!is.null(info$shape_properties$vertex_colors)) {
        info$shape_properties$vertex_colors = FALSE
      }
      if (!is.null(info$shape_properties$override_material)) {
        info$shape_properties$override_material = TRUE
      }
    }
    original = info$shape_properties$original_scene
    if (!is.null(original)) {
      info$shape_properties$original_scene[[1]] = prepare_subsurface(original[[
        1
      ]])
    }
    scene$shape_info[[i]] = info
  }
  scene
}
