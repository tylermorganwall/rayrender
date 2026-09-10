#' Procedural Cloud Object
#' @md
#'
#' @description
#' Create a cloud volume with billowing Perlin-noise density. Add it to a scene
#' with [add_object()] and render with `integrator_type = "nee"`. Requires the
#' suggested package `ambient` to generate the density field.
#'
#' @param x Default `0`. x-coordinate of the center of the cloud's bounding box.
#' @param y Default `0`. y-coordinate of the center of the cloud's bounding box.
#' @param z Default `0`. z-coordinate of the center of the cloud's bounding box.
#' @param width Default `100`. Width of the cloud volume along its local x-axis.
#' @param height Default `25`. Height of the cloud volume along its local y-axis.
#' @param depth Default `75`. Depth of the cloud volume along its local z-axis.
#' @param style Default `c("cumulus", "stratus")`. Cloud shape. `"cumulus"` creates
#'   rounded bodies with billowing tops; `"stratus"` creates a shallow cloud bank.
#' @param seed Default `42`. Nonnegative integer seed for the cloud shape.
#'   Generating a cloud preserves the caller's random-number state.
#' @param resolution Default `128`. Integer number of density cells along the
#'   longest dimension, at least 24. Other dimensions follow the aspect ratio,
#'   with at least eight cells each. Higher values add detail and use more memory.
#' @param coverage Default `0.5`. Number between zero and one controlling the
#'   size and connection of cloud bodies. Zero does not make the volume empty;
#'   use `optical_depth = 0` for a non-scattering cloud.
#' @param detail Default `0.35`. Number between zero and one controlling the
#'   strength of small-scale Perlin detail.
#' @param optical_depth Default `8`. Nonnegative extinction through a fully
#'   dense column of length `height`. Actual optical depth depends on the density
#'   along the ray. Extinction is 99.9% scattering and 0.1% absorption.
#' @param g Default `0.65`. Henyey-Greenstein scattering asymmetry, strictly
#'   between -1 and 1. Positive values scatter forward along the light direction.
#' @param angle Default `c(0, 0, 0)`. Rotation in degrees around the x, y, and z
#'   axes, applied in the order specified by `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. Order of rotations, referring to
#'   x, y, and z. Must be a permutation of `c(1, 2, 3)`.
#' @param scale Default `c(1, 1, 1)`. Nonzero scale factors along x, y, and z.
#'   A single value scales uniformly. Scales the density field and its boundary
#'   together, retaining extinction per world-space unit.
#'
#' @details The local bounding box is centered at zero before object transforms,
#' spanning `-c(width, height, depth) / 2` to `c(width, height, depth) / 2`.
#' The density fades to vacuum at all six faces; the box has no visible surface.
#' An unrotated, unscaled cloud with base altitude `b` has `y = b + height / 2`.
#' The center describes the box, not the irregular density's center of mass.
#'
#' Positions and dimensions use scene units. Setting the dimensions rebuilds
#' the field and normalizes extinction by `height`. Applying `scale` stretches
#' the existing volume without renormalizing extinction, so stretching it along
#' a ray increases that ray's optical depth. Standard [group_objects()],
#' [animate_objects()], and [create_instances()] operations transform the cloud
#' using the same object machinery as other closed shapes.
#'
#' Cloud scattering is separate from [sky_light()]'s clear-air atmospheric haze.
#' Avoid intersecting, non-nested cloud boxes, including their empty edge cells:
#' these are separate medium boundaries. Use one larger field for a connected
#' bank. The attached [grid_medium()] is stored in
#' `object$shape_info[[1]]$medium` for further density or scattering adjustments.
#'
#' @return A single-row `ray_scene` containing an invisible box with an attached
#'   cloud density grid.
#' @seealso [sky_light()], [grid_medium()], [set_medium()]
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' if (requireNamespace("ambient", quietly = TRUE)) {
#'   # The default cloud is centered at the origin. Raise a cloud by its center
#'   # to put its base above the ground, then rotate the entire density field.
#'   puff = cloud(y = 20, width = 60, height = 20, depth = 40,
#'                angle = c(0, 25, 0), resolution = 64, optical_depth = 4)
#'   scene = generate_ground(material = diffuse("#699447")) |>
#'     add_object(puff) |>
#'     add_object(sphere(x = -50, y = 80, z = -30, radius = 15,
#'                       material = light(intensity = 40)))
#'   render_scene(scene, lookfrom = c(80, 35, -100), lookat = c(0, 20, 0),
#'                fov = 35, integrator_type = "nee", samples = 64,
#'                clamp_value = Inf, aperture = 0)
#'
#'   # Reuse the shape at another position, or change the style and its detail.
#'   bank = cloud(x = 100, y = 30, z = 20, style = "stratus", seed = 17,
#'                width = 80, height = 10, depth = 50,
#'                coverage = 0.7, detail = 0.2, optical_depth = 6, g = 0.6)
#'   scaled = cloud(scale = c(1.5, 1, 0.75), angle = c(10, 30, 0),
#'                  order_rotation = c(2, 1, 3), resolution = 64)
#' }
cloud = function(
  x = 0,
  y = 0,
  z = 0,
  width = 100,
  height = 25,
  depth = 75,
  style = c("cumulus", "stratus"),
  seed = 42,
  resolution = 128,
  coverage = 0.5,
  detail = 0.35,
  optical_depth = 8,
  g = 0.65,
  angle = c(0, 0, 0),
  order_rotation = c(1, 2, 3),
  scale = c(1, 1, 1)
) {
  for (field in c("x", "y", "z")) {
    validate_cloud_scalar(get(field), field)
  }
  for (field in c("width", "height", "depth")) {
    value = get(field)
    validate_cloud_scalar(value, field)
    if (value <= 0) {
      stop(field, " must be greater than zero.", call. = FALSE)
    }
  }
  validate_cloud_scalar(optical_depth, "optical_depth")
  if (optical_depth < 0 || !is.finite(optical_depth / height)) {
    stop(
      "optical_depth must be nonnegative, with finite optical_depth / height.",
      call. = FALSE
    )
  }
  validate_cloud_scalar(g, "g")
  if (abs(g) >= 1) {
    stop("g must be strictly between -1 and 1.", call. = FALSE)
  }
  if (!is.numeric(angle) || length(angle) != 3L || any(!is.finite(angle))) {
    stop("angle must contain three finite rotation angles.", call. = FALSE)
  }
  if (
    !is.numeric(order_rotation) ||
      length(order_rotation) != 3L ||
      anyNA(order_rotation) ||
      !all(sort(order_rotation) == 1:3)
  ) {
    stop("order_rotation must be a permutation of c(1, 2, 3).", call. = FALSE)
  }
  if (
    !is.numeric(scale) ||
      !length(scale) %in% c(1L, 3L) ||
      any(!is.finite(scale)) ||
      any(scale == 0)
  ) {
    stop(
      "scale must contain one or three finite nonzero values.",
      call. = FALSE
    )
  }

  # Generate in local coordinates so object placement never changes the field.
  # The same transform then moves both the density and its invisible boundary.
  size = c(width, height, depth)
  density = perlin_cloud_density(
    size,
    resolution,
    style,
    seed,
    coverage,
    detail
  )
  extinction = optical_depth / height
  medium = grid_medium(
    density,
    bounds = rbind(-size / 2, size / 2),
    sigma_s = extinction * 0.999,
    sigma_a = extinction * 0.001,
    g = g
  )
  set_medium(
    cube(
      x = x,
      y = y,
      z = z,
      xwidth = width,
      ywidth = height,
      zwidth = depth,
      angle = angle,
      order_rotation = order_rotation,
      scale = scale
    ),
    medium
  )
}

#' @keywords internal
validate_cloud_scalar = function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop(name, " must be one finite number.", call. = FALSE)
  }
}

#' @keywords internal
cloud_smoothstep = function(x) {
  x = pmin(1, pmax(0, x))
  x * x * (3 - 2 * x)
}

# All dimensions and cloud positions use rayrender world units; y is up.
# Increasing resolution adds detail while keeping the large cloud forms similar.
#' @keywords internal
perlin_cloud_density = function(
  size = c(100, 25, 75),
  resolution = 128,
  style = c("cumulus", "stratus"),
  seed = 42,
  coverage = 0.5,
  detail = 0.35
) {
  if (!requireNamespace("ambient", quietly = TRUE)) {
    stop("Install the ambient package to generate Perlin cloud density.")
  }
  style = match.arg(style)
  if (
    !is.numeric(size) ||
      length(size) != 3L ||
      any(!is.finite(size)) ||
      any(size <= 0)
  ) {
    stop("size must give positive finite x, y, z dimensions.", call. = FALSE)
  }
  validate_cloud_scalar(resolution, "resolution")
  if (
    resolution < 24 ||
      resolution > .Machine$integer.max ||
      resolution != floor(resolution)
  ) {
    stop(
      "resolution must be an integer from 24 to .Machine$integer.max.",
      call. = FALSE
    )
  }
  validate_cloud_scalar(coverage, "coverage")
  validate_cloud_scalar(detail, "detail")
  if (coverage < 0 || coverage > 1 || detail < 0 || detail > 1) {
    stop("coverage and detail must be between zero and one.", call. = FALSE)
  }
  validate_cloud_scalar(seed, "seed")
  if (
    seed < 0 ||
      seed > .Machine$integer.max - 1 ||
      seed != floor(seed)
  ) {
    stop(
      "seed must be a nonnegative integer smaller than .Machine$integer.max.",
      call. = FALSE
    )
  }
  dims = pmax(8L, as.integer(ceiling(resolution * (size / max(size)))))
  # noise_perlin() draws its seed from R's RNG. Restore the caller's RNG state.
  had_seed = exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed = get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit(
    {
      if (had_seed) {
        assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    },
    add = TRUE
  )
  set.seed(seed)
  broad = ambient::noise_perlin(
    dim = dims,
    frequency = 5 / max(dims),
    fractal = "fbm",
    octaves = 4,
    gain = 0.5,
    lacunarity = 2
  )
  fine = ambient::noise_perlin(
    dim = dims,
    frequency = 11 / max(dims),
    fractal = "fbm",
    octaves = 3,
    gain = 0.5,
    lacunarity = 2
  )
  # Cell centers, in the same x/y/z order as grid_medium().
  x = (seq_len(dims[1]) - 0.5) / dims[1]
  y = (seq_len(dims[2]) - 0.5) / dims[2]
  z = (seq_len(dims[3]) - 0.5) / dims[3]
  density = array(0, dims)
  if (style == "cumulus") {
    # Broad domes establish the cloud bodies; Perlin erodes and wrinkles them.
    # A shared base profile below makes flatter bottoms and billowing tops.
    centers = rbind(
      c(0.22, 0.30),
      c(0.50, 0.27),
      c(0.76, 0.45),
      c(0.32, 0.72),
      c(0.65, 0.77)
    )
    centers = centers + matrix(runif(length(centers), -0.045, 0.045), ncol = 2)
    radius = 0.16 + 0.12 * coverage
    for (j in seq_along(y)) {
      envelope = matrix(-Inf, dims[1], dims[3])
      for (i in seq_len(nrow(centers))) {
        puff_radius = radius * (0.84 + 0.14 * (i %% 3))
        horizontal = outer(
          (x - centers[i, 1]) / puff_radius,
          (z - centers[i, 2]) / (puff_radius * 1.15),
          function(a, b) a^2 + b^2
        )
        center_y = 0.30 + 0.09 * ((2 * i) %% 3)
        dome = 1 - horizontal - ((y[j] - center_y) / (0.36 + 0.07 * (i %% 3)))^2
        envelope = pmax(envelope, dome)
      }
      shape = envelope + 2.2 * broad[, j, ] + detail * 1.8 * fine[, j, ]
      density[, j, ] = cloud_smoothstep(shape / 0.65)
    }
  } else {
    for (j in seq_along(y)) {
      # A shallow, uneven deck, with larger gaps at lower coverage.
      shape = 0.55 +
        0.55 * (coverage - 0.5) +
        1.3 * broad[, j, ] +
        detail * fine[, j, ] -
        ((y[j] - 0.4) / 0.40)^2
      footprint = outer((x - 0.5) / 0.67, (z - 0.5) / 0.70, function(a, b) {
        cloud_smoothstep(1.2 - a^2 - b^2)
      })
      density[, j, ] = cloud_smoothstep(shape / 0.55) * footprint
    }
  }
  # Fade all six edges to vacuum; otherwise a density grid looks cut from a box.
  edge_x = cloud_smoothstep(pmin(x, 1 - x) / 0.10)
  edge_z = cloud_smoothstep(pmin(z, 1 - z) / 0.10)
  edge_y = cloud_smoothstep(y / 0.12) * cloud_smoothstep((1 - y) / 0.15)
  for (j in seq_along(y)) {
    density[, j, ] = density[, j, ] * outer(edge_x, edge_z) * edge_y[j]
  }
  density[c(1, dims[1]), , ] = 0
  density[, c(1, dims[2]), ] = 0
  density[,, c(1, dims[3])] = 0
  density
}
