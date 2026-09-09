#' Homogeneous Participating Medium
#'
#' Describe a medium independently of the surface that contains it. Attach it to
#' closed objects with [set_medium()] and render with `integrator_type = "nee"`.
#'
#' @param sigma_a Default `0`. Absorption coefficient, a nonnegative number or RGB
#'   vector, per world-space distance unit.
#' @param sigma_s Default `1`. Scattering coefficient, a nonnegative number or RGB
#'   vector, per world-space distance unit.
#' @param density_scale Default `1`. Nonnegative multiplier for both coefficients.
#' @param g Default `0`. Henyey-Greenstein asymmetry, strictly between -1 and 1.
#'   Positive values scatter forward along the incident light direction.
#' @param emission Default `0`. Nonnegative scalar or RGB emitted radiance `Le`.
#'   Color names are also accepted. The volume source is `sigma_a * Le`.
#' @param temperature Default `NULL`. Blackbody temperature in kelvin. Mutually
#'   exclusive with nonzero `emission`.
#' @param emission_scale Default `1`. Nonnegative multiplier for emitted radiance.
#' @param temperature_scale Default `1`. Nonnegative multiplier applied after
#'   subtracting `temperature_offset` from the temperature.
#' @param temperature_offset Default `0`. Offset subtracted from temperatures.
#' @param medium_transform Default `diag(4)`. Invertible affine matrix mapping
#'   medium coordinates into the containing object's local coordinates.
#' @return A reusable `ray_medium` description.
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' fog = homogeneous_medium(sigma_s = 5, g = 0.65)
#' scene = set_medium(cube(width=1.2), fog) |>
#' add_object(sphere(y=1,x=1,radius=0.3,material=light(intensity=20))) |>
#' add_object(generate_studio(material=diffuse(color="dodgerblue")	))
#' render_scene(scene, integrator_type = "nee", sample_method="sobol")
#'
#' # Foggy cornell box
#' fog = homogeneous_medium(sigma_s = 0.1, sigma_a = 0.1, g = 0, density_scale = 0.01, temperature = 0)
#'  scene = generate_cornell(
#'  lightwidth = 10,
#'  lightdepth = 10,
#'  lightintensity = 3000
#') |>
#'  add_object(set_medium(
#'    cube(x = 555 / 2, y = 555 / 2, z = 555 / 2, width = 554),
#'    fog
#'  ))
#'render_scene(
#'  scene,
#'  integrator_type = "nee",
#'  fov = 40,
#'  samples = 256
#')
homogeneous_medium = function(
  sigma_a = 0,
  sigma_s = 1,
  density_scale = 1,
  g = 0,
  emission = 0,
  temperature = NULL,
  emission_scale = 1,
  temperature_scale = 1,
  temperature_offset = 0,
  medium_transform = diag(4)
) {
  new_medium(
    "homogeneous",
    sigma_a,
    sigma_s,
    density_scale,
    g,
    emission,
    temperature,
    emission_scale,
    temperature_scale,
    temperature_offset,
    medium_transform
  )
}

#' Dense Grid Participating Medium
#'
#' Density samples are cell centered and trilinearly interpolated. Array indices
#' correspond to x, y, z, with x varying fastest. Outside `bounds` the medium is
#' vacuum. Inside the bounds, interpolation clamps to the outermost samples.
#'
#' @inheritParams homogeneous_medium
#' @param density A nonnegative numeric array with dimensions `c(nx, ny, nz)`.
#' @param bounds Default `rbind(c(-0.5, -0.5, -0.5), c(0.5, 0.5, 0.5))`.
#'   Two rows giving minimum and maximum medium-space coordinates.
#' @param emission Default `0`. Scalar/RGB radiance or an array with dimensions
#'   `c(nx, ny, nz, 3)`. The volume source is `sigma_a * Le`.
#' @param temperature Default `NULL`. Scalar kelvin temperature or an array with
#'   the same dimensions as `density`. Mutually exclusive with nonzero emission.
#' @return A reusable `ray_medium` description.
#' @export
#' @examples
#' noise_field = ambient::noise_perlin(dim = c(300, 300, 300), frequency=0.02,octaves=8, gain=0.1)
#' noise_field[noise_field < 0] = 0
#' smoke = grid_medium(noise_field, g=0.2, sigma_a = 10)
#' scene = set_medium(cube(), smoke) |>
#'   add_object(sphere(x=3, material=light(intensity=10)))
#' render_scene(scene, integrator_type = "nee")
grid_medium = function(
  density,
  sigma_a = 0,
  sigma_s = 1,
  density_scale = 1,
  g = 0,
  bounds = rbind(c(-0.5, -0.5, -0.5), c(0.5, 0.5, 0.5)),
  emission = 0,
  temperature = NULL,
  emission_scale = 1,
  temperature_scale = 1,
  temperature_offset = 0,
  medium_transform = diag(4)
) {
  validate_medium_array(density, "density", 3L)
  if (
    !is.matrix(bounds) ||
      !identical(dim(bounds), c(2L, 3L)) ||
      !is.numeric(bounds) ||
      any(!is.finite(bounds)) ||
      any(bounds[2, ] <= bounds[1, ])
  ) {
    stop(
      "`bounds` must be a finite 2 by 3 matrix with minimum < maximum.",
      call. = FALSE
    )
  }
  if (
    !is.null(dim(emission)) && !identical(dim(emission), c(dim(density), 3L))
  ) {
    stop("Emission array dimensions must be c(dim(density), 3).", call. = FALSE)
  }
  if (
    !is.null(temperature) &&
      length(temperature) != 1L &&
      !identical(dim(temperature), dim(density))
  ) {
    stop(
      "Temperature must be a scalar or have the same dimensions as density.",
      call. = FALSE
    )
  }
  out = new_medium(
    "grid",
    sigma_a,
    sigma_s,
    density_scale,
    g,
    emission,
    temperature,
    emission_scale,
    temperature_scale,
    temperature_offset,
    medium_transform
  )
  out$density = density
  out$bounds = bounds
  out
}

#' NanoVDB Participating Medium
#'
#' Read uncompressed float grids from a NanoVDB file. Native grid coordinates and
#' transforms are retained; use `medium_transform` to place them inside the
#' boundary. Convert OpenVDB or compressed NanoVDB files externally before use.
#'
#' @inheritParams homogeneous_medium
#' @param filename Path to an uncompressed `.nvdb` file.
#' @param density_grid Default `"density"`. Name of the float density grid.
#' @param temperature_grid Default `NULL`. Optional float temperature grid name.
#' @return A reusable `ray_medium` description. Files are loaded during scene construction.
#' @export
nanovdb_medium = function(
  filename,
  sigma_a = 0,
  sigma_s = 1,
  density_scale = 1,
  g = 0,
  density_grid = "density",
  temperature_grid = NULL,
  emission = 0,
  emission_scale = 1,
  temperature_scale = 1,
  temperature_offset = 0,
  medium_transform = diag(4)
) {
  if (
    !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !file.exists(path.expand(filename))
  ) {
    stop("`filename` must name an existing NanoVDB file.", call. = FALSE)
  }
  for (name in list(density_grid, temperature_grid)) {
    if (
      !is.null(name) &&
        (!is.character(name) ||
          length(name) != 1L ||
          is.na(name) ||
          !nzchar(name))
    ) {
      stop("Grid names must be nonempty strings.", call. = FALSE)
    }
  }
  if (is.null(density_grid)) {
    stop("`density_grid` cannot be NULL.", call. = FALSE)
  }
  out = new_medium(
    "nanovdb",
    sigma_a,
    sigma_s,
    density_scale,
    g,
    emission,
    if (is.null(temperature_grid)) NULL else 0,
    emission_scale,
    temperature_scale,
    temperature_offset,
    medium_transform
  )
  out$filename = normalizePath(path.expand(filename), mustWork = TRUE)
  out$density_grid = density_grid
  out["temperature_grid"] = list(temperature_grid)
  out
}

#' Attach a Medium to Closed Objects
#'
#' Media nest: the innermost medium replaces the surrounding one. Cameras inside
#' media are detected automatically. A medium with zero absorption and scattering
#' creates a vacuum cavity. Intersecting, non-nested volume boundaries are unsupported.
#'
#' @param scene A `ray_scene` containing closed spheres, cubes, ellipsoids, or
#'   watertight consistently oriented triangle meshes.
#' @param medium A description from [homogeneous_medium()], [grid_medium()], or
#'   [nanovdb_medium()]. Use `NULL` to remove an attachment.
#' @param keep_surface Default `FALSE`. Keep the object's surface material when
#'   `TRUE`, for example to place a medium inside glass. Dielectric attenuation
#'   adds absorption to the medium; it does not add emission.
#' @return The scene with medium attachments stored independently of materials.
#' @export
#' @examples
#' water = homogeneous_medium(sigma_a = c(0.3, 0.05, 0.02), sigma_s = 0.01)
#' glass = sphere(material = dielectric())
#' scene = set_medium(glass, water, keep_surface = TRUE)
set_medium = function(scene, medium, keep_surface = FALSE) {
  if (!inherits(scene, "ray_scene")) {
    stop("`scene` must be a ray_scene.", call. = FALSE)
  }
  if (!is.null(medium) && !inherits(medium, "ray_medium")) {
    stop("`medium` must be created by a medium constructor.", call. = FALSE)
  }
  if (
    !is.logical(keep_surface) ||
      length(keep_surface) != 1L ||
      is.na(keep_surface)
  ) {
    stop("`keep_surface` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.null(medium) &&
      any(
        !scene$shape %in%
          c("sphere", "box", "ellipsoid", "obj", "ply", "mesh3d", "raymesh")
      )
  ) {
    stop(
      "Medium boundaries require spheres, cubes, ellipsoids, or closed triangle meshes. Attach media before creating instances.",
      call. = FALSE
    )
  }
  for (i in seq_len(nrow(scene))) {
    scene$shape_info[[i]]$medium = medium
    scene$shape_info[[i]]$medium_keep_surface = if (is.null(medium)) {
      NULL
    } else {
      keep_surface
    }
  }
  scene
}

#' @keywords internal
medium_rgb = function(x, name) {
  if (is.character(x) && name == "emission") {
    x = convert_color(x)
  }
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < 0)) {
    stop(
      sprintf("`%s` must contain finite nonnegative numbers.", name),
      call. = FALSE
    )
  }
  if (is.null(dim(x))) {
    if (!length(x) %in% c(1L, 3L)) {
      stop(sprintf("`%s` must be scalar or RGB.", name), call. = FALSE)
    }
    return(rep(x, length.out = 3L))
  }
  if (name != "emission" || length(dim(x)) != 4L || dim(x)[4] != 3L) {
    stop(sprintf("Invalid array dimensions for `%s`.", name), call. = FALSE)
  }
  x
}

#' @keywords internal
validate_medium_array = function(x, name, dimensions) {
  if (
    !is.numeric(x) ||
      length(dim(x)) != dimensions ||
      any(dim(x) < 1L) ||
      any(!is.finite(x)) ||
      any(x < 0)
  ) {
    stop(
      sprintf(
        "`%s` must be a finite nonnegative %sD numeric array.",
        name,
        dimensions
      ),
      call. = FALSE
    )
  }
}

#' @keywords internal
new_medium = function(
  type,
  sigma_a,
  sigma_s,
  density_scale,
  g,
  emission,
  temperature,
  emission_scale,
  temperature_scale,
  temperature_offset,
  medium_transform
) {
  sigma_a = medium_rgb(sigma_a, "sigma_a")
  sigma_s = medium_rgb(sigma_s, "sigma_s")
  emission = medium_rgb(emission, "emission")
  if (type != "grid" && !is.null(dim(emission))) {
    stop("Only grid_medium() accepts an emission array.", call. = FALSE)
  }
  controls = list(
    density_scale = density_scale,
    emission_scale = emission_scale,
    temperature_scale = temperature_scale,
    temperature_offset = temperature_offset,
    g = g
  )
  for (name in names(controls)) {
    x = controls[[name]]
    if (
      !is.numeric(x) ||
        length(x) != 1L ||
        !is.finite(x) ||
        (name %in%
          c("density_scale", "emission_scale", "temperature_scale") &&
          x < 0)
    ) {
      stop(sprintf("Invalid `%s`.", name), call. = FALSE)
    }
  }
  if (abs(g) >= 1) {
    stop("`g` must be strictly between -1 and 1.", call. = FALSE)
  }
  if (!is.null(temperature)) {
    if (
      !is.numeric(temperature) ||
        !length(temperature) ||
        any(!is.finite(temperature)) ||
        any(temperature < 0)
    ) {
      stop(
        "`temperature` must contain finite nonnegative kelvin values.",
        call. = FALSE
      )
    }
    if (any(emission != 0)) {
      stop("Specify emission or temperature, not both.", call. = FALSE)
    }
    if (type == "homogeneous" && length(temperature) != 1L) {
      stop("A homogeneous temperature must be scalar.", call. = FALSE)
    }
  }
  if (
    !is.matrix(medium_transform) ||
      !identical(dim(medium_transform), c(4L, 4L)) ||
      !is.numeric(medium_transform) ||
      any(!is.finite(medium_transform)) ||
      any(medium_transform[4, ] != c(0, 0, 0, 1)) ||
      det(medium_transform) == 0
  ) {
    stop(
      "`medium_transform` must be a finite invertible affine 4 by 4 matrix.",
      call. = FALSE
    )
  }
  structure(
    c(
      list(
        type = type,
        sigma_a = sigma_a,
        sigma_s = sigma_s,
        emission = emission,
        temperature = temperature,
        medium_transform = medium_transform
      ),
      controls
    ),
    class = "ray_medium"
  )
}

#' Convert Accumulated Volume Foreground to Straight RGB
#' @param rgb_mat Renderer output containing premultiplied RGB and scalar opacity.
#' @keywords internal
straight_volume_rgb = function(rgb_mat) {
  visible = rgb_mat$a > 0
  for (channel in c("r", "g", "b")) {
    rgb_mat[[channel]][visible] = rgb_mat[[channel]][visible] /
      rgb_mat$a[visible]
    rgb_mat[[channel]][!visible] = 0
  }
  rgb_mat$premultiplied = FALSE
  rgb_mat
}

#' Inspect Medium Attachments Before Scene Construction
#' @param scene A scene, including processed scenes stored inside instances.
#' @keywords internal
scene_medium_features = function(scene) {
  attached = FALSE
  emissive = FALSE
  for (info in scene$shape_info) {
    medium = info[["medium"]]
    if (!is.null(medium)) {
      attached = TRUE
      emissive = emissive ||
        (medium$density_scale > 0 &&
          medium$emission_scale > 0 &&
          any(medium$sigma_a > 0) &&
          (!is.null(medium$temperature) || any(medium$emission > 0)))
    }
    original = info$shape_properties$original_scene
    if (!is.null(original)) {
      child = scene_medium_features(original[[1]])
      attached = attached || child$attached
      emissive = emissive || child$emissive
    }
  }
  list(attached = attached, emissive = emissive)
}
