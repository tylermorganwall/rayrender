ray_schema_version = 2L

schema_choice = function(value, choices, arg) {
  value = match.arg(value, choices)
  value
}

schema_stop = function(path, message) {
  stop(sprintf("%s: %s", path, message), call. = FALSE)
}

is_finite_numeric = function(value, length = NULL) {
  is.numeric(value) &&
    (is.null(length) || base::length(value) == length) &&
    all(is.finite(value))
}

check_finite_numeric = function(value, path, length = NULL) {
  if (!is_finite_numeric(value, length = length)) {
    schema_stop(path, "must be finite numeric")
  }
  invisible(value)
}

check_nonnegative_numeric = function(value, path, length = NULL) {
  check_finite_numeric(value, path, length = length)
  if (any(value < 0)) {
    schema_stop(path, "must be non-negative")
  }
  invisible(value)
}

check_positive_numeric = function(value, path, length = NULL) {
  check_finite_numeric(value, path, length = length)
  if (any(value <= 0)) {
    schema_stop(path, "must be positive")
  }
  invisible(value)
}

check_scalar_logical = function(value, path) {
  if (!is.logical(value) || length(value) != 1 || is.na(value)) {
    schema_stop(path, "must be TRUE or FALSE")
  }
  invisible(value)
}

check_scalar_character = function(value, path, allow_empty = FALSE) {
  if (!is.character(value) || length(value) != 1 || is.na(value)) {
    schema_stop(path, "must be a string")
  }
  if (!allow_empty && !nzchar(value)) {
    schema_stop(path, "must not be empty")
  }
  invisible(value)
}

check_color_space = function(value, path) {
  check_scalar_character(value, path)
  if (!identical(value, "sRGB")) {
    schema_stop(path, "only sRGB is supported in schema v2 PR6")
  }
  invisible(value)
}

new_ray_descriptor = function(family, type, fields = list(), class = NULL) {
  check_scalar_character(family, "family")
  check_scalar_character(type, "type")
  descriptor = c(
    list(
      schema_version = ray_schema_version,
      type = type
    ),
    fields
  )
  class(descriptor) = unique(c(
    class %||% paste0("ray_", family, "_", type),
    paste0("ray_", family),
    "ray_descriptor"
  ))
  descriptor
}

`%||%` = function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

descriptor_family = function(x) {
  family = grep("^ray_[a-z_]+$", class(x), value = TRUE)
  family = family[!family %in% c("ray_descriptor", "ray_material")]
  if (length(family) == 0) {
    return("ray_descriptor")
  }
  family[[length(family)]]
}

#' @export
print.ray_descriptor = function(x, ...) {
  cat(sprintf(
    "<%s:%s schema_version=%i>\n",
    descriptor_family(x),
    x$type,
    x$schema_version
  ))
  fields = setdiff(names(x), c("schema_version", "type"))
  if (length(fields) > 0) {
    cat(paste0("  ", paste(fields, collapse = ", "), "\n"))
  }
  invisible(x)
}

normalize_rgb_input = function(value, role, path) {
  if (is.character(value)) {
    value = convert_color(value)
  }
  if (length(value) == 1 && is.numeric(value)) {
    value = rep(value, 3)
  }
  check_nonnegative_numeric(value, path, length = 3)
  if (role %in% c("albedo", "illuminant") && any(value > 1)) {
    schema_stop(path, "albedo and illuminant RGB values must be in [0, 1]")
  }
  as.numeric(value)
}

ensure_spectrum_descriptor = function(
  value,
  role = "albedo",
  path = "spectrum"
) {
  if (inherits(value, "ray_spectrum")) {
    return(value)
  }
  if (is.numeric(value) && length(value) == 1) {
    return(spectrum_constant(value))
  }
  spectrum_rgb(value, role = role)
}

ensure_texture_descriptor = function(
  value,
  value_type = "auto",
  role = "albedo",
  path = "texture"
) {
  if (inherits(value, "ray_texture")) {
    return(value)
  }
  texture_constant(value, value_type = value_type)
}

ensure_float_texture_descriptor = function(value, path = "texture") {
  if (inherits(value, "ray_texture")) {
    if (!is.null(value$value_type) && !identical(value$value_type, "float")) {
      schema_stop(path, "must be a float texture descriptor")
    }
    return(value)
  }
  texture_constant(value, value_type = "float")
}

check_named_list = function(value, path) {
  if (!is.list(value)) {
    schema_stop(path, "must be a list")
  }
  if (
    length(value) > 0 && (is.null(names(value)) || any(!nzchar(names(value))))
  ) {
    schema_stop(path, "must use non-empty names")
  }
  invisible(value)
}

#' Constant Spectrum Descriptor
#'
#' @param value Spectrum value.
#'
#' @return A schema-v2 spectrum descriptor.
#' @export
spectrum_constant = function(value) {
  check_nonnegative_numeric(value, "spectrum_constant(value)", length = 1)
  new_ray_descriptor("spectrum", "constant", list(value = as.numeric(value)))
}

#' RGB Spectrum Descriptor
#'
#' @param value RGB value as a color string or numeric vector.
#' @param role Spectrum role: albedo, illuminant, or unbounded.
#' @param color_space Default `"sRGB"`. Input RGB color space.
#' @param encoding Default `NULL`. Input encoding.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 spectrum descriptor.
#' @export
spectrum_rgb = function(
  value,
  role = c("albedo", "illuminant", "unbounded"),
  color_space = "sRGB",
  encoding = NULL,
  scale = 1
) {
  role = schema_choice(role, c("albedo", "illuminant", "unbounded"), "role")
  check_color_space(color_space, "spectrum_rgb(color_space)")
  if (!is.null(encoding)) {
    encoding = schema_choice(
      encoding,
      c("linear", "srgb", "legacy", "auto"),
      "encoding"
    )
  }
  check_nonnegative_numeric(scale, "spectrum_rgb(scale)", length = 1)
  new_ray_descriptor(
    "spectrum",
    "rgb",
    list(
      value = normalize_rgb_input(value, role, "spectrum_rgb(value)"),
      role = role,
      color_space = color_space,
      encoding = encoding,
      scale = as.numeric(scale)
    )
  )
}

#' Sampled Spectrum Descriptor
#'
#' @param wavelength_nm Wavelength samples in nanometers.
#' @param value Spectrum values.
#' @param role Default `NULL`. Optional spectrum role.
#' @param interpolation Default `"linear"`. Interpolation policy.
#' @param extrapolation Default `"zero"`. Extrapolation policy.
#' @param normalize Default `"none"`. Normalization policy.
#' @param units Default `NULL`. Value units.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 spectrum descriptor.
#' @export
spectrum_sampled = function(
  wavelength_nm,
  value,
  role = NULL,
  interpolation = "linear",
  extrapolation = "zero",
  normalize = "none",
  units = NULL,
  scale = 1
) {
  check_finite_numeric(wavelength_nm, "spectrum_sampled(wavelength_nm)")
  check_nonnegative_numeric(value, "spectrum_sampled(value)")
  if (length(wavelength_nm) != length(value)) {
    schema_stop(
      "spectrum_sampled(value)",
      "must have the same length as wavelength_nm"
    )
  }
  if (length(wavelength_nm) < 2) {
    schema_stop(
      "spectrum_sampled(wavelength_nm)",
      "must contain at least two samples"
    )
  }
  if (any(diff(wavelength_nm) <= 0)) {
    schema_stop(
      "spectrum_sampled(wavelength_nm)",
      "must be strictly increasing"
    )
  }
  if (!is.null(role)) {
    role = schema_choice(role, c("albedo", "illuminant", "unbounded"), "role")
  }
  interpolation = schema_choice(interpolation, "linear", "interpolation")
  extrapolation = schema_choice(
    extrapolation,
    c("zero", "constant", "error"),
    "extrapolation"
  )
  normalize = schema_choice(normalize, c("none", "y", "max"), "normalize")
  check_nonnegative_numeric(scale, "spectrum_sampled(scale)", length = 1)
  new_ray_descriptor(
    "spectrum",
    "sampled",
    list(
      wavelength_nm = as.numeric(wavelength_nm),
      value = as.numeric(value),
      role = role,
      interpolation = interpolation,
      extrapolation = extrapolation,
      normalize = normalize,
      units = units,
      scale = as.numeric(scale)
    )
  )
}

#' Blackbody Spectrum Descriptor
#'
#' @param temperature Blackbody temperature in kelvin.
#' @param scale Default `1`. Multiplicative scale.
#' @param normalize Default `TRUE`. Whether to normalize the blackbody spectrum.
#'
#' @return A schema-v2 spectrum descriptor.
#' @export
spectrum_blackbody = function(temperature, scale = 1, normalize = TRUE) {
  check_positive_numeric(
    temperature,
    "spectrum_blackbody(temperature)",
    length = 1
  )
  check_nonnegative_numeric(scale, "spectrum_blackbody(scale)", length = 1)
  check_scalar_logical(normalize, "spectrum_blackbody(normalize)")
  new_ray_descriptor(
    "spectrum",
    "blackbody",
    list(
      temperature = as.numeric(temperature),
      scale = as.numeric(scale),
      normalize = normalize
    )
  )
}

#' Named Spectrum Descriptor
#'
#' @param name Named spectrum identifier.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 spectrum descriptor.
#' @export
spectrum_named = function(name, scale = 1) {
  check_scalar_character(name, "spectrum_named(name)")
  check_nonnegative_numeric(scale, "spectrum_named(scale)", length = 1)
  new_ray_descriptor(
    "spectrum",
    "named",
    list(name = name, scale = as.numeric(scale))
  )
}

#' Cauchy IOR Spectrum Descriptor
#'
#' @param A Cauchy A coefficient.
#' @param B Default `0`. Cauchy B coefficient.
#' @param C Default `0`. Cauchy C coefficient.
#'
#' @return A schema-v2 spectrum descriptor.
#' @export
spectrum_cauchy_ior = function(A, B = 0, C = 0) {
  check_positive_numeric(A, "spectrum_cauchy_ior(A)", length = 1)
  check_finite_numeric(B, "spectrum_cauchy_ior(B)", length = 1)
  check_finite_numeric(C, "spectrum_cauchy_ior(C)", length = 1)
  new_ray_descriptor("spectrum", "cauchy_ior", list(A = A, B = B, C = C))
}

#' Sellmeier IOR Spectrum Descriptor
#'
#' @param B Sellmeier B coefficients.
#' @param C Sellmeier C coefficients.
#'
#' @return A schema-v2 spectrum descriptor.
#' @export
spectrum_sellmeier_ior = function(B, C) {
  check_finite_numeric(B, "spectrum_sellmeier_ior(B)")
  check_finite_numeric(C, "spectrum_sellmeier_ior(C)")
  if (length(B) != length(C)) {
    schema_stop("spectrum_sellmeier_ior(C)", "must have the same length as B")
  }
  new_ray_descriptor(
    "spectrum",
    "sellmeier_ior",
    list(B = as.numeric(B), C = as.numeric(C))
  )
}

new_texture_descriptor = function(type, fields) {
  new_ray_descriptor("texture", type, fields)
}

#' Constant Texture Descriptor
#'
#' @param value Texture value.
#' @param value_type Default `"auto"`. Texture value type.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_constant = function(
  value,
  value_type = c("auto", "float", "spectrum")
) {
  value_type = schema_choice(
    value_type,
    c("auto", "float", "spectrum"),
    "value_type"
  )
  if (value_type == "auto") {
    value_type = if (inherits(value, "ray_spectrum") || length(value) != 1) {
      "spectrum"
    } else {
      "float"
    }
  }
  if (value_type == "float") {
    check_finite_numeric(value, "texture_constant(value)", length = 1)
    value = as.numeric(value)
  } else {
    value = ensure_spectrum_descriptor(value, path = "texture_constant(value)")
  }
  new_texture_descriptor(
    "constant",
    list(value = value, value_type = value_type)
  )
}

validate_image_texture_fields = function(
  filename,
  encoding,
  wrap,
  filter,
  scale,
  invert,
  path
) {
  check_scalar_character(filename, paste0(path, "(filename)"))
  encoding = schema_choice(encoding, c("auto", "srgb", "linear"), "encoding")
  wrap = schema_choice(wrap, c("repeat", "black", "clamp"), "wrap")
  filter = schema_choice(filter, c("ewa", "bilinear", "nearest"), "filter")
  check_nonnegative_numeric(scale, paste0(path, "(scale)"), length = 1)
  check_scalar_logical(invert, paste0(path, "(invert)"))
  list(encoding = encoding, wrap = wrap, filter = filter)
}

#' Image Texture Descriptor
#'
#' @param filename Image filename.
#' @param value_type Default `"spectrum"`. Texture value type.
#' @param role Default `NULL`. Spectrum role for color images.
#' @param color_space Default `"sRGB"`. Input color space.
#' @param encoding Default `"auto"`. Input encoding.
#' @param wrap Default `"repeat"`. Wrap policy.
#' @param filter Default `"ewa"`. Filter policy.
#' @param mapping Default `NULL`. Texture mapping descriptor.
#' @param scale Default `1`. Multiplicative scale.
#' @param invert Default `FALSE`. Whether to invert scalar values.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_image = function(
  filename,
  value_type = c("spectrum", "float"),
  role = NULL,
  color_space = "sRGB",
  encoding = "auto",
  wrap = "repeat",
  filter = "ewa",
  mapping = NULL,
  scale = 1,
  invert = FALSE
) {
  value_type = schema_choice(value_type, c("spectrum", "float"), "value_type")
  fields = validate_image_texture_fields(
    filename,
    encoding,
    wrap,
    filter,
    scale,
    invert,
    "texture_image"
  )
  if (value_type == "spectrum") {
    if (is.null(role)) {
      schema_stop(
        "texture_image(role)",
        "must be supplied for spectrum image textures"
      )
    }
    role = schema_choice(role, c("albedo", "illuminant", "unbounded"), "role")
    check_color_space(color_space, "texture_image(color_space)")
  } else {
    role = NULL
    color_space = NULL
  }
  new_texture_descriptor(
    "image",
    c(
      list(
        filename = filename,
        value_type = value_type,
        role = role,
        color_space = color_space,
        mapping = mapping,
        scale = as.numeric(scale),
        invert = invert
      ),
      fields
    )
  )
}

#' Color Image Texture Descriptor
#'
#' @inheritParams texture_image
#' @param role Spectrum role.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_image_color = function(
  filename,
  color_space = "sRGB",
  encoding = "srgb",
  role = c("albedo", "illuminant", "unbounded"),
  wrap = "repeat",
  filter = "ewa",
  mapping = NULL,
  scale = 1
) {
  role = schema_choice(role, c("albedo", "illuminant", "unbounded"), "role")
  texture_image(
    filename = filename,
    value_type = "spectrum",
    role = role,
    color_space = color_space,
    encoding = encoding,
    wrap = wrap,
    filter = filter,
    mapping = mapping,
    scale = scale,
    invert = FALSE
  )
}

#' Scalar Image Texture Descriptor
#'
#' @inheritParams texture_image
#' @param channel Default `"luminance"`. Scalar channel.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_image_scalar = function(
  filename,
  channel = "luminance",
  encoding = "linear",
  wrap = "repeat",
  filter = "ewa",
  mapping = NULL,
  scale = 1
) {
  channel = schema_choice(
    channel,
    c("luminance", "red", "green", "blue", "alpha"),
    "channel"
  )
  texture = texture_image(
    filename = filename,
    value_type = "float",
    encoding = encoding,
    wrap = wrap,
    filter = filter,
    mapping = mapping,
    scale = scale,
    invert = FALSE
  )
  texture$channel = channel
  class(texture) = c(
    "ray_texture_image_scalar",
    setdiff(class(texture), "ray_texture_image")
  )
  texture
}

#' Checker Texture Descriptor
#'
#' @param tex1 First texture.
#' @param tex2 Second texture.
#' @param dimension Default `2`. Checker dimension.
#' @param uscale Default `1`. U scale.
#' @param vscale Default `1`. V scale.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_checker = function(tex1, tex2, dimension = 2, uscale = 1, vscale = 1) {
  if (!dimension %in% c(2, 3)) {
    schema_stop("texture_checker(dimension)", "must be 2 or 3")
  }
  check_positive_numeric(uscale, "texture_checker(uscale)", length = 1)
  check_positive_numeric(vscale, "texture_checker(vscale)", length = 1)
  new_texture_descriptor(
    "checker",
    list(
      tex1 = ensure_texture_descriptor(tex1),
      tex2 = ensure_texture_descriptor(tex2),
      dimension = as.integer(dimension),
      uscale = as.numeric(uscale),
      vscale = as.numeric(vscale)
    )
  )
}

#' Mix Texture Descriptor
#'
#' @param tex1 First texture.
#' @param tex2 Second texture.
#' @param amount Default `0.5`. Mix amount.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_mix = function(tex1, tex2, amount = 0.5) {
  if (inherits(amount, "ray_texture")) {
    amount_descriptor = amount
  } else {
    check_finite_numeric(amount, "texture_mix(amount)", length = 1)
    amount_descriptor = as.numeric(amount)
  }
  new_texture_descriptor(
    "mix",
    list(
      tex1 = ensure_texture_descriptor(tex1),
      tex2 = ensure_texture_descriptor(tex2),
      amount = amount_descriptor
    )
  )
}

#' Scale Texture Descriptor
#'
#' @param tex Texture descriptor.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_scale = function(tex, scale = 1) {
  check_finite_numeric(scale, "texture_scale(scale)", length = 1)
  new_texture_descriptor(
    "scale",
    list(texture = ensure_texture_descriptor(tex), scale = as.numeric(scale))
  )
}

#' Float Texture Descriptor
#'
#' @param value Float value.
#' @param range Default `NULL`. Optional allowed range.
#' @param encoding Default `"linear"`. Encoding.
#'
#' @return A schema-v2 texture descriptor.
#' @export
texture_float = function(value, range = NULL, encoding = "linear") {
  check_finite_numeric(value, "texture_float(value)", length = 1)
  encoding = schema_choice(encoding, c("linear", "srgb"), "encoding")
  if (!is.null(range)) {
    check_finite_numeric(range, "texture_float(range)", length = 2)
    if (range[[1]] > range[[2]]) {
      schema_stop("texture_float(range)", "minimum must be <= maximum")
    }
  }
  new_texture_descriptor(
    "float",
    list(value = as.numeric(value), range = range, encoding = encoding)
  )
}

#' Triplanar Texture Mapping Descriptor
#'
#' @param scale Default `1`. Mapping scale.
#' @param blend Default `0.2`. Blend width.
#'
#' @return A schema-v2 mapping descriptor.
#' @export
triplanar_mapping = function(scale = 1, blend = 0.2) {
  check_positive_numeric(scale, "triplanar_mapping(scale)", length = 1)
  check_nonnegative_numeric(blend, "triplanar_mapping(blend)", length = 1)
  new_ray_descriptor(
    "texture_mapping",
    "triplanar",
    list(scale = as.numeric(scale), blend = as.numeric(blend))
  )
}

#' Object Texture Mapping Descriptor
#'
#' @param transform Default `NULL`. Object mapping transform descriptor.
#'
#' @return A schema-v2 mapping descriptor.
#' @export
object_mapping = function(transform = NULL) {
  new_ray_descriptor("texture_mapping", "object", list(transform = transform))
}

#' Uniform Light Sampler Descriptor
#'
#' @return A schema-v2 light sampler descriptor.
#' @export
uniform_light_sampler = function() {
  new_ray_descriptor("light_sampler", "uniform", list())
}

#' Path Integrator Descriptor
#'
#' @param max_depth Default `50`. Maximum path depth.
#' @param regularize Default `FALSE`. Whether to regularize paths.
#' @param light_sampler Default `uniform_light_sampler()`. Light sampler.
#' @param russian_roulette Default `TRUE`. Whether to use Russian roulette.
#'
#' @return A schema-v2 integrator descriptor.
#' @export
path_integrator = function(
  max_depth = 50,
  regularize = FALSE,
  light_sampler = uniform_light_sampler(),
  russian_roulette = TRUE
) {
  check_positive_numeric(max_depth, "path_integrator(max_depth)", length = 1)
  check_scalar_logical(regularize, "path_integrator(regularize)")
  check_scalar_logical(russian_roulette, "path_integrator(russian_roulette)")
  new_ray_descriptor(
    "integrator",
    "path",
    list(
      max_depth = as.integer(max_depth),
      regularize = regularize,
      light_sampler = light_sampler,
      russian_roulette = russian_roulette
    )
  )
}

#' Volumetric Path Integrator Descriptor
#'
#' @inheritParams path_integrator
#'
#' @return A schema-v2 integrator descriptor.
#' @export
volpath_integrator = function(
  max_depth = 50,
  regularize = FALSE,
  light_sampler = uniform_light_sampler()
) {
  check_positive_numeric(max_depth, "volpath_integrator(max_depth)", length = 1)
  check_scalar_logical(regularize, "volpath_integrator(regularize)")
  new_ray_descriptor(
    "integrator",
    "volpath",
    list(
      max_depth = as.integer(max_depth),
      regularize = regularize,
      light_sampler = light_sampler
    )
  )
}

#' Random Walk Integrator Descriptor
#'
#' @param max_depth Default `50`. Maximum path depth.
#'
#' @return A schema-v2 integrator descriptor.
#' @export
random_walk_integrator = function(max_depth = 50) {
  check_positive_numeric(
    max_depth,
    "random_walk_integrator(max_depth)",
    length = 1
  )
  new_ray_descriptor(
    "integrator",
    "random_walk",
    list(max_depth = as.integer(max_depth))
  )
}

#' Sobol Sampler Descriptor
#'
#' @param pixel_samples Default `128`. Samples per pixel.
#' @param randomize Default `TRUE`. Whether to randomize samples.
#'
#' @return A schema-v2 sampler descriptor.
#' @export
sobol_sampler = function(pixel_samples = 128, randomize = TRUE) {
  check_positive_numeric(
    pixel_samples,
    "sobol_sampler(pixel_samples)",
    length = 1
  )
  check_scalar_logical(randomize, "sobol_sampler(randomize)")
  new_ray_descriptor(
    "sampler",
    "sobol",
    list(pixel_samples = as.integer(pixel_samples), randomize = randomize)
  )
}

#' Stratified Sampler Descriptor
#'
#' @param x_samples Default `4`. Horizontal samples.
#' @param y_samples Default `4`. Vertical samples.
#' @param jitter Default `TRUE`. Whether to jitter samples.
#'
#' @return A schema-v2 sampler descriptor.
#' @export
stratified_sampler = function(x_samples = 4, y_samples = 4, jitter = TRUE) {
  check_positive_numeric(x_samples, "stratified_sampler(x_samples)", length = 1)
  check_positive_numeric(y_samples, "stratified_sampler(y_samples)", length = 1)
  check_scalar_logical(jitter, "stratified_sampler(jitter)")
  new_ray_descriptor(
    "sampler",
    "stratified",
    list(
      x_samples = as.integer(x_samples),
      y_samples = as.integer(y_samples),
      jitter = jitter
    )
  )
}

#' Independent Sampler Descriptor
#'
#' @param pixel_samples Default `128`. Samples per pixel.
#'
#' @return A schema-v2 sampler descriptor.
#' @export
independent_sampler = function(pixel_samples = 128) {
  check_positive_numeric(
    pixel_samples,
    "independent_sampler(pixel_samples)",
    length = 1
  )
  new_ray_descriptor(
    "sampler",
    "independent",
    list(pixel_samples = as.integer(pixel_samples))
  )
}

#' Perspective Camera Descriptor
#'
#' @param lookfrom Default `c(0, 1, -10)`. Camera position.
#' @param lookat Default `c(0, 0, 0)`. Target point.
#' @param up Default `c(0, 1, 0)`. Up direction.
#' @param fov Default `20`. Field of view.
#' @param aperture Default `0`. Aperture.
#' @param focal_distance Default `NULL`. Focal distance.
#' @param shutteropen Default `0`. Shutter open time.
#' @param shutterclose Default `1`. Shutter close time.
#' @param initial_regions Default `character()`. Initial optical regions.
#'
#' @return A schema-v2 camera descriptor.
#' @export
perspective_camera = function(
  lookfrom = c(0, 1, -10),
  lookat = c(0, 0, 0),
  up = c(0, 1, 0),
  fov = 20,
  aperture = 0,
  focal_distance = NULL,
  shutteropen = 0,
  shutterclose = 1,
  initial_regions = character()
) {
  check_finite_numeric(lookfrom, "perspective_camera(lookfrom)", length = 3)
  check_finite_numeric(lookat, "perspective_camera(lookat)", length = 3)
  check_finite_numeric(up, "perspective_camera(up)", length = 3)
  check_positive_numeric(fov, "perspective_camera(fov)", length = 1)
  check_nonnegative_numeric(
    aperture,
    "perspective_camera(aperture)",
    length = 1
  )
  if (!is.null(focal_distance)) {
    check_positive_numeric(
      focal_distance,
      "perspective_camera(focal_distance)",
      length = 1
    )
  }
  check_finite_numeric(
    c(shutteropen, shutterclose),
    "perspective_camera(shutter)"
  )
  new_ray_descriptor(
    "camera",
    "perspective",
    list(
      lookfrom = as.numeric(lookfrom),
      lookat = as.numeric(lookat),
      up = as.numeric(up),
      fov = as.numeric(fov),
      aperture = as.numeric(aperture),
      focal_distance = focal_distance,
      shutteropen = as.numeric(shutteropen),
      shutterclose = as.numeric(shutterclose),
      initial_regions = as.character(initial_regions)
    )
  )
}

#' Orthographic Camera Descriptor
#'
#' @param ... Named orthographic camera fields.
#'
#' @return A schema-v2 camera descriptor.
#' @export
orthographic_camera = function(...) {
  new_ray_descriptor("camera", "orthographic", list(params = list(...)))
}

#' Realistic Camera Descriptor
#'
#' @param ... Named realistic camera fields.
#'
#' @return A schema-v2 camera descriptor.
#' @export
realistic_camera = function(...) {
  new_ray_descriptor("camera", "realistic", list(params = list(...)))
}

#' CIE 1931 Sensor Descriptor
#'
#' @param output_color_space Default `"sRGB"`. Output color space.
#' @param white_balance Default `NULL`. White balance.
#'
#' @return A schema-v2 sensor descriptor.
#' @export
cie1931_sensor = function(output_color_space = "sRGB", white_balance = NULL) {
  check_color_space(output_color_space, "cie1931_sensor(output_color_space)")
  new_ray_descriptor(
    "sensor",
    "cie1931",
    list(
      output_color_space = output_color_space,
      white_balance = white_balance
    )
  )
}

#' RGB Sensor Descriptor
#'
#' @param color_space Default `"sRGB"`. Sensor color space.
#' @param white_balance Default `NULL`. White balance.
#'
#' @return A schema-v2 sensor descriptor.
#' @export
rgb_sensor = function(color_space = "sRGB", white_balance = NULL) {
  check_color_space(color_space, "rgb_sensor(color_space)")
  new_ray_descriptor(
    "sensor",
    "rgb",
    list(color_space = color_space, white_balance = white_balance)
  )
}

#' Measured Sensor Descriptor
#'
#' @param red Red channel response spectrum.
#' @param green Green channel response spectrum.
#' @param blue Blue channel response spectrum.
#' @param imaging_ratio Default `1`. Imaging ratio.
#' @param white_balance Default `NULL`. White balance.
#'
#' @return A schema-v2 sensor descriptor.
#' @export
measured_sensor = function(
  red,
  green,
  blue,
  imaging_ratio = 1,
  white_balance = NULL
) {
  check_positive_numeric(
    imaging_ratio,
    "measured_sensor(imaging_ratio)",
    length = 1
  )
  new_ray_descriptor(
    "sensor",
    "measured",
    list(
      red = ensure_spectrum_descriptor(red, role = "unbounded"),
      green = ensure_spectrum_descriptor(green, role = "unbounded"),
      blue = ensure_spectrum_descriptor(blue, role = "unbounded"),
      imaging_ratio = as.numeric(imaging_ratio),
      white_balance = white_balance
    )
  )
}

#' RGB Film Descriptor
#'
#' @param width Default `NULL`. Film width.
#' @param height Default `NULL`. Film height.
#' @param sensor Default `cie1931_sensor()`. Pixel sensor descriptor.
#' @param output_color_space Default `"sRGB"`. Output color space.
#' @param white_balance Default `NULL`. White balance.
#' @param filename Default `NULL`. Output filename.
#' @param save_linear Default `TRUE`. Whether to save linear output.
#' @param alpha Default `TRUE`. Whether to store alpha.
#'
#' @return A schema-v2 film descriptor.
#' @export
rgb_film = function(
  width = NULL,
  height = NULL,
  sensor = cie1931_sensor(),
  output_color_space = "sRGB",
  white_balance = NULL,
  filename = NULL,
  save_linear = TRUE,
  alpha = TRUE
) {
  if (!is.null(width)) {
    check_positive_numeric(width, "rgb_film(width)", length = 1)
  }
  if (!is.null(height)) {
    check_positive_numeric(height, "rgb_film(height)", length = 1)
  }
  if (!inherits(sensor, "ray_sensor")) {
    schema_stop("rgb_film(sensor)", "must be a ray_sensor descriptor")
  }
  check_color_space(output_color_space, "rgb_film(output_color_space)")
  check_scalar_logical(save_linear, "rgb_film(save_linear)")
  check_scalar_logical(alpha, "rgb_film(alpha)")
  new_ray_descriptor(
    "film",
    "rgb",
    list(
      width = width,
      height = height,
      sensor = sensor,
      output_color_space = output_color_space,
      white_balance = white_balance,
      filename = filename,
      save_linear = save_linear,
      alpha = alpha
    )
  )
}

#' Spectral Options Descriptor
#'
#' @param sampling Default `"pbrt_visible_4"`. Wavelength sampling mode.
#' @param wavelength_range Default `c(360, 830)`. Wavelength interval.
#' @param input_color_space Default `"sRGB"`. Input color space.
#' @param rgb_numeric_encoding Default `"legacy"`. Numeric RGB encoding.
#' @param rgb_table_resolution Default `64`. RGB table resolution.
#' @param table_cache Default `TRUE`. Whether to cache tables.
#'
#' @return A schema-v2 spectral-options descriptor.
#' @export
spectral_options = function(
  sampling = "pbrt_visible_4",
  wavelength_range = c(360, 830),
  input_color_space = "sRGB",
  rgb_numeric_encoding = "legacy",
  rgb_table_resolution = 64,
  table_cache = TRUE
) {
  sampling = schema_choice(sampling, c("pbrt_visible_4"), "sampling")
  check_finite_numeric(
    wavelength_range,
    "spectral_options(wavelength_range)",
    length = 2
  )
  if (wavelength_range[[1]] >= wavelength_range[[2]]) {
    schema_stop(
      "spectral_options(wavelength_range)",
      "minimum must be less than maximum"
    )
  }
  check_color_space(input_color_space, "spectral_options(input_color_space)")
  rgb_numeric_encoding = schema_choice(
    rgb_numeric_encoding,
    c("legacy", "linear", "srgb"),
    "rgb_numeric_encoding"
  )
  if (!identical(as.integer(rgb_table_resolution), 64L)) {
    schema_stop(
      "spectral_options(rgb_table_resolution)",
      "must be 64 in schema v2 PR6"
    )
  }
  check_scalar_logical(table_cache, "spectral_options(table_cache)")
  new_ray_descriptor(
    "spectral_options",
    "spectral_options",
    list(
      sampling = sampling,
      wavelength_range = as.numeric(wavelength_range),
      input_color_space = input_color_space,
      rgb_numeric_encoding = rgb_numeric_encoding,
      rgb_table_resolution = as.integer(rgb_table_resolution),
      table_cache = table_cache
    ),
    class = "ray_spectral_options"
  )
}

#' Scene Validation Descriptor
#'
#' @param mode Default `"warn"`. Validation mode.
#' @param check_watertight Default `TRUE`. Whether to check watertightness.
#' @param check_region_consistency Default `TRUE`. Whether to check regions.
#' @param check_spectrum_bounds Default `TRUE`. Whether to check spectrum bounds.
#' @param check_light_geometry Default `TRUE`. Whether to check light geometry.
#' @param allow_open_regions Default `FALSE`. Whether open regions are allowed.
#' @param allow_unsampled_emitters Default `FALSE`. Whether unsampled emitters are allowed.
#'
#' @return A schema-v2 validation descriptor.
#' @export
scene_validation = function(
  mode = c("warn", "strict", "advanced", "none"),
  check_watertight = TRUE,
  check_region_consistency = TRUE,
  check_spectrum_bounds = TRUE,
  check_light_geometry = TRUE,
  allow_open_regions = FALSE,
  allow_unsampled_emitters = FALSE
) {
  mode = schema_choice(mode, c("warn", "strict", "advanced", "none"), "mode")
  for (field in c(
    "check_watertight",
    "check_region_consistency",
    "check_spectrum_bounds",
    "check_light_geometry",
    "allow_open_regions",
    "allow_unsampled_emitters"
  )) {
    check_scalar_logical(get(field), paste0("scene_validation(", field, ")"))
  }
  new_ray_descriptor(
    "validation",
    "scene_validation",
    list(
      mode = mode,
      check_watertight = check_watertight,
      check_region_consistency = check_region_consistency,
      check_spectrum_bounds = check_spectrum_bounds,
      check_light_geometry = check_light_geometry,
      allow_open_regions = allow_open_regions,
      allow_unsampled_emitters = allow_unsampled_emitters
    ),
    class = "ray_scene_validation"
  )
}

#' Constructor for schema-v2 ray materials
#'
#' @param type Material type string.
#' @param params Default `list()`. Named pbrt-style parameters.
#' @param textures Default `list()`. Named texture parameters.
#' @param normal_map Default `NULL`. Normal map descriptor.
#' @param displacement Default `NULL`. Displacement descriptor.
#' @param flags Default `list()`. Material flags.
#' @param legacy Default `list()`. Legacy provenance.
#'
#' @return A schema-v2 material descriptor.
#' @export
new_ray_material = function(
  type,
  params = list(),
  textures = list(),
  normal_map = NULL,
  displacement = NULL,
  flags = list(),
  legacy = list()
) {
  check_scalar_character(type, "new_ray_material(type)")
  check_named_list(params, "new_ray_material(params)")
  check_named_list(textures, "new_ray_material(textures)")
  check_named_list(flags, "new_ray_material(flags)")
  check_named_list(legacy, "new_ray_material(legacy)")
  vctrs::new_vctr(
    list(list(
      schema_version = ray_schema_version,
      type = type,
      params = params,
      textures = textures,
      normal_map = normal_map,
      displacement = displacement,
      flags = flags,
      legacy = legacy
    )),
    class = c("ray_material_v2", "ray_material")
  )
}

#' @export
print.ray_material_v2 = function(x, ...) {
  cat(cli::col_grey(sprintf("ray_material_v2 list <%i>\n", length(x))))
  for (i in seq_along(x)) {
    material = x[[i]]
    cat(sprintf("[[%i]] %s\n", i, cli::col_blue(material$type)))
  }
  invisible(x)
}

#' @export
pillar_shaft.ray_material_v2 = function(x, ...) {
  labels = vapply(x, function(item) sprintf("<%s>", item$type), character(1))
  pillar::new_pillar_shaft_simple(labels, width = 18)
}

#' Interface Material Descriptor
#'
#' @return A schema-v2 material descriptor.
#' @export
interface_material = function() {
  new_ray_material("interface")
}

#' @rdname interface_material
#' @export
null_material = interface_material

#' Conductor Material Descriptor
#'
#' @param eta Default `spectrum_named("metal-Ag-eta")`. Eta spectrum.
#' @param k Default `spectrum_named("metal-Ag-k")`. K spectrum.
#' @param reflectance Default `NULL`. Compatibility reflectance spectrum used only when `eta` and `k` are omitted.
#' @param roughness Default `0`. Isotropic roughness.
#' @param u_roughness Default `NULL`. U roughness.
#' @param v_roughness Default `NULL`. V roughness.
#' @param remap_roughness Default `TRUE`. Whether to remap roughness.
#' @param displacement Default `NULL`. Displacement texture.
#' @param normal_map Default `NULL`. Normal map texture.
#'
#' @return A schema-v2 material descriptor.
#' @export
conductor = function(
  eta = spectrum_named("metal-Ag-eta"),
  k = spectrum_named("metal-Ag-k"),
  reflectance = NULL,
  roughness = 0,
  u_roughness = NULL,
  v_roughness = NULL,
  remap_roughness = TRUE,
  displacement = NULL,
  normal_map = NULL
) {
  eta_missing = missing(eta)
  k_missing = missing(k)
  check_nonnegative_numeric(roughness, "conductor(roughness)", length = 1)
  check_scalar_logical(remap_roughness, "conductor(remap_roughness)")
  if (!is.null(reflectance) && (!eta_missing || !k_missing)) {
    schema_stop(
      "conductor(reflectance)",
      "cannot be combined with eta or k"
    )
  }
  if (!is.null(u_roughness)) {
    u_roughness = ensure_float_texture_descriptor(
      u_roughness,
      "conductor(u_roughness)"
    )
  }
  if (!is.null(v_roughness)) {
    v_roughness = ensure_float_texture_descriptor(
      v_roughness,
      "conductor(v_roughness)"
    )
  }
  if (is.null(u_roughness)) {
    u_roughness = ensure_float_texture_descriptor(
      roughness,
      "conductor(roughness)"
    )
  }
  if (is.null(v_roughness)) {
    v_roughness = ensure_float_texture_descriptor(
      roughness,
      "conductor(roughness)"
    )
  }
  conductor_params = list(
    u_roughness = u_roughness,
    v_roughness = v_roughness,
    remap_roughness = remap_roughness
  )
  if (is.null(reflectance)) {
    if (eta_missing) {
      eta = spectrum_named("metal-Ag-eta")
    }
    if (k_missing) {
      k = spectrum_named("metal-Ag-k")
    }
    conductor_params$eta = ensure_spectrum_descriptor(
      eta,
      role = "unbounded",
      path = "conductor(eta)"
    )
    conductor_params$k = ensure_spectrum_descriptor(
      k,
      role = "unbounded",
      path = "conductor(k)"
    )
  } else {
    conductor_params$reflectance = ensure_spectrum_descriptor(
      reflectance,
      role = "albedo",
      path = "conductor(reflectance)"
    )
  }
  new_ray_material(
    "conductor",
    params = conductor_params,
    normal_map = normal_map,
    displacement = displacement
  )
}

#' Dielectric Interface Material Descriptor
#'
#' @inheritParams conductor
#'
#' @return A schema-v2 material descriptor.
#' @export
dielectric_interface = function(
  roughness = 0,
  u_roughness = NULL,
  v_roughness = NULL,
  remap_roughness = TRUE,
  displacement = NULL,
  normal_map = NULL
) {
  check_nonnegative_numeric(
    roughness,
    "dielectric_interface(roughness)",
    length = 1
  )
  check_scalar_logical(remap_roughness, "dielectric_interface(remap_roughness)")
  new_ray_material(
    "dielectric_interface",
    params = list(
      roughness = roughness,
      u_roughness = u_roughness,
      v_roughness = v_roughness,
      remap_roughness = remap_roughness
    ),
    normal_map = normal_map,
    displacement = displacement
  )
}

#' Thin Dielectric Material Descriptor
#'
#' @param eta Default `spectrum_constant(1.5)`. Eta spectrum.
#' @param displacement Default `NULL`. Displacement texture.
#' @param normal_map Default `NULL`. Normal map texture.
#'
#' @return A schema-v2 material descriptor.
#' @export
thin_dielectric = function(
  eta = spectrum_constant(1.5),
  displacement = NULL,
  normal_map = NULL
) {
  new_ray_material(
    "thin_dielectric",
    params = list(
      eta = ensure_spectrum_descriptor(eta, role = "unbounded")
    ),
    normal_map = normal_map,
    displacement = displacement
  )
}

#' Generic schema-v2 material descriptors
#'
#' @param ... Named material fields.
#'
#' @return A schema-v2 material descriptor.
#' @name schema_v2_materials
NULL

#' @rdname schema_v2_materials
#' @export
coated_diffuse = function(...) {
  new_ray_material("coated_diffuse", params = list(...))
}

#' @rdname schema_v2_materials
#' @export
coated_conductor = function(...) {
  new_ray_material("coated_conductor", params = list(...))
}

#' @rdname schema_v2_materials
#' @export
diffuse_transmission = function(...) {
  new_ray_material("diffuse_transmission", params = list(...))
}

#' @rdname schema_v2_materials
#' @export
measured_material = function(filename) {
  check_scalar_character(filename, "measured_material(filename)")
  new_ray_material("measured", params = list(filename = filename))
}

#' @rdname schema_v2_materials
#' @export
mix_material = function(mat1, mat2, amount = 0.5) {
  check_finite_numeric(amount, "mix_material(amount)", length = 1)
  new_ray_material(
    "mix",
    params = list(
      mat1 = material_to_spectral(mat1),
      mat2 = material_to_spectral(mat2),
      amount = amount
    )
  )
}

is_legacy_material_payload = function(material) {
  is.list(material) &&
    !is.null(material$type) &&
    !is.null(material$properties) &&
    is.null(material$schema_version)
}

material_payload = function(material, path) {
  if (inherits(material, "ray_material")) {
    return(material[[1]])
  }
  if (is_legacy_material_payload(material)) {
    return(material)
  }
  schema_stop(path, "must be a ray_material object")
}

legacy_material_color = function(legacy, path) {
  color = legacy$properties[[1]]
  if (!is.numeric(color) || length(color) < 3) {
    schema_stop(path, "must contain at least three RGB color values")
  }
  as.numeric(color[seq_len(3)])
}

legacy_dielectric_properties = function(legacy, path) {
  properties = legacy$properties[[1]]
  if (!is.numeric(properties) || length(properties) < 8) {
    schema_stop(path, "legacy dielectric material is missing optical fields")
  }
  list(
    color = as.numeric(properties[seq_len(3)]),
    refraction = as.numeric(properties[[4]]),
    attenuation = as.numeric(properties[5:7]),
    priority = as.integer(properties[[8]])
  )
}

legacy_glossy_properties = function(legacy, path) {
  glossyinfo = legacy$glossyinfo[[1]]
  if (!is.numeric(glossyinfo) || length(glossyinfo) < 6) {
    schema_stop(path, "legacy glossy material is missing microfacet fields")
  }
  f0 = mean(as.numeric(glossyinfo[4:6]))
  f0 = max(0, min(0.99, f0))
  sqrt_f0 = sqrt(f0)
  eta = if (sqrt_f0 >= 1) {
    1
  } else {
    (1 + sqrt_f0) / (1 - sqrt_f0)
  }
  list(
    alpha_x = max(0, as.numeric(glossyinfo[[2]])),
    alpha_y = max(0, as.numeric(glossyinfo[[3]])),
    eta = eta
  )
}

#' Convert a material to a schema-v2 spectral material
#'
#' @param material Material descriptor.
#'
#' @return A schema-v2 material descriptor.
#' @export
material_to_spectral = function(material) {
  if (inherits(material, "ray_material_v2")) {
    return(material)
  }
  if (
    is.list(material) &&
      identical(material$schema_version, ray_schema_version) &&
      !is.null(material$params)
  ) {
    return(new_ray_material(
      material$type,
      params = material$params,
      textures = material$textures %||% list(),
      normal_map = material$normal_map,
      displacement = material$displacement,
      flags = material$flags %||% list(),
      legacy = material$legacy %||% list()
    ))
  }
  legacy = material_payload(material, "material_to_spectral(material)")
  type = get_material_name(legacy$type)
  warning = NULL
  spectral_type = paste0("legacy_", type)
  params = list()
  if (type %in% c("diffuse", "oren-nayar")) {
    params$reflectance = spectrum_rgb(
      legacy_material_color(legacy, "material_to_spectral(material)"),
      role = "albedo",
      encoding = "legacy"
    )
  } else if (type == "metal") {
    properties = legacy$properties[[1]]
    fuzz = if (is.numeric(properties) && length(properties) >= 4) {
      properties[[4]]
    } else {
      0
    }
    fuzz = max(0, min(1, as.numeric(fuzz)))
    params$reflectance = spectrum_rgb(
      legacy_material_color(legacy, "material_to_spectral(material)"),
      role = "albedo",
      encoding = "legacy"
    )
    params$u_roughness = texture_constant(fuzz, value_type = "float")
    params$v_roughness = texture_constant(fuzz, value_type = "float")
    params$remap_roughness = FALSE
    spectral_type = "compat_rgb_metal_conductor"
    warning = "`metal()` is adapted to a spectral compatibility conductor from legacy RGB reflectance and fuzz. Prefer `conductor(eta = ..., k = ...)` for measured spectral metals."
  } else if (type == "glossy") {
    glossy = legacy_glossy_properties(
      legacy,
      "material_to_spectral(material)"
    )
    params$reflectance = spectrum_rgb(
      legacy_material_color(legacy, "material_to_spectral(material)"),
      role = "albedo",
      encoding = "legacy"
    )
    params$u_roughness = texture_constant(glossy$alpha_x, value_type = "float")
    params$v_roughness = texture_constant(glossy$alpha_y, value_type = "float")
    params$eta = spectrum_constant(glossy$eta)
    params$remap_roughness = FALSE
    params$thickness = texture_constant(0.01, value_type = "float")
    params$albedo = spectrum_constant(0)
    params$g = texture_constant(0, value_type = "float")
    spectral_type = "coated_diffuse"
    warning = "`glossy()` is adapted to a spectral coated diffuse material using its legacy RGB base color, microfacet alpha, and normal-incidence reflectance."
  } else if (type %in% c("light", "spotlight")) {
    params$emission = spectrum_rgb(
      legacy_material_color(legacy, "material_to_spectral(material)"),
      role = "illuminant",
      encoding = "legacy"
    )
    params$scale = legacy$lightintensity
    warning = "`light()` as a material is a legacy shorthand. In spectral mode, prefer `with_light(area_light(...))`."
  } else if (type == "dielectric") {
    properties = legacy_dielectric_properties(
      legacy,
      "material_to_spectral(material)"
    )
    params$eta = spectrum_constant(properties$refraction)
    params$priority = properties$priority
    spectral_type = "legacy_dielectric_interface"
    warning = "`dielectric()` is adapted to a spectral dielectric interface plus optical-region metadata during scene conversion."
  }
  new_ray_material(
    spectral_type,
    params = params,
    legacy = list(payload = legacy, warning = warning)
  )
}

#' Convert a material to a legacy material
#'
#' @param material Material descriptor.
#'
#' @return A legacy material descriptor.
#' @export
material_to_legacy = function(material) {
  if (
    inherits(material, "ray_material") && !inherits(material, "ray_material_v2")
  ) {
    return(material)
  }
  schema_stop(
    "material_to_legacy(material)",
    "schema-v2 materials cannot be rendered by the legacy RGB renderer"
  )
}

#' Phase Function Descriptor
#'
#' @param g Default `0`. Henyey-Greenstein asymmetry.
#'
#' @return A schema-v2 phase descriptor.
#' @export
henyey_greenstein_phase = function(g = 0) {
  check_finite_numeric(g, "henyey_greenstein_phase(g)", length = 1)
  if (g <= -1 || g >= 1) {
    schema_stop(
      "henyey_greenstein_phase(g)",
      "must be greater than -1 and less than 1"
    )
  }
  new_ray_descriptor("phase", "henyey_greenstein", list(g = as.numeric(g)))
}

#' Homogeneous Medium Descriptor
#'
#' @param sigma_a Default `spectrum_constant(0)`. Absorption spectrum. Must remain zero until spectral participating media are implemented.
#' @param sigma_s Default `spectrum_constant(0)`. Scattering spectrum.
#' @param scale Default `1`. Density scale.
#' @param phase Default `henyey_greenstein_phase(0)`. Phase function.
#' @param emission Default `spectrum_constant(0)`. Emission spectrum.
#'
#' @return A schema-v2 medium descriptor.
#' @export
homogeneous_medium = function(
  sigma_a = spectrum_constant(0),
  sigma_s = spectrum_constant(0),
  scale = 1,
  phase = henyey_greenstein_phase(0),
  emission = spectrum_constant(0)
) {
  check_nonnegative_numeric(scale, "homogeneous_medium(scale)", length = 1)
  new_ray_descriptor(
    "medium",
    "homogeneous",
    list(
      sigma_a = ensure_spectrum_descriptor(sigma_a, role = "unbounded"),
      sigma_s = ensure_spectrum_descriptor(sigma_s, role = "unbounded"),
      scale = as.numeric(scale),
      phase = phase,
      emission = ensure_spectrum_descriptor(emission, role = "illuminant")
    )
  )
}

#' Optical Region Descriptor
#'
#' @param id Region identifier.
#' @param eta Default `spectrum_constant(1)`. Eta spectrum.
#' @param priority Default `0L`. Region priority.
#' @param medium Default `NULL`. Optional medium descriptor.
#'
#' @return A schema-v2 optical-region descriptor.
#' @export
optical_region = function(
  id,
  eta = spectrum_constant(1),
  priority = 0L,
  medium = NULL
) {
  check_scalar_character(id, "optical_region(id)")
  check_finite_numeric(priority, "optical_region(priority)", length = 1)
  new_ray_descriptor(
    "region",
    "optical_region",
    list(
      id = id,
      eta = ensure_spectrum_descriptor(eta, role = "unbounded"),
      priority = as.integer(priority),
      medium = medium
    )
  )
}

#' Region Boundary Descriptor
#'
#' @param region Region descriptor or region id.
#' @param side Boundary side.
#'
#' @return A schema-v2 region-boundary descriptor.
#' @export
region_boundary = function(
  region,
  side = c("negative_normal", "positive_normal", "inside")
) {
  side = schema_choice(
    side,
    c("negative_normal", "positive_normal", "inside"),
    "side"
  )
  region_id = if (inherits(region, "ray_region")) {
    region$id
  } else {
    check_scalar_character(region, "region_boundary(region)")
    region
  }
  new_ray_descriptor(
    "region_boundary",
    "region_boundary",
    list(region = region_id, side = side)
  )
}

#' Dielectric Region Convenience Descriptor
#'
#' @param eta Default `spectrum_constant(1.5)`. Eta spectrum.
#' @param priority Default `0L`. Region priority.
#' @param sigma_a Default `spectrum_constant(0)`. Absorption spectrum.
#' @param medium Default `NULL`. Optional medium descriptor.
#' @param id Default `NULL`. Region id.
#'
#' @return A schema-v2 optical-region descriptor.
#'
#' Nonzero absorption is rejected in the constant-IOR spectral region stage.
#' @export
dielectric_region = function(
  eta = spectrum_constant(1.5),
  priority = 0L,
  sigma_a = spectrum_constant(0),
  medium = NULL,
  id = NULL
) {
  if (!is.null(medium) && !identical(sigma_a, spectrum_constant(0))) {
    schema_stop("dielectric_region(sigma_a)", "cannot be supplied with medium")
  }
  if (is.null(medium) && !identical(sigma_a, spectrum_constant(0))) {
    schema_stop(
      "dielectric_region(sigma_a)",
      "nonzero absorption is deferred until spectral participating media"
    )
  }
  if (is.null(id)) {
    id = "dielectric_region"
  }
  optical_region(id = id, eta = eta, priority = priority, medium = medium)
}

new_light_descriptor = function(type, fields, area = FALSE, infinite = FALSE) {
  descriptor = new_ray_descriptor("light", type, fields)
  attr(descriptor, "area") = area
  attr(descriptor, "infinite") = infinite
  attr(descriptor, "free") = !area && !infinite
  descriptor
}

#' Point Light Descriptor
#'
#' @param position Light position.
#' @param intensity Default `NULL`. Spectral intensity.
#' @param power Default `NULL`. Spectral power.
#' @param emission Default `spectrum_rgb("white", role = "illuminant")`. Emission spectrum.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 light descriptor.
#' @export
point_light = function(
  position,
  intensity = NULL,
  power = NULL,
  emission = spectrum_rgb("white", role = "illuminant"),
  scale = 1
) {
  check_finite_numeric(position, "point_light(position)", length = 3)
  if (!is.null(intensity) && !is.null(power)) {
    schema_stop("point_light(power)", "cannot be supplied with intensity")
  }
  check_nonnegative_numeric(scale, "point_light(scale)", length = 1)
  new_light_descriptor(
    "point",
    list(
      position = as.numeric(position),
      intensity = intensity,
      power = power,
      emission = ensure_spectrum_descriptor(emission, role = "illuminant"),
      scale = as.numeric(scale)
    )
  )
}

#' Spot Light Descriptor
#'
#' @param from Light position.
#' @param to Light target.
#' @param intensity Default `NULL`. Spectral intensity.
#' @param emission Default `spectrum_rgb("white", role = "illuminant")`. Emission spectrum.
#' @param cone_angle Default `30`. Cone angle in degrees.
#' @param cone_delta_angle Default `5`. Cone falloff angle in degrees.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 light descriptor.
#' @export
spot_light = function(
  from,
  to,
  intensity = NULL,
  emission = spectrum_rgb("white", role = "illuminant"),
  cone_angle = 30,
  cone_delta_angle = 5,
  scale = 1
) {
  check_finite_numeric(from, "spot_light(from)", length = 3)
  check_finite_numeric(to, "spot_light(to)", length = 3)
  check_positive_numeric(cone_angle, "spot_light(cone_angle)", length = 1)
  check_nonnegative_numeric(
    cone_delta_angle,
    "spot_light(cone_delta_angle)",
    length = 1
  )
  check_nonnegative_numeric(scale, "spot_light(scale)", length = 1)
  new_light_descriptor(
    "spot",
    list(
      from = as.numeric(from),
      to = as.numeric(to),
      intensity = intensity,
      emission = ensure_spectrum_descriptor(emission, role = "illuminant"),
      cone_angle = as.numeric(cone_angle),
      cone_delta_angle = as.numeric(cone_delta_angle),
      scale = as.numeric(scale)
    )
  )
}

#' Distant Light Descriptor
#'
#' @param direction Light direction.
#' @param radiance Default `spectrum_rgb("white", role = "illuminant")`. Radiance spectrum.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 light descriptor.
#' @export
distant_light = function(
  direction,
  radiance = spectrum_rgb("white", role = "illuminant"),
  scale = 1
) {
  check_finite_numeric(direction, "distant_light(direction)", length = 3)
  check_nonnegative_numeric(scale, "distant_light(scale)", length = 1)
  new_light_descriptor(
    "distant",
    list(
      direction = as.numeric(direction),
      radiance = ensure_spectrum_descriptor(radiance, role = "illuminant"),
      scale = as.numeric(scale)
    )
  )
}

#' Uniform Infinite Light Descriptor
#'
#' @param radiance Default `spectrum_rgb("white", role = "illuminant")`. Radiance spectrum.
#' @param scale Default `1`. Multiplicative scale.
#'
#' @return A schema-v2 infinite-light descriptor.
#' @export
uniform_infinite_light = function(
  radiance = spectrum_rgb("white", role = "illuminant"),
  scale = 1
) {
  check_nonnegative_numeric(scale, "uniform_infinite_light(scale)", length = 1)
  new_light_descriptor(
    "uniform_infinite",
    list(
      radiance = ensure_spectrum_descriptor(radiance, role = "illuminant"),
      scale = as.numeric(scale)
    ),
    infinite = TRUE
  )
}

#' Image Infinite Light Descriptor
#'
#' @param filename Environment image filename.
#' @param scale Default `1`. Multiplicative scale.
#' @param rotation Default `0`. Rotation in degrees.
#' @param encoding Default `"auto"`. Input encoding.
#' @param color_space Default `"sRGB"`. Input color space.
#' @param importance_sample Default `TRUE`. Whether to importance sample.
#'
#' @return A schema-v2 infinite-light descriptor.
#' @export
image_infinite_light = function(
  filename,
  scale = 1,
  rotation = 0,
  encoding = "auto",
  color_space = "sRGB",
  importance_sample = TRUE
) {
  check_scalar_character(filename, "image_infinite_light(filename)")
  check_nonnegative_numeric(scale, "image_infinite_light(scale)", length = 1)
  check_finite_numeric(rotation, "image_infinite_light(rotation)", length = 1)
  encoding = schema_choice(encoding, c("auto", "srgb", "linear"), "encoding")
  check_color_space(color_space, "image_infinite_light(color_space)")
  check_scalar_logical(
    importance_sample,
    "image_infinite_light(importance_sample)"
  )
  new_light_descriptor(
    "image_infinite",
    list(
      filename = filename,
      scale = as.numeric(scale),
      rotation = as.numeric(rotation),
      encoding = encoding,
      color_space = color_space,
      importance_sample = importance_sample
    ),
    infinite = TRUE
  )
}

#' Area Light Descriptor
#'
#' @param emission Default `spectrum_rgb("white", role = "illuminant")`. Emission spectrum.
#' @param radiance Default `NULL`. Explicit radiance spectrum.
#' @param scale Default `1`. Multiplicative scale.
#' @param power Default `NULL`. Optional power target.
#' @param filename Default `NULL`. Optional emission texture filename.
#' @param two_sided Default `FALSE`. Whether the light emits from both sides.
#' @param visible Default `TRUE`. Whether the emitting primitive is visible.
#' @param sampling Default `"auto"`. Sampling policy.
#'
#' @return A schema-v2 area-light descriptor.
#' @export
area_light = function(
  emission = spectrum_rgb("white", role = "illuminant"),
  radiance = NULL,
  scale = 1,
  power = NULL,
  filename = NULL,
  two_sided = FALSE,
  visible = TRUE,
  sampling = c("auto", "sampled", "none")
) {
  sampling = schema_choice(sampling, c("auto", "sampled", "none"), "sampling")
  check_nonnegative_numeric(scale, "area_light(scale)", length = 1)
  check_scalar_logical(two_sided, "area_light(two_sided)")
  check_scalar_logical(visible, "area_light(visible)")
  if (!is.null(radiance)) {
    radiance = ensure_spectrum_descriptor(radiance, role = "illuminant")
  }
  new_light_descriptor(
    "area",
    list(
      emission = ensure_spectrum_descriptor(emission, role = "illuminant"),
      radiance = radiance,
      scale = as.numeric(scale),
      power = power,
      filename = filename,
      two_sided = two_sided,
      visible = visible,
      sampling = sampling
    ),
    area = TRUE
  )
}

empty_ray_scene_v2 = function(validation = scene_validation()) {
  scene = data.frame(
    x = numeric(),
    y = numeric(),
    z = numeric(),
    shape = character(),
    stringsAsFactors = FALSE
  )
  scene$material = list()
  scene$shape_info = list()
  scene$transforms = list()
  scene$animation_info = list()
  scene$light = list()
  scene$region_boundaries = list()
  scene$medium_interface = list()
  scene$shape_capabilities = list()
  scene$object_id = integer()
  scene$object_name = character()
  scene$visibility = list()
  scene$user_data = list()
  class(scene) = c("ray_scene_v2", "ray_scene", "tbl_df", "tbl", "data.frame")
  restore_ray_scene_attrs(scene, list(validation = validation))
}

as_v2_list_column = function(column, n) {
  if (n == 0) {
    return(list())
  }
  lapply(seq_len(n), function(i) column[i])
}

normalize_v2_row_columns = function(scene) {
  n = nrow(scene)
  vctrs_columns = list(
    material = "ray_material",
    shape_info = "ray_shape_info",
    transforms = "ray_transform",
    animation_info = "ray_animated_transform"
  )
  for (name in names(vctrs_columns)) {
    if (
      name %in% names(scene) && inherits(scene[[name]], vctrs_columns[[name]])
    ) {
      scene[[name]] = as_v2_list_column(scene[[name]], n)
    }
  }
  scene
}

default_shape_capabilities = function(shape) {
  is_csg = identical(shape, "csg_object") || grepl("^csg", shape)
  list(
    uv = !is_csg,
    generated_mapping = TRUE,
    area_light = !is_csg,
    region_boundary = TRUE,
    coherent_inside = is_csg
  )
}

default_visibility = function() {
  list(
    camera = TRUE,
    shadow = TRUE,
    diffuse = TRUE,
    glossy = TRUE,
    transmission = TRUE
  )
}

ray_scene_attrs = function(scene) {
  list(
    ray_schema_version = attr(scene, "ray_schema_version") %||%
      ray_schema_version,
    regions = attr(scene, "regions") %||% list(),
    lights = attr(scene, "lights") %||% list(),
    environment = attr(scene, "environment"),
    named_textures = attr(scene, "named_textures") %||% list(),
    named_materials = attr(scene, "named_materials") %||% list(),
    render_defaults = attr(scene, "render_defaults") %||% list(),
    validation = attr(scene, "validation") %||% scene_validation(),
    conversion_report = attr(scene, "conversion_report")
  )
}

restore_ray_scene_attrs = function(scene, attrs) {
  attr(scene, "ray_schema_version") = attrs$ray_schema_version %||%
    ray_schema_version
  attr(scene, "regions") = attrs$regions %||% list()
  attr(scene, "lights") = attrs$lights %||% list()
  attr(scene, "environment") = attrs$environment
  attr(scene, "named_textures") = attrs$named_textures %||% list()
  attr(scene, "named_materials") = attrs$named_materials %||% list()
  attr(scene, "render_defaults") = attrs$render_defaults %||% list()
  attr(scene, "validation") = attrs$validation %||% scene_validation()
  attr(scene, "conversion_report") = attrs$conversion_report
  scene
}

merge_named_registry = function(lhs, rhs, label) {
  duplicate = intersect(names(lhs), names(rhs))
  if (length(duplicate) > 0) {
    schema_stop(label, sprintf("duplicate name `%s`", duplicate[[1]]))
  }
  c(lhs, rhs)
}

merge_ray_scene_attrs = function(lhs, rhs) {
  list(
    ray_schema_version = ray_schema_version,
    regions = merge_named_registry(
      lhs$regions %||% list(),
      rhs$regions %||% list(),
      "regions"
    ),
    lights = c(lhs$lights %||% list(), rhs$lights %||% list()),
    environment = rhs$environment %||% lhs$environment,
    named_textures = merge_named_registry(
      lhs$named_textures %||% list(),
      rhs$named_textures %||% list(),
      "named_textures"
    ),
    named_materials = merge_named_registry(
      lhs$named_materials %||% list(),
      rhs$named_materials %||% list(),
      "named_materials"
    ),
    render_defaults = c(
      lhs$render_defaults %||% list(),
      rhs$render_defaults %||% list()
    ),
    validation = lhs$validation %||% rhs$validation %||% scene_validation(),
    conversion_report = lhs$conversion_report %||% rhs$conversion_report
  )
}

assign_missing_object_ids = function(scene) {
  if (!"object_id" %in% names(scene)) {
    scene$object_id = rep(NA_integer_, nrow(scene))
  }
  missing_id = is.na(scene$object_id)
  if (any(missing_id)) {
    max_id = suppressWarnings(max(
      scene$object_id[!missing_id],
      0L,
      na.rm = TRUE
    ))
    scene$object_id[missing_id] = seq.int(
      max_id + 1L,
      length.out = sum(missing_id)
    )
  }
  scene
}

ensure_ray_scene_v2 = function(scene) {
  if (is.null(scene)) {
    return(empty_ray_scene_v2())
  }
  if (!inherits(scene, "data.frame")) {
    schema_stop("scene", "must be a ray scene data frame")
  }
  scene = normalize_v2_row_columns(scene)
  n = nrow(scene)
  add_list_col = function(scene, name, value_fun) {
    if (!name %in% names(scene)) {
      values = rep(list(NULL), n)
      for (i in seq_len(n)) {
        values[i] = list(value_fun(i))
      }
      scene[[name]] = values
    }
    scene
  }
  scene = add_list_col(scene, "light", function(i) NULL)
  scene = add_list_col(scene, "region_boundaries", function(i) list())
  scene = add_list_col(scene, "medium_interface", function(i) NULL)
  scene = add_list_col(scene, "shape_capabilities", function(i) {
    default_shape_capabilities(scene$shape[[i]])
  })
  scene = add_list_col(scene, "visibility", function(i) default_visibility())
  scene = add_list_col(scene, "user_data", function(i) list())
  if (!"object_id" %in% names(scene)) {
    scene$object_id = rep(NA_integer_, n)
  }
  if (!"object_name" %in% names(scene)) {
    scene$object_name = rep(NA_character_, n)
  }
  class(scene) = unique(c(
    "ray_scene_v2",
    "ray_scene",
    "tbl_df",
    "tbl",
    "data.frame"
  ))
  scene = restore_ray_scene_attrs(scene, ray_scene_attrs(scene))
  assign_missing_object_ids(scene)
}

#' Schema-v2 Scene Constructor
#'
#' @param objects Default `NULL`. Optional object rows.
#' @param validation Default `scene_validation()`. Scene validation descriptor.
#'
#' @return A schema-v2 ray scene.
#' @export
ray_scene_v2 = function(objects = NULL, validation = scene_validation()) {
  scene = empty_ray_scene_v2(validation = validation)
  if (!is.null(objects)) {
    scene = add_object(scene, objects)
  }
  scene
}

#' Attach an area light to object rows
#'
#' @param object Object row or scene.
#' @param light Area light descriptor.
#' @param replace Default `TRUE`. Whether to replace an existing area light.
#'
#' @return Object rows with schema-v2 light metadata.
#' @export
with_light = function(object, light, replace = TRUE) {
  if (!inherits(light, "ray_light")) {
    schema_stop("with_light(light)", "must be a ray_light object")
  }
  if (!isTRUE(attr(light, "area"))) {
    schema_stop(
      "with_light(light)",
      "only accepts area lights; use add_light() for free lights"
    )
  }
  object = ensure_ray_scene_v2(object)
  has_light = !vapply(object$light, is.null, logical(1))
  if (any(has_light) && !replace) {
    schema_stop(
      "with_light(object)",
      "object already has an area light; use replace = TRUE"
    )
  }
  object$light = rep(list(light), nrow(object))
  object
}

#' Attach optical-region boundaries to object rows
#'
#' @param object Object row or scene.
#' @param boundary Region-boundary descriptor or list of descriptors.
#' @param append Default `TRUE`. Whether to append to existing boundaries.
#'
#' @return Object rows with schema-v2 region metadata.
#' @export
with_region_boundary = function(object, boundary, append = TRUE) {
  object = ensure_ray_scene_v2(object)
  if (inherits(boundary, "ray_region_boundary")) {
    boundary = list(boundary)
  }
  if (
    !is.list(boundary) ||
      !all(vapply(boundary, inherits, logical(1), "ray_region_boundary"))
  ) {
    schema_stop(
      "with_region_boundary(boundary)",
      "must be a ray_region_boundary or list of boundaries"
    )
  }
  check_scalar_logical(append, "with_region_boundary(append)")
  if (!append) {
    object$region_boundaries = rep(list(boundary), nrow(object))
  } else {
    object$region_boundaries = Map(
      function(old) c(old, boundary),
      object$region_boundaries
    )
  }
  object
}

#' Add a free light to a schema-v2 scene
#'
#' @param scene Scene.
#' @param light Free or infinite light descriptor.
#'
#' @return A schema-v2 scene.
#' @export
add_light = function(scene, light) {
  if (!inherits(light, "ray_light")) {
    schema_stop("add_light(light)", "must be a ray_light object")
  }
  if (isTRUE(attr(light, "area"))) {
    schema_stop(
      "add_light(light)",
      "area lights must be attached to geometry with with_light()"
    )
  }
  if (isTRUE(attr(light, "infinite"))) {
    return(set_environment(scene, light))
  }
  scene = ensure_ray_scene_v2(scene)
  attrs = ray_scene_attrs(scene)
  attrs$lights = c(attrs$lights, list(light))
  restore_ray_scene_attrs(scene, attrs)
}

#' Set the scene environment light
#'
#' @param scene Scene.
#' @param environment Infinite light descriptor or `NULL`.
#'
#' @return A schema-v2 scene.
#' @export
set_environment = function(scene, environment) {
  if (
    !is.null(environment) &&
      (!inherits(environment, "ray_light") ||
        !isTRUE(attr(environment, "infinite")))
  ) {
    schema_stop(
      "set_environment(environment)",
      "must be an infinite ray_light object or NULL"
    )
  }
  scene = ensure_ray_scene_v2(scene)
  attrs = ray_scene_attrs(scene)
  attrs$environment = environment
  restore_ray_scene_attrs(scene, attrs)
}

#' Add an optical region to a schema-v2 scene
#'
#' @param scene Scene.
#' @param region Region descriptor.
#'
#' @return A schema-v2 scene.
#' @export
add_region = function(scene, region) {
  if (!inherits(region, "ray_region")) {
    schema_stop("add_region(region)", "must be a ray_region object")
  }
  scene = ensure_ray_scene_v2(scene)
  attrs = ray_scene_attrs(scene)
  if (region$id %in% names(attrs$regions)) {
    schema_stop(sprintf("region[%s]", region$id), "duplicate region id")
  }
  attrs$regions[[region$id]] = region
  restore_ray_scene_attrs(scene, attrs)
}

#' Add a named texture to a schema-v2 scene
#'
#' @param scene Scene.
#' @param name Texture name.
#' @param texture Texture descriptor.
#'
#' @return A schema-v2 scene.
#' @export
add_named_texture = function(scene, name, texture) {
  check_scalar_character(name, "add_named_texture(name)")
  if (!inherits(texture, "ray_texture")) {
    schema_stop(
      "add_named_texture(texture)",
      "must be a ray_texture descriptor"
    )
  }
  scene = ensure_ray_scene_v2(scene)
  attrs = ray_scene_attrs(scene)
  if (name %in% names(attrs$named_textures)) {
    schema_stop(sprintf("named_textures[%s]", name), "duplicate named texture")
  }
  attrs$named_textures[[name]] = texture
  restore_ray_scene_attrs(scene, attrs)
}

#' Add a named material to a schema-v2 scene
#'
#' @param scene Scene.
#' @param name Material name.
#' @param material Material descriptor.
#'
#' @return A schema-v2 scene.
#' @export
add_named_material = function(scene, name, material) {
  check_scalar_character(name, "add_named_material(name)")
  material = material_to_spectral(material)
  scene = ensure_ray_scene_v2(scene)
  attrs = ray_scene_attrs(scene)
  if (name %in% names(attrs$named_materials)) {
    schema_stop(
      sprintf("named_materials[%s]", name),
      "duplicate named material"
    )
  }
  attrs$named_materials[[name]] = material
  restore_ray_scene_attrs(scene, attrs)
}

#' Reference a named texture
#'
#' @param name Texture name.
#'
#' @return A named texture reference descriptor.
#' @export
named_texture = function(name) {
  check_scalar_character(name, "named_texture(name)")
  new_ray_descriptor(
    "texture_reference",
    "named_texture",
    list(name = name),
    class = "ray_named_texture"
  )
}

#' Reference a named material
#'
#' @param name Material name.
#'
#' @return A named material reference descriptor.
#' @export
named_material = function(name) {
  check_scalar_character(name, "named_material(name)")
  new_ray_descriptor(
    "material_reference",
    "named_material",
    list(name = name),
    class = "ray_named_material"
  )
}

validation_mode = function(validation) {
  if (inherits(validation, "ray_scene_validation")) {
    validation$mode
  } else {
    "warn"
  }
}

emit_validation = function(mode, path, message) {
  if (mode == "none") {
    return(invisible(FALSE))
  }
  text = sprintf("%s: %s", path, message)
  if (mode == "strict") {
    stop(text, call. = FALSE)
  }
  warning(text, call. = FALSE)
  invisible(FALSE)
}

#' Validate a schema-v2 scene
#'
#' @param scene Scene.
#' @param validation Default `attr(scene, "validation") %||% scene_validation()`. Validation descriptor.
#'
#' @return The input scene, invisibly.
#' @export
validate_ray_scene_v2 = function(
  scene,
  validation = attr(scene, "validation") %||% scene_validation()
) {
  scene = ensure_ray_scene_v2(scene)
  mode = validation_mode(validation)
  if (mode == "none") {
    return(invisible(scene))
  }
  regions = attr(scene, "regions") %||% list()
  region_names = names(regions)
  for (i in seq_len(nrow(scene))) {
    for (boundary in scene$region_boundaries[[i]]) {
      if (!boundary$region %in% region_names) {
        emit_validation(
          mode,
          sprintf("object[%i].region[%s]", i, boundary$region),
          "region is not registered"
        )
      }
    }
    if (!is.null(scene$light[[i]]) && isTRUE(attr(scene$light[[i]], "area"))) {
      caps = scene$shape_capabilities[[i]]
      if (!isTRUE(caps$area_light)) {
        emit_validation(
          mode,
          sprintf("object[%i].light", i),
          "area light is not supported by this shape capability"
        )
      }
    }
  }
  invisible(scene)
}

new_conversion_report = function(warnings = character(), entries = list()) {
  structure(
    list(warnings = unique(warnings), entries = entries),
    class = "ray_schema_conversion_report"
  )
}

#' @export
print.ray_schema_conversion_report = function(x, ...) {
  cat(sprintf(
    "<ray_schema_conversion_report warnings=%i entries=%i>\n",
    length(x$warnings),
    length(x$entries)
  ))
  if (length(x$warnings) > 0) {
    cat(paste0("  - ", x$warnings, collapse = "\n"), "\n")
  }
  invisible(x)
}

emit_conversion_warnings = function(report, validation) {
  mode = validation_mode(validation)
  if (mode %in% c("warn", "advanced") && length(report$warnings) > 0) {
    for (warning_text in report$warnings) {
      warning(warning_text, call. = FALSE)
    }
  }
}

#' Convert a legacy scene to schema v2
#'
#' @param scene Legacy or schema-v2 scene.
#' @param validation Default `scene_validation()`. Validation descriptor.
#'
#' @return A schema-v2 scene with a conversion report.
#' @export
legacy_scene_to_schema_v2 = function(scene, validation = scene_validation()) {
  scene = ensure_ray_scene_v2(scene)
  warnings = character()
  entries = vector("list", nrow(scene))
  attrs = ray_scene_attrs(scene)
  for (i in seq_len(nrow(scene))) {
    material = scene$material[[i]]
    if (
      (inherits(material, "ray_material") ||
        is_legacy_material_payload(material)) &&
        !inherits(material, "ray_material_v2")
    ) {
      converted = material_to_spectral(material)
      legacy = material_payload(material, "legacy_scene_to_schema_v2(scene)")
      legacy_warning = converted[[1]]$legacy$warning
      if (!is.null(legacy_warning)) {
        warnings = c(warnings, legacy_warning)
      }
      legacy_type = get_material_name(legacy$type)
      if (legacy_type %in% c("light", "spotlight")) {
        scene$light[[i]] = area_light(
          emission = converted[[1]]$params$emission,
          scale = converted[[1]]$params$scale %||% 1
        )
        scene$material[[i]] = interface_material()
      } else if (identical(legacy_type, "dielectric")) {
        properties = legacy_dielectric_properties(
          legacy,
          "legacy_scene_to_schema_v2(scene)"
        )
        if (any(properties$attenuation != 0)) {
          warnings = c(
            warnings,
            "Legacy `dielectric(attenuation = ...)` absorption is deferred in spectral mode; conversion preserves IOR and priority only."
          )
        }
        region_id = sprintf(
          "legacy_dielectric_%s",
          scene$object_id[[i]] %||% i
        )
        attrs$regions[[region_id]] = dielectric_region(
          eta = spectrum_constant(properties$refraction),
          priority = properties$priority,
          id = region_id
        )
        scene$region_boundaries[[i]] = c(
          scene$region_boundaries[[i]],
          list(region_boundary(region_id, side = "negative_normal"))
        )
        scene$material[[i]] = dielectric_interface()
      } else {
        scene$material[[i]] = converted
      }
      entries[[i]] = list(
        object_id = scene$object_id[[i]],
        material_type = converted[[1]]$type
      )
    }
  }
  report = new_conversion_report(warnings = warnings, entries = entries)
  scene = restore_ray_scene_attrs(scene, attrs)
  attr(scene, "conversion_report") = report
  attr(scene, "validation") = validation
  emit_conversion_warnings(report, validation)
  scene
}

#' Serialize a schema-v2 descriptor
#'
#' @param x Descriptor or scene object.
#'
#' @return A raw vector.
#' @export
ray_schema_serialize = function(x) {
  serialize(x, NULL, version = 3)
}

#' Unserialize a schema-v2 descriptor
#'
#' @param x Raw serialized descriptor.
#'
#' @return The unserialized object.
#' @export
ray_schema_unserialize = function(x) {
  unserialize(x)
}

#' Round-trip a schema-v2 descriptor through serialization
#'
#' @param x Descriptor or scene object.
#'
#' @return The unserialized object.
#' @export
ray_schema_roundtrip = function(x) {
  ray_schema_unserialize(ray_schema_serialize(x))
}

scene_row_to_compiler_list = function(scene, i) {
  result = list()
  for (name in names(scene)) {
    value = scene[[name]]
    result[[name]] = if (is.list(value)) value[[i]] else value[[i]]
  }
  result
}

#' Convert a schema-v2 scene to compiler input
#'
#' @param scene Scene.
#' @param validation Default `attr(scene, "validation") %||% scene_validation()`. Validation descriptor.
#'
#' @return A schema-v2 compiler input list.
#' @export
as_scene_compiler_input = function(
  scene,
  validation = attr(scene, "validation") %||% scene_validation()
) {
  scene = ensure_ray_scene_v2(scene)
  validate_ray_scene_v2(scene, validation = validation)
  structure(
    list(
      schema_version = ray_schema_version,
      objects = lapply(seq_len(nrow(scene)), function(i) {
        scene_row_to_compiler_list(scene, i)
      }),
      regions = attr(scene, "regions") %||% list(),
      lights = attr(scene, "lights") %||% list(),
      environment = attr(scene, "environment"),
      named_textures = attr(scene, "named_textures") %||% list(),
      named_materials = attr(scene, "named_materials") %||% list(),
      render_defaults = attr(scene, "render_defaults") %||% list(),
      validation = validation
    ),
    class = "ray_scene_compiler_input"
  )
}
