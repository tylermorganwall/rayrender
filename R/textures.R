#' Composable Surface Textures
#'
#' Build scalar and color texture graphs for `diffuse(color)`,
#' `microfacet(color, roughness)`, and `openpbr(base_color, specular_roughness)`.
#' Constants are accepted wherever a child texture is expected. Scalar values
#' broadcast to colors; use `texture_channel()` to convert a color to a scalar.
#' Graphs are evaluated in C++ at surface hits, without callbacks to R.
#'
#' @param value A finite scalar, linear RGB triple, R color, or texture.
#' @param space Default `"object"`. Coordinate space: `"object"`, `"world"`, or `"uv"`.
#' @param scale Default `1`. Scalar or three-component coordinate multiplier;
#' UV mappings also accept two components.
#' @param offset Default `c(0, 0, 0)`. Coordinate translation with three components,
#' or two for UV mappings.
#' @param rotation Default `0`. Rotation in degrees about the texture z axis,
#' applied after scaling and before translation (UV rotation is in the UV plane).
#' @param coordinates Default `texture_coordinates()`. Independent texture mapping.
#' @param a First value or texture.
#' @param b Second value or texture.
#' @param weight Default `0.5`. Scalar blend weight, clamped to [0, 1].
#' @param direction Default `c(0, 1, 0)`. Nonzero direction in the selected space.
#' @param absolute Default `TRUE`. Use the absolute geometric-normal dot product.
#' If false, use the positive part, so opposite-facing surfaces receive zero weight.
#' @param octaves Default `4L`. Number of noise octaves, between one and sixteen.
#' @param seed Default `0L`. Integer seed; independent of the renderer's sampling seed.
#' @param filename Image filename.
#' @param type Default `"color"`. Output type: `"color"` or `"scalar"`.
#' @param encoding Default `"auto"`. `"auto"`, `"linear"`, or `"srgb"`.
#' Auto treats scalar images and HDR/EXR images as linear, other color images as sRGB.
#' @param wrap Default `"repeat"`. Image addressing: `"repeat"` or `"clamp"`.
#' @param channel Default `"luminance"`. Color reduction: `"r"`, `"g"`, `"b"`,
#' `"average"`, or linear RGB `"luminance"`.
#' @param factor Scalar multiplier or scalar texture.
#' @param axis Default `"y"`. Coordinate axis for a clamped zero-to-one gradient.
#' @details
#' `texture_mix(a, b, weight)` returns `(1-weight)*a + weight*b`.
#' `texture_direction_mix(a, b)` follows PBRT: the absolute dot product of the
#' geometric normal and normalized direction weights **a**, not b. It defaults
#' to world space. Bump mapping does not change this orientation weight.
#'
#' Numeric colors are linear RGB, following rayrender's material convention.
#' Image color decoding happens before filtering and mixing. Scalar images use
#' the average of RGB channels. Image filtering is bilinear; procedural textures
#' currently evaluate a point sample, without ray-footprint antialiasing.
#' Object coordinates refer to the primitive's original coordinates, including
#' inside instances; world coordinates include the instance placement.
#' Noise is seeded smooth value noise with normalized, decreasing octave weights.
#' Roughness textures are clamped to [0, 1] at the material input. The microfacet
#' material uses the same conversion as its numeric roughness argument.
#' Spatial microfacet roughness uses a minimum alpha of `1e-6`; a constant zero
#' roughness uses the existing perfectly smooth material. Transmitting microfacet
#' color graphs require nonzero roughness.
#' Graph inputs cannot be combined with the legacy image/procedural settings for
#' the same input. Bump, alpha, emission, and volume inputs do not yet accept graphs.
#' @return A serializable `ray_texture` or `ray_texture_coordinates` descriptor.
#' @name textures
#' @export
texture_constant = function(value) {
  if (inherits(value, "ray_texture")) {
    return(value)
  }
  if (is.character(value)) {
    value = convert_color(value)
  }
  if (
    !is.numeric(value) ||
      !length(value) %in% c(1L, 3L) ||
      !is.null(dim(value)) ||
      any(!is.finite(value))
  ) {
    stop(
      "A texture constant must be a finite scalar or RGB triple.",
      call. = FALSE
    )
  }
  texture_node(
    "constant",
    if (length(value) == 1L) "scalar" else "color",
    value = as.numeric(value)
  )
}

#' @rdname textures
#' @export
texture_coordinates = function(
  space = "object",
  scale = 1,
  offset = c(0, 0, 0),
  rotation = 0
) {
  space = match.arg(space, c("object", "world", "uv"))
  if (space == "uv" && is.numeric(scale) && length(scale) == 2L) {
    scale = c(scale, 1)
  }
  if (space == "uv" && is.numeric(offset) && length(offset) == 2L) {
    offset = c(offset, 0)
  }
  if (
    !is.numeric(scale) ||
      !length(scale) %in% c(1L, 3L) ||
      any(!is.finite(scale)) ||
      !is.numeric(offset) ||
      length(offset) != 3L ||
      any(!is.finite(offset)) ||
      !is.numeric(rotation) ||
      length(rotation) != 1L ||
      !is.finite(rotation)
  ) {
    stop(
      "Invalid texture coordinate scale, offset, or rotation.",
      call. = FALSE
    )
  }
  structure(
    list(
      space = space,
      scale = rep(scale, length.out = 3L),
      offset = offset,
      rotation = rotation * pi / 180
    ),
    class = "ray_texture_coordinates"
  )
}

#' @rdname textures
#' @export
texture_mix = function(a, b, weight = 0.5) {
  a = texture_constant(a)
  b = texture_constant(b)
  weight = texture_scalar(weight)
  type = if (a$type == "color" || b$type == "color") "color" else "scalar"
  if (weight$op == "constant") {
    w = pmin(1, pmax(0, weight$value))
    if (a$op == "constant" && b$op == "constant") {
      return(texture_constant((1 - w) * a$value + w * b$value))
    }
  }
  texture_node("mix", type, a = a, b = b, weight = weight)
}

#' @rdname textures
#' @export
texture_direction_mix = function(
  a,
  b,
  direction = c(0, 1, 0),
  absolute = TRUE,
  space = "world"
) {
  space = match.arg(space, c("world", "object"))
  if (
    !is.numeric(direction) ||
      length(direction) != 3L ||
      any(!is.finite(direction)) ||
      !any(direction != 0)
  ) {
    stop("`direction` must be a finite nonzero three-vector.", call. = FALSE)
  }
  if (!is.logical(absolute) || length(absolute) != 1L || is.na(absolute)) {
    stop("`absolute` must be TRUE or FALSE.", call. = FALSE)
  }
  direction = direction / max(abs(direction))
  direction = direction / sqrt(sum(direction^2))
  weight = texture_node(
    "direction",
    "scalar",
    direction = direction,
    absolute = absolute,
    space = space
  )
  texture_mix(b, a, weight)
}

#' @rdname textures
#' @export
texture_noise = function(
  scale = 1,
  coordinates = texture_coordinates(),
  octaves = 4L,
  seed = 0L
) {
  coordinates = texture_mapping(coordinates)
  if (
    !is.numeric(scale) ||
      length(scale) != 1L ||
      !is.finite(scale) ||
      scale <= 0 ||
      !is.numeric(octaves) ||
      length(octaves) != 1L ||
      !is.finite(octaves) ||
      octaves < 1 ||
      octaves > 16 ||
      octaves != floor(octaves) ||
      !is.numeric(seed) ||
      length(seed) != 1L ||
      !is.finite(seed) ||
      abs(seed) > .Machine$integer.max ||
      seed != floor(seed)
  ) {
    stop("Invalid noise scale, octave count, or integer seed.", call. = FALSE)
  }
  coordinates$scale = coordinates$scale * scale
  texture_node(
    "noise",
    "scalar",
    coordinates = coordinates,
    octaves = as.integer(octaves),
    seed = as.integer(seed)
  )
}

#' @rdname textures
#' @export
texture_checker = function(a, b, coordinates = texture_coordinates()) {
  a = texture_constant(a)
  b = texture_constant(b)
  texture_node(
    "checker",
    if (a$type == "color" || b$type == "color") "color" else "scalar",
    a = a,
    b = b,
    coordinates = texture_mapping(coordinates)
  )
}

#' @rdname textures
#' @export
texture_gradient = function(
  a,
  b,
  coordinates = texture_coordinates(),
  axis = "y"
) {
  axis = match.arg(axis, c("x", "y", "z"))
  weight = texture_node(
    "gradient",
    "scalar",
    coordinates = texture_mapping(coordinates),
    axis = match(axis, c("x", "y", "z")) - 1L
  )
  texture_mix(a, b, weight)
}

#' @rdname textures
#' @export
texture_image = function(
  filename,
  type = "color",
  encoding = "auto",
  coordinates = texture_coordinates("uv"),
  wrap = "repeat"
) {
  type = match.arg(type, c("color", "scalar"))
  encoding = match.arg(encoding, c("auto", "linear", "srgb"))
  wrap = match.arg(wrap, c("repeat", "clamp"))
  if (
    !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !file.exists(path.expand(filename))
  ) {
    stop("`filename` must name an existing texture image.", call. = FALSE)
  }
  if (encoding == "auto") {
    encoding = if (
      type == "scalar" ||
        tolower(tools::file_ext(filename)) %in% c("exr", "hdr")
    ) {
      "linear"
    } else {
      "srgb"
    }
  }
  texture_node(
    "image",
    type,
    filename = normalizePath(path.expand(filename)),
    encoding = encoding,
    coordinates = texture_mapping(coordinates),
    wrap = wrap
  )
}

#' @rdname textures
#' @export
texture_channel = function(value, channel = "luminance") {
  value = texture_constant(value)
  channel = match.arg(channel, c("r", "g", "b", "average", "luminance"))
  if (value$type == "scalar") {
    return(value)
  }
  texture_node(
    "channel",
    "scalar",
    child = value,
    channel = match(channel, c("r", "g", "b", "average", "luminance")) - 1L
  )
}

#' @rdname textures
#' @export
texture_scale = function(value, factor) {
  value = texture_constant(value)
  factor = texture_scalar(factor)
  if (value$op == "constant" && factor$op == "constant") {
    return(texture_constant(value$value * factor$value))
  }
  texture_node("scale", value$type, child = value, factor = factor)
}

#' @param op Operation name.
#' @param type Output type.
#' @param ... Node fields.
#' @return Typed texture descriptor.
#' @keywords internal
#' @noRd
texture_node = function(op, type, ...) {
  structure(c(list(op = op, type = type), list(...)), class = "ray_texture")
}

#' @param value Scalar value or texture.
#' @return Scalar texture, rejecting implicit color reduction.
#' @keywords internal
#' @noRd
texture_scalar = function(value) {
  value = texture_constant(value)
  if (value$type != "scalar") {
    stop(
      "Expected a scalar texture; use texture_channel() to reduce colors.",
      call. = FALSE
    )
  }
  value
}

#' @param coordinates Coordinate descriptor.
#' @return Validated coordinate descriptor.
#' @keywords internal
#' @noRd
texture_mapping = function(coordinates) {
  if (!inherits(coordinates, "ray_texture_coordinates")) {
    stop("Expected texture_coordinates().", call. = FALSE)
  }
  coordinates
}

#' @param image Legacy image input.
#' @param checker Legacy checker color.
#' @param noise Legacy noise scale.
#' @param gradient Legacy gradient color.
#' @return Nothing; rejects ambiguous color input combinations.
#' @keywords internal
#' @noRd
texture_check_legacy = function(image, checker, noise, gradient) {
  if (
    !identical(image, "") ||
      any(!is.na(checker)) ||
      noise != 0 ||
      any(!is.na(gradient))
  ) {
    stop(
      "Choose a color texture graph or legacy image/checker/noise/gradient arguments, not both.",
      call. = FALSE
    )
  }
}
