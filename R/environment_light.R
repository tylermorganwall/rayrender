#' Prepare environment light white balance
#'
#' @param environment_light Environment light filename.
#' @param environment_light_bake_white Default `FALSE`. Whether to bake the
#' environment light white point into a temporary copy.
#' @param environment_light_bake_white_target Default `"D65"`. Target white
#' point.
#' @return A list with `environment_light` and `cleanup` entries.
#' @keywords internal
prepare_environment_light_white_balance = function(
  environment_light,
  environment_light_bake_white = FALSE,
  environment_light_bake_white_target = "D65"
) {
  if (
    !is.logical(environment_light_bake_white) ||
      length(environment_light_bake_white) != 1L ||
      is.na(environment_light_bake_white)
  ) {
    stop("`environment_light_bake_white` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!environment_light_bake_white || is.null(environment_light)) {
    return(list(environment_light = environment_light, cleanup = character()))
  }

  if (
    !is.character(environment_light) ||
      length(environment_light) != 1L ||
      is.na(environment_light)
  ) {
    stop(
      "`environment_light_bake_white` requires `environment_light` to be a single filename.",
      call. = FALSE
    )
  }

  environment_light_path = path.expand(environment_light)
  if (tolower(tools::file_ext(environment_light_path)) != "exr") {
    warning(
      "`environment_light_bake_white` currently only supports EXR environment lights; using the original environment light.",
      call. = FALSE
    )
    return(list(environment_light = environment_light, cleanup = character()))
  }

  if (
    !file.exists(environment_light_path) || dir.exists(environment_light_path)
  ) {
    return(list(environment_light = environment_light, cleanup = character()))
  }

  if (!requireNamespace("rayimage", quietly = TRUE)) {
    stop(
      "`environment_light_bake_white` requires the 'rayimage' package.",
      call. = FALSE
    )
  }

  environment_image = rayimage::ray_read_image(
    environment_light_path,
    normalize = FALSE
  )
  reference_white = attr(environment_image, "white_current", exact = TRUE)
  if (is.null(reference_white)) {
    warning(
      "`environment_light_bake_white` could not find `white_current` on the EXR; using the original environment light.",
      call. = FALSE
    )
    return(list(environment_light = environment_light, cleanup = character()))
  }

  target_white = environment_light_white_xyz(
    environment_light_bake_white_target
  )
  if (isTRUE(all.equal(reference_white, target_white, tolerance = 1e-8))) {
    return(list(environment_light = environment_light, cleanup = character()))
  }

  baked_environment = rayimage::render_white_balance(
    environment_image,
    reference_white = reference_white,
    target_white = target_white,
    bake = TRUE
  )
  baked_environment = clamp_negative_environment_rgb(baked_environment)
  attr(baked_environment, "white_current") = target_white

  baked_environment_path = tempfile(
    pattern = "rayrender-environment-light-white-",
    fileext = ".exr"
  )
  rayimage::ray_write_image(
    baked_environment,
    baked_environment_path,
    clamp = FALSE,
    write_linear = TRUE
  )

  list(
    environment_light = baked_environment_path,
    cleanup = baked_environment_path
  )
}

#' Resolve environment light white point
#'
#' @param white_point Named white point or XYZ vector.
#' @return Numeric XYZ white point with Y = 1.
#' @keywords internal
environment_light_white_xyz = function(white_point) {
  white_points = list(
    D65 = c(0.95047, 1, 1.08883),
    D60 = c(0.95264, 1, 1.00827),
    D55 = c(0.95560, 1, 0.92149),
    D50 = c(0.96422, 1, 0.82521),
    D75 = c(0.94972, 1, 1.22638),
    E = c(1, 1, 1)
  )

  if (is.character(white_point) && length(white_point) == 1L) {
    resolved_white = white_points[[toupper(white_point)]]
    if (is.null(resolved_white)) {
      stop(
        "Unknown `environment_light_bake_white_target`: ",
        white_point,
        call. = FALSE
      )
    }
    return(resolved_white)
  }

  if (
    !is.numeric(white_point) ||
      length(white_point) != 3L ||
      any(!is.finite(white_point))
  ) {
    stop(
      "`environment_light_bake_white_target` must be a named white point or a length-3 numeric XYZ vector.",
      call. = FALSE
    )
  }

  if (white_point[2] <= 0) {
    stop(
      "`environment_light_bake_white_target` must have a positive Y value.",
      call. = FALSE
    )
  }

  white_point / white_point[2]
}

#' Clamp negative environment RGB values
#'
#' @param image Environment image.
#' @return Environment image with non-negative RGB channels.
#' @keywords internal
clamp_negative_environment_rgb = function(image) {
  image_dims = dim(image)
  if (length(image_dims) != 3L || image_dims[3] < 3L) {
    return(image)
  }

  rgb_channels = seq_len(3L)
  image[,, rgb_channels] = pmax(image[,, rgb_channels], 0)
  image
}
