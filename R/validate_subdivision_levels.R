#' @param subdivision_levels Number of subdivision passes requested.
#' @return A validated nonnegative integer count stored as a double.
#' @keywords internal
#' @noRd
validate_subdivision_levels = function(subdivision_levels) {
  if (
    !is.numeric(subdivision_levels) ||
      length(subdivision_levels) != 1L ||
      !is.finite(subdivision_levels) ||
      subdivision_levels < 0 ||
      subdivision_levels > .Machine$integer.max ||
      subdivision_levels != floor(subdivision_levels)
  ) {
    stop(
      'subdivision_levels must be a single nonnegative integer.',
      call. = FALSE
    )
  }
  # Canonical storage also lets PBRT batching recognize an explicit integer 0L.
  as.double(subdivision_levels)
}
