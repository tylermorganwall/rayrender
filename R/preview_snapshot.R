#' Find the next available preview snapshot filename
#'
#' @param filename Default `NA_character_`. Source filename supplied to the renderer.
#'
#' @return An available filename for a preview snapshot.
#'
#' @keywords internal
next_preview_snapshot_filename = function(filename = NA_character_) {
  has_source_filename = length(filename) == 1 &&
    !is.na(filename) &&
    nzchar(filename) &&
    nzchar(tools::file_ext(filename))

  if (has_source_filename) {
    snapshot_base = tools::file_path_sans_ext(filename)
    snapshot_extension = tools::file_ext(filename)
  } else {
    snapshot_base = "rayrender_snapshot"
    snapshot_extension = "png"
  }

  snapshot_number = 1L
  repeat {
    candidate = paste0(
      snapshot_base,
      snapshot_number,
      ".",
      snapshot_extension
    )
    if (!file.exists(candidate)) {
      return(candidate)
    }
    snapshot_number = snapshot_number + 1L
  }
}

#' Save an interactive preview snapshot
#'
#' @param image Preview image array.
#' @param filename Default `NA_character_`. Source filename supplied to the renderer.
#'
#' @return Invisibly returns the saved snapshot filename.
#'
#' @keywords internal
save_preview_snapshot = function(image, filename = NA_character_) {
  snapshot_filename = next_preview_snapshot_filename(filename)
  # The native preview supplies display-encoded RGB pixels. Identify their
  # primaries and decode once, rather than treating an untagged array as ACEScg.
  image = rayimage::ray_read_image(
    image,
    source_linear = FALSE,
    assume_colorspace = rayimage::CS_SRGB
  )
  rayimage::ray_write_image(
    image,
    snapshot_filename
  )
  message("Saved preview snapshot: ", snapshot_filename)
  invisible(snapshot_filename)
}
