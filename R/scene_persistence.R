#' @param scene Scene with optional saved defaults.
#' @param supplied Names of arguments explicitly supplied to the renderer.
#' @param frame Renderer evaluation environment containing argument promises.
#' @param renderer Rendering function.
#' @return NULL, or a merged argument list for one renderer re-entry.
#' @keywords internal
saved_scene_render_call = function(scene, supplied, frame, renderer) {
  settings = attr(scene, "ray_render_settings", exact = TRUE)
  pending = setdiff(
    intersect(names(settings), names(formals(renderer))),
    supplied
  )
  if (!length(pending)) {
    return(NULL)
  }
  # Re-enter with the merged arguments before validation or camera resolution.
  # This preserves R's missing() semantics; explicit NULL is still an override.
  c(mget(supplied, envir = frame, inherits = FALSE), settings[pending])
}
