#' Set Saved Render Settings
#'
#' Update a scene's default rendering arguments, including sample count and image
#' dimensions. The defaults are used by both renderers and retained in USD exports.
#' @param scene A rayrender scene, including one loaded with [read_usd()].
#' @param ... Named arguments accepted by [render_scene()] or [render_animation()],
#'   such as `samples`, `width`, and `height`. Only the supplied settings change.
#'   An explicit `NULL` is retained as a saved value.
#' @return A rayrender scene with updated render defaults. Assign the result or
#'   pass it onward with a pipe to use the changes.
#' @export
#' @md
#' @details
#' Explicit arguments passed to a renderer override the saved defaults. Settings
#' supported by only one renderer are ignored by the other. Names must match
#' renderer arguments exactly; argument values are validated when rendering.
#' Changing these settings does not rewrite the source USD file. Use [write_usd()]
#' to save the modified scene.
#' @examples
#' # Save image quality and dimensions with the scene.
#' scene = sphere() |>
#'   set_render_settings(samples = 512, width = 1920, height = 1080)
#'
#' # Change sample count while keeping the saved dimensions.
#' scene = set_render_settings(scene, samples = 128)
set_render_settings = function(scene, ...) {
  if (!inherits(scene, "ray_scene")) {
    stop("scene must be a ray_scene.", call. = FALSE)
  }
  updates = list(...)
  if (!length(updates)) {
    return(scene)
  }
  setting_names = names(updates)
  if (
    is.null(setting_names) ||
      anyNA(setting_names) ||
      any(!nzchar(setting_names))
  ) {
    stop("All render settings must be named.", call. = FALSE)
  }
  if (anyDuplicated(setting_names)) {
    stop("Render setting names must be unique.", call. = FALSE)
  }
  supported = setdiff(
    union(names(formals(render_scene)), names(formals(render_animation))),
    c("scene", "...")
  )
  unknown = setdiff(setting_names, supported)
  if (length(unknown)) {
    stop(
      "Unknown render settings: ",
      paste(unknown, collapse = ", "),
      call. = FALSE
    )
  }
  settings = attr(scene, "ray_render_settings", exact = TRUE)
  if (is.null(settings)) {
    settings = list()
  }
  # List subassignment retains NULL entries so they can override saved defaults.
  settings[setting_names] = updates
  attr(scene, "ray_render_settings") = settings
  scene
}


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
