#' Prepare an isolated native editor scene
#'
#' @param scene A processed scene.
#' @keywords internal
native_editor_scene = function(scene) {
  groups = vector("list", 0L)
  scene$preview_groups = integer(nrow(scene))
  for (i in seq_len(nrow(scene))) {
    marker = attr(scene$transforms[[i]], "rayrender_preview_group")
    if (!is.null(marker)) {
      same = which(vapply(groups, identical, logical(1), marker))
      if (!length(same)) {
        groups[[length(groups) + 1L]] = marker
        same = length(groups)
      }
      scene$preview_groups[i] = same[1L]
    }
    # Editing a selected object must not modify another object's shared material.
    scene$shape_info[[i]]$material_id = NA_integer_
    if (scene$shape[i] == 15L) {
      properties = scene$shape_info[[i]]$shape_properties
      properties$original_scene[[
        1L
      ]] = native_editor_scene(properties$original_scene[[1L]])
      scene$shape_info[[i]]$shape_properties = properties
    }
  }
  scene
}
