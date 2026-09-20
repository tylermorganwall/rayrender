#' Prepare an isolated native editor scene
#'
#' @param scene A processed scene.
#' @keywords internal
native_editor_scene = function(scene) {
  groups = list()
  scene$preview_groups = integer(nrow(scene))
  scene$preview_paths = vector("list", nrow(scene))
  group_id = function(marker) {
    same = which(vapply(groups, identical, logical(1), marker))
    if (!length(same)) {
      groups[[length(groups) + 1L]] <<- marker
      return(length(groups))
    }
    same[1L]
  }
  for (i in seq_len(nrow(scene))) {
    path = attr(scene$transforms[[i]], "rayrender_preview_groups")
    marker = attr(scene$transforms[[i]], "rayrender_preview_group")
    if (is.null(path) && !is.null(marker)) {
      path = list(marker)
    }
    scene$preview_paths[[i]] = vapply(path, group_id, integer(1))
    if (length(path)) {
      scene$preview_groups[i] = scene$preview_paths[[i]][1L]
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
