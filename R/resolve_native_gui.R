# Resolve the optional native provider. Keep gui out of prepare_scene_list's renderer data.
#' @keywords internal
#' @param gui One of auto, imgui, legacy or none.
#' @param preview Existing effective preview choice.
resolve_native_gui = function(gui, preview) {
  gui = match.arg(gui, c("auto", "imgui", "legacy", "none"))
  if (gui == "none" || (gui == "auto" && !isTRUE(preview))) {
    return(list(mode = "none", api = NULL, fallback = FALSE))
  }
  if (gui == "legacy") {
    return(list(mode = "legacy", api = NULL, fallback = FALSE))
  }
  unavailable = function(reason) {
    if (gui == "imgui") {
      stop("Native editor unavailable: ", reason, call. = FALSE)
    }
    message("Native editor unavailable (", reason, "); using existing preview.")
    list(mode = "legacy", api = NULL, fallback = FALSE)
  }
  if (!requireNamespace("rayimgui", quietly = TRUE)) {
    return(unavailable("provider_absent"))
  }
  api = tryCatch(
    rayimgui::acquire_api(major = 1L),
    rayimgui_unavailable = identity
  )
  if (inherits(api, "rayimgui_unavailable")) {
    return(unavailable(api$reason))
  }
  list(mode = "imgui", api = api, fallback = gui == "auto")
}
