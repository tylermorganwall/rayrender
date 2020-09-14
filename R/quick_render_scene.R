#' Render Preview
#' 
#' Takes the scene description and renders an image, either to the device or to a filename. 
#'
#' @param ... All arguments that would be passed to `render_scene()`.
#' @export
#' @importFrom  grDevices col2rgb
#' @return Raytraced plot to current device, or an image saved to a file. 
#'
#' @examples
render_preview = function(..., light_direction = c(0,-1,0), exponent = 6) {
  stopifnot(length(light_direction) == 3 && is.numeric(light_direction))
  light_direction = light_direction/sqrt(sum(light_direction*light_direction))
  screen_colors = render_scene(..., debug_channel = c(-light_direction, exponent))
  invisible(screen_colors)
}