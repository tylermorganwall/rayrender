#' Render Preview
#' 
#' Takes the scene description and renders an image, either to the device or to a filename. 
#'
#' @param ... All arguments that would be passed to `render_scene()`.
#' @param light_direction Default `c(0,-1,0)`. Vector specifying the orientation for the global light using for phong shading.
#' @param exponent Default `6`. Phong exponent.  
#' @export
#' @importFrom  grDevices col2rgb
#' @return Raytraced plot to current device, or an image saved to a file. 
#'
#' @examples
#' \donttest{
#' generate_ground(material=diffuse(color="darkgreen")) %>% 
#'   add_object(sphere(material=diffuse(checkercolor="red"))) %>% 
#'   render_preview()
#'   
#' #Change the light direction
#' generate_ground(material=diffuse(color="darkgreen")) %>% 
#'   add_object(sphere(material=diffuse(checkercolor="red"))) %>% 
#'   render_preview(light_direction = c(-1,-1,0))
#'   
#' #Change the Phong exponent
#' generate_ground(material=diffuse(color="darkgreen")) %>% 
#'   add_object(sphere(material=diffuse(checkercolor="red"))) %>% 
#'   render_preview(light_direction = c(-1,-1,0), exponent=100)
#' }
render_preview = function(..., light_direction = c(0,-1,0), exponent = 6) {
  stopifnot(length(light_direction) == 3 && is.numeric(light_direction))
  light_direction = light_direction/sqrt(sum(light_direction*light_direction))
  screen_colors = render_scene(..., debug_channel = c(-light_direction, exponent))
  invisible(screen_colors)
}