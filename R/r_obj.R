#' R 3D Model
#' 
#' 3D obj model of R logo (created from the R SVG logo with the `raybevel` package), 
#' to be used with `obj_model()`
#' 
#' @param simple_r Default `FALSE`. If `TRUE`, this will return a 3D R (instead of the R logo).
#' @return File location of the 3d_r_logo.obj file (saved with a .txt extension)
#' @export
#'
#' @examples
#' #Load and render the included example R object file.
#' if(run_documentation()) {
#' generate_ground(material = diffuse(noise = TRUE, noisecolor = "grey20")) %>%
#'   add_object(sphere(x = 2, y = 3, z = 2, radius = 1,
#'                     material = light(intensity = 10))) %>%
#'   add_object(obj_model(r_obj(), y = -1, material = diffuse(color="red"))) %>%
#'   render_scene(parallel=TRUE, lookfrom = c(0, 1, 10), clamp_value = 5, samples = 200)
#' }
r_obj = function(simple_r = FALSE) {
  if(!simple_r) {
    system.file("extdata", "3d_r_logo.txt", package="rayrender")
  } else {
    system.file("extdata", "r_obj.txt", package="rayrender")
  }
}
