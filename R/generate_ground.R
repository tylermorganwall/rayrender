#' Title
#'
#' @param depth Default `-1`. Depth of the surface.
#' @param spheresize Default `1000`. Size of the sphere representing the surface.
#' @param color Default `#ccff00`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' #Generate the ground and add a single red sphere
#' scene = generate_ground() %>%
#'   add_sphere(lambertian(color = "#ff0000"))
#'   
#' render_scene(scene)
generate_ground = function(depth = -1, spheresize = 1000, color = "#ccff00", checkercolor=NA, noise=0) {
  lambertian(x = 0,y = -spheresize + depth, z = 0, radius = spheresize, color = color, checkercolor=checkercolor,noise=noise)
}