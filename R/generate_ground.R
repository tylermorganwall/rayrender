#' Generate Ground
#' 
#' Generates a large sphere that can be used as the ground for a scene.
#'
#' @param depth Default `-1`. Depth of the surface.
#' @param spheresize Default `10000`. Size of the sphere representing the surface.
#' @param color Default `#ccff00`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param material Default  \code{\link{lambertian(color= "#ccff00"}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#'
#' @return Single row of a tibble describing the ground.
#' @export
#'
#' @examples
#' #Generate the ground and add some objects
#' scene = generate_ground(depth=-0.5,material = lambertian(noise=1,noisephase=10)) %>%
#'   add_object(cube(z=0.7,material=metal(color="red",fuzz=0.025),angle=c(0,-15,0))) %>%
#'   add_object(sphere(z=-0.7,radius=0.5,material=dielectric(color="green")))
#' render_scene(scene)
generate_ground = function(depth = -1, spheresize = 10000, material = lambertian(color = "#ccff00")) {
  sphere(x = 0,y = -spheresize + depth, z = 0, radius = spheresize, material=material)
}