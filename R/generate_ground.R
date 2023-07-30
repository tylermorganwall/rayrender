#' Generate Ground
#' 
#' Generates a large sphere that can be used as the ground for a scene.
#'
#' @param depth Default `-1`. Depth of the surface.
#' @param spheresize Default `1000`. Radius of the sphere representing the surface.
#' @param color Default `#ccff00`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param material Default  \code{\link{diffuse}} with `color= "#ccff00"`.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#'
#' @return Single row of a tibble describing the ground.
#' @export
#'
#' @examples
#' #Generate the ground and add some objects
#' scene = generate_ground(depth=-0.5,
#'                         material = diffuse(noise=1,noisecolor="blue",noisephase=10)) %>%
#'   add_object(cube(x=0.7,material=diffuse(color="red"),angle=c(0,-15,0))) %>%
#'   add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))
#' if(run_documentation()) {
#' render_scene(scene, parallel=TRUE,lookfrom=c(0,2,10))
#' }
#' 
#' # Make the sphere representing the ground larger and make it a checkered surface.
#' scene = generate_ground(depth=-0.5, spheresize=10000,
#'                         material = diffuse(checkercolor="grey50")) %>%
#'   add_object(cube(x=0.7,material=diffuse(color="red"),angle=c(0,-15,0))) %>%
#'   add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))
#' if(run_documentation()) {
#' render_scene(scene, parallel=TRUE,lookfrom=c(0,1,10))
#' }
generate_ground = function(depth = -1, spheresize = 1000, material = diffuse(color = "#ccff00")) {
  sphere(x = 0, y = -spheresize + depth, z = 0, radius = spheresize, material=material)
}
