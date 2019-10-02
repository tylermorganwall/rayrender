#' Generate Ground
#' 
#' Generates a large sphere that can be used as the ground for a scene.
#'
#' @param depth Default `-1`. Depth of the surface.
#' @param spheresize Default `1000`. Radius of the sphere representing the surface.
#' @param color Default `#ccff00`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param material Default  \code{\link{lambertian}} with `color= "#ccff00"`.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#'
#' @return Single row of a tibble describing the ground.
#' @export
#'
#' @examples
#' #Generate the ground and add some objects
#' scene = generate_ground(depth=-0.5,
#'                         material = lambertian(noise=1,noisecolor="blue",noisephase=10)) %>%
#'   add_object(cube(x=0.7,material=lambertian(color="red"),angle=c(0,-15,0))) %>%
#'   add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))
#' \donttest{
#' render_scene(scene, parallel=TRUE,lookfrom=c(0,2,10))
#' }
#' 
#' # Make the sphere representing the ground larger and make it a checkered surface.
#' scene = generate_ground(depth=-0.5, spheresize=10000,
#'                         material = lambertian(checkercolor="grey50")) %>%
#'   add_object(cube(x=0.7,material=lambertian(color="red"),angle=c(0,-15,0))) %>%
#'   add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))
#' \donttest{
#' render_scene(scene, parallel=TRUE,lookfrom=c(0,1,10))
#' }
generate_ground = function(depth = -1, spheresize = 1000, material = lambertian(color = "#ccff00")) {
  sphere(x = 0, y = -spheresize + depth, z = 0, radius = spheresize, material=material)
}
