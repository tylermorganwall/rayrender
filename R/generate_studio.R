#' Generate Studio
#' 
#' Generates a curved studio backdrop.
#'
#' @param depth Default `-1`. Depth of the ground in the scene.
#' @param distance Default `-10`. Distance to the backdrop in the scene from the origin, on the z-axis.
#' @param width Default `100`. Width of the backdrop.
#' @param height Default `100`. height of the backdrop.
#' @param curvature Default `2`. Radius of the curvature connecting the bottom plane to the vertical
#' backdrop.
#' @param material Default  \code{\link{diffuse}} with `color= "#ccff00"`.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#'
#' @return Tibble representing the scene.
#' @export
#'
#' @examples
#' #Generate the ground and add some objects
#' scene = generate_ground(depth=-0.5,
#'                         material = diffuse(noise=1,noisecolor="blue",noisephase=10)) %>%
#'   add_object(cube(x=0.7,material=diffuse(color="red"),angle=c(0,-15,0))) %>%
#'   add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))
#' \donttest{
#' render_scene(scene, parallel=TRUE,lookfrom=c(0,2,10))
#' }
#' 
#' # Make the sphere representing the ground larger and make it a checkered surface.
#' scene = generate_ground(depth=-0.5, spheresize=10000,
#'                         material = diffuse(checkercolor="grey50")) %>%
#'   add_object(cube(x=0.7,material=diffuse(color="red"),angle=c(0,-15,0))) %>%
#'   add_object(sphere(x=-0.7,radius=0.5,material=dielectric(color="white")))
#' \donttest{
#' render_scene(scene, parallel=TRUE,lookfrom=c(0,1,10))
#' }
generate_studio = function(depth = -1, distance = -10, width = 100, height = 100,
                           curvature = 8, material = diffuse()) {
  xz_rect(y=depth, xwidth=width, zwidth=width,z = distance + width/2 + curvature, 
          material=material) %>%
    add_object(xy_rect(z=distance, y = depth + height/2 + curvature, xwidth = width,
                       ywidth= height,
                       material=material)) %>%
    add_object(cylinder(z=distance + curvature, y = depth + curvature, radius=curvature,
                        phi_min = 270, phi_max = 360,
                        angle=c(90,90,0),
                        material=material, length = width))
  
}
