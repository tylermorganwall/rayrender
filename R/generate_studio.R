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
#' scene = generate_studio(depth=-1, material = diffuse(color="white")) %>%
#'    add_object(obj_model(r_obj(),y=-0.5,x=0.5, scale=1.2,
#'                         material=glossy(color="darkred"),angle=c(0,-20,0))) %>%
#'    add_object(sphere(x=-0.5,radius=0.5,material=dielectric())) %>% 
#'    add_object(sphere(y=3,x=-2,z=20,material=light(intensity=600)))
#' if(run_documentation()) {
#' render_scene(scene, parallel = TRUE, lookfrom = c(0,2,10), lookat=c(0,-0.25,0),
#'              fov = 14, clamp_value = 10, samples = 16)
#' }
#' 
#' #Zooming out to show the full default scene
#' if(run_documentation()) {
#' render_scene(scene, parallel=TRUE,lookfrom=c(0,200,400),clamp_value=10,samples=16)
#' }
generate_studio = function(depth = -1, distance = -10, width = 100, height = 100,
                           curvature = 8, material = diffuse()) {
  xz_rect(y=depth, xwidth=width, zwidth=width,z = distance + width/2 + curvature, 
          material=material) %>%
    add_object(xy_rect(z=distance, y = depth + height/2 + curvature, xwidth = width,
                       ywidth= height,
                       material=material)) %>%
    add_object(cylinder(z=distance + curvature, y = depth + curvature, radius=curvature,
                        phi_min = 0, phi_max = 90,  flipped = T,
                        angle=c(90,90,0), capped = FALSE,
                        material=material, length = width))
  
}
