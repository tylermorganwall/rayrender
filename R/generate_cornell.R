#' Generate Cornell Box
#'
#' @param light Default `TRUE`. Whether to include a light on the ceiling of the box.
#' @param lightintensity Default `5`. The intensity of the light.
#' @param lightcolor Default `white`. The color the of the light.
#' @param lightwidth Default `332`. Width (z) of the light.
#' @param lightdepth Default `343`. Depth (x) of the light.
#'
#' @return Tibble containing the scene description of the Cornell box.
#' @export
#'
#' @examples
#' #Generate and render the default Cornell box.
#' scene = generate_cornell()
#' render_scene(scene, samples=200,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' 
#' #Make a much smaller light in the center of the room.
#' scene = generate_cornell(lightwidth=200,lightdepth=200)
#' render_scene(scene, samples=200,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' 
#' #Place a sphere in the middle of the box.
#' scene = scene %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/4))
#' render_scene(scene, samples=200,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
generate_cornell = function(light = TRUE, lightintensity = 5,lightcolor = "white",lightwidth = 332, lightdepth=343) {
  scene = yz_rect(x=555,y=555/2,z=555/2,555,555,
      material = lambertian(color = "#1f7326"),flipped = TRUE) %>%
    add_object(yz_rect(x=0,y=555/2,z=555/2,555,555,
      material = lambertian(color = "#a60d0d"))) %>%
    add_object(xz_rect(x=555/2,y=555,z=555/2,555,555,
      material = lambertian(color="#bababa"),flipped = TRUE)) %>%
    add_object(xz_rect(x=555/2,y=0,z=555/2,555,555,
      material = lambertian(color="#bababa"))) %>%
    add_object(xy_rect(x=555/2,y=555/2,z=555,555,555,
      material = lambertian(color = "#bababa"),flipped = TRUE))
  if(light) {
    scene = scene %>%
      add_object(xz_rect(x=555/2,y=554,z=555/2,lightdepth,lightwidth,
                       material = lambertian(color=lightcolor,
                                             lightintensity=lightintensity,implicit_sample = TRUE),
                       flipped=TRUE)) 
  }
  attr(scene,"cornell") = TRUE
  scene
}