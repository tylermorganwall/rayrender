#' Generate Cornell Box
#'
#' @param light Default `TRUE`. Whether to include a light on the ceiling of the box.
#' @param lightintensity Default `5`. The intensity of the light.
#' @param lightcolor Default `white`. The color the of the light.
#' @param lightwidth Default `332`. Width (z) of the light.
#' @param lightdepth Default `343`. Depth (x) of the light.
#' @param light_position Default `c(555/2,554,555/2)`. Position of the light.
#' @param sigma Default `0`. Oren-Nayar microfacet angle.
#' @param leftcolor Default `#1f7326` (green).
#' @param rightcolor Default `#a60d0d` (red).
#' @param roomcolor Default `#bababa` (light grey).
#' @param importance_sample Default `TRUE`. Importance sample the light in the room.
#'
#' @return Tibble containing the scene description of the Cornell box.
#' @export
#'
#' @examples
#' #Generate and render the default Cornell box.
#' if(run_documentation()) {
#' render_scene(generate_cornell(), 
#'              samples=16,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' if(run_documentation()) {
#' #Make a much smaller light in the center of the room.
#' render_scene(generate_cornell(lightwidth=200,lightdepth=200), 
#'              samples=16,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' if(run_documentation()) {
#' #Place a sphere in the middle of the box.
#' scene = generate_cornell(lightwidth=200,lightdepth=200) %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/4))
#' render_scene(scene, samples=16,aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' if(run_documentation()) {
#' #Reduce "fireflies" by setting a clamp_value in render_scene()
#' render_scene(scene, samples=16,aperture=0, fov=40, ambient_light=FALSE, 
#'              parallel=TRUE,clamp_value=3)
#' }
#' if(run_documentation()) {
#' # Change the color scheme of the cornell box
#' generate_cornell(leftcolor="purple", rightcolor="yellow") |> 
#'   render_scene(samples=16,aperture=0, fov=40, ambient_light=FALSE, 
#'                parallel=TRUE,clamp_value=3)
#' }
generate_cornell = function(light = TRUE, lightintensity = 5, 
                            lightcolor = "white",lightwidth = 332, lightdepth=343, 
                            light_position = c(555/2,554,555/2),
                            sigma=0,
                            leftcolor = "#1f7326", rightcolor = "#a60d0d", roomcolor = "#bababa",
                            importance_sample = TRUE) {
  scene = cube(x=555+5,y=555/2,z=555/2,scale = c(10,555,555),
      material = diffuse(color = leftcolor, sigma = sigma)) %>%
    add_object(cube(x=0-5,y=555/2,z=555/2,scale = c(10,555,555),
      material = diffuse(color = rightcolor, sigma = sigma))) %>%
    add_object(cube(x=555/2,y=555+5,z=555/2,scale = c(575,10,555),
      material = diffuse(color=roomcolor, sigma = sigma),flipped=TRUE)) %>%
    add_object(cube(x=555/2,y=0-5,z=555/2,scale = c(575,10,555),
      material = diffuse(color=roomcolor, sigma = sigma))) %>%
    add_object(cube(x=555/2,y=555/2,z=555+5,scale = c(555,555,10),
      material = diffuse(color = roomcolor, sigma = sigma),flipped=TRUE))
  if(light) {
    scene = scene %>%
      add_object(xz_rect(x=light_position[1],y=light_position[2],z=light_position[3],
        lightdepth,lightwidth, flipped = TRUE,
                       material = light(color=lightcolor,intensity=lightintensity,
                                        importance_sample = importance_sample))) 
  }
  attr(scene,"cornell") = TRUE
  scene
}
