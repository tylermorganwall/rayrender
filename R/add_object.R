#' Add Object
#'
#' @param scene Tibble of pre-existing object locations and properties.
#' @param objects A tibble row or collection of rows representing each object.
#'
#' @return Tibble of object locations and properties.
#' @export
#'
#' @examples
#' #Generate the ground and add some objects
#' scene = generate_ground(depth=-0.5,material = lambertian(checkercolor="blue")) %>%
#'   add_object(cube(z=0.7,
#'                   material=lambertian(noise=5,noisecolor="purple",color="black",noisephase=45),
#'                   angle=c(0,30,0))) %>%
#'   add_object(sphere(z=-0.7,radius=0.5,material=metal(color="gold")))
#' \dontrun{
#' render_scene(scene,parallel=TRUE)
#' }
add_object = function(scene, objects) {
  rbind(scene,objects)
}