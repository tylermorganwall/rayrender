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
#' scene = generate_ground(depth=-0.5,material = diffuse(checkercolor="blue")) %>%
#'   add_object(cube(x=0.7,
#'                   material=diffuse(noise=5,noisecolor="purple",color="black",noisephase=45),
#'                   angle=c(0,-30,0))) %>%
#'   add_object(sphere(x=-0.7,radius=0.5,material=metal(color="gold")))
#' \donttest{
#' render_scene(scene,parallel=TRUE)
#' }
add_object = function(scene, objects) {
  newscene = rbind(scene,objects)
  if(!is.null(attr(objects,"cornell")) || !is.null(attr(scene,"cornell"))) {
    attr(newscene,"cornell") = TRUE
  }
  return(newscene)
}
