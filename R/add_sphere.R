#' Add Object
#'
#' @param scene Tibble of object locations and properties. Generate each row using one of the
#' material functions: [lambertian()], [metal()], and [dielectric()].
#' @param spheres A tibble row or collection of rows representing each object.
#'
#' @return Tibble of object locations and properties.
#' @export
#'
#' @examples
#' #Generate the ground and add a single metallic red sphere
#' scene = generate_ground() %>%
#'   add_object(metal(color = "#ff0000",fuzz=0))
#'   
#' render_scene(scene)
#' 
#' #add multiple glass spheres in front
#' scene = scene %>%
#'   add_object(dielectric(x=1.5,y=-0.5,z=0.75,radius=0.5)) %>%
#'   add_object(dielectric(x=1.5,y=-0.5,z=-0.75,radius=0.5)) %>%
#'   add_object(dielectric(x=1.5,y=0.75,z=0.75,radius=0.5)) %>%
#'   add_object(dielectric(x=1.5,y=0.75,z=-0.75,radius=0.5))
#'   
#' render_scene(scene)
add_object = function(scene, objects) {
  rbind(scene,objects)
}