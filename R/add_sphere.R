#' Add Sphere
#'
#' @param scene Tibble of sphere locations and properties. Generate each row using one of the
#' material functions: [lambertian()], [metal()], and [dielectric()].
#' @param spheres A tibble row or collection of rows representing each sphere.
#'
#' @return Tibble of sphere locations and properties.
#' @export
#'
#' @examples
#' #Generate the ground and add a single metallic red sphere
#' scene = generate_ground() %>%
#'   add_sphere(metal(color = "#ff0000",fuzz=0))
#'   
#' render_scene(scene)
#' 
#' #add multiple glass spheres in front
#' scene = scene %>%
#'   add_sphere(dielectric(x=1.5,y=-0.5,z=0.75,radius=0.5)) %>%
#'   add_sphere(dielectric(x=1.5,y=-0.5,z=-0.75,radius=0.5)) %>%
#'   add_sphere(dielectric(x=1.5,y=0.75,z=0.75,radius=0.5)) %>%
#'   add_sphere(dielectric(x=1.5,y=0.75,z=-0.75,radius=0.5))
#'   
#' render_scene(scene)
add_sphere = function(scene, spheres) {
  rbind(scene,spheres)
}