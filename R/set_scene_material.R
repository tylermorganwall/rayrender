#' Set Material for All Objects
#'
#' @param scene A ray_scene object.
#' @param material A material specification created by diffuse(), metal(), dielectric(), etc.
#'
#' @return A modified ray_scene with the new material applied to all objects
#' @export
#'
#' @examples
#' # Create a scene with different materials
#' scene = generate_cornell() %>%
#'   add_object(sphere(x=555/2, y=555/2, z=555/2, radius=100))
#'   
#' # Set all objects to be metallic
#' scene = set_scene_material(scene, metal(color="gold"))
#'
#' # Set all objects to be glass
#' scene = set_scene_material(scene, dielectric())
set_scene_material = function(scene, material) {
  if (!inherits(scene, "ray_scene")) {
    stop("Input must be a ray_scene object")
  }
  if (!inherits(material, "ray_material")) {
    stop("Material must be created using material functions (diffuse(), metal(), etc.)")
  }
  
  scene$material = do.call(c, replicate(nrow(scene), material, simplify = FALSE))
  
  return(scene)
}