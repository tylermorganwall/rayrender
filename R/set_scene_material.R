#' Set Material for All Objects
#'
#' @param scene A ray_scene object.
#' @param material A material specification created by diffuse(), metal(), dielectric(), etc.
#'
#' @return A modified ray_scene with the new material applied to all objects
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # Create a scene with different materials
#' scene = generate_cornell() |>
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
    stop(
      "Material must be created using material functions (diffuse(), metal(), etc.)"
    )
  }

  for (i in seq_len(nrow(scene))) {
    if (identical(scene$shape_info[[i]]$medium_owner, "subsurface")) {
      scene$shape_info[[i]]$medium = NULL
      scene$shape_info[[i]]$medium_keep_surface = NULL
      scene$shape_info[[i]]$medium_owner = NULL
    }
    original = scene$shape_info[[i]]$shape_properties$original_scene
    if (!is.null(original)) {
      scene$shape_info[[i]]$shape_properties$original_scene[[1]] =
        set_scene_material(original[[1]], material)
    }
  }
  scene$material = do.call(
    c,
    replicate(nrow(scene), material, simplify = FALSE)
  )

  # An instance row is a placement, not a closed material boundary.
  for (i in seq_len(nrow(scene))) {
    if (!is.null(scene$shape_info[[i]]$shape_properties$original_scene)) {
      scene$material[[i]]$subsurface = NULL
      scene$shape_info[[i]]$shape_properties$any_light =
        material[[1]]$type %in% c(5L, 8L)
    }
  }

  return(scene)
}
