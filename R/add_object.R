#' Add Object
#'
#' @param scene Tibble of pre-existing object locations and properties.
#' @param objects A tibble row or collection of rows representing each object.
#'
#' @return Tibble of object locations and properties.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' #Generate the ground and add some objects
#' scene = generate_ground(depth=-0.5,material = diffuse(checkercolor="blue")) |>
#'   add_object(cube(x=0.7,
#'                   material=diffuse(noise=5,noisecolor="purple",color="black",noisephase=45),
#'                   angle=c(0,-30,0))) |>
#'   add_object(sphere(x=-0.7,radius=0.5,material=metal(color="gold")))
#' render_scene(scene,parallel=TRUE)
add_object = function(scene, objects = NULL) {
  if (is.null(objects)) {
    return(scene)
  }
  if (inherits(scene, "ray_scene_v2") || inherits(objects, "ray_scene_v2")) {
    scene = ensure_ray_scene_v2(scene)
    objects = ensure_ray_scene_v2(objects)
    scene_attrs = ray_scene_attrs(scene)
    object_attrs = ray_scene_attrs(objects)
    newscene = rbind(as.data.frame(scene), as.data.frame(objects))
    class(newscene) = class(scene)
    newscene = restore_ray_scene_attrs(
      newscene,
      merge_ray_scene_attrs(scene_attrs, object_attrs)
    )
    if (
      !is.null(attr(objects, "cornell")) || !is.null(attr(scene, "cornell"))
    ) {
      attr(newscene, "cornell") = TRUE
    }
    return(assign_missing_object_ids(newscene))
  }
  newscene = rbind(scene, objects)
  if (!is.null(attr(objects, "cornell")) || !is.null(attr(scene, "cornell"))) {
    attr(newscene, "cornell") = TRUE
  }
  return(newscene)
}
