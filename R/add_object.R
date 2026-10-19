#' Add Object
#'
#' @param scene Tibble of pre-existing object locations and properties.
#' @param objects Default `NULL`. A tibble row or collection of object rows to add.
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
  objects = remap_model_instances(objects, scene)
  newscene = rbind.data.frame(scene, objects)
  preserve_ray_scene_attrs(newscene, scene, objects)
}

#' @param objects Scene fragment being combined.
#' @param scene Existing scene whose instance IDs must remain unique.
#' @return The fragment with colliding instance IDs and nested references remapped.
#' @keywords internal
#' @noRd
remap_model_instances = function(objects, scene) {
  used = model_instance_ids(scene)
  # One identity set spans the entire added fragment, including models inside
  # ordinary create_instances() wrappers. Assets themselves remain immutable.
  remap = function(fragment) {
    for (i in seq_len(nrow(fragment))) {
      binding = fragment$animation_info[[i]]$model_instance
      if (!is.null(binding)) {
        if (binding$id %in% used) {
          binding$id = basename(tempfile("instance-"))
          fragment$animation_info[[i]]$model_instance = binding
          prototype = fragment$shape_info[[i]]$shape_properties$original_scene
          fragment$shape_info[[i]]$shape_properties$original_scene = list(
            rebind_model_nodes(prototype[[1]], binding)
          )
        }
        used <<- c(used, binding$id)
        next
      }
      prototype = fragment$shape_info[[i]]$shape_properties$original_scene
      if (!is.null(prototype) && length(model_instance_ids(prototype[[1]]))) {
        fragment$shape_info[[i]]$shape_properties$original_scene = list(
          remap(prototype[[1]])
        )
      }
    }
    fragment
  }
  remap(objects)
}

#' @param scene Scene fragment, possibly containing nested model instances.
#' @return Stable model instance IDs found in the fragment and its descendants.
#' @keywords internal
#' @noRd
model_instance_ids = function(scene) {
  ids = character()
  for (i in seq_len(nrow(scene))) {
    binding = scene$animation_info[[i]]$model_instance
    if (!is.null(binding)) {
      ids = c(ids, binding$id)
      # A model asset contains static geometry, never another model instance.
      next
    }
    prototype = scene$shape_info[[i]]$shape_properties$original_scene
    if (!is.null(prototype)) {
      ids = c(ids, model_instance_ids(prototype[[1]]))
    }
  }
  ids
}

#' @param scene Processed nested model hierarchy.
#' @param binding Replacement instance identity and playback configuration.
#' @return Hierarchy with all node references pointing to the new binding.
#' @keywords internal
#' @noRd
rebind_model_nodes = function(scene, binding) {
  for (i in seq_len(nrow(scene))) {
    node = scene$animation_info[[i]]$model_node
    if (!is.null(node)) {
      node$instance = binding
      scene$animation_info[[i]]$model_node = node
    }
    prototype = scene$shape_info[[i]]$shape_properties$original_scene
    # Only hierarchy wrappers carry animation; immutable geometry prototypes
    # need no changes and must retain their shared storage identity.
    if (!is.null(prototype) && scene_has_keyframes(prototype[[1]])) {
      scene$shape_info[[
        i
      ]]$shape_properties$original_scene = list(rebind_model_nodes(
        prototype[[1]],
        binding
      ))
    }
  }
  scene
}

#' @param ... Scene fragments to combine.
#' @param deparse.level Default `1`. Retained for compatibility with base rbind().
#' @return A combined scene with preserved metadata and independent instance IDs.
#' @keywords internal
#' @noRd
#' @export
rbind.ray_scene = function(..., deparse.level = 1) {
  scenes = Filter(Negate(is.null), list(...))
  Reduce(add_object, scenes)
}
