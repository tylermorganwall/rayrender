#' Author a model hierarchy
#'
#' Share geometry and animation clips across independently placed model instances.
#'
#' @param nodes List of animation_model_node() objects with unique IDs.
#' @param clips Default `list()`. List of animation_model_clip() objects with unique IDs.
#' @param name Default `NULL`. Optional descriptive asset label.
#' @return A locked, serializable ray_model_asset environment containing source
#'   geometry, rest transforms, and clips.
#' @export
#' @details
#' Nodes use column-vector transforms: global = parent_global * local. Clip
#' properties replace the corresponding rest property; missing properties retain
#' rest values. Node IDs are unique within the asset; names are only labels.
#' Geometry retains its existing static object and group transforms inside the
#' node transform. Existing group_objects() behavior is unchanged.
#' Assets own ordinary R data and locked bindings, with no serialized native
#' pointer. Copies and instances share the asset. Rendering compiles owned native
#' data once per asset and shares source geometry and property tracks. Individual
#' instances own clocks and transform descriptors. Removing an instance removes
#' its complete hierarchy; subsetting scene rows cannot leave dangling children.
#' @seealso \code{vignette("model-animation", package = "rayrender")}
#' @examples
#' # A wing rotates around its local origin under the body node.
#' body = animation_model_node("body", sphere(radius = 0.3))
#' wing = animation_model_node("wing", cube(x = 0.5, xwidth = 1, ywidth = 0.05),
#'   parent = "body")
#' flap = animation_model_clip("flap", list(wing = animation_transform_track(c(0, 0.5, 1),
#'   angle = rbind(c(0, 0, -30), c(0, 0, 30), c(0, 0, -30)))))
#' bird = animation_model_asset(list(body, wing), list(flap))
#' animation_model_instance(bird, animation = animation_clip_player("flap", loop = "repeat"))
animation_model_asset = function(nodes, clips = list(), name = NULL) {
  if (
    !is.list(nodes) ||
      !length(nodes) ||
      !all(vapply(nodes, inherits, logical(1), "ray_model_node"))
  ) {
    stop(
      "nodes must be a nonempty list of animation_model_node() objects.",
      call. = FALSE
    )
  }
  if (
    !is.list(clips) ||
      !all(vapply(clips, inherits, logical(1), "ray_model_clip"))
  ) {
    stop(
      "clips must be a list of animation_model_clip() objects.",
      call. = FALSE
    )
  }
  if (!is.null(name)) {
    validate_model_label(name, "name")
  }
  asset = new.env(parent = emptyenv())
  asset$version = if (
    any(vapply(nodes, function(node) !is.null(node$deformation), logical(1)))
  ) {
    2L
  } else {
    1L
  }
  asset$id = basename(tempfile("asset-"))
  asset$name = name
  asset$nodes = nodes
  asset$clips = clips
  class(asset) = "ray_model_asset"
  # Validate the entire graph and all channels before preparing any geometry.
  animation_model_clips(asset)
  prototypes = vector("list", length(nodes))
  node_ids = vapply(nodes, `[[`, character(1), "id")
  geometry_nodes = setNames(node_ids, node_ids)
  metadata = sphere()[FALSE, ]
  for (i in seq_along(nodes)) {
    scene = nodes[[i]]$scene
    if (is.null(scene) || !nrow(scene)) {
      next
    }
    if (model_scene_has_animation(scene)) {
      stop(
        "Model geometry must be static; author local animation with animation_model_clip().",
        call. = FALSE
      )
    }
    metadata = preserve_ray_scene_attrs(metadata, metadata, scene)
    matching = which(vapply(
      seq_len(i - 1L),
      function(j) {
        identical(scene, nodes[[j]]$scene)
      },
      logical(1)
    ))
    prototypes[[i]] = if (length(matching)) {
      geometry_nodes[i] = geometry_nodes[matching[1]]
      prototypes[[matching[1]]]
    } else {
      create_instances(scene)
    }
  }
  asset$prototypes = prototypes
  asset$geometry_nodes = geometry_nodes
  asset$metadata = metadata
  lockEnvironment(asset, bindings = TRUE)
  asset
}

#' @param value Label or ID to validate.
#' @param argument Name used in the error message.
#' @return Invisibly the validated string.
#' @keywords internal
#' @noRd
validate_model_label = function(value, argument) {
  if (
    !is.character(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !nzchar(value)
  ) {
    stop(argument, " must be one nonempty string.", call. = FALSE)
  }
  invisible(value)
}

#' @param scene Source geometry, including nested processed prototypes.
#' @return Whether any row or descendant already owns animation.
#' @keywords internal
#' @noRd
model_scene_has_animation = function(scene) {
  for (info in scene$animation_info) {
    if (
      !is.null(info$keyframes) ||
        !is.null(info$model_node) ||
        !is.null(info$model_instance) ||
        nrow(info$start_transform_animation[[1]]) == 4L
    ) {
      return(TRUE)
    }
  }
  for (info in scene$shape_info) {
    prototype = info$shape_properties$original_scene
    if (is.null(prototype)) {
      prototype = info$original_scene
    }
    if (!is.null(prototype) && model_scene_has_animation(prototype[[1]])) {
      return(TRUE)
    }
  }
  FALSE
}

#' Define a model part and its parent
#'
#' Attach static geometry to a local coordinate system with a stable target ID.
#'
#' @param id Unique node ID within the model asset.
#' @param scene Default `NULL`. Static rayrender geometry in this node's coordinates.
#' @param parent Default `NULL`. Parent node ID; NULL denotes a top-level node.
#' @param transform Default `animation_transform_track()`. Constant rest transform, with
#'   at most one key per property. Its pivot also applies to animated channels.
#' @param name Default `id`. Display label, which need not be unique.
#' @return A serializable ray_model_node authoring object.
#' @export
#' @seealso \code{vignette("model-animation", package = "rayrender")}
#' @examples
#' # Offset a wing's hinge while retaining its geometry relative to that hinge.
#' animation_model_node("left", cube(x = -0.5, xwidth = 1, ywidth = 0.05),
#'   parent = "body", transform = animation_transform_track(0, position = c(-0.2, 0, 0)))
animation_model_node = function(
  id,
  scene = NULL,
  parent = NULL,
  transform = animation_transform_track(),
  name = id
) {
  validate_model_label(id, "id")
  validate_model_label(name, "name")
  if (!is.null(parent)) {
    validate_model_label(parent, "parent")
  }
  if (!is.null(scene) && !inherits(scene, "data.frame")) {
    stop("scene must be a rayrender scene or NULL.", call. = FALSE)
  }
  if (!inherits(transform, "ray_transform_track")) {
    stop("transform must be an animation_transform_track().", call. = FALSE)
  }
  animation_evaluate_transform(transform, numeric())
  for (channel in transform[c("position", "rotation", "scale")]) {
    if (!is.null(channel) && length(channel$time) != 1L) {
      stop("Rest properties require one key.", call. = FALSE)
    }
  }
  if (!is.null(transform$player)) {
    stop("Rest transforms cannot have a player.", call. = FALSE)
  }
  structure(
    list(
      id = id,
      name = name,
      parent = parent,
      scene = scene,
      transform = transform
    ),
    class = "ray_model_node"
  )
}

#' Author a named model clip
#'
#' Animate independently timed properties of explicitly identified model nodes.
#'
#' @param id Stable clip ID, unique within the asset.
#' @param channels Named list of animation_transform_track() objects, named by target node
#'   IDs. Each property replaces that node's corresponding rest property.
#' @param name Default `id`. Display name; select by ID when names are ambiguous.
#' @param weights Default `list()`. Named list of morph channels by node ID.
#'   Each entry is a list of scalar animation_keyframe_track() objects, one per target.
#' @return A serializable ray_model_clip object. Its interval spans all channel keys.
#' @export
#' @details
#' Each node occurs at most once in channels; combine its independent properties
#' in one animation_transform_track(). Set pivots on animation_model_node(), not on clip channels.
#' A property's stable identity is the tuple (asset ID, clip ID, node ID, property).
#' All channels share the clip interval and instance clock, even when individual
#' properties have shorter key ranges. Values hold outside each property's keys.
#' @seealso \code{vignette("model-animation", package = "rayrender")}
#' @examples
#' # A short wing cycle does not wrap independently of the complete clip.
#' animation_model_clip("motion", list(
#'   body = animation_transform_track(c(0, 2), position = rbind(c(0, 0, 0), c(1, 0, 0))),
#'   wing = animation_transform_track(c(0, 1), angle = rbind(c(0, 0, -30), c(0, 0, 30)))))
animation_model_clip = function(id, channels, name = id, weights = list()) {
  validate_model_label(id, "id")
  validate_model_label(name, "name")
  if (
    !is.list(channels) ||
      (length(channels) &&
        (is.null(names(channels)) ||
          anyNA(names(channels)) ||
          any(!nzchar(names(channels))) ||
          anyDuplicated(names(channels))))
  ) {
    stop("channels must be a list with unique target node IDs.", call. = FALSE)
  }
  if (!all(vapply(channels, inherits, logical(1), "ray_transform_track"))) {
    stop(
      "Each clip channel must be an animation_transform_track().",
      call. = FALSE
    )
  }
  if (
    !is.list(weights) ||
      (length(weights) &&
        (is.null(names(weights)) ||
          anyNA(names(weights)) ||
          any(!nzchar(names(weights))) ||
          anyDuplicated(names(weights))))
  ) {
    stop("weights must be a list with unique target node IDs.", call. = FALSE)
  }
  for (node in weights) {
    if (
      !is.list(node) ||
        !length(node) ||
        !all(vapply(
          node,
          function(track) {
            inherits(track, "ray_keyframe_track") &&
              identical(track$property, "scalar")
          },
          logical(1)
        ))
    ) {
      stop("Morph weights require lists of scalar tracks.", call. = FALSE)
    }
  }
  structure(
    list(id = id, name = name, channels = channels, weights = weights),
    class = "ray_model_clip"
  )
}

#' Inspect an asset's animation clips
#'
#' List stable clip IDs, display names, and local time intervals in seconds.
#'
#' @param asset An asset created by animation_model_asset().
#' @return A data frame with id, name, start, end, duration, and targets columns.
#' @export
#' @examples
#' # A model with no clips is still a reusable static hierarchy.
#' animation_list_clips(animation_model_asset(list(animation_model_node("body", sphere()))))
animation_list_clips = function(asset) {
  if (!inherits(asset, "ray_model_asset")) {
    stop("asset must be a model asset.", call. = FALSE)
  }
  result = animation_model_clips(asset)
  result$targets = lapply(asset$clips, function(clip) {
    union(names(clip$channels), names(clip$weights))
  })
  result
}


#' @rdname animation_list_clips
#' @export
list_animations = animation_list_clips
