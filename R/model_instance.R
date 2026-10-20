#' Place and animate a shared model
#'
#' Give a shared model its own outer path and independently timed internal clip.
#'
#' @param asset An animation_model_asset() containing immutable geometry and clips.
#' @param transform Default `NULL`. Outer animation_transform_track() in scene seconds.
#'   NULL retains identity placement and allows a later animation_keyframes() call.
#' @param animation Default `animation_clip_player()`. Internal clip selection and clock.
#'   A NULL clip selects the rest pose, without guessing the first clip.
#' @return One scene row owning a complete model instance and stable instance ID.
#' @export
#' @details
#' World transforms are outer placement * parent globals * local pose. Each
#' instance shares the asset's geometry and tracks, but owns its playback settings.
#' In explicit-time rendering, the default `animation_motion_blur = "full"`
#' evaluates internal nodes, skins, morphs, and outer placement at a shared
#' shutter time per image sample. Geometry and bounds are updated between samples
#' while BVH topology is retained. Select `animation_motion_blur = "placement"`
#' in [render_scene()] or [render_animation()] to freeze internal poses at frame
#' time and blur only placement. A zero-width interval samples the complete model
#' at time plus its shutter offset. Seeking never depends on the preceding pose.
#' Static group_objects() changes stay inside the outer animated placement.
#' add_object() remaps colliding instance IDs and every nested reference when a
#' copied instance is combined. Node IDs remain asset-local and are qualified by
#' instance ID. Serialization preserves shared assets without native pointers.
#' Rendered node clocks use the outer instance's playback settings. Reloaded
#' scenes resolve unchanged geometry through the shared asset; explicit
#' set_scene_material() edits retain private geometry without changing the asset.
#' @seealso \code{vignette("model-animation", package = "rayrender")}
#' @examples
#' # Two paths can share an asset while using different phases and speeds.
#' asset = animation_model_asset(list(animation_model_node("body", cube())))
#' first = animation_model_instance(asset, animation_transform_track(0, position = c(-1, 0, 0)))
#' second = animation_model_instance(asset, animation_transform_track(0, position = c(1, 0, 0)))
#' scene = add_object(first, second)
animation_model_instance = function(
  asset,
  transform = NULL,
  animation = animation_clip_player()
) {
  if (!inherits(asset, "ray_model_asset")) {
    stop("asset must be a model asset.", call. = FALSE)
  }
  if (!inherits(animation, "ray_clip_player")) {
    stop("animation must be an animation_clip_player().", call. = FALSE)
  }
  if (!is.null(transform) && !inherits(transform, "ray_transform_track")) {
    stop(
      "transform must be an animation_transform_track() or NULL.",
      call. = FALSE
    )
  }
  model_pose_native(
    asset,
    animation,
    if (is.null(transform)) animation_transform_track() else transform,
    numeric()
  )
  binding = list(
    id = basename(tempfile("instance-")),
    asset = asset,
    player = animation
  )
  result = build_model_hierarchy(asset, binding)
  result$animation_info[[1]]$model_instance = binding
  if (!is.null(transform)) {
    result = animation_keyframes(result, transform)
  }
  preserve_ray_scene_attrs(result, asset$metadata)
}

#' @param asset Validated immutable model data.
#' @param binding Instance identity, asset reference, and clock.
#' @return A nested ordinary instance scene, sharing source geometry prototypes.
#' @keywords internal
#' @noRd
build_model_hierarchy = function(asset, binding) {
  nodes = asset$nodes
  ids = vapply(nodes, `[[`, character(1), "id")
  parents = vapply(
    nodes,
    function(node) if (is.null(node$parent)) NA_character_ else node$parent,
    character(1)
  )
  parent_index = match(parents, ids)
  children = lapply(seq_along(nodes), function(i) which(parent_index == i))
  remaining = lengths(children)
  queue = which(remaining == 0L)
  built = vector("list", length(nodes))
  cursor = 1L
  # Build leaves first. Each wrapper is an actual local transform node, not a
  # baked world-space key sequence. Geometry prototypes retain shared R storage.
  while (cursor <= length(queue)) {
    i = queue[cursor]
    contents = asset$prototypes[[i]]
    if (!is.null(contents)) {
      # Resolve geometry through the shared asset after serialization, which
      # preserves environment identity but can duplicate ordinary R lists.
      contents$animation_info[[1]]$model_geometry = ids[i]
    }
    for (child in children[[i]]) {
      if (is.null(contents)) {
        contents = built[[child]]
      } else if (!is.null(built[[child]])) {
        contents = add_object(contents, built[[child]])
      }
    }
    if (!is.null(contents)) {
      entry = create_instances(contents)
      entry$animation_info[[1]]$model_node = list(
        id = ids[i],
        instance = binding
      )
      built[[i]] = entry
    }
    parent = parent_index[i]
    if (!is.na(parent)) {
      remaining[parent] = remaining[parent] - 1L
      if (remaining[parent] == 0L) queue = c(queue, parent)
    }
    cursor = cursor + 1L
  }
  roots = Filter(Negate(is.null), built[is.na(parent_index)])
  if (!length(roots)) {
    stop("A rendered model instance requires geometry.", call. = FALSE)
  }
  create_instances(Reduce(add_object, roots))
}

#' Inspect exact model poses
#'
#' Evaluate a model instance at arbitrary scene times without opening a renderer.
#'
#' @param instance One scene row returned by animation_model_instance().
#' @param time Finite vector of scene times in seconds, in any order.
#' @return A list per requested time with clip_time, named local/global/world
#'   matrices, morph weights, and separate placement_changed and pose_changed flags relative to
#'   the preceding requested result. Matrices never depend on that result.
#' @export
#' @details
#' Global matrices exclude placement. World matrices include the instance's
#' static grouping and animated placement. The flags describe changed matrices;
#' a different scene time need not change a held or repeated pose.
#' Legacy animate_objects() matrix attachments remain renderable, but this
#' inspector requires an animation_transform_track() placement and rejects legacy attachments.
#' @seealso \code{vignette("model-animation", package = "rayrender")}
#' @examples
#' # Repeated and backward queries produce the same matrices.
#' asset = animation_model_asset(list(animation_model_node("body", sphere())))
#' animation_model_pose(animation_model_instance(asset), c(2, 0, 2))
animation_model_pose = function(instance, time) {
  if (
    !inherits(instance, "data.frame") ||
      nrow(instance) != 1L ||
      is.null(instance$animation_info[[1]]$model_instance)
  ) {
    stop(
      "instance must contain one animation_model_instance() row.",
      call. = FALSE
    )
  }
  if (!is.numeric(time) || any(!is.finite(time))) {
    stop("time must contain finite scene seconds.", call. = FALSE)
  }
  info = instance$animation_info[[1]]
  if (nrow(info$start_transform_animation[[1]]) == 4L) {
    stop(
      "animation_model_pose() does not inspect legacy matrix attachments; use an animation_transform_track() placement.",
      call. = FALSE
    )
  }
  binding = info$model_instance
  track = if (is.null(info$keyframes)) {
    animation_transform_track()
  } else {
    info$keyframes
  }
  result = model_pose_native(binding$asset, binding$player, track, time)
  group = instance$transforms[[1]]$group_transform[[1]]
  if (nrow(group) == 4L) {
    placements = animation_evaluate_transform(track, time)
    for (i in seq_along(result)) {
      result[[i]]$world = lapply(result[[i]]$global, function(global) {
        placements[[i]] %*% group %*% global
      })
    }
  }
  result
}
