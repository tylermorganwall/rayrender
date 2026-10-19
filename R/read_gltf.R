#' Load a local glTF model and its animation clips
#'
#' Import node, skeletal, and morph animation into shared scene-time model assets.
#'
#' @param filename Path to a local .gltf or .glb file.
#' @param material Default `diffuse()`. Rayrender material applied to every
#'   primitive. glTF materials and textures are not imported.
#' @param scene Default `NULL`. One-based glTF scene index. NULL selects the
#'   document's default scene, then its first scene, or all roots if none exist.
#' @param recalculate_normals Default `FALSE`. Recompute smooth shading normals
#'   from the posed triangles at each requested time instead of deforming authored
#'   normals. Useful for skins with large differences in neighboring weights.
#' @return A serializable model asset for animation_model_instance() and animation_list_clips().
#' @export
#' @details
#' Local relative buffers, embedded base64 buffers, and GLB binary buffers are
#' supported. Network URIs and required extensions are rejected. Geometry uses
#' triangle primitives, float positions/normals, and optional texture coordinates.
#' Normalized integer texture coordinates and sparse/interleaved accessors are
#' decoded with their authored offsets and strides. Scene cameras and lighting
#' remain under the caller's control; imported cameras/lights are diagnosed.
#'
#' Coordinates retain glTF's right-handed, +Y-up convention, with column vectors
#' and x/y/z/w quaternions. No axis flip is applied. Authored node matrices must
#' be decomposable into nonsingular TRS and cannot also carry TRS animation.
#' STEP, linear/SLERP, and explicit cubic Hermite channels retain their independent
#' key schedules. Cubic quaternion signs and per-second tangents are preserved.
#' Node IDs are node-0, node-1, and so on; clip IDs are clip-0, clip-1, and so on.
#' Names are display labels and may be duplicated; select ambiguous clips by ID.
#'
#' Skinning and morphs evaluate from immutable rest geometry at the requested
#' frame time; rendering requires explicit time/times. Morphs precede skinning, including all joint-weight sets. Normals
#' use the inverse transpose of the blended skin map; singular blends fail.
#' With recalculate_normals = TRUE, deformed meshes use corner-angle-weighted
#' face normals, preserving splits in the imported vertex indexing. Authored
#' normals and morph normal deltas are replaced; tangents are reorthogonalized.
#' Vertices without a nondegenerate incident surface fail clearly. Static meshes
#' use mesh3d_model()'s existing normal recalculation.
#' Nonzero shutter intervals freeze the complete internal pose at frame time
#' while retaining motion blur on the outer instance placement. Zero-width
#' intervals evaluate the exact time plus shutter offset. Deformation blur is
#' not implemented. Unsupported compression and non-triangle primitives fail
#' clearly. read_gltf() remains a compatibility alias.
#' Source arrays are copied into ordinary R data before parser storage is freed;
#' no native pointer or source-buffer file is required after a successful load.
#' @examples
#' # The shipped original fixture has two asynchronous rigid clips.
#' file = system.file("animation", "fixtures", "rigid.gltf", package = "rayrender")
#' if (nzchar(file)) {
#'   asset = animation_read_gltf(file)
#'   animation_list_clips(asset)
#'   bird = animation_model_instance(asset, animation = animation_clip_player("clip-0", loop = "repeat"))
#'   animation_model_pose(bird, c(1, 0, 1))
#' }
animation_read_gltf = function(
  filename,
  material = diffuse(),
  scene = NULL,
  recalculate_normals = FALSE
) {
  if (
    !is.logical(recalculate_normals) ||
      length(recalculate_normals) != 1L ||
      is.na(recalculate_normals)
  ) {
    stop("recalculate_normals must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !file.exists(filename)
  ) {
    stop(
      "filename must identify one existing local glTF or GLB file.",
      call. = FALSE
    )
  }
  if (!inherits(material, "ray_material")) {
    stop("material must be a rayrender material.", call. = FALSE)
  }
  source = read_gltf_native(normalizePath(filename, mustWork = TRUE))
  if (!length(source$nodes)) {
    stop("glTF contains no nodes.", call. = FALSE)
  }
  if (source$cameras || source$lights) {
    warning(
      "glTF cameras and lights are not imported; use rayrender scene cameras and lights.",
      call. = FALSE
    )
  }
  if (
    !is.null(scene) &&
      (!is.numeric(scene) ||
        length(scene) != 1L ||
        !is.finite(scene) ||
        scene != trunc(scene))
  ) {
    stop(
      "scene must be NULL or one valid one-based scene index.",
      call. = FALSE
    )
  }
  ids = vapply(source$nodes, `[[`, character(1), "id")
  if (length(source$scenes)) {
    if (is.null(scene)) {
      scene = if (is.na(source$default_scene)) 1L else source$default_scene
    }
    if (scene < 1L || scene > length(source$scenes)) {
      stop("glTF scene index is out of range.", call. = FALSE)
    }
    selected = source$scenes[[scene]]
  } else {
    if (!is.null(scene)) {
      stop("glTF has no indexed scenes.", call. = FALSE)
    }
    selected = ids[vapply(
      source$nodes,
      function(node) is.null(node$parent),
      logical(1)
    )]
  }
  # Expand only selected scene roots. Keep unselected nodes as empty transform
  # nodes so clip identities and full clip intervals remain document-wide.
  cursor = 1L
  while (cursor <= length(selected)) {
    parent = selected[cursor]
    selected = c(
      selected,
      ids[vapply(
        source$nodes,
        function(node) identical(node$parent, parent),
        logical(1)
      )]
    )
    cursor = cursor + 1L
  }
  geometry = lapply(source$meshes, function(mesh) {
    Reduce(
      add_object,
      lapply(
        mesh$primitives,
        gltf_primitive_scene,
        material = material,
        recalculate_normals = recalculate_normals
      )
    )
  })
  nodes = lapply(source$nodes, function(node) {
    part = if (node$id %in% selected && !is.na(node$mesh)) {
      geometry[[node$mesh]]
    } else {
      NULL
    }
    rest = do.call(transform_track, c(list(time = 0), node$rest))
    result = animation_model_node(
      node$id,
      scene = part,
      parent = node$parent,
      transform = rest,
      name = node$name
    )
    if (!is.na(node$mesh)) {
      mesh = source$meshes[[node$mesh]]
      targets = length(mesh$primitives[[1]]$targets)
      if (!is.na(node$skin) || targets > 0L) {
        weights = if (length(node$weights)) node$weights else mesh$weights
        if (!length(weights)) {
          weights = numeric(targets)
        }
        result$deformation = list(
          primitives = mesh$primitives,
          recalculate_normals = recalculate_normals,
          weights = weights,
          skin = if (!is.na(node$skin)) source$skins[[node$skin]] else NULL
        )
        if (!is.null(part)) {
          for (i in seq_len(nrow(part))) {
            part$animation_info[[i]]$model_deformation = list(
              node = node$id,
              primitive = i
            )
          }
          result$scene = part
        }
      }
    }
    result
  })
  clips = lapply(source$clips, gltf_model_clip)
  animation_model_asset(nodes, clips, name = basename(filename))
}

#' @param primitive Decoded triangle primitive with owned attribute matrices.
#' @param material Rayrender material overriding source materials.
#' @param recalculate_normals Whether to replace authored normals.
#' @return A static mesh scene in glTF coordinates.
#' @keywords internal
#' @noRd
gltf_primitive_scene = function(primitive, material, recalculate_normals) {
  attributes = primitive$attributes
  positions = attributes$POSITION
  normals = attributes$NORMAL
  uv = attributes$TEXCOORD_0
  if (!is.null(normals) && ncol(normals) != 3L) {
    stop("glTF normals must be VEC3.", call. = FALSE)
  }
  if (!is.null(uv) && ncol(uv) != 2L) {
    stop("glTF texture coordinates must be VEC2.", call. = FALSE)
  }
  mesh = structure(
    list(
      vb = t(positions),
      it = matrix(primitive$indices, nrow = 3L),
      normals = if (!is.null(normals)) t(normals),
      texcoords = if (!is.null(uv)) t(uv)
    ),
    class = "mesh3d"
  )
  mesh3d_model(
    mesh,
    material = material,
    override_material = TRUE,
    recalculate_normals = recalculate_normals
  )
}

#' @param clip Decoded channels retaining glTF sampler semantics.
#' @return A model clip using the shared native track evaluator.
#' @keywords internal
#' @noRd
gltf_model_clip = function(clip) {
  channels = weights = list()
  for (channel in clip$channels) {
    if (channel$property == "weights") {
      rows = length(channel$time) *
        if (channel$interpolation == "cubic") 3L else 1L
      width = length(channel$value) / rows
      if (width < 1L || width != trunc(width)) {
        stop(
          "glTF morph sampler output does not match its key count.",
          call. = FALSE
        )
      }
      value = matrix(as.numeric(channel$value), nrow = rows, byrow = TRUE)
      weights[[channel$target]] = lapply(seq_len(width), function(i) {
        incoming = outgoing = NULL
        values = value[, i]
        if (channel$interpolation == "cubic") {
          keys = seq_along(channel$time) * 3L
          incoming = matrix(values[keys - 2L], ncol = 1L)
          outgoing = matrix(values[keys], ncol = 1L)
          values = values[keys - 1L]
        }
        animation_keyframe_track(
          channel$time,
          values,
          property = "scalar",
          interpolation = channel$interpolation,
          incoming = incoming,
          outgoing = outgoing
        )
      })
      next
    }
    incoming = outgoing = NULL
    value = channel$value
    if (channel$interpolation == "cubic") {
      rows = seq_len(length(channel$time)) * 3L
      incoming = value[rows - 2L, , drop = FALSE]
      outgoing = value[rows, , drop = FALSE]
      value = value[rows - 1L, , drop = FALSE]
    }
    property = animation_keyframe_track(
      channel$time,
      value,
      channel$property,
      interpolation = channel$interpolation,
      incoming = incoming,
      outgoing = outgoing
    )
    if (is.null(channels[[channel$target]])) {
      channels[[channel$target]] = animation_transform_track()
    }
    channels[[channel$target]][[channel$property]] = property
  }
  animation_model_clip(clip$id, channels, name = clip$name, weights = weights)
}

#' @rdname animation_read_gltf
#' @export
read_gltf = animation_read_gltf
