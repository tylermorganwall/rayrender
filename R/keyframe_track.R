#' Author a property animation track
#'
#' Describe an independently timed property using keys in scene seconds.
#'
#' @param time Strictly increasing finite key times, in seconds.
#' @param value Numeric matrix with one row per key. A vector is accepted for
#'   a single key or a scalar track.
#' @param property Default `"position"`. One of position, rotation, scale, scalar.
#'   Rotation values are quaternions in x, y, z, w order.
#' @param interpolation Default `"linear"`. One of step, linear, cubic, smooth.
#'   Linear rotation uses shortest-arc SLERP; smooth is native cubic rotation.
#' @param incoming Default `NULL`. Matrix of incoming cubic derivatives per second.
#' @param outgoing Default `NULL`. Matrix of outgoing cubic derivatives per second.
#' @return A serializable `ray_keyframe_track` authoring object.
#' @export
#' @details
#' Values are held outside the key range. STEP selects the new value at an exact
#' key. One-key constant tracks are valid. Cubic vector tracks without tangents
#' use nonuniform centered differences weighted by opposite interval lengths;
#' endpoints use one-sided differences. Smooth rotation unwraps quaternion signs
#' and projects these shared derivatives onto the unit-quaternion tangent plane.
#' Normalizing the resulting Hermite curve preserves angular velocity continuity.
#' Explicit cubic quaternion tracks retain authored signs and require unit keys
#' and both tangent matrices. Singular or uncertifiably near-singular curves fail.
#' Continuous scale curves may not cross zero; negative scales are supported.
#' Scalar/vector trajectories and pivots support magnitudes up to 1e12. Scale
#' magnitudes must stay at least 1e-8 away from zero.
#' Quaternion endpoints do not encode revolution counts: supply intermediate
#' orientations less than 180 degrees apart to author multiple turns.
#' @examples
#' # Position and rotation can use different key times.
#' position = keyframe_track(c(0, 1, 3), rbind(c(0, 0, 0), c(1, 2, 0), c(4, 0, 0)),
#'   interpolation = "cubic")
#' rotation = keyframe_track(c(0, 2), rbind(c(0, 0, 0, 1), c(0, 1, 0, 0)),
#'   property = "rotation")
keyframe_track = function(
  time,
  value,
  property = "position",
  interpolation = "linear",
  incoming = NULL,
  outgoing = NULL
) {
  property = match.arg(property, c("position", "rotation", "scale", "scalar"))
  interpolation = match.arg(
    interpolation,
    c("step", "linear", "cubic", "smooth")
  )
  if (!is.numeric(time) || !is.numeric(value)) {
    stop("Track times and values must be numeric.", call. = FALSE)
  }
  dimension = switch(property, rotation = 4L, scalar = 1L, 3L)
  if (is.null(dim(value))) {
    if (dimension == 1L) {
      value = matrix(value, ncol = 1)
    } else if (length(time) == 1L && length(value) == dimension) {
      value = matrix(value, nrow = 1)
    } else {
      stop(
        "Supply a matrix with one row per key; values are never recycled.",
        call. = FALSE
      )
    }
  }
  result = structure(
    list(
      version = 1L,
      property = property,
      time = time,
      value = value,
      interpolation = interpolation,
      incoming = incoming,
      outgoing = outgoing
    ),
    class = "ray_keyframe_track"
  )
  animation_evaluate_track(result, numeric())
  result
}

#' Author an outer transform animation
#'
#' Combine position, orientation, and scale with independently timed channels.
#'
#' @param time Default `NULL`. Shared key times for matrix-valued channels, in seconds.
#' @param position Default `NULL`. Position matrix or a position keyframe track.
#' @param angle Default `NULL`. Euler angle matrix in degrees, applied in x, y, z
#'   order using rayrender's existing rotation convention.
#' @param rotation Default `NULL`. Quaternion matrix in x, y, z, w order, or a
#'   rotation keyframe track. Mutually exclusive with angle.
#' @param scale Default `NULL`. Scale matrix or a scale keyframe track.
#' @param pivot Default `c(0, 0, 0)`. Fixed local pivot for rotation and scale.
#' @param position_interpolation Default `"linear"`. Position interpolation.
#' @param rotation_interpolation Default `"linear"`. Rotation interpolation.
#' @param scale_interpolation Default `"linear"`. Scale interpolation.
#' @return A serializable `ray_transform_track` authoring object.
#' @export
#' @details
#' The column-vector composition is T(position) T(pivot) R S T(-pivot), outside
#' existing static object and group transforms. Absent channels retain identity.
#' A supplied keyframe track retains its own key times and interpolation.
#' @examples
#' # Four keys make a curved path without constraining output frame times.
#' flight = transform_track(time = c(0, 2, 5, 8),
#'   position = rbind(c(0, 1, 0), c(2, 3, 0), c(6, 4, 2), c(9, 2, 4)),
#'   angle = rbind(c(0, 0, 0), c(0, 15, 0), c(0, 40, 0), c(0, 60, 0)),
#'   position_interpolation = "cubic", rotation_interpolation = "smooth")
transform_track = function(
  time = NULL,
  position = NULL,
  angle = NULL,
  rotation = NULL,
  scale = NULL,
  pivot = c(0, 0, 0),
  position_interpolation = "linear",
  rotation_interpolation = "linear",
  scale_interpolation = "linear"
) {
  if (!is.null(angle) && !is.null(rotation)) {
    stop("Supply only one of angle or rotation.", call. = FALSE)
  }
  if (!is.null(angle)) {
    if (
      is.numeric(angle) &&
        is.null(dim(angle)) &&
        length(time) == 1L &&
        length(angle) == 3L
    ) {
      angle = matrix(angle, nrow = 1)
    }
    if (
      !is.numeric(angle) ||
        !is.matrix(angle) ||
        ncol(angle) != 3L ||
        nrow(angle) != length(time) ||
        any(!is.finite(angle))
    ) {
      stop(
        "angle requires a finite three-column matrix with one row per time.",
        call. = FALSE
      )
    }
    rotation = animation_euler_quaternions(angle)
  }
  channels = list(position = position, rotation = rotation, scale = scale)
  methods = c(
    position_interpolation,
    rotation_interpolation,
    scale_interpolation
  )
  for (i in seq_along(channels)) {
    value = channels[[i]]
    if (is.null(value)) {
      next
    }
    if (!inherits(value, "ray_keyframe_track")) {
      channels[i] = list(keyframe_track(
        time,
        value,
        names(channels)[i],
        methods[i]
      ))
    }
  }
  result = structure(
    c(list(version = 1L), channels, list(pivot = pivot)),
    class = "ray_transform_track"
  )
  animation_evaluate_transform(result, numeric())
  result
}

#' Attach a keyframed outer transform
#'
#' Move a scene fragment along a deterministic path while retaining its static transforms.
#'
#' @param scene Scene or object fragment to animate.
#' @param track A transform track created by transform_track().
#' @param player Default `NULL`. A clip_player() mapping scene seconds to the
#'   track's local key interval. NULL evaluates keys directly in scene time.
#' @return The scene with an animation attachment in each selected row.
#' @export
#' @details
#' All rows share the authored transform about its fixed pivot. Bind, copy, and
#' subset operations retain attachments in the existing animation_info column.
#' Repeated attachments, including a legacy animate_objects() attachment, are
#' rejected. Static group_objects() transforms stay inside the animated transform.
#' Use render_scene(time = ...) for snapshots or render_animation(times = ...)
#' for an explicit output sequence. Device key times use validated frame-relative
#' floats; an unrepresentable key schedule produces an explicit CPU fallback.
#' @examples
#' # Keep the sphere's static offset inside the animated placement.
#' path = transform_track(c(0, 1, 2, 3),
#'   position = rbind(c(0, 0, 0), c(1, 1, 0), c(2, 0, 0), c(3, 1, 0)))
#' scene = animate_keyframes(sphere(x = 0.2), path)
animate_keyframes = function(scene, track, player = NULL) {
  if (!inherits(track, "ray_transform_track")) {
    stop("track must be a transform_track().", call. = FALSE)
  }
  if (!is.null(player) && !inherits(player, "ray_clip_player")) {
    stop("player must be a clip_player().", call. = FALSE)
  }
  track$player = player
  animation_evaluate_transform(track, numeric())
  for (i in seq_len(nrow(scene))) {
    info = scene$animation_info[[i]]
    if (
      !is.null(info$keyframes) ||
        nrow(info$start_transform_animation[[1]]) == 4L
    ) {
      stop("The fragment already has an animation attachment.", call. = FALSE)
    }
    scene$animation_info[[i]]$keyframes = track
  }
  scene
}
