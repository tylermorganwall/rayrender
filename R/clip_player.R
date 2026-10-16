#' Play animation on an independent clock
#'
#' Map absolute scene seconds to a clip using an offset, playback speed, and
#' loop policy. The same scene time always selects the same local pose.
#'
#' @param clip Default `NULL`. Clip name or identifier for an animated asset.
#'   Use NULL when attaching a player to a bare transform track.
#' @param start_time Default `0`. Scene seconds anchoring the local clock.
#' @param offset Default `0`. Local seconds from the beginning of the clip at
#'   start_time. Negative speeds do not automatically change this offset.
#' @param speed Default `1`. Local seconds per scene second; negative and zero
#'   values select reverse and frozen playback.
#' @param loop Default `"hold"`. One of hold, repeat, or pingpong.
#' @param root_motion Default `"preserve"`. Root-motion policy. Only preserve
#'   is currently supported; removal requires an explicit locomotion root.
#' @return A serializable ray_clip_player configuration.
#' @export
#' @details
#' For an interval [a,b], local displacement is
#' speed * (scene_time - start_time) + offset. Hold clamps to [a,b], repeat wraps
#' to [a,b), and pingpong reflects at both ends. Exact repeat boundaries select
#' a; pingpong and hold can select b. A zero-duration clip stays at a. Start time
#' is an anchor, not a visibility gate, so negative scene times are defined.
#' For animate_keyframes(), the interval spans all channel keys. The player maps
#' each ray's time, including exposures crossing a loop boundary. Authored
#' endpoints are retained, so repeating different endpoints creates a jump.
#' animation_loop_seam() checks endpoint pose agreement without changing keys.
#' @examples
#' # Reverse from the end of a two-second clip, then reflect at each end.
#' path = transform_track(c(0, 1, 2),
#'   position = rbind(c(-1, 0, 0), c(0, 1, 0), c(1, 0, 0)))
#' player = clip_player(offset = 2, speed = -1, loop = "pingpong")
#' scene = animate_keyframes(sphere(), path, player = player)
clip_player = function(
  clip = NULL,
  start_time = 0,
  offset = 0,
  speed = 1,
  loop = "hold",
  root_motion = "preserve"
) {
  loop = match.arg(loop, c("hold", "repeat", "pingpong"))
  for (value in list(start_time, offset, speed)) {
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
      stop("Clip clock fields must be finite scalars.", call. = FALSE)
    }
  }
  if (
    !is.null(clip) &&
      (!is.character(clip) ||
        length(clip) != 1L ||
        is.na(clip) ||
        !nzchar(clip))
  ) {
    stop("clip must be NULL or one nonempty name or identifier.", call. = FALSE)
  }
  if (!identical(root_motion, "preserve")) {
    stop(
      "Only preserve root motion is currently supported; removal requires an explicit locomotion root.",
      call. = FALSE
    )
  }
  structure(
    list(
      version = 1L,
      clip = clip,
      start_time = start_time,
      offset = offset,
      speed = speed,
      loop = loop,
      root_motion = root_motion
    ),
    class = "ray_clip_player"
  )
}

#' Check the pose discontinuity at a transform loop seam
#'
#' Compare a transform's first and last local poses before configuring repetition.
#'
#' @param track A transform_track() to inspect.
#' @param tolerance Default `1e-6`. Maximum absolute difference between matrix
#'   entries considered continuous; translation entries use scene distance units.
#' @return A list with the local interval, maximum matrix-entry gap, and a
#'   pose_continuous flag. This does not certify velocity continuity.
#' @export
#' @examples
#' # This path returns to its initial position.
#' path = transform_track(c(0, 1, 2),
#'   position = rbind(c(0, 0, 0), c(1, 1, 0), c(0, 0, 0)))
#' animation_loop_seam(path)
animation_loop_seam = function(track, tolerance = 1e-6) {
  if (!inherits(track, "ray_transform_track")) {
    stop("track must be a transform_track().", call. = FALSE)
  }
  if (
    !is.numeric(tolerance) ||
      length(tolerance) != 1L ||
      !is.finite(tolerance) ||
      tolerance < 0
  ) {
    stop("tolerance must be one finite nonnegative number.", call. = FALSE)
  }
  times = unlist(lapply(track[c("position", "rotation", "scale")], function(x) {
    x$time
  }))
  interval = if (length(times)) range(times) else c(0, 0)
  track$player = NULL
  poses = animation_evaluate_transform(track, interval)
  gap = max(abs(poses[[1]] - poses[[2]]))
  list(
    interval = interval,
    matrix_gap = gap,
    pose_continuous = gap <= tolerance
  )
}
