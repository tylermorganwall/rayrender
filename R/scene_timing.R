#' Resolve scene timing
#' @param time Requested scene time or NULL for legacy timing.
#' @param timeline Requested playback range or NULL.
#' @param fps Positive output/display frame rate.
#' @param shutter_offsets Exposure offsets in seconds or NULL.
#' @param legacy_supplied Whether competing legacy shutter arguments were supplied.
#' @return Validated explicit timing or NULL for legacy timing.
#' @keywords internal
resolve_scene_timing = function(
  time,
  timeline,
  fps,
  shutter_offsets,
  legacy_supplied
) {
  if (!is.numeric(fps) || length(fps) != 1L || !is.finite(fps) || fps <= 0) {
    stop("fps must be a positive finite number.", call. = FALSE)
  }
  if (is.null(time)) {
    if (!is.null(timeline) || !is.null(shutter_offsets)) {
      stop(
        "timeline and shutter_offsets require an explicit time.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  if (!is.numeric(time) || length(time) != 1L || !is.finite(time)) {
    stop("time must be one finite scene time in seconds.", call. = FALSE)
  }
  if (legacy_supplied) {
    stop(
      "Explicit time cannot be combined with legacy shutter controls.",
      call. = FALSE
    )
  }
  if (is.null(shutter_offsets)) {
    shutter_offsets = c(0, 0)
  }
  if (
    !is.numeric(shutter_offsets) ||
      length(shutter_offsets) != 2L ||
      any(!is.finite(shutter_offsets)) ||
      shutter_offsets[1] > shutter_offsets[2] ||
      any(abs(shutter_offsets) > 1e6) ||
      any(!is.finite(time + shutter_offsets))
  ) {
    stop(
      "shutter_offsets must contain two finite ordered seconds within +/-1e6.",
      call. = FALSE
    )
  }
  if (
    !is.null(timeline) &&
      (!is.numeric(timeline) ||
        length(timeline) != 2L ||
        any(!is.finite(timeline)) ||
        timeline[1] > timeline[2] ||
        time < timeline[1] ||
        time > timeline[2])
  ) {
    stop(
      "timeline must be a finite ordered range containing time.",
      call. = FALSE
    )
  }
  validate_relative_scene_times(shutter_offsets)
  list(time = time, offsets = shutter_offsets, timeline = timeline, fps = fps)
}

#' Validate frame-relative float precision
#' @param values Finite scene-relative times in seconds.
#' @return Invisibly the float-representable values in input order.
#' @keywords internal
validate_relative_scene_times = function(values) {
  packed = readBin(
    writeBin(as.double(values), raw(), size = 4),
    "double",
    n = length(values),
    size = 4
  )
  ordered = order(values)
  if (
    any(!is.finite(packed)) ||
      any(abs(packed - values) > 1e-5) ||
      any(diff(values[ordered]) > 0 & diff(packed[ordered]) <= 0)
  ) {
    stop(
      "Scene-relative times exceed float precision; use a shorter time interval.",
      call. = FALSE
    )
  }
  invisible(packed)
}

#' Detect keyframe attachments
#' @param scene Scene or processed instance prototype.
#' @return Whether the scene contains a keyframe attachment.
#' @keywords internal
scene_has_keyframes = function(scene) {
  if (
    any(vapply(
      scene$animation_info,
      function(info) !is.null(info$keyframes),
      logical(1)
    ))
  ) {
    return(TRUE)
  }
  for (info in scene$shape_info) {
    if (
      !is.null(info$original_scene) &&
        scene_has_keyframes(info$original_scene[[1]])
    ) {
      return(TRUE)
    }
  }
  FALSE
}

#' Sample a camera at explicit scene times
#' @param camera Authored camera, with optional time column in its motion.
#' @param times Finite requested output times in seconds.
#' @param fps Rate used to timestamp camera rows without a time column.
#' @param shutter_offsets Ordered exposure offsets in seconds.
#' @return Camera sampled at the requested times with a shared double epoch.
#' @keywords internal
camera_at_scene_times = function(camera, times, fps, shutter_offsets) {
  motion = as.data.frame(
    if (!is.null(camera$authored_motion)) {
      camera$authored_motion
    } else {
      camera$motion
    }
  )
  # A lens descriptor overrides the row's projection in both render APIs.
  # Preserve that choice in the authored track as well as sampled rows.
  if (!is.na(camera$camera_description_file)) {
    motion$fov = -1
  }
  key_times = if ("time" %in% names(motion)) {
    motion$time
  } else {
    (seq_len(nrow(motion)) - 1) / fps
  }
  if (
    !is.numeric(key_times) ||
      any(!is.finite(key_times)) ||
      any(diff(key_times) <= 0)
  ) {
    stop(
      "Camera key times must be finite and strictly increasing.",
      call. = FALSE
    )
  }
  if (
    nrow(motion) > 1L &&
      any(shutter_offsets != 0) &&
      any(motion$fov < 0) &&
      (any(motion$aperture != motion$aperture[1]) ||
        any(motion$focal != motion$focal[1]))
  ) {
    stop(
      "Physical-camera exposure requires constant aperture and focus; use offline zero-width snapshots for varying optics.",
      call. = FALSE
    )
  }
  # Keep the authored rows for native timeline sampling. Copies own ordinary R
  # data; the sampled output rows must not replace the source of later seeks.
  camera$authored_motion = motion
  camera$authored_motion$time = key_times
  # Projection changes are cuts. Within one projection, interpolate camera
  # parameters against authored time, never against a previously shown pose.
  index = pmax(1L, pmin(nrow(motion), findInterval(times, key_times)))
  sampled = motion[index, , drop = FALSE]
  next_index = pmin(index + 1L, nrow(motion))
  fraction = ifelse(
    index == next_index,
    0,
    (times - key_times[index]) / (key_times[next_index] - key_times[index])
  )
  fraction = pmax(0, pmin(1, fraction))
  type = function(fov) {
    ifelse(fov < 0, -1L, ifelse(fov == 0, 0L, ifelse(fov == 360, 360L, 1L)))
  }
  fraction[type(motion$fov[index]) != type(motion$fov[next_index])] = 0
  for (name in ray_camera_motion_columns()) {
    active = which(fraction > 0)
    sampled[[name]][active] = motion[[name]][index[active]] *
      (1 - fraction[active]) +
      motion[[name]][next_index[active]] * fraction[active]
  }
  epoch = min(times)
  offsets = times - epoch
  # The native camera stores relative seconds as floats. Refuse a requested
  # output sequence if this representation changes a time by over ten us or
  # collapses distinct requested times; the authored track times remain doubles.
  validate_relative_scene_times(offsets)
  open = validate_relative_scene_times(offsets + shutter_offsets[1])
  close = validate_relative_scene_times(offsets + shutter_offsets[2])
  if (shutter_offsets[1] < shutter_offsets[2] && any(open >= close)) {
    stop(
      "The exposure collapses at frame-relative float precision.",
      call. = FALSE
    )
  }
  sampled$scene_offset = offsets
  sampled$time = times
  camera$motion = as_ray_camera_motion(sampled)
  camera$scene_time = epoch
  camera$shutter_offsets = shutter_offsets
  camera$shutteropen = min(offsets) + shutter_offsets[1]
  camera$shutterclose = max(offsets) + shutter_offsets[2]
  camera$shutter_speed = 1
  camera$camera_motion_blur = FALSE
  camera
}
