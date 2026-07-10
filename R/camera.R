#' Camera
#'
#' Creates a rayrender camera that can be attached to a `ray_scene` with
#' `add_camera()`.
#'
#' @param lookfrom Default `c(0, 1, -10)`. Location of the camera.
#' @param lookat Default `c(0, 0, 0)`. Location where the camera is pointed.
#' @param camera_up Default `c(0, 1, 0)`. Vector indicating the up direction of the camera.
#' @param fov Default `20`. Field of view, in degrees.
#' @param aperture Default `0.1`. Aperture of the camera.
#' @param focal_distance Default `NULL`. Focal distance. If `NULL`, this is the distance between `lookfrom` and `lookat`.
#' @param ortho_dimensions Default `c(1, 1)`. Width and height of the orthographic camera when `fov = 0`.
#' @param motion Default `NULL`. Camera motion data frame from `generate_camera_motion()`.
#' @param keyframe_motion_args Default `list()`. Named list of additional arguments passed to
#' `generate_camera_motion()` when pressing `M` in interactive preview to preview the saved keyframes.
#' The saved keyframes always supply the camera positions. Defaults are `type = "linear"`,
#' 30 frames per saved keyframe, and `damp_motion = TRUE`.
#' @param name Default `"camera"`. Camera name.
#' @param filename Default `NA_character_`. Optional output filename or animation filename pattern.
#' @param camera_description_file Default `NA`. Filename of a realistic camera description file.
#' @param camera_scale Default `1`. Amount to scale a realistic camera.
#' @param iso Default `100`. Camera exposure.
#' @param film_size Default `22`. Film size in millimeters for realistic cameras.
#' @param shutteropen Default `0`. Time at which the shutter opens.
#' @param shutterclose Default `1`. Time at which the shutter closes.
#'
#' @return A `ray_camera` object.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # Static camera attached to the scene, equivalent to passing lookfrom/lookat
#' # directly to render_scene().
#' scene = generate_ground(depth = -0.5, material = diffuse(checkercolor = "blue")) |>
#'   add_object(sphere(y = 0.5, radius = 0.5, material = diffuse(color = "red"))) |>
#'   add_camera(camera(
#'     name = "main",
#'     lookfrom = c(7, 1.5, 10),
#'     lookat = c(0, 0.5, 0),
#'     fov = 15,
#'     filename = NA_character_
#'   ))
#' render_scene(scene, samples = 16, parallel = TRUE)
#'
#' # Animated camera attached to the scene, equivalent to passing camera_motion
#' # directly to render_animation().
#' camera_pos = list(c(0, 1, 15), c(5, -5, 5), c(-5, 5, -5), c(0, 1, -15))
#' camera_motion = generate_camera_motion(
#'   positions = camera_pos,
#'   lookats = camera_pos,
#'   offset_lookat = 1,
#'   fovs = 80,
#'   frames = 12,
#'   type = "bezier"
#' )
#' animated_scene = generate_ground(material = diffuse(checkercolor = "grey20"), depth = -10) |>
#'   add_object(sphere(y = 50, radius = 10, material = light(intensity = 30))) |>
#'   add_object(path(camera_pos, y = -0.2, material = diffuse(color = "red"))) |>
#'   add_camera(camera(
#'     name = "flythrough",
#'     motion = camera_motion,
#'     filename = NA_character_
#'   ))
#' render_scene(
#'   animated_scene,
#'   camera = "flythrough",
#'   mode = "animation",
#'   samples = 16,
#'   sample_method = "sobol_blue",
#'   clamp_value = 10,
#'   width = 400,
#'   height = 400
#' )
camera = function(
  lookfrom = c(0, 1, -10),
  lookat = c(0, 0, 0),
  camera_up = c(0, 1, 0),
  fov = 20,
  aperture = 0.1,
  focal_distance = NULL,
  ortho_dimensions = c(1, 1),
  motion = NULL,
  keyframe_motion_args = list(),
  name = "camera",
  filename = NA_character_,
  camera_description_file = NA,
  camera_scale = 1,
  iso = 100,
  film_size = 22,
  shutteropen = 0,
  shutterclose = 1
) {
  validate_camera_name(name)
  keyframe_motion_args = normalize_keyframe_motion_args(keyframe_motion_args)

  if (is.null(motion)) {
    validate_static_camera_inputs(
      lookfrom = lookfrom,
      lookat = lookat,
      camera_up = camera_up,
      fov = fov,
      aperture = aperture,
      focal_distance = focal_distance,
      ortho_dimensions = ortho_dimensions
    )
    if (is.null(focal_distance)) {
      focal_distance = sqrt(sum((lookfrom - lookat)^2))
    }
    motion = data.frame(
      x = lookfrom[1],
      y = lookfrom[2],
      z = lookfrom[3],
      dx = lookat[1],
      dy = lookat[2],
      dz = lookat[3],
      aperture = aperture,
      fov = fov,
      focal = focal_distance,
      orthox = ortho_dimensions[1],
      orthoy = ortho_dimensions[2],
      upx = camera_up[1],
      upy = camera_up[2],
      upz = camera_up[3]
    )
  } else {
    motion = as.data.frame(motion)
    validate_camera_motion(motion)
  }

  motion = as_ray_camera_motion(motion)
  validate_camera_filename(filename, nrow(motion))

  structure(
    list(
      name = name,
      motion = motion,
      keyframe_motion_args = keyframe_motion_args,
      filename = filename,
      camera_description_file = camera_description_file,
      camera_scale = camera_scale,
      iso = iso,
      film_size = film_size,
      shutteropen = shutteropen,
      shutterclose = shutterclose
    ),
    class = "ray_camera"
  )
}

#' @export
print.ray_camera = function(x, ...) {
  format_vec = function(values) {
    formatted = format(
      round(as.numeric(values), 3),
      trim = TRUE,
      scientific = FALSE
    )
    sprintf("c(%s)", paste(formatted, collapse = ", "))
  }
  format_range = function(values) {
    values = as.numeric(values)
    if (length(unique(values)) == 1) {
      return(format(values[1], trim = TRUE, scientific = FALSE))
    }
    sprintf(
      "%s to %s",
      format(min(values), trim = TRUE, scientific = FALSE),
      format(max(values), trim = TRUE, scientific = FALSE)
    )
  }
  projection_type = function(fov) {
    types = ifelse(
      fov < 0,
      "realistic",
      ifelse(
        fov == 0,
        "orthographic",
        ifelse(fov == 360, "environment", "perspective")
      )
    )
    types = unique(types)
    if (length(types) == 1) {
      return(types)
    }
    paste("mixed", paste(types, collapse = ", "), sep = ": ")
  }
  format_filename = function(filename) {
    if (length(filename) == 1 && is.na(filename)) {
      return("not set")
    }
    if (length(filename) == 1) {
      return(filename)
    }
    sprintf("%d filenames", length(filename))
  }

  motion = as.data.frame(x$motion)
  frames = nrow(motion)
  first_frame = motion[1, , drop = FALSE]
  last_frame = motion[frames, , drop = FALSE]
  animated = frames > 1

  cat(sprintf("ray_camera <%s>\n", x$name))
  cat(sprintf(
    "  type: %s %s\n",
    projection_type(motion$fov),
    if (animated) "animation" else "static"
  ))
  cat(sprintf("  frames: %d\n", frames))
  cat(sprintf(
    "  lookfrom: %s\n",
    format_vec(c(first_frame$x, first_frame$y, first_frame$z))
  ))
  cat(sprintf(
    "  lookat: %s\n",
    format_vec(c(first_frame$dx, first_frame$dy, first_frame$dz))
  ))
  if (animated) {
    cat(sprintf(
      "  last lookfrom: %s\n",
      format_vec(c(last_frame$x, last_frame$y, last_frame$z))
    ))
    cat(sprintf(
      "  last lookat: %s\n",
      format_vec(c(last_frame$dx, last_frame$dy, last_frame$dz))
    ))
  }
  cat(sprintf("  fov: %s\n", format_range(motion$fov)))
  cat(sprintf("  aperture: %s\n", format_range(motion$aperture)))
  cat(sprintf("  focal distance: %s\n", format_range(motion$focal)))
  if (!is.na(x$camera_description_file)) {
    cat(sprintf("  camera description: %s\n", x$camera_description_file))
  }
  cat(sprintf("  output: %s\n", format_filename(x$filename)))

  invisible(x)
}

#' Add Camera
#'
#' Adds a camera to a `ray_scene`.
#'
#' @param scene Scene to modify.
#' @param camera Camera created with `camera()`.
#' @param name Default `NULL`. Optional name that overrides `camera$name`.
#' @param active Default `TRUE`. Whether to set this camera as the active scene camera.
#' @param replace Default `FALSE`. Whether to replace an existing camera with the same name.
#'
#' @return A modified `ray_scene`.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' scene = generate_ground(material=diffuse(color="grey20")) |>
#'   add_object(sphere()) |>
#'   add_camera(camera(
#'     name = "main",
#'     lookfrom = c(0, 1, -10),
#'     lookat = c(0, 0, 0),
#'     fov = 35
#'   ))
#'
#' render_scene(scene, samples = 16)
add_camera = function(
  scene,
  camera,
  name = NULL,
  active = TRUE,
  replace = FALSE
) {
  if (!inherits(camera, "ray_camera")) {
    stop("camera must inherit from class 'ray_camera'")
  }
  if (!is.null(name)) {
    validate_camera_name(name)
    camera$name = name
  } else {
    validate_camera_name(camera$name)
  }

  cameras = ray_scene_cameras(scene)
  if (!replace && camera$name %in% names(cameras)) {
    stop(
      "Camera name '",
      camera$name,
      "' already exists. Use replace = TRUE to replace it."
    )
  }
  cameras[[camera$name]] = camera
  names(cameras)[names(cameras) == ""] = camera$name

  attr(scene, "ray_cameras") = cameras
  if (isTRUE(active)) {
    attr(scene, "active_camera") = camera$name
  }
  scene
}

#' Get Camera
#'
#' Gets a camera from a `ray_scene`.
#'
#' @param scene Scene containing cameras.
#' @param camera Default `NULL`. Camera name, `ray_camera` object, or `NULL` to use the active camera.
#'
#' @return A `ray_camera` object.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' scene = generate_ground(material=diffuse(color="grey20")) |>
#'   add_camera(camera(name = "wide", fov = 55), active = FALSE) |>
#'   add_camera(camera(name = "main", fov = 35), active = TRUE)
#'
#' # Returns the active camera.
#' get_camera(scene)
#'
#' # Returns a named camera.
#' get_camera(scene, "wide")
get_camera = function(scene, camera = NULL) {
  if (inherits(camera, "ray_camera")) {
    return(camera)
  }

  cameras = ray_scene_cameras(scene)
  camera_names = names(cameras)
  if (length(cameras) == 0) {
    stop("No cameras are attached to this scene.")
  }

  if (is.null(camera)) {
    active_camera = attr(scene, "active_camera")
    if (!is.null(active_camera) && active_camera %in% camera_names) {
      return(cameras[[active_camera]])
    }
    if (length(cameras) == 1) {
      return(cameras[[1]])
    }
    stop(
      "Multiple cameras are attached and no active camera is set. Available cameras: ",
      paste(camera_names, collapse = ", ")
    )
  }

  validate_camera_name(camera)
  if (!camera %in% camera_names) {
    stop(
      "Camera '",
      camera,
      "' was not found. Available cameras: ",
      paste(camera_names, collapse = ", ")
    )
  }
  cameras[[camera]]
}

#' List Cameras
#'
#' Lists cameras attached to a `ray_scene`.
#'
#' @param scene Scene containing cameras.
#'
#' @return A data frame with one row per camera.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' scene = generate_ground(material=diffuse(color="grey20")) |>
#'   add_camera(camera(name = "wide", fov = 55, filename = "wide.png")) |>
#'   add_camera(camera(name = "detail", fov = 20, filename = "detail.png"))
#'
#' list_cameras(scene)
list_cameras = function(scene) {
  cameras = ray_scene_cameras(scene)
  active_camera = attr(scene, "active_camera")
  if (length(cameras) == 0) {
    return(data.frame(
      name = character(),
      frames = integer(),
      animated = logical(),
      filename = character(),
      active = logical()
    ))
  }
  data.frame(
    name = names(cameras),
    frames = vapply(cameras, function(cam) nrow(cam$motion), integer(1)),
    animated = vapply(cameras, function(cam) nrow(cam$motion) > 1, logical(1)),
    filename = vapply(
      cameras,
      function(cam) paste(cam$filename, collapse = ","),
      character(1)
    ),
    active = names(cameras) == active_camera,
    row.names = NULL
  )
}

#' Remove Camera
#'
#' Removes a camera from a `ray_scene`.
#'
#' @param scene Scene containing cameras.
#' @param camera Camera name.
#'
#' @return A modified `ray_scene`.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' scene = generate_ground(material=diffuse(color="grey20")) |>
#'   add_camera(camera(name = "wide"), active = FALSE) |>
#'   add_camera(camera(name = "detail"), active = TRUE)
#'
#' scene = remove_camera(scene, "detail")
#' list_cameras(scene)
remove_camera = function(scene, camera) {
  validate_camera_name(camera)
  cameras = ray_scene_cameras(scene)
  if (!camera %in% names(cameras)) {
    stop("Camera '", camera, "' was not found.")
  }
  cameras[[camera]] = NULL
  attr(scene, "ray_cameras") = cameras

  active_camera = attr(scene, "active_camera")
  if (!is.null(active_camera) && identical(active_camera, camera)) {
    if (length(cameras) == 1) {
      attr(scene, "active_camera") = names(cameras)
    } else {
      attr(scene, "active_camera") = NULL
    }
  }
  scene
}

#' Set Active Camera
#'
#' Sets the active camera for a `ray_scene`.
#'
#' @param scene Scene containing cameras.
#' @param camera Camera name.
#'
#' @return A modified `ray_scene`.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' scene = generate_ground(material=diffuse(color="grey20")) |>
#'   add_camera(camera(name = "wide", fov = 55), active = FALSE) |>
#'   add_camera(camera(name = "detail", fov = 20), active = FALSE)
#'
#' scene = set_active_camera(scene, "detail")
#' get_camera(scene)
set_active_camera = function(scene, camera) {
  validate_camera_name(camera)
  cameras = ray_scene_cameras(scene)
  if (!camera %in% names(cameras)) {
    stop("Camera '", camera, "' was not found.")
  }
  attr(scene, "active_camera") = camera
  scene
}

#' Preview Camera
#'
#' Previews a scene-attached camera without writing image files.
#'
#' @param scene Scene containing cameras.
#' @param camera Default `NULL`. Camera name or `ray_camera` object to preview.
#' @param width Default `800`. Preview width, in pixels.
#' @param height Default `800`. Preview height, in pixels.
#' @param fps Default `24`. Intended preview frame rate for animated cameras.
#' @param samples Default `1`. Number of samples per pixel for the preview.
#' @param ... Additional arguments passed to `render_scene()`.
#'
#' @return Invisibly returns the rendered preview result.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' #Zooming over the R logo
#' motion = generate_camera_motion(
#'   positions = list(c(0, 1, -10), c(0, 1, -2), c(0, 1, 10)),
#'   lookats = list(c(0, 0, 0), c(0, 0.5, 0), c(0, 0, 0)),
#'   fovs = c(60, 45, 90),
#'   frames = 24,
#'   type = "linear",
#'   damp_motion = TRUE,
#'   closed = TRUE
#' )
#'
#' scene = generate_ground(material=diffuse(color="grey20")) |>
#'   add_object(obj_model(r_obj())) |>
#'   add_object(sphere(y=10,x=-10,z=-5,material=light(intensity=100))) |>
#'   add_camera(camera(name = "flythrough", motion = motion))
#'
#' preview_camera(scene, camera = "flythrough", width = 400, height = 400)
preview_camera = function(
  scene,
  camera = NULL,
  width = 800,
  height = 800,
  fps = 24,
  samples = 1,
  ...
) {
  invisible(fps)
  render_scene(
    scene = scene,
    width = width,
    height = height,
    camera = camera,
    mode = "preview",
    samples = samples,
    ...
  )
}

#' @keywords internal
ray_camera_motion_columns = function() {
  c(
    "x",
    "y",
    "z",
    "dx",
    "dy",
    "dz",
    "aperture",
    "fov",
    "focal",
    "orthox",
    "orthoy",
    "upx",
    "upy",
    "upz"
  )
}

#' @keywords internal
as_ray_camera_motion = function(motion) {
  motion = as.data.frame(motion)
  class(motion) = unique(c("ray_camera_motion", "data.frame"))
  motion
}

#' @keywords internal
validate_camera_motion = function(motion) {
  missing_columns = setdiff(ray_camera_motion_columns(), colnames(motion))
  if (length(missing_columns) > 0) {
    stop(
      "motion must contain camera-motion columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  invisible(TRUE)
}

#' @keywords internal
normalize_keyframe_motion_args = function(args = list()) {
  if (is.null(args)) {
    args = list()
  }
  if (!is.list(args)) {
    stop("keyframe_motion_args must be a named list.")
  }
  if (length(args) > 0 && (is.null(names(args)) || any(!nzchar(names(args))))) {
    stop("keyframe_motion_args must be a named list.")
  }

  keyframe_supplied_args = c(
    "positions",
    "lookats",
    "apertures",
    "fovs",
    "focal_distances",
    "ortho_dims",
    "camera_ups"
  )
  conflicting_args = intersect(names(args), keyframe_supplied_args)
  if (length(conflicting_args) > 0) {
    stop(
      "keyframe_motion_args cannot include arguments supplied by saved keyframes: ",
      paste(conflicting_args, collapse = ", ")
    )
  }

  utils::modifyList(
    list(type = "linear", damp_motion = TRUE),
    args
  )
}

#' @keywords internal
validate_camera_name = function(name) {
  if (
    !is.character(name) ||
      length(name) != 1 ||
      is.na(name) ||
      !nzchar(name)
  ) {
    stop("Camera names must be nonempty scalar character strings.")
  }
  invisible(TRUE)
}

#' @keywords internal
validate_static_camera_inputs = function(
  lookfrom,
  lookat,
  camera_up,
  fov,
  aperture,
  focal_distance,
  ortho_dimensions
) {
  validate_numeric_vector(lookfrom, 3, "lookfrom")
  validate_numeric_vector(lookat, 3, "lookat")
  validate_numeric_vector(camera_up, 3, "camera_up")
  if (all(lookfrom == lookat)) {
    stop("lookfrom and lookat must not be identical.")
  }
  lookvec = lookat - lookfrom
  camera_cross = c(
    lookvec[2] * camera_up[3] - lookvec[3] * camera_up[2],
    lookvec[3] * camera_up[1] - lookvec[1] * camera_up[3],
    lookvec[1] * camera_up[2] - lookvec[2] * camera_up[1]
  )
  if (all(camera_cross == 0)) {
    stop("camera_up must not be exactly collinear with lookat - lookfrom.")
  }
  validate_numeric_scalar(aperture, "aperture")
  if (aperture < 0) {
    stop("aperture must be nonnegative.")
  }
  validate_numeric_scalar(fov, "fov")
  if (!is.null(focal_distance)) {
    validate_numeric_scalar(focal_distance, "focal_distance")
    if (focal_distance <= 0) {
      stop("focal_distance must be positive if supplied.")
    }
  }
  validate_numeric_vector(ortho_dimensions, 2, "ortho_dimensions")
  if (any(ortho_dimensions <= 0)) {
    stop("ortho_dimensions must contain positive values.")
  }
  invisible(TRUE)
}

#' @keywords internal
validate_numeric_vector = function(x, len, arg) {
  if (
    !is.numeric(x) ||
      length(x) != len ||
      any(is.na(x)) ||
      any(!is.finite(x))
  ) {
    stop(arg, " must be a numeric length-", len, " vector.")
  }
  invisible(TRUE)
}

#' @keywords internal
validate_numeric_scalar = function(x, arg) {
  if (
    !is.numeric(x) ||
      length(x) != 1 ||
      is.na(x) ||
      !is.finite(x)
  ) {
    stop(arg, " must be a finite numeric scalar.")
  }
  invisible(TRUE)
}

#' @keywords internal
validate_camera_filename = function(filename, n_frames) {
  if (length(filename) == 1 && is.na(filename)) {
    return(invisible(TRUE))
  }
  if (!is.character(filename)) {
    stop("filename must be NA, a scalar character, or a character vector.")
  }
  if (!length(filename) %in% c(1, n_frames)) {
    stop(
      "filename must have length 1 or the number of camera frames (",
      n_frames,
      ")."
    )
  }
  if (any(is.na(filename))) {
    stop("filename values must not be NA unless filename is a scalar NA.")
  }
  invisible(TRUE)
}

#' @keywords internal
ray_scene_cameras = function(scene) {
  cameras = attr(scene, "ray_cameras")
  if (is.null(cameras)) {
    return(list())
  }
  if (!is.list(cameras)) {
    stop("ray_cameras scene attribute must be a named list.")
  }
  cameras
}

#' @keywords internal
preserve_ray_scene_attrs = function(newscene, scene, objects = NULL) {
  if (!is.null(attr(scene, "cornell")) || !is.null(attr(objects, "cornell"))) {
    attr(newscene, "cornell") = TRUE
  }

  scene_cameras = ray_scene_cameras(scene)
  object_cameras = ray_scene_cameras(objects)
  duplicate_cameras = intersect(names(scene_cameras), names(object_cameras))
  if (length(duplicate_cameras) > 0) {
    stop(
      "Duplicate camera names found while combining scenes: ",
      paste(duplicate_cameras, collapse = ", ")
    )
  }
  merged_cameras = c(scene_cameras, object_cameras)
  if (length(merged_cameras) > 0) {
    attr(newscene, "ray_cameras") = merged_cameras
  }

  scene_active = attr(scene, "active_camera")
  object_active = attr(objects, "active_camera")
  if (!is.null(scene_active) && scene_active %in% names(merged_cameras)) {
    attr(newscene, "active_camera") = scene_active
  } else if (
    !is.null(object_active) &&
      object_active %in% names(merged_cameras)
  ) {
    attr(newscene, "active_camera") = object_active
  }

  newscene
}

#' @keywords internal
camera_has_frame_marker = function(filename) {
  grepl("(^|[^%])%[0-9]*[di]", filename, perl = TRUE)
}

#' @keywords internal
append_png_if_missing_extension = function(filename) {
  missing_extension = tools::file_ext(filename) == ""
  filename[missing_extension] = paste0(filename[missing_extension], ".png")
  filename
}

#' @keywords internal
camera_frame_filenames = function(
  filename,
  n_frames,
  frame_numbers = seq_len(n_frames)
) {
  validate_camera_filename(filename, n_frames)

  if (length(filename) == 1 && is.na(filename)) {
    return(list(
      filenames = rep("", n_frames),
      write_image = FALSE
    ))
  }

  if (length(filename) == n_frames && n_frames != 1) {
    return(list(
      filenames = append_png_if_missing_extension(filename),
      write_image = TRUE
    ))
  }

  if (camera_has_frame_marker(filename)) {
    filenames = tryCatch(
      sprintf(filename, frame_numbers),
      error = function(e) {
        stop("Could not expand filename with frame numbers: ", e$message)
      }
    )
    return(list(
      filenames = append_png_if_missing_extension(filenames),
      write_image = TRUE
    ))
  }

  extension = tools::file_ext(filename)
  if (extension == "") {
    filenames = paste0(filename, frame_numbers, ".png")
  } else {
    base = tools::file_path_sans_ext(filename)
    filenames = paste0(base, frame_numbers, ".", extension)
  }
  list(
    filenames = filenames,
    write_image = TRUE
  )
}

#' @keywords internal
camera_static_filename = function(filename) {
  validate_camera_filename(filename, 1)
  if (length(filename) == 1 && is.na(filename)) {
    return(NA_character_)
  }
  if (tools::file_ext(filename) == "") {
    return(paste0(filename, ".png"))
  }
  filename
}

#' @keywords internal
camera_frame_range = function(n_frames, start_frame = 1, end_frame = NA) {
  validate_numeric_scalar(start_frame, "start_frame")
  if (start_frame != as.integer(start_frame)) {
    stop("start_frame must be an integer.")
  }
  start_frame = as.integer(start_frame)
  if (start_frame < 1 || start_frame > n_frames) {
    stop("start_frame must be between 1 and the number of camera frames.")
  }

  if (length(end_frame) == 1 && is.na(end_frame)) {
    end_frame = n_frames
  } else {
    validate_numeric_scalar(end_frame, "end_frame")
    if (end_frame != as.integer(end_frame)) {
      stop("end_frame must be an integer.")
    }
    end_frame = as.integer(end_frame)
  }
  if (end_frame < start_frame || end_frame > n_frames) {
    stop(
      "end_frame must be between start_frame and the number of camera frames."
    )
  }
  start_frame:end_frame
}

#' @keywords internal
camera_render_mode = function(camera, mode = "auto") {
  if (mode == "auto") {
    if (nrow(camera$motion) == 1) {
      return("image")
    }
    return("animation")
  }
  mode
}

#' @keywords internal
camera_image_filename = function(
  camera,
  start_frame = 1,
  filename_override = NULL,
  filename_supplied = FALSE
) {
  filename = if (isTRUE(filename_supplied)) {
    filename_override
  } else {
    camera$filename
  }
  if (nrow(camera$motion) == 1) {
    return(camera_static_filename(filename))
  }
  frame_files = camera_frame_filenames(
    filename,
    nrow(camera$motion),
    seq_len(nrow(camera$motion))
  )
  if (!frame_files$write_image) {
    return(NA_character_)
  }
  frame_files$filenames[start_frame]
}

#' @keywords internal
camera_written_filenames = function(
  camera,
  mode = "auto",
  start_frame = 1,
  end_frame = NA,
  filename_override = NULL,
  filename_supplied = FALSE
) {
  mode = camera_render_mode(camera, mode)
  if (mode == "preview") {
    return(character())
  }

  filename = if (isTRUE(filename_supplied)) {
    filename_override
  } else {
    camera$filename
  }

  if (mode == "image") {
    frame_range = camera_frame_range(
      nrow(camera$motion),
      start_frame,
      start_frame
    )
    filename = camera_image_filename(
      camera,
      start_frame = frame_range[1],
      filename_override = filename_override,
      filename_supplied = filename_supplied
    )
    if (is.na(filename)) {
      return(character())
    }
    return(filename)
  }

  frame_range = camera_frame_range(nrow(camera$motion), start_frame, end_frame)
  filenames = camera_frame_filenames(
    filename,
    nrow(camera$motion),
    seq_len(nrow(camera$motion))
  )
  if (!filenames$write_image) {
    return(character())
  }
  filenames$filenames[frame_range]
}

#' @keywords internal
validate_camera_output_filenames = function(
  cameras,
  mode = "auto",
  start_frame = 1,
  end_frame = NA,
  filename_override = NULL,
  filename_supplied = FALSE
) {
  filenames = unlist(lapply(
    cameras,
    camera_written_filenames,
    mode = mode,
    start_frame = start_frame,
    end_frame = end_frame,
    filename_override = filename_override,
    filename_supplied = filename_supplied
  ))
  filenames = filenames[nzchar(filenames)]
  duplicate_filenames = unique(filenames[duplicated(filenames)])
  if (length(duplicate_filenames) > 0) {
    stop(
      "Duplicate camera output filenames detected: ",
      paste(duplicate_filenames, collapse = ", ")
    )
  }
  invisible(filenames)
}

#' @keywords internal
is_camera_arg_supplied = function(supplied, arg) {
  arg %in% names(supplied) && isTRUE(supplied[[arg]])
}

#' @keywords internal
apply_camera_overrides = function(camera, overrides, supplied) {
  override_names = intersect(names(overrides), names(supplied))
  for (override_name in override_names) {
    if (isTRUE(supplied[[override_name]])) {
      camera[[override_name]] = overrides[[override_name]]
    }
  }
  validate_camera_filename(camera$filename, nrow(camera$motion))
  camera
}

#' @keywords internal
resolve_scene_camera = function(
  scene,
  camera = NULL,
  legacy_camera = NULL,
  legacy_camera_supplied = FALSE,
  allow_all = FALSE,
  default_camera = NULL,
  warn_legacy_override = TRUE
) {
  cameras = ray_scene_cameras(scene)

  if (isTRUE(legacy_camera_supplied)) {
    if (
      isTRUE(warn_legacy_override) &&
        (length(cameras) > 0 || !is.null(camera))
    ) {
      warning(
        "Explicit render-time camera arguments override scene-attached cameras."
      )
    }
    return(list(legacy_camera))
  }

  if (!is.null(camera)) {
    if (inherits(camera, "ray_camera")) {
      return(list(camera))
    }
    validate_camera_name(camera)
    if (identical(camera, "all")) {
      if (!isTRUE(allow_all)) {
        stop("camera = 'all' is not supported here.")
      }
      if (length(cameras) == 0) {
        stop("camera = 'all' requires at least one scene-attached camera.")
      }
      return(cameras)
    }
    return(list(get_camera(scene, camera)))
  }

  if (length(cameras) > 0) {
    return(list(get_camera(scene)))
  }

  if (!is.null(default_camera)) {
    return(list(default_camera))
  }

  stop("No camera could be resolved for this scene.")
}

#' @keywords internal
render_scene_legacy_camera = function(
  scene,
  supplied,
  lookfrom,
  lookat,
  camera_up,
  fov,
  aperture,
  focal_distance,
  ortho_dimensions,
  filename,
  camera_description_file,
  camera_scale,
  iso,
  film_size,
  shutteropen,
  shutterclose,
  message_cornell = TRUE
) {
  if (!is.null(attr(scene, "cornell"))) {
    corn_message = "Setting default values for Cornell box: "
    missing_corn = FALSE
    if (!is_camera_arg_supplied(supplied, "lookfrom")) {
      lookfrom = c(278, 278, -800)
      corn_message = paste0(corn_message, "lookfrom `c(278,278,-800)` ")
      missing_corn = TRUE
    }
    if (!is_camera_arg_supplied(supplied, "lookat")) {
      lookat = c(278, 278, 555 / 2)
      corn_message = paste0(corn_message, "lookat `c(278,278,555/2)` ")
      missing_corn = TRUE
    }
    if (
      !is_camera_arg_supplied(supplied, "fov") &&
        is.na(camera_description_file)
    ) {
      fov = 40
      corn_message = paste0(corn_message, "fov `40` ")
      missing_corn = TRUE
    }
    if (
      fov == 0 &&
        !is_camera_arg_supplied(supplied, "ortho_dimensions") &&
        is.na(camera_description_file)
    ) {
      ortho_dimensions = c(580, 580)
      corn_message = paste0(corn_message, "ortho_dimensions `c(580, 580)` ")
      missing_corn = TRUE
    }
    corn_message = paste0(corn_message, ".")
    if (missing_corn && isTRUE(message_cornell)) {
      message(corn_message)
    }
  }

  camera(
    lookfrom = lookfrom,
    lookat = lookat,
    camera_up = camera_up,
    fov = fov,
    aperture = aperture,
    focal_distance = focal_distance,
    ortho_dimensions = ortho_dimensions,
    filename = filename,
    camera_description_file = camera_description_file,
    camera_scale = camera_scale,
    iso = iso,
    film_size = film_size,
    shutteropen = shutteropen,
    shutterclose = shutterclose
  )
}

#' @keywords internal
camera_frame_args = function(camera, frame = 1) {
  frame = camera_frame_range(nrow(camera$motion), frame, frame)
  motion = camera$motion[frame, , drop = FALSE]
  list(
    lookfrom = c(motion$x, motion$y, motion$z),
    lookat = c(motion$dx, motion$dy, motion$dz),
    camera_up = c(motion$upx, motion$upy, motion$upz),
    fov = motion$fov,
    aperture = motion$aperture,
    focal_distance = motion$focal,
    ortho_dimensions = c(motion$orthox, motion$orthoy),
    camera_description_file = camera$camera_description_file,
    camera_scale = camera$camera_scale,
    iso = camera$iso,
    film_size = camera$film_size,
    shutteropen = camera$shutteropen,
    shutterclose = camera$shutterclose
  )
}

#' @keywords internal
camera_batch_metadata = function(camera) {
  camera_description_file = camera$camera_description_file
  if (length(camera_description_file) == 1 && is.na(camera_description_file)) {
    camera_description_file = NA_character_
  } else {
    camera_description_file = as.character(camera_description_file)
  }

  list(
    camera_description_file = camera_description_file,
    camera_scale = as.numeric(camera$camera_scale),
    iso = as.numeric(camera$iso),
    film_size = as.numeric(camera$film_size),
    shutteropen = as.numeric(camera$shutteropen),
    shutterclose = as.numeric(camera$shutterclose)
  )
}

#' @keywords internal
camera_batch_metadata_compatible = function(cameras) {
  if (length(cameras) < 2) {
    return(TRUE)
  }
  reference = camera_batch_metadata(cameras[[1]])
  all(vapply(
    cameras[-1],
    function(camera) identical(camera_batch_metadata(camera), reference),
    logical(1)
  ))
}

#' @keywords internal
camera_mode_frame_range = function(
  camera,
  mode = "auto",
  start_frame = 1,
  end_frame = NA
) {
  render_mode = camera_render_mode(camera, mode)
  if (render_mode == "animation") {
    return(camera_frame_range(nrow(camera$motion), start_frame, end_frame))
  }
  if (render_mode == "preview" && nrow(camera$motion) > 1) {
    return(camera_frame_range(nrow(camera$motion), start_frame, end_frame))
  }
  camera_frame_range(nrow(camera$motion), start_frame, start_frame)
}

#' @keywords internal
camera_batch_frame_filenames = function(
  camera,
  mode = "auto",
  frame_range,
  force_no_write = FALSE
) {
  render_mode = camera_render_mode(camera, mode)
  if (isTRUE(force_no_write) || render_mode == "preview") {
    return(rep("", length(frame_range)))
  }

  if (render_mode == "image") {
    filename = camera_image_filename(camera, start_frame = frame_range[1])
    if (is.na(filename)) {
      return("")
    }
    return(filename)
  }

  filename_info = camera_frame_filenames(
    camera$filename,
    nrow(camera$motion),
    seq_len(nrow(camera$motion))
  )
  if (!filename_info$write_image) {
    return(rep("", length(frame_range)))
  }
  filename_info$filenames[frame_range]
}

#' @keywords internal
camera_batch_plan = function(
  cameras,
  mode = "auto",
  start_frame = 1,
  end_frame = NA,
  force_no_write = FALSE
) {
  if (!is.list(cameras) || length(cameras) == 0) {
    stop("cameras must be a nonempty list of ray_camera objects.")
  }

  motions = vector("list", length(cameras))
  filename_list = vector("list", length(cameras))
  camera_index = integer()
  frame_numbers = integer()
  camera_names = character(length(cameras))

  for (i in seq_along(cameras)) {
    camera = cameras[[i]]
    if (!inherits(camera, "ray_camera")) {
      stop("Each camera must inherit from class 'ray_camera'.")
    }
    validate_camera_name(camera$name)

    frame_range = camera_mode_frame_range(
      camera,
      mode = mode,
      start_frame = start_frame,
      end_frame = end_frame
    )
    motions[[i]] = as.data.frame(camera$motion)[frame_range, , drop = FALSE]
    filename_list[[i]] = camera_batch_frame_filenames(
      camera,
      mode = mode,
      frame_range = frame_range,
      force_no_write = force_no_write
    )
    camera_index = c(camera_index, rep(i, length(frame_range)))
    frame_numbers = c(frame_numbers, frame_range)
    camera_names[i] = camera$name
  }

  motion = do.call(rbind, motions)
  rownames(motion) = NULL
  filenames = unlist(filename_list, use.names = FALSE)

  list(
    motion = as_ray_camera_motion(motion),
    filenames = filenames,
    write_image = any(nzchar(filenames)),
    camera_index = camera_index,
    frame_numbers = frame_numbers,
    camera_names = camera_names
  )
}
