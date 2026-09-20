#' Apply saved native editor changes
#'
#' Restore committed transforms and material settings, including overrides within
#' individual instances. Changes are applied when the scene is rendered, leaving
#' the original geometry and materials available for subsequent editing.
#'
#' @param scene A rayrender scene.
#' @param edits Object changes returned in an image's `scene_edits` attribute or
#'   written by the native editor's Export R code button.
#' @param sky Default `NULL`. Committed sky settings written by the editor.
#' @return A scene carrying the supplied changes. Pass it to [render_scene()].
#' @export
apply_scene_edits = function(scene, edits, sky = NULL) {
  if (
    !is.data.frame(scene) || !is.list(edits) || (!is.null(sky) && !is.list(sky))
  ) {
    stop(
      "Supply a rayrender scene and a list of editor changes.",
      call. = FALSE
    )
  }
  # Native replay validates rows, instance placements, invertible transforms,
  # material layouts and values before tracing. Store value-only data here.
  attr(scene, "rayrender_scene_edits") = edits
  attr(scene, "rayrender_sky_edit") = sky
  scene
}

#' @keywords internal
native_export_render_args = function(args, state) {
  camera = state$camera
  args$lookfrom = unname(unlist(camera[c("x", "y", "z")]))
  args$lookat = unname(unlist(camera[c("dx", "dy", "dz")]))
  args$camera_up = unname(unlist(camera[c("upx", "upy", "upz")]))
  args$fov = camera$fov
  args$aperture = camera$aperture
  args$focal_distance = camera$focal
  args$ortho_dimensions = unname(unlist(camera[c("orthox", "orthoy")]))
  rotation = if (is.null(args$rotate_env)) 0 else args$rotate_env
  args$rotate_env = rotation + camera$env_rotation
  args$camera_motion_blur = state$camera_motion_blur
  args$shutter_speed = state$shutter_speed
  args$denoise = state$denoise
  # Freeze the calibrated viewport brightness. Re-running auto exposure against
  # a different number of samples could otherwise change the exported view.
  args$exposure = state$exposure
  args$auto_exposure = FALSE
  args$bloom = FALSE # Bloom is final-image processing, absent from the viewport.
  args$camera = NULL
  args$mode = "image"
  args$gui = "none"
  args$preview = FALSE
  args$interactive = FALSE
  args$deferred_render = FALSE
  args$filename = NA_character_
  args$plot_scene = TRUE
  args$integrator_type = state$integrator_type
  args
}

#' @keywords internal
native_export_assets = function(value, directory, relative, source_directory) {
  copied = new.env(parent = emptyenv())
  copy_asset = function(path) {
    path = normalizePath(path, mustWork = TRUE)
    if (exists(path, envir = copied, inherits = FALSE)) {
      return(get(path, envir = copied))
    }
    if (!dir.exists(directory)) {
      dir.create(directory)
    }
    name = paste0(
      length(ls(copied, all.names = TRUE)) + 1L,
      "-",
      basename(path)
    )
    target = file.path(directory, name)
    result = file.path(relative, name)
    assign(path, result, envir = copied)
    extension = tolower(tools::file_ext(path))
    if (extension %in% c("obj", "mtl")) {
      lines = readLines(path, warn = FALSE)
      resolve = function(value) {
        value = gsub('^["\']|["\']$', "", value)
        if (!grepl("^(/|[A-Za-z]:[/\\\\])", value)) {
          value = file.path(dirname(path), value)
        }
        value
      }
      # Rewrite companion references into this export's asset directory. Merely
      # copying a temporary OBJ would lose its MTL and texture files on restart.
      for (i in seq_along(lines)) {
        if (extension == "obj" && grepl("^\\s*mtllib\\s+", lines[i])) {
          value = sub("^\\s*mtllib\\s+", "", trimws(lines[i]))
          references = if (file.exists(resolve(value))) {
            value
          } else {
            strsplit(value, "\\s+")[[1L]]
          }
          paths = vapply(
            references,
            function(value) basename(copy_asset(resolve(value))),
            character(1)
          )
          lines[i] = paste("mtllib", paste(paths, collapse = " "))
        } else if (
          extension == "mtl" &&
            grepl("^\\s*(map_\\w+|bump|disp|decal|norm)\\s+", lines[i])
        ) {
          tokens = strsplit(trimws(lines[i]), "\\s+")[[1L]]
          matched = FALSE
          for (first in seq.int(2L, length(tokens))) {
            candidate = resolve(paste(
              tokens[first:length(tokens)],
              collapse = " "
            ))
            if (file.exists(candidate) && !dir.exists(candidate)) {
              lines[i] = paste(
                c(tokens[seq_len(first - 1L)], basename(copy_asset(candidate))),
                collapse = " "
              )
              matched = TRUE
              break
            }
          }
          if (!matched) {
            stop(
              "Missing texture referenced by temporary material file: ",
              path,
              call. = FALSE
            )
          }
        }
      }
      writeLines(lines, target, useBytes = TRUE)
    } else if (!file.copy(path, target)) {
      stop("Unable to preserve export asset: ", path, call. = FALSE)
    }
    result
  }
  visit = function(value, field = "") {
    if (is.character(value)) {
      for (i in seq_along(value)) {
        path = value[i]
        if (is.na(path) || !nzchar(path)) {
          next
        }
        # Do not mistake a color or label for an unrelated file in the working
        # directory. Named file fields also allow filenames without extensions.
        if (
          !grepl("[/\\\\]", path) &&
            !nzchar(tools::file_ext(path)) &&
            !grepl("file|path|image|texture", field, ignore.case = TRUE)
        ) {
          next
        }
        candidate = path.expand(path)
        if (!grepl("^(/|[A-Za-z]:[/\\\\])", candidate)) {
          candidate = file.path(source_directory, candidate)
        }
        if (!file.exists(candidate) || dir.exists(candidate)) {
          next
        }
        path = normalizePath(candidate, mustWork = TRUE)
        # Permanent model paths retain their directories, including OBJ/MTL
        # dependencies. Temporary generated maps must survive this R session.
        temporary = startsWith(path, normalizePath(tempdir())) ||
          grepl("^/(private/)?(tmp/|var/folders/)", path)
        if (temporary) {
          path = copy_asset(path)
        }
        value[i] = path
      }
    } else if (is.list(value)) {
      for (i in seq_along(value)) {
        child_field = if (is.null(names(value))) field else names(value)[i]
        value[i] = list(visit(value[[i]], child_field))
      }
    }
    # Scene lights and grouping identity live in attributes. RDS preserves the
    # shared group tokens; visit only data attributes, never their environments.
    metadata = attributes(value)
    if (length(metadata)) {
      for (name in setdiff(
        names(metadata),
        c("names", "row.names", "class", "dim", "dimnames")
      )) {
        if (!is.environment(metadata[[name]])) {
          metadata[[name]] = visit(metadata[[name]], name)
        }
      }
      attributes(value) = metadata
    }
    value
  }
  visit(value)
}

#' @keywords internal
write_native_scene_export = function(
  scene,
  args,
  state,
  filename,
  source_directory
) {
  if (length(filename) != 1L || is.na(filename) || !nzchar(trimws(filename))) {
    stop("Enter a filename for the R script.", call. = FALSE)
  }
  filename = path.expand(filename)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", filename)) {
    filename = file.path(source_directory, filename)
  }
  if (tolower(tools::file_ext(filename)) != "r") {
    filename = paste0(filename, ".R")
  }
  parent = normalizePath(dirname(filename), mustWork = TRUE)
  stem = tools::file_path_sans_ext(basename(filename))
  number = 0L
  repeat {
    name = paste0(stem, if (number) paste0("-", number) else "")
    script = file.path(parent, paste0(name, ".R"))
    data = file.path(parent, paste0(name, "-scene.rds"))
    assets = file.path(parent, paste0(name, "-assets"))
    if (!any(file.exists(c(script, data, assets)))) {
      break
    }
    number = number + 1L
  }
  success = FALSE
  on.exit(
    if (!success) unlink(c(script, data, assets), recursive = TRUE),
    add = TRUE
  )
  args = native_export_render_args(args, state)
  # Export the base scene once, and write the cumulative changes separately as
  # ordinary R values. Reopening an exported scene never applies them twice.
  attr(scene, "rayrender_scene_edits") = NULL
  attr(scene, "rayrender_sky_edit") = NULL
  values = native_export_assets(
    list(scene = scene, edits = state$objects, sky = state$sky, args = args),
    assets,
    basename(assets),
    source_directory
  )
  saveRDS(values$scene, data, version = 3)
  literal = function(value) {
    capture.output(dput(
      value,
      control = c("keepNA", "keepInteger", "niceNames", "showAttributes")
    ))
  }
  render = as.call(c(
    list(quote(rayrender::render_scene), quote(scene)),
    values$args
  ))
  lines = c(
    "# Exported from the rayrender native editor.",
    "# Keep this script with its companion -scene.rds and any -assets directory.",
    "# Existing external model/texture files are referenced by absolute path.",
    "# Object edits and sky settings below are cumulative, committed changes.",
    "local({",
    "  files = Filter(Negate(is.null), lapply(sys.frames(), function(frame) frame$ofile))",
    "  command = grep('^--file=', commandArgs(FALSE), value = TRUE)",
    paste0(
      "  script = if (length(files)) utils::tail(files, 1L)[[1L]] else if (length(command)) sub('^--file=', '', command[1L]) else ",
      literal(script)
    ),
    "  directory = dirname(normalizePath(script, mustWork = TRUE))",
    "  previous_directory = setwd(directory)",
    "  on.exit(setwd(previous_directory), add = TRUE)",
    paste0("  scene = readRDS(", literal(basename(data)), ")"),
    "",
    "  # Transforms are 4 x 4 matrices in the enclosing scene's coordinates.",
    paste0("  ", c("scene_edits =", literal(values$edits))),
    paste0("  ", c("sky =", literal(values$sky))),
    "  scene = rayrender::apply_scene_edits(scene, scene_edits, sky)",
    "",
    "  # Current camera and viewport appearance; increase samples as needed.",
    paste0("  ", deparse(render, width.cutoff = 90L)),
    "})"
  )
  writeLines(lines, script, useBytes = TRUE)
  # Parsing verifies the exact file written before reporting a successful export.
  parse(script)
  success = TRUE
  normalizePath(script, mustWork = TRUE)
}

#' @keywords internal
native_scene_exporter = function(scene, args) {
  force(scene)
  force(args)
  source_directory = getwd()
  function(state, filename) {
    write_native_scene_export(scene, args, state, filename, source_directory)
  }
}
