#' Export a Scene to OpenUSD
#'
#' Write a USD stage and collect its referenced files for later rendering.
#' Asset collection happens only during explicit export.
#' @param scene A rayrender scene, optionally with attached cameras and render settings.
#' @param filename Destination `.usd`, `.usdc`, or `.usda` file. USD and USDC use
#'   binary storage; USDA is readable text.
#' @param copy_assets Default `TRUE`. Copy referenced images, models, volume files,
#'   and camera lens files beside the stage. Otherwise retain absolute references.
#' @param overwrite Default `FALSE`. Replace an existing stage.
#' @return Invisibly, the normalized path of the written stage.
#' @export
#' @md
#' @details
#' Uses the registered native interface from `libopenusd` and the optional
#' `digest` package for asset hashes. No Python interpreter is used.
#' Standard USD geometry, cameras, instancers, materials and volume references
#' provide an interchange representation. Typed, versioned rayrender metadata
#' preserves renderer descriptors for [read_usd()]; no RDS payload is stored.
#' Features without a standard USD representation are retained in that metadata
#' and reported as warnings. Other renderers need not reproduce procedural
#' shaders, atmospheric sky, CSG or NanoVDB support.
#' Installed sky datasets remain external dependencies. OBJ files collect their
#' MTL files and texture dependencies. Unsupported R environments, functions and
#' external pointers fail explicitly instead of producing incomplete archives.
#' Export streams numeric data and metadata to a temporary USDA file beside the
#' destination, then uses native OpenUSD to convert binary output. Temporary disk
#' space is required for that text file and the output; they are cleaned up if
#' export fails. Existing destination files are replaced only after validation.
#' @examplesIf requireNamespace("libopenusd", quietly = TRUE) && requireNamespace("digest", quietly = TRUE)
#' # Export explicitly, then restore the scene from a binary USD file.
#' scene = sphere() |> add_camera(camera(lookfrom = c(3, 2, 5)))
#' path = tempfile(fileext = ".usdc")
#' write_usd(scene, path)
#' restored = read_usd(path)
#' unlink(path)
write_usd = function(scene, filename, copy_assets = TRUE, overwrite = FALSE) {
  if (!inherits(scene, "ray_scene")) {
    stop("scene must be a ray_scene.")
  }
  if (
    !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !tolower(tools::file_ext(filename)) %in% c("usd", "usdc", "usda")
  ) {
    stop("filename must end in .usd, .usdc, or .usda.")
  }
  if (
    !is.logical(copy_assets) ||
      length(copy_assets) != 1L ||
      is.na(copy_assets) ||
      !is.logical(overwrite) ||
      length(overwrite) != 1L ||
      is.na(overwrite)
  ) {
    stop("copy_assets and overwrite must be TRUE or FALSE.")
  }
  filename = path.expand(filename)
  if (file.exists(filename) && !overwrite) {
    stop("Destination exists: ", filename)
  }
  if (!dir.exists(dirname(filename))) {
    stop("Destination directory does not exist.")
  }
  filename = file.path(
    normalizePath(dirname(filename), winslash = "/"),
    basename(filename)
  )
  usd_require_native()
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("USD export requires digest.")
  }
  assets = new.env(hash = TRUE, parent = emptyenv())
  usd_collect_assets(scene, assets)
  collected = usd_copy_assets(as.list(assets), filename, copy_assets)
  scene = usd_remap_assets(scene, collected$paths, collected$directories)
  document = usd_stage_document(
    scene,
    collected$dependencies,
    dirname(filename)
  )
  extension = if (tolower(tools::file_ext(filename)) == "usda") {
    ".usda"
  } else {
    ".usdc"
  }
  temporary = tempfile(
    ".rayrender-usd-",
    tmpdir = dirname(filename),
    fileext = extension
  )
  source = if (extension == ".usda") {
    temporary
  } else {
    tempfile(
      ".rayrender-usd-source-",
      tmpdir = dirname(filename),
      fileext = ".usda"
    )
  }
  on.exit(unlink(unique(c(source, temporary))), add = TRUE)
  # Stream directly to disk, including each numeric array and metadata record.
  # OpenUSD reads that file natively; no stage-sized R or C++ string is created.
  usd_stream_native(document, source)
  usd_export_file_native(source, temporary)
  if (source != temporary) {
    unlink(source)
  }
  if (file.exists(filename) && !overwrite) {
    stop("Destination exists: ", filename)
  }
  if (!file.rename(temporary, filename)) {
    stop("Cannot publish USD stage; previous file retained.")
  }
  invisible(filename)
}

#' Import a Rayrender OpenUSD Scene
#'
#' Reconstruct a scene exported by [write_usd()], resolving external assets
#' relative to the stage. The returned scene can be rendered or edited in R.
#' @param filename USD stage exported by `write_usd()`.
#' @return A rayrender scene with cameras, animation and render settings.
#' @export
#' @md
#' @details
#' Uses the registered native `libopenusd` API to read schema version 1 archives.
#' This importer reconstructs rayrender's typed descriptors; it does not infer
#' shaders from arbitrary third-party stages. Edits to standard interchange
#' prims do not replace those descriptors. Unsupported schema versions fail.
#' Geometry and settings use typed USD values, not serialized R objects.
#' @examplesIf interactive() && requireNamespace("libopenusd", quietly = TRUE)
#' # Restore an exported scene after restarting R.
#' # scene = read_usd("sphere.usdc")
#' # render_scene(scene, samples = 256)
read_usd = function(filename) {
  if (
    !is.character(filename) ||
      length(filename) != 1L ||
      is.na(filename) ||
      !file.exists(path.expand(filename))
  ) {
    stop("filename must name an existing USD stage.")
  }
  usd_require_native()
  filename = normalizePath(path.expand(filename), winslash = "/")
  # Select the default root and its direct object records, without traversing
  # large instance prototypes or translating standard geometry back through text.
  root = usd_inspect_native(filename, "", "rayrender:assets")
  if (
    !length(root) || !identical(root[[1]]$custom_data$rayrender$version, 1L)
  ) {
    stop(
      "This importer requires a rayrender USD stage with schema version 1; arbitrary third-party stages are not supported."
    )
  }
  root = root[[1]]
  metadata = root$custom_data$rayrender
  attrs = lapply(metadata$attributes, usd_unpack_value)
  columns = lapply(
    metadata$columns,
    usd_unpack_value,
    restore_attributes = FALSE
  )
  column_attributes = lapply(metadata$columns, function(node) {
    lapply(node$attributes, usd_unpack_value)
  })
  objects = usd_inspect_native(
    filename,
    paste0(root$path, "/Objects/*"),
    character()
  )
  if (length(objects)) {
    objects = objects[order(vapply(objects, `[[`, "", "path"))]
  }
  records = lapply(objects, function(object) {
    record = object$custom_data$rayrender$record
    if (is.null(record)) {
      stop("Missing rayrender record on ", object$path)
    }
    lapply(record, usd_unpack_value)
  })
  for (name in attrs$names) {
    template = columns[[name]]
    values = lapply(records, `[[`, name)
    if (is.list(template)) {
      attributes(values) = column_attributes[[name]]
      columns[name] = list(values)
    } else {
      combined = do.call(c, values)
      if (!length(values)) {
        combined = template
      }
      attributes(combined) = column_attributes[[name]]
      columns[name] = list(combined)
    }
  }
  scene = columns[attrs$names]
  attributes(scene) = attrs
  if (!inherits(scene, "ray_scene")) {
    stop("USD stage does not contain a rayrender scene.")
  }
  paths = directories = character()
  for (asset in root$properties[["rayrender:assets"]]) {
    resolved = asset$resolved_path
    if (!nzchar(resolved)) {
      resolved = if (grepl("^(/|[A-Za-z]:[/\\\\])", asset$path)) {
        asset$path
      } else {
        file.path(dirname(filename), asset$path)
      }
    }
    if (!file.exists(resolved) || dir.exists(resolved)) {
      stop("Missing USD asset: ", resolved)
    }
    resolved = normalizePath(resolved, winslash = "/")
    paths[asset$path] = resolved
    directories[dirname(asset$path)] = dirname(resolved)
  }
  usd_remap_assets(scene, paths, directories)
}

#' @return NULL, invisibly, after loading the native provider.
#' @keywords internal
#' @noRd
usd_require_native = function() {
  if (
    !requireNamespace("libopenusd", quietly = TRUE) ||
      utils::packageVersion("libopenusd") < "26.8.2"
  ) {
    stop("USD import/export requires libopenusd >= 26.8-2.")
  }
  invisible(NULL)
}

#' @param assets Discovered source files and their original aliases.
#' @param filename Destination stage.
#' @param copy_assets Whether to collect files beside the stage.
#' @return Asset aliases, directory aliases, and dependency paths.
#' @keywords internal
#' @noRd
usd_copy_assets = function(assets, filename, copy_assets) {
  directory = paste0(tools::file_path_sans_ext(basename(filename)), "_assets")
  root = dirname(filename)
  paths = new.env(hash = TRUE, parent = emptyenv())
  active = new.env(hash = TRUE, parent = emptyenv())
  dependencies = character()
  # Model hashes include rewritten dependency names, so changed textures cannot
  # mutate the contents referenced by an older archive with the same model text.
  collect = function(source) {
    source = normalizePath(source, winslash = "/", mustWork = TRUE)
    if (exists(source, paths, inherits = FALSE)) {
      return(paths[[source]])
    }
    if (exists(source, active, inherits = FALSE)) {
      stop("Cyclic USD asset dependency: ", source)
    }
    active[[source]] = TRUE
    on.exit(rm(list = source, envir = active), add = TRUE)
    suffix = tolower(tools::file_ext(source))
    rewritten = NULL
    if (suffix %in% c("obj", "mtl", "txt")) {
      lines = readLines(source, warn = FALSE)
      changed = lines
      for (i in seq_along(lines)) {
        words = scan(
          text = lines[i],
          what = "",
          quote = "\"'",
          comment.char = "#",
          quiet = TRUE
        )
        if (!length(words)) {
          next
        }
        if (words[1] == "mtllib") {
          names = words[-1]
          if (
            file.exists(file.path(
              dirname(source),
              paste(names, collapse = " ")
            ))
          ) {
            names = paste(names, collapse = " ")
          }
          targets = vapply(file.path(dirname(source), names), collect, "")
          if (copy_assets) {
            changed[i] = paste(
              "mtllib",
              paste(basename(targets), collapse = " ")
            )
          }
        } else if (
          suffix == "mtl" &&
            (startsWith(tolower(words[1]), "map_") ||
              tolower(words[1]) %in% c("bump", "disp", "decal", "refl"))
        ) {
          if (length(words) < 2L) {
            stop("Missing material texture in ", source)
          }
          candidates = seq.int(2L, length(words))
          candidates = candidates[vapply(
            candidates,
            function(j) {
              file.exists(file.path(
                dirname(source),
                paste(words[j:length(words)], collapse = " ")
              ))
            },
            TRUE
          )]
          if (!length(candidates)) {
            stop("Missing material texture in ", source, ": ", lines[i])
          }
          start = candidates[1]
          target = collect(file.path(
            dirname(source),
            paste(words[start:length(words)], collapse = " ")
          ))
          if (copy_assets) {
            changed[i] = paste(
              c(words[seq_len(start - 1L)], basename(target)),
              collapse = " "
            )
          }
        }
      }
      if (!identical(changed, lines)) {
        rewritten = paste0(paste(changed, collapse = "\n"), "\n")
      }
    }
    target = source
    if (copy_assets) {
      hash = if (is.null(rewritten)) {
        digest::digest(file = source, algo = "sha256", serialize = FALSE)
      } else {
        digest::digest(rewritten, algo = "sha256", serialize = FALSE)
      }
      target = paste0(
        directory,
        "/",
        hash,
        if (nzchar(suffix)) paste0(".", suffix) else ""
      )
      destination = file.path(root, target)
      dir.create(file.path(root, directory), showWarnings = FALSE)
      if (!file.exists(destination)) {
        if (is.null(rewritten)) {
          if (!file.copy(source, destination)) {
            stop("Cannot copy USD asset: ", source)
          }
        } else {
          writeChar(rewritten, destination, eos = NULL, useBytes = TRUE)
        }
      }
    }
    paths[[source]] = target
    dependencies <<- c(dependencies, target)
    target
  }
  aliases = directories = character()
  for (asset in assets) {
    target = collect(asset$source)
    aliases[c(asset$original, asset$source)] = target
    directories[c(dirname(asset$original), dirname(asset$source))] = dirname(
      target
    )
  }
  list(
    paths = aliases,
    directories = directories,
    dependencies = sort(unique(dependencies))
  )
}

#' @param node Typed USD metadata node.
#' @param restore_attributes Default `TRUE`. Restore attributes after decoding.
#' @return Reconstructed R value with attributes and missing values.
#' @keywords internal
#' @noRd
usd_unpack_value = function(node, restore_attributes = TRUE) {
  type = node$type
  result = switch(
    type,
    "NULL" = NULL,
    list = lapply(node$items[order(names(node$items))], usd_unpack_value),
    logical = as.logical(node$data),
    integer = as.integer(node$data),
    double = as.double(node$data),
    character = as.character(node$data),
    raw = as.raw(node$data),
    stop("Unsupported USD R value type: ", type)
  )
  if (type == "list") {
    names(result) = NULL
  }
  if (length(node$missingIndices)) {
    if (
      anyNA(node$missingIndices) ||
        any(node$missingIndices < 1 | node$missingIndices > length(result))
    ) {
      stop("Invalid USD missing-value indices.")
    }
    for (i in seq_along(node$missingIndices)) {
      result[node$missingIndices[i]] = switch(
        node$missingCodes[i],
        "NA" = NA,
        "NaN" = NaN,
        "Inf" = Inf,
        "-Inf" = -Inf,
        stop("Invalid USD missing-value code.")
      )
    }
  }
  if (restore_attributes && !is.null(node$attributes)) {
    attributes(result) = lapply(node$attributes, usd_unpack_value)
  }
  result
}

#' @param x Scene data.
#' @param paths Source-to-destination aliases.
#' @param directories Directory aliases for model basenames.
#' @param field Default `""`. Enclosing field name.
#' @return Scene data with remapped file references.
#' @keywords internal
#' @noRd
usd_remap_assets = function(x, paths, directories, field = "") {
  attrs = attributes(x)
  if (is.character(x)) {
    matches = match(x, names(paths))
    x[!is.na(matches)] = unname(paths[matches[!is.na(matches)]])
    if (field == "basename") {
      matches = match(x, names(directories))
      x[!is.na(matches)] = unname(directories[matches[!is.na(matches)]])
    }
  } else if (is.list(x)) {
    fields = names(x)
    values = unclass(x)
    if (is.null(fields)) {
      fields = rep(field, length(values))
    }
    x = lapply(seq_along(values), function(i) {
      usd_remap_assets(values[[i]], paths, directories, fields[i])
    })
  }
  if (length(attrs)) {
    attributes(x) = lapply(
      attrs,
      usd_remap_assets,
      paths = paths,
      directories = directories
    )
  }
  x
}
#' @param x Scene data.
#' @param assets Environment accumulating referenced files by source path.
#' @param field Default `""`. Enclosing field name.
#' @return NULL, invisibly; referenced files are recorded in assets.
#' @keywords internal
#' @noRd
usd_collect_assets = function(x, assets, field = "") {
  image_field = endsWith(field, "_texture") ||
    endsWith(field, "_texname") ||
    field %in% c("image", "alphaimage", "texture", "environment_light")
  if (is.array(x) && image_field) {
    # Retain the source array in typed USD data for exact rendering semantics.
    # These are not filesystem dependencies and are never copied by discovery.
    return(invisible(NULL))
  }
  file_field = image_field ||
    field %in% c("filename", "fileinfo", "camera_description_file")
  if (is.character(x) && file_field) {
    for (path in x[!is.na(x) & nzchar(x)]) {
      if (exists(path, envir = assets, inherits = FALSE)) {
        next
      }
      if (!file.exists(path.expand(path)) || dir.exists(path.expand(path))) {
        stop("USD asset is missing or is not a file: ", path)
      }
      assets[[path]] = list(
        original = path,
        source = normalizePath(path.expand(path))
      )
    }
    return(invisible(NULL))
  }
  if (!is.list(x)) {
    return(invisible(NULL))
  }
  attrs = attributes(x)
  values = unclass(x)
  fields = names(values)
  if (is.null(fields)) {
    fields = rep(field, length(values))
  }
  # Camera output filenames are not input dependencies.
  if (inherits(x, "ray_camera")) {
    fields[fields == "filename"] = "output_filename"
  }
  for (i in seq_along(values)) {
    usd_collect_assets(values[[i]], assets, fields[i])
  }
  for (key in c(
    "ray_cameras",
    "ray_render_settings",
    "ray_infinite_lights",
    "ray_lights"
  )) {
    if (!is.null(attrs[[key]])) usd_collect_assets(attrs[[key]], assets, key)
  }
  invisible(NULL)
}
