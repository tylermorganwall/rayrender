#' Read a PBRT Scene
#'
#' Convert a PBRT text file to an editable rayrender scene. Parsing and scene
#' conversion use a PEGTL grammar and R scene translation.
#'
#' @md
#' @param filename Path to a PBRT scene, optionally compressed with gzip.
#' @param strict Default `TRUE`. Stop on unsupported scene features. With `FALSE`,
#'   warn and return a diagnostic for each approximation or omitted feature.
#'   Invalid syntax, missing assets, and undefined references always stop.
#' @param asset_dir Default `tempfile("rayrender-pbrt-")`. Directory for generated
#'   texture, environment, and mesh assets. Choose a permanent directory before saving
#'   the returned scene for use in another R session. Existing assets are not
#'   overwritten.
#' @param mesh_storage Default `c("memory", "ply")`. With `"ply"`, store compatible
#'   inline surface meshes as binary PLY assets under `asset_dir`, retaining only
#'   lazy file references in R. Native geometry is still built at render time.
#'   Supplied normals and UVs are preserved; attributes use PLY float32 precision.
#'   Other meshes keep their in-memory representation. Keep the generated assets
#'   for as long as the returned scene is needed.
#' @return A list with `scene`, `render_args`, `diagnostics`, `source_files`, and
#'   `assets`. Render with `do.call(render_scene, c(list(scene = x$scene),
#'   x$render_args))`. Camera and render arguments can be edited before rendering.
#' @details
#' Supports static and animated affine transforms, attribute scopes, named coordinate systems,
#' includes and serial imports (relative to the main file, as in PBRT), object instances, perspective
#' orthographic and realistic cameras, spheres, disks, cylinders, triangle/PLY meshes,
#' bilinear patches (triangulated), cubic Bezier and uniform B-spline curves, diffuse area lights,
#' point/spot lights, uniform and equal-area image infinite lights,
#' RGB diffuse/translucent/conductor/dielectric/hair materials, OpenPBR
#' approximations of coated and Disney materials, constant/image/Ptex/UV checker/direction-mix textures, and
#' homogeneous, uniform-grid, and NanoVDB interior media. Material displacement
#' becomes a bump map; PLY displacement uses rayvertex with simple subdivision.
#' PBRT v4 is the target; common v3 spellings are
#' accepted, including `WorldEnd`, `TransformBegin`, and `matte`/`glass`/`metal`.
#' Environment images always use the v4 equal-area convention; legacy v3
#' latitude-longitude maps must first be converted to equal-area maps.
#' Imported files inherit graphics state without exporting subsequent changes.
#' Object references can resolve after later files are read, including forward and
#' sibling-import references; repeated placements share prepared geometry.
#' Large adjacent compatible diffuse meshes may be batched into multi-shape
#' meshes, and repeated restoring placement-only files use nested instances.
#' Large adjacent sets of compatible nonemissive curves use packed numeric
#' descriptors, retaining each span's control points, widths, normals, parameter
#' interval and subdivision depth. Emissive and medium-bearing curves remain separate.
#' These optimizations preserve geometry but may reduce the number of scene rows.
#' Large static populations are stored as packed affine transforms on instance
#' rows; small populations and animated placements retain individual rows.
#' Parsed-file caching is bounded, but individual files are still parsed in full
#' and the converted scene must fit in memory.
#' Image textures are decoded to linear EXR assets. UV checkerboards are baked
#' at 256 by 256 pixels. Image environments retain the input height and use
#' twice that width; their light transforms are baked into the resampling.
#' EXR environments use their stored RGB primaries, converted to linear sRGB
#' without white adaptation; untagged images default to sRGB, as in PBRT.
#' Ptex textures use the \pkg{ptexr} runtime package. Source face IDs
#' survive batching, instances, and lazy PLY storage. See [texture_ptex()] for
#' encoding, filtering, cache controls, and the required face-local UV layout.
#'
#' This is a scene converter, not a PBRT renderer. RGB values are linear;
#' rayrender's sampler, integrator, reconstruction filter, and tone mapping are
#' used. PBRT output filenames are not passed to `render_scene()`. Unsupported
#' parameters are diagnosed, including spectral data, animated light transforms,
#' measured materials, and procedural media. Approximate conversions (including
#' layered materials, mesh refinement, and finite substitutes for distant lights) require
#' `strict = FALSE`. Missing geometry is never silently substituted.
#' Disney materials map to [openpbr()], including color/bump/roughness maps,
#' metallic, anisotropic specular, clearcoat, sheen, transmission, and subsurface
#' controls. Disney's additive lobes and diffusion profile are approximated by
#' OpenPBR's layered BSDF and random-walk interior. Thin Disney diffuse transmission
#' uses OpenPBR's thin-wall subsurface lobe; fake-subsurface flatness is diagnosed.
#' Opaque solid surfaces with `eta >= 1` use fitted diffuse, partial-metal,
#' specular-tint, sheen, and coat mappings, calibrated against PBRT's directional
#' response. This adjusts the effective specular IOR for metal/coat blending.
#' Thin, transmissive, and subsurface materials retain their requested interface
#' IOR and use the direct lobe approximation instead. The surface fit cannot
#' exactly reproduce Disney retroreflection or its additive GTR1 coat and sheen.
#' Solid Disney subsurface color and distance use the paired Hyperion (2018)
#' conversion to volume albedo and extinction mean free path, then re-encode
#' color for OpenPBR's van de Hulst parameterization. This fit assumes isotropic
#' scattering and an internal-reflection IOR of 1.4; the requested surface IOR
#' is retained. Different IORs, roughnesses, and thin geometry remain approximate.
#' Random walks may need a larger `max_depth` than PBRT's diffusion BSSRDF.
#' Imported render arguments start Russian roulette at bounce 5; override
#' `roulette_active_depth` to change this without reducing `max_depth`.
#' Curves preserve PBRT's `splitdepth` (default 3) as native BVH subdivisions.
#' Shape alpha masks retain their own UV scale/offset and use rayrender's native
#' 8-bit coverage lookup rather than PBRT's mipmapped filtering. Curve alpha
#' masks remain unsupported and are diagnosed.
#' Constant-zero-alpha area emitters use [light()] with `invisible = TRUE`,
#' preserving illumination and camera/specular invisibility. Unlike PBRT's
#' fully invisible emitters, diffuse rays may still intersect these lights;
#' this approximation requires `strict = FALSE`.
#' Direction-mix graphs retain their declaration-time direction transform and
#' geometric-normal weighting on supported color inputs and dielectric/conductor
#' roughness. Roughness conversion follows texture evaluation. Procedural roughness
#' on layered or subsurface materials is still diagnosed as unsupported.
#'
#' PBRT's world coordinates are preserved. Its camera transform maps world to
#' camera and is inverted for rayrender. PBRT's field of view applies to the
#' shorter image dimension; rayrender's applies to image height. This conversion
#' adjusts portrait fields of view and retains the original framing.
#' Horizontal camera handedness, including reflected camera transforms, is
#' retained using `camera_flip_x`. Ordinary rayrender cameras keep their defaults.
#' Film ISO is preserved. Diffuse transmission uses [translucent()] with
#' independent reflection/transmission colors and textures; bump maps on this
#' material are diagnosed and omitted. Combined energy above one is normalized.
#' Animated cameras return a two-pose [camera()] with `mode = "image"` and
#' motion blur enabled. Camera transform times differing from the shutter
#' interval require manual adjustment. Spatial media retain their definition
#' transform at the start of motion and move with an animated enclosing boundary.
#' Light intensity is mapped in rayrender's relative linear RGB units, without
#' PBRT's spectral-to-photometric normalization. NanoVDB files must contain
#' float grids (uncompressed, ZIP or BLOSC); blackbody volume emission uses
#' rayrender's RGB model with PBRT's CIE film normalization accounted for.
#' Constant spectral absorption/scattering tables covering 360--830 nm are
#' imported as neutral RGB coefficients. Other spectral tables remain unsupported.
#'
#' @references [PBRT v4 file format](https://pbrt.org/fileformat-v4)
#' @export
#' @examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # read_pbrt(): import a tiny Newton's cradle, with one pendulum pulled aside.
#' # PBRT object instances share the metal ball and its two suspension wires.
#' # Camera, geometry, and materials come from the file; add rayrender disk lights afterward.
#' file = tempfile(fileext = ".pbrt")
#' lines = c(
#'   'LookAt 4 3 -7  0 1.2 0  0 1 0',
#'   'Camera "perspective" "float fov" 32',
#'   'WorldBegin',
#'   'ObjectBegin "pendulum"',
#'   'Material "conductor" "rgb eta" [.2 .9 1.1] "rgb k" [3.9 2.5 2.1] "float roughness" .01',
#'   'AttributeBegin Translate 0 .76 0 Shape "sphere" "float radius" .24 AttributeEnd',
#'   'Material "diffuse" "rgb reflectance" [.12 .15 .2]'
#' )
#' for (side in c(-1, 1)) {
#'   p = rbind(
#'     c(0, 2.4, side * .4),
#'     c(0, 1.933, side * .267),
#'     c(0, 1.467, side * .133),
#'     c(0, 1, 0)
#'   )
#'   lines = c(
#'     lines,
#'     sprintf(
#'       'Shape "curve" "string type" "cylinder" "point3 P" [%s] "float width" .012',
#'       paste(t(p), collapse = ' ')
#'     )
#'   )
#' }
#' lines = c(lines, 'ObjectEnd')
#' for (x in c(-.72, -.24, .24, .72)) {
#'   lines = c(
#'     lines,
#'     sprintf(
#'       'AttributeBegin Translate %g 0 0 ObjectInstance "pendulum" AttributeEnd',
#'       x
#'     )
#'   )
#' }
#' lines = c(
#'   lines,
#'   'AttributeBegin Translate -1.2 2.4 0 Rotate -28 0 0 1',
#'   'Translate 0 -2.4 0 ObjectInstance "pendulum" AttributeEnd',
#'   'Material "diffuse" "rgb reflectance" [.16 .36 .42]'
#' )
#' for (z in c(-.4, .4)) {
#'   lines = c(
#'     lines,
#'     sprintf(
#'       'AttributeBegin Translate -2 2.4 %g Rotate 90 0 1 0',
#'       z
#'     ),
#'     'Shape "cylinder" "float radius" .055 "float zmax" 3.4 AttributeEnd'
#'   )
#'   for (x in c(-2, 1.4)) {
#'     lines = c(
#'       lines,
#'       sprintf(
#'         'AttributeBegin Translate %g 0 %g Rotate -90 1 0 0',
#'         x,
#'         z
#'       ),
#'       'Shape "cylinder" "float radius" .055 "float zmax" 2.4 AttributeEnd'
#'     )
#'   }
#' }
#' lines = c(
#'   lines,
#'   'Material "diffuse" "rgb reflectance" [.3 .25 .18]',
#'   'Shape "trianglemesh" "point3 P" [-100 0 4 100 0 4 100 100 4 -100 100 4]',
#'   '"integer indices" [0 2 1 0 3 2]',
#'   'Shape "trianglemesh" "point3 P" [-100 0 -100 100 0 -100 100 0 100 -100 0 100]',
#'   '"integer indices" [0 2 1 0 3 2]'
#' )
#' writeLines(lines, file)
#' imported = read_pbrt(file)
#' scene = imported$scene |>
#'   add_infinite_light(disk_light(
#'     direction = c(-1, 1, -1),
#'     angular_diameter = 35,
#'     intensity = 12,
#'     name = "key"
#'   )) |>
#'   add_infinite_light(disk_light(
#'     direction = c(1, .5, 1),
#'     angular_diameter = 30,
#'     intensity = 2,
#'     color = "#c3dfff",
#'     name = "rim"
#'   ))
#' settings = imported$render_args
#' settings$width = 500
#' settings$height = 400
#' settings$samples = 16
#' settings$denoise = TRUE
#' do.call(render_scene, c(list(scene = scene), settings))
#' unlink(file)
read_pbrt = function(
  filename,
  strict = TRUE,
  asset_dir = tempfile("rayrender-pbrt-"),
  mesh_storage = c("memory", "ply")
) {
  if (!is.character(filename) || length(filename) != 1L || is.na(filename)) {
    stop("`filename` must be a single file path.", call. = FALSE)
  }
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  if (
    !is.character(asset_dir) ||
      length(asset_dir) != 1L ||
      is.na(asset_dir) ||
      !nzchar(asset_dir)
  ) {
    stop("`asset_dir` must be a nonempty directory path.", call. = FALSE)
  }
  mesh_storage = match.arg(mesh_storage)
  filename = normalizePath(path.expand(filename), mustWork = TRUE)
  # This environment owns only import-wide resources. Graphics state is a value
  # copied by the attribute stack, so scopes cannot mutate their parents.
  context = new.env(parent = emptyenv())
  context$root = dirname(filename)
  context$strict = strict
  context$asset_dir = path.expand(asset_dir)
  context$assets = character()
  context$mesh_cache = if (mesh_storage == 'ply') {
    cache = new.env(parent = emptyenv())
    cache$directory = tempfile('meshes-', tmpdir = context$asset_dir)
    cache
  } else {
    NULL
  }
  on.exit(
    {
      if (!complete && !is.null(context$mesh_cache)) {
        unlink(context$mesh_cache$directory, recursive = TRUE)
      }
    },
    add = TRUE
  )
  complete = FALSE
  on.exit(
    if (!complete && length(context$assets)) unlink(context$assets),
    add = TRUE
  )
  context$diagnostics = pbrt_row_buffer()
  context$parsed = new.env(parent = emptyenv())
  context$parse_cache = new.env(parent = emptyenv())
  context$parse_cache$bytes = 0
  context$parse_cache$order = character()
  context$source_files = character()
  context$include_stack = character()
  context$objects = new.env(hash = TRUE, parent = emptyenv())
  context$placement_files = new.env(hash = TRUE, parent = emptyenv())
  # Inline meshes use an outer rayrender material. Keep their valid default
  # rayvertex material shared instead of rebuilding it for every small mesh.
  context$mesh_material = rayvertex::material_list()
  context$rows = pbrt_row_buffer(context$mesh_cache)
  context$lights = pbrt_row_buffer()
  context$point_lights = pbrt_row_buffer()
  context$stack = list()
  context$object = NULL
  context$world = FALSE
  context$ended = FALSE
  context$coordinates = list()
  context$materials = list()
  context$textures = list()
  context$alpha_images = new.env(hash = TRUE, parent = emptyenv())
  context$media = list()
  context$settings = list()
  state = list(
    transform = diag(4),
    end_transform = diag(4),
    material = diffuse(color = rep(.5, 3)),
    interface = FALSE,
    area = NULL,
    reverse = FALSE,
    inside = "",
    outside = "",
    defaults = list(),
    active = "All"
  )
  state = pbrt_execute_file(filename, state, context)
  if (!context$world) {
    stop("PBRT scene has no WorldBegin: ", filename, call. = FALSE)
  }
  if (length(context$stack) || !is.null(context$object)) {
    stop("Unclosed PBRT attribute, transform, or object scope.", call. = FALSE)
  }
  rows = pbrt_resolve_instances(pbrt_buffer_rows(context$rows), context)
  scene = if (length(rows)) {
    vctrs::list_unchop(rows, ptype = rows[[1]][FALSE, ])
  } else {
    sphere()[FALSE, ]
  }
  # Constructors already validate these descriptors. Attach each complete list
  # once: add_light() in a loop would revalidate all preceding lights each time.
  lights = pbrt_buffer_rows(context$lights)
  if (length(lights)) {
    names(lights) = paste0("pbrt-light-", seq_along(lights))
    for (i in seq_along(lights)) {
      lights[[i]]$name = names(lights)[i]
    }
    attr(scene, "ray_infinite_lights") = lights
  }
  lights = pbrt_buffer_rows(context$point_lights)
  if (length(lights)) {
    names(lights) = paste0("pbrt-point-", seq_along(lights))
    for (i in seq_along(lights)) {
      lights[[i]]$name = names(lights)[i]
    }
    attr(scene, "ray_lights") = lights
  }
  render_args = pbrt_render_arguments(context)
  notes = pbrt_buffer_rows(context$diagnostics)
  diagnostics = if (length(notes)) {
    unique(vctrs::list_unchop(notes, ptype = notes[[1]][FALSE, ]))
  } else {
    data.frame(
      file = character(),
      line = integer(),
      directive = character(),
      message = character()
    )
  }
  if (nrow(diagnostics)) {
    warning(
      sprintf(
        "PBRT conversion recorded %d diagnostic(s); inspect $diagnostics. First: %s:%d: %s",
        nrow(diagnostics),
        diagnostics$file[1],
        diagnostics$line[1],
        diagnostics$message[1]
      ),
      call. = FALSE
    )
  }
  complete = TRUE
  list(
    scene = scene,
    render_args = render_args,
    diagnostics = diagnostics,
    source_files = context$source_files,
    assets = c(
      context$assets,
      if (!is.null(context$mesh_cache)) {
        list.files(context$mesh_cache$directory, full.names = TRUE)
      }
    )
  )
}

#' @param mesh_cache Default `NULL`. Shared lazy-mesh asset directory.
#' @return An import-local buffer of bounded lists, preserving insertion order.
#' @keywords internal
#' @noRd
pbrt_row_buffer = function(mesh_cache = NULL) {
  buffer = new.env(hash = TRUE, parent = emptyenv())
  buffer$mesh_cache = mesh_cache
  buffer$pending = list()
  buffer$chunks = 0L
  buffer$new_instances = 0L
  buffer
}

#' @param filename Absolute input path.
#' @param state Current graphics state.
#' @param context Import resources and output accumulator.
#' @param continuation Default `list()`. Parameters continuing the included file's final directive.
#' @return Updated graphics state after textual inclusion.
#' @keywords internal
#' @noRd
pbrt_execute_file = function(filename, state, context, continuation = list()) {
  if (filename %in% context$include_stack) {
    stop(
      "Cyclic PBRT Include: ",
      paste(c(context$include_stack, filename), collapse = " -> "),
      call. = FALSE
    )
  }
  context$include_stack = c(context$include_stack, filename)
  on.exit(context$include_stack <- utils::head(context$include_stack, -1L))
  share_placements = context$world &&
    is.null(context$object) &&
    !length(continuation) &&
    identical(state$active, "All") &&
    identical(state$transform, state$end_transform)
  if (
    share_placements &&
      pbrt_share_placement_file(filename, NULL, state, context)
  ) {
    context$source_files = unique(c(context$source_files, filename))
    return(state)
  }
  if (!exists(filename, context$parsed, inherits = FALSE)) {
    commands = pbrt_parse_file(filename)
    context$source_files = unique(c(context$source_files, filename))
    # Repeated small includes benefit from caching, but retaining the parsed
    # commands of production geometry duplicates much of the converted scene.
    # Keep an import-wide 32 MiB LRU cache, shared by nested Import contexts.
    bytes = as.numeric(utils::object.size(commands))
    budget = 32 * 1024^2
    if (bytes <= budget) {
      while (context$parse_cache$bytes + bytes > budget) {
        oldest = context$parse_cache$order[1]
        context$parse_cache$bytes = context$parse_cache$bytes -
          get(oldest, context$parsed, inherits = FALSE)$bytes
        rm(list = oldest, envir = context$parsed)
        context$parse_cache$order = context$parse_cache$order[-1]
      }
      assign(filename, list(commands = commands, bytes = bytes), context$parsed)
      context$parse_cache$bytes = context$parse_cache$bytes + bytes
      context$parse_cache$order = c(context$parse_cache$order, filename)
    }
  } else {
    commands = get(filename, context$parsed, inherits = FALSE)$commands
    context$parse_cache$order = c(
      setdiff(context$parse_cache$order, filename),
      filename
    )
  }
  if (length(continuation)) {
    # Include is textual substitution. A parent may finish the final included
    # directive's parameter list. Modify this invocation, not the cached parse:
    # the same fragment can be included again with different parameter values.
    last = length(commands)
    if (!last || !isTRUE(commands[[last]]$accepts_parameters)) {
      stop(
        "Included file does not end in a parameterized directive: ",
        filename,
        call. = FALSE
      )
    }
    duplicates = intersect(names(commands[[last]]$params), names(continuation))
    for (name in duplicates) {
      if (!identical(commands[[last]]$params[[name]], continuation[[name]])) {
        pbrt_error(commands[[last]], paste("Duplicate parameter:", name))
      }
    }
    commands[[last]]$params = c(
      commands[[last]]$params,
      continuation[setdiff(names(continuation), duplicates)]
    )
  }
  if (
    share_placements &&
      pbrt_share_placement_file(filename, commands, state, context)
  ) {
    return(state)
  }
  if (
    share_placements && pbrt_pack_placement_blocks(commands, state, context)
  ) {
    return(state)
  }
  for (command in commands) {
    state = pbrt_execute(command, state, context)
  }
  if (is.null(context$object) && context$rows$new_instances >= 65536L) {
    # Lower completed populations at file boundaries, rather than retaining
    # tens of millions of reference records until the root file finishes.
    # Unknown forward references stay deferred, with their source locations.
    rows = pbrt_resolve_instances(
      pbrt_buffer_rows(context$rows),
      context,
      partial = TRUE
    )
    rm(list = ls(context$rows, all.names = TRUE), envir = context$rows)
    context$rows$chunks = 0L
    context$rows$pending = list()
    context$rows$new_instances = 0L
    for (row in rows) {
      pbrt_buffer_append(context$rows, row)
    }
    # Retained forward references were just checked. Revisit them after new
    # instances accumulate, not after every unrelated mesh/curve include.
    context$rows$new_instances = 0L
  }
  state
}

#' @param commands Parsed file directives.
#' @param state Static world graphics state, restored after every block.
#' @param context Import resources and output buffer.
#' @return Whether the complete file was handled as independent placements.
#' @keywords internal
#' @noRd
pbrt_pack_placement_blocks = function(commands, state, context) {
  count = length(commands)
  pattern = c(
    'AttributeBegin',
    'ConcatTransform',
    'ObjectInstance',
    'AttributeEnd'
  )
  if (
    count < 256L ||
      count %% 4L ||
      !identical(vapply(commands[1:4], `[[`, character(1), 'name'), pattern)
  ) {
    return(FALSE)
  }
  names = matrix(vapply(commands, `[[`, character(1), 'name'), nrow = 4L)
  if (!all(names == pattern)) {
    return(FALSE)
  }
  instances = seq.int(3L, count, 4L)
  objects = vapply(
    commands[instances],
    function(command) command$args[1],
    character(1)
  )
  if (
    !all(vapply(
      unique(objects),
      pbrt_instance_ready,
      logical(1),
      context = context
    ))
  ) {
    return(FALSE)
  }
  # Validate every local matrix through the ordinary transform implementation.
  # Compute precisely the same product once, instead of maintaining two equal
  # graphics-state matrices and a deferred reference for each static placement.
  transforms = vapply(
    commands[instances - 1L],
    function(command) {
      as.vector(state$transform %*% pbrt_transform(command))
    },
    numeric(16)
  )
  groups = split(seq_along(instances), objects)
  rows = vector('list', length(instances))
  for (indices in groups) {
    prototype = pbrt_instance_prototype(
      commands[[instances[indices[1]]]],
      context
    )
    if (!nrow(prototype)) {
      next
    }
    if (length(indices) >= 64L) {
      prototype$shape_info[[
        1
      ]]$shape_properties$instance_transforms = transforms[,
        indices,
        drop = FALSE
      ]
      rows[[indices[1]]] = prototype
    } else {
      for (i in indices) {
        local = state
        local$transform = local$end_transform = matrix(transforms[, i], 4L, 4L)
        rows[[i]] = pbrt_place(prototype, local, diag(4), context)
      }
    }
  }
  for (row in rows) {
    if (!is.null(row)) {
      pbrt_buffer_append(context$rows, row)
    }
  }
  TRUE
}

#' @param filename Resolved input filename.
#' @param commands Parsed commands, or NULL for a cached-prototype lookup.
#' @param state Current static graphics state.
#' @param context Import resources and output buffer.
#' @return Whether a shared placement prototype was appended.
#' @keywords internal
#' @noRd
pbrt_share_placement_file = function(filename, commands, state, context) {
  entry = context$placement_files[[filename]]
  if (is.null(commands)) {
    if (is.null(entry) || is.null(entry$prototype)) {
      return(FALSE)
    }
    pbrt_buffer_append(
      context$rows,
      pbrt_place(entry$prototype, state, diag(4), context)
    )
    return(TRUE)
  }
  if (is.null(entry)) {
    context$placement_files[[filename]] = list(eligible = NA, prototype = NULL)
    return(FALSE)
  }
  if (identical(entry$eligible, FALSE)) {
    return(FALSE)
  }
  # Only restoring relative-transform scopes can be compiled independently of
  # the caller. Absolute transforms, resource definitions and graphics changes
  # retain normal textual execution. Never cache unresolved forward references.
  scopes = character()
  names = character()
  count = 0L
  eligible = TRUE
  for (command in commands) {
    name = command$name
    if (name %in% c('AttributeBegin', 'TransformBegin')) {
      scopes = c(scopes, name)
    } else if (name %in% c('AttributeEnd', 'TransformEnd')) {
      expected = if (name == 'AttributeEnd') {
        'AttributeBegin'
      } else {
        'TransformBegin'
      }
      if (!length(scopes) || utils::tail(scopes, 1) != expected) {
        eligible = FALSE
        break
      }
      scopes = utils::head(scopes, -1L)
    } else if (name == 'ObjectInstance') {
      count = count + 1L
      names = union(names, command$args[1])
    } else if (
      !(name %in% c('ConcatTransform', 'Translate', 'Scale', 'Rotate')) ||
        !length(scopes)
    ) {
      eligible = FALSE
      break
    }
  }
  if (!eligible || length(scopes) || count < 64L) {
    context$placement_files[[filename]] = list(
      eligible = FALSE,
      prototype = NULL
    )
    return(FALSE)
  }
  if (!all(vapply(names, pbrt_instance_ready, logical(1), context = context))) {
    return(FALSE)
  }
  child = list2env(as.list(context, all.names = TRUE), parent = emptyenv())
  child$rows = pbrt_row_buffer(child$mesh_cache)
  child$stack = list()
  local = state
  local$transform = local$end_transform = diag(4)
  for (command in commands) {
    local = pbrt_execute(command, local, child)
  }
  rows = pbrt_resolve_instances(pbrt_buffer_rows(child$rows), child)
  prototype = if (length(rows)) {
    create_instances(vctrs::list_unchop(rows, ptype = rows[[1]][FALSE, ]))
  } else {
    sphere()[FALSE, ]
  }
  context$placement_files[[filename]] = list(
    eligible = TRUE,
    prototype = prototype
  )
  pbrt_buffer_append(
    context$rows,
    pbrt_place(prototype, state, diag(4), context)
  )
  TRUE
}

#' @param rows Geometry rows and deferred instance references.
#' @param context Import resources with the complete object registry.
#' @param partial Default `FALSE`. Leave unresolved forward references for a later file.
#' @return Geometry rows with shared, prepared instance prototypes.
#' @keywords internal
#' @noRd
pbrt_resolve_instances = function(rows, context, partial = FALSE) {
  ready = new.env(parent = emptyenv())
  references = which(vapply(
    rows,
    inherits,
    logical(1),
    "pbrt_instance_reference"
  ))
  for (i in references) {
    name = rows[[i]]$command$args[1]
    if (is.null(ready[[name]])) {
      ready[[name]] = !partial || pbrt_instance_ready(name, context)
    }
  }
  static = which(vapply(
    rows,
    function(row) {
      inherits(row, "pbrt_instance_reference") &&
        ready[[row$command$args[1]]] &&
        identical(row$transform, row$end_transform)
    },
    logical(1)
  ))
  groups = split(
    static,
    vapply(rows[static], function(row) row$command$args[1], character(1))
  )
  for (indices in groups) {
    # Keep small scenes individually editable. Large repeated populations need
    # one matrix per placement, not a full set of scene columns per placement.
    if (length(indices) < 64L) {
      next
    }
    prototype = pbrt_instance_prototype(rows[[indices[1]]]$command, context)
    if (nrow(prototype)) {
      transforms = vapply(
        rows[indices],
        function(row) as.vector(row$transform),
        numeric(16)
      )
      prototype$shape_info[[
        1
      ]]$shape_properties$instance_transforms = transforms
    }
    rows[indices] = list(NULL)
    rows[[indices[1]]] = prototype
  }
  for (i in seq_along(rows)) {
    reference = rows[[i]]
    if (!inherits(reference, "pbrt_instance_reference")) {
      next
    }
    if (!ready[[reference$command$args[1]]]) {
      next
    }
    prototype = pbrt_instance_prototype(reference$command, context)
    placement_context = list(settings = list(TransformTimes = reference$times))
    rows[[i]] = pbrt_place(
      prototype,
      reference,
      diag(4),
      placement_context
    )
  }
  Filter(
    function(row) {
      inherits(row, "pbrt_instance_reference") ||
        (!is.null(row) && nrow(row) > 0L)
    },
    rows
  )
}

#' @param name Object name to check before incremental instance lowering.
#' @param context Import resources and object definitions.
#' @param ancestors Default `character()`. Definitions currently being checked.
#' @return Whether this definition and its nested instances can be resolved now.
#' @keywords internal
#' @noRd
pbrt_instance_ready = function(name, context, ancestors = character()) {
  if (!exists(name, context$objects, inherits = FALSE)) {
    return(FALSE)
  }
  # Let the strict resolver report cycles with the original directive location.
  if (name %in% ancestors) {
    return(TRUE)
  }
  definition = get(name, context$objects, inherits = FALSE)
  if (!is.null(definition$prototype)) {
    return(TRUE)
  }
  references = Filter(
    function(row) inherits(row, "pbrt_instance_reference"),
    definition$rows
  )
  all(vapply(
    references,
    function(row) {
      pbrt_instance_ready(row$command$args[1], context, c(ancestors, name))
    },
    logical(1)
  ))
}

#' @param command Deferred ObjectInstance directive, for name and error location.
#' @param context Import resources with the complete object registry.
#' @return One prepared instance row, or an empty scene for an empty definition.
#' @keywords internal
#' @noRd
pbrt_instance_prototype = function(command, context) {
  name = command$args[1]
  if (!exists(name, context$objects, inherits = FALSE)) {
    pbrt_error(command, paste("Undefined object:", name))
  }
  definition = get(name, context$objects, inherits = FALSE)
  if (definition$resolving) {
    pbrt_error(command, paste("Cyclic object instance:", name))
  }
  if (is.null(definition$prototype)) {
    definition$resolving = TRUE
    geometry = pbrt_resolve_instances(definition$rows, context)
    definition$prototype = if (length(geometry)) {
      create_instances(
        vctrs::list_unchop(geometry, ptype = geometry[[1]][FALSE, ]),
        x = 0
      )
    } else {
      sphere()[FALSE, ]
    }
    definition$resolving = FALSE
    definition$rows = NULL
  }
  definition$prototype
}

#' @param buffer An import-local row buffer.
#' @return The accumulated rows in insertion order, as a list.
#' @keywords internal
#' @noRd
pbrt_buffer_rows = function(buffer) {
  chunks = mget(
    as.character(seq_len(buffer$chunks)),
    envir = buffer,
    inherits = FALSE
  )
  chunks = c(
    chunks,
    list(pbrt_store_mesh_rows(
      pbrt_batch_curve_rows(buffer$pending),
      buffer$mesh_cache
    ))
  )
  unlist(chunks, recursive = FALSE, use.names = FALSE)
}

#' @param context Import resources and camera/settings descriptions.
#' @return Named arguments for render_scene().
#' @keywords internal
#' @noRd
pbrt_render_arguments = function(context) {
  settings = context$settings
  result = list(
    width = 1280L,
    height = 720L,
    samples = 16L,
    shutter_speed = 1,
    aperture = 0,
    ambient_light = FALSE,
    backgroundhigh = "black",
    backgroundlow = "black",
    integrator_type = "nee",
    # PBRT path integrators use roulette. Avoid carrying negligible paths all
    # the way to an imported high maxdepth; callers can override this argument.
    roulette_active_depth = 5L
  )
  if (!is.null(settings$Film)) {
    p = pbrt_parameters(settings$Film, context)
    result$width = pbrt_get(p, "xresolution", 1280L, 1L)
    result$height = pbrt_get(p, "yresolution", 720L, 1L)
    result$film_size = pbrt_get(p, "diagonal", 35, 1L)
    result$iso = pbrt_get(p, "iso", 100, 1L)
    if (!is.finite(result$iso) || result$iso < 0) {
      pbrt_error(settings$Film, "Film iso must be finite and nonnegative.")
    }
    pbrt_get(p, "filename") # Never adopt a source file's output destination.
    if (!settings$Film$args[1] %in% c("rgb", "image")) {
      pbrt_note(
        context,
        settings$Film,
        "Film is rendered with rayrender's RGB film."
      )
    }
    pbrt_unused(p)
  }
  if (
    any(
      c(result$width, result$height) < 1 |
        c(result$width, result$height) %% 1 != 0
    )
  ) {
    stop("PBRT film resolution must be positive integers.", call. = FALSE)
  }
  if (!is.null(settings$Sampler)) {
    p = pbrt_parameters(settings$Sampler, context)
    result$samples = pbrt_get(p, "pixelsamples", 16L, 1L)
    pbrt_unused(p)
  }
  if (result$samples < 1 || result$samples %% 1 != 0) {
    stop("PBRT pixelsamples must be a positive integer.", call. = FALSE)
  }
  if (!is.null(settings$Integrator)) {
    p = pbrt_parameters(settings$Integrator, context)
    result$max_depth = pbrt_get(p, "maxdepth", 5L, 1L)
    pbrt_get(p, "lightsampler")
    pbrt_unused(p)
  }
  camera = settings$Camera
  if (is.null(camera)) {
    camera = list(
      name = "Camera",
      args = "perspective",
      params = list(),
      file = context$source_files[1],
      line = 1L,
      transform = diag(4),
      end_transform = diag(4)
    )
  }
  p = pbrt_parameters(camera, context)
  camera_to_world = solve(camera$transform)
  basis = camera_to_world[1:3, 1:3, drop = FALSE]
  if (max(abs(crossprod(basis) - diag(3))) > 1e-6) {
    pbrt_note(
      context,
      camera,
      "Scaled or sheared cameras are reduced to a LookAt frame."
    )
  }
  result$lookfrom = camera_to_world[1:3, 4]
  result$lookat = result$lookfrom + basis[, 3]
  result$camera_up = basis[, 2]
  # Our conventional cameras and PBRT use the same right/up/forward frame,
  # but rayrender reverses x when displaying its native raster. Its realistic
  # camera already reverses film x in the integrator. Retain reflected PBRT
  # frames (e.g. Scale -1 1 1) by reversing the film rather than the world.
  result$camera_flip_x = if (camera$args[1] == "realistic") {
    det(basis) < 0
  } else {
    det(basis) > 0
  }
  aspect = result$width / result$height
  frame_aspect = pbrt_get(p, "frameaspectratio", aspect, 1L)
  if (frame_aspect <= 0) {
    pbrt_error(camera, "frameaspectratio must be positive.")
  }
  screen = if (frame_aspect > 1) {
    c(-frame_aspect, frame_aspect, -1, 1)
  } else {
    c(-1, 1, -1 / frame_aspect, 1 / frame_aspect)
  }
  screen = pbrt_get(p, "screenwindow", screen, 4L)
  if (screen[2] <= screen[1] || screen[4] <= screen[3]) {
    pbrt_error(camera, "Invalid camera screenwindow.")
  }
  if (
    abs(sum(screen[1:2])) > 1e-8 ||
      abs(sum(screen[3:4])) > 1e-8 ||
      abs(diff(screen[1:2]) / diff(screen[3:4]) - aspect) > 1e-6
  ) {
    pbrt_note(
      context,
      camera,
      "Off-center or non-square-pixel screenwindow is approximated by centered framing."
    )
  }
  if (camera$args[1] == "realistic") {
    result$camera_description_file = pbrt_asset(
      pbrt_get(p, "lensfile", "", 1L),
      context,
      camera
    )
    result$fov = -1
  } else if (camera$args[1] == "orthographic") {
    result$fov = 0
    result$ortho_dimensions = c(diff(screen[1:2]), diff(screen[3:4]))
  } else {
    if (camera$args[1] != "perspective") {
      pbrt_note(
        context,
        camera,
        "Unsupported camera replaced by a perspective camera."
      )
    }
    fov = pbrt_get(p, "fov", 90, 1L)
    if (fov <= 0 || fov >= 180) {
      pbrt_error(camera, "Perspective fov must lie between 0 and 180 degrees.")
    }
    result$fov = 360 / pi * atan(tan(fov * pi / 360) * diff(screen[3:4]) / 2)
  }
  result$aperture = if (camera$args[1] == "realistic") {
    pbrt_get(p, "aperturediameter", 1, 1L)
  } else {
    2 * pbrt_get(p, "lensradius", 0, 1L)
  }
  result$focal_distance = pbrt_get(
    p,
    if (camera$args[1] == "realistic") "focusdistance" else "focaldistance",
    if (result$aperture > 0) 10 else 1,
    1L
  )
  if (result$aperture < 0 || result$focal_distance <= 0) {
    pbrt_error(
      camera,
      "Lens radius must be nonnegative and focal distance positive."
    )
  }
  result$shutteropen = pbrt_get(p, "shutteropen", 0, 1L)
  result$shutterclose = pbrt_get(p, "shutterclose", 1, 1L)
  if (result$shutterclose < result$shutteropen) {
    pbrt_error(camera, "shutterclose must not precede shutteropen.")
  }
  pbrt_unused(p)
  if (!isTRUE(all.equal(camera$transform, camera$end_transform))) {
    times = if (is.null(settings$TransformTimes)) {
      c(0, 1)
    } else {
      settings$TransformTimes$args
    }
    if (
      !isTRUE(all.equal(
        as.numeric(times),
        c(result$shutteropen, result$shutterclose)
      ))
    ) {
      pbrt_note(
        context,
        camera,
        "Camera transform endpoints are mapped to the shutter interval; differing TransformTimes require manual adjustment."
      )
    }
    constructor = get("camera", mode = "function")
    fields = intersect(names(result), names(formals(constructor)))
    moving_camera = do.call(constructor, result[fields])
    end = solve(camera$end_transform)
    endpoint = moving_camera$motion[1, ]
    endpoint[1, c("x", "y", "z")] = end[1:3, 4]
    endpoint[1, c("dx", "dy", "dz")] = end[1:3, 4] + end[1:3, 3]
    endpoint[1, c("upx", "upy", "upz")] = end[1:3, 2]
    moving_camera$motion = rbind(moving_camera$motion, endpoint)
    moving_camera$camera_motion_blur = TRUE
    moving_camera$shutter_speed = 1
    result[fields] = NULL
    result$camera = moving_camera
    result$mode = "image"
  }
  result
}

#' @param filename Absolute input path.
#' @return Parsed directives with typed parameters and source locations.
#' @keywords internal
#' @noRd
pbrt_parse_file = function(filename) {
  lexed = pbrt_file_tokens(filename)
  tokens = lexed$tokens
  line = lexed$line
  arrays = lexed$arrays
  quoted = startsWith(tokens, "\"")
  starts = which(
    !quoted &
      grepl("^[A-Z][A-Za-z]+$", tokens) &
      !tokens %in% c("All", "StartTime", "EndTime", "Inf", "NaN")
  )
  if (!length(tokens)) {
    return(list())
  }
  if (!length(starts) || starts[1] != 1L) {
    stop(filename, ":", line[1], ": expected a PBRT directive.", call. = FALSE)
  }
  ends = c(starts[-1] - 1L, length(tokens))
  arity = c(
    Include = 1,
    Import = 1,
    Camera = 1,
    Film = 1,
    Sampler = 1,
    Integrator = 1,
    PixelFilter = 1,
    Accelerator = 1,
    Material = 1,
    MakeNamedMaterial = 1,
    NamedMaterial = 1,
    Shape = 1,
    LightSource = 1,
    AreaLightSource = 1,
    MakeNamedMedium = 1,
    MediumInterface = 2,
    Texture = 3,
    Attribute = 1,
    ObjectBegin = 1,
    ObjectInstance = 1,
    CoordinateSystem = 1,
    CoordSysTransform = 1,
    ColorSpace = 1,
    ActiveTransform = 1,
    Translate = 3,
    Scale = 3,
    Rotate = 4,
    LookAt = 9,
    TransformTimes = 2,
    Transform = 16,
    ConcatTransform = 16,
    Identity = 0,
    WorldBegin = 0,
    WorldEnd = 0,
    AttributeBegin = 0,
    AttributeEnd = 0,
    TransformBegin = 0,
    TransformEnd = 0,
    ObjectEnd = 0,
    ReverseOrientation = 0,
    Option = 2
  )
  commands = vector("list", length(starts))
  # Parameter declarations repeat across millions of primitives. Decode each
  # spelling once; values and directive-local validation still run separately.
  declarations = new.env(hash = TRUE, parent = emptyenv())
  parameter_directives = c(
    "Camera",
    "Film",
    "Sampler",
    "Integrator",
    "PixelFilter",
    "Accelerator",
    "Material",
    "MakeNamedMaterial",
    "Shape",
    "LightSource",
    "AreaLightSource",
    "MakeNamedMedium",
    "Texture",
    "Attribute",
    "Include"
  )
  for (i in seq_along(starts)) {
    from = starts[i]
    command = list(
      name = tokens[from],
      args = character(),
      params = list(),
      file = filename,
      line = line[from]
    )
    command$accepts_parameters = command$name %in% parameter_directives
    t = if (ends[i] > from) {
      tokens[seq.int(from + 1L, ends[i])]
    } else {
      character()
    }
    if (
      command$name %in%
        c("Transform", "ConcatTransform") &&
        length(t) >= 3L &&
        t[1] == "[" &&
        startsWith(t[2], "@array:")
    ) {
      values = arrays[[as.integer(substring(t[2], 8L))]]
      if (is.numeric(values)) {
        # Placement-heavy files contain millions of these matrices. Keep the
        # lexer's doubles rather than allocating sixteen strings per matrix.
        # pbrt_transform() performs the same finite/affine checks at execution.
        if (length(values) != 16L) {
          pbrt_error(command, "Transform requires 16 bracketed numbers.")
        }
        if (length(t) != 3L) {
          pbrt_error(command, "Unexpected extra operands or parameters.")
        }
        command$args = values
        commands[[i]] = command
        next
      }
      t = c(
        "[",
        as.character(values),
        "]",
        utils::tail(t, -3L)
      )
    }
    n = unname(arity[command$name])
    if (is.na(n)) {
      packed = which(startsWith(t, "@array:"))
      if (length(packed)) {
        pieces = as.list(t)
        for (j in packed) {
          pieces[[j]] = as.character(arrays[[as.integer(substring(t[j], 8L))]])
        }
        t = unlist(pieces, use.names = FALSE)
      }
      command$args = t
      commands[[i]] = command
      next
    }
    if (
      command$name %in%
        c("Transform", "ConcatTransform") &&
        length(t) &&
        t[1] == "["
    ) {
      if (length(t) < 18L || t[18L] != "]") {
        pbrt_error(command, "Transform requires 16 bracketed numbers.")
      }
      t = c(t[2:17], t[-seq_len(18)])
    }
    if (length(t) < n) {
      pbrt_error(command, paste("Expected", n, "operands."))
    }
    if (n) {
      command$args = pbrt_unquote(t[seq_len(n)], command)
      t = t[-seq_len(n)]
    }
    if (length(t) && !command$accepts_parameters) {
      pbrt_error(command, "Unexpected extra operands or parameters.")
    }
    at = 1L
    while (at <= length(t)) {
      if (!startsWith(t[at], "\"")) {
        pbrt_error(command, "Expected a quoted parameter declaration.")
      }
      declaration = declarations[[t[at]]]
      if (is.null(declaration)) {
        declaration = strsplit(
          trimws(pbrt_unquote(t[at], command)),
          "[[:space:]]+"
        )[[1]]
        if (length(declaration) != 2L) {
          pbrt_error(command, "Parameters must be declared as \"type name\".")
        }
        declarations[[t[at]]] = declaration
      }
      type = declaration[1]
      name = declaration[2]
      at = at + 1L
      if (at > length(t)) {
        pbrt_error(command, paste("Missing value for", name))
      }
      if (t[at] == "[") {
        stop_at = which(t == "]" & seq_along(t) > at)[1]
        if (is.na(stop_at)) {
          pbrt_error(command, paste("Unclosed array for", name))
        }
        value = if (stop_at == at + 1L) {
          character()
        } else {
          t[seq.int(at + 1L, stop_at - 1L)]
        }
        at = stop_at + 1L
      } else {
        value = t[at]
        at = at + 1L
      }
      if (length(value) == 1L && startsWith(value, "@array:")) {
        value = arrays[[as.integer(substring(value, 8L))]]
      }
      string_value = if (is.numeric(value)) {
        FALSE
      } else {
        startsWith(value, "\"")
      }
      if (!is.numeric(value)) {
        value = pbrt_unquote(value, command)
      }
      if (
        type %in%
          c("string", "texture") ||
          (type == "spectrum" && all(string_value))
      ) {
        if (!all(string_value)) {
          pbrt_error(command, paste(name, "requires quoted string values."))
        }
      } else if (type == "bool") {
        if (!all(value %in% c("true", "false"))) {
          pbrt_error(command, paste(name, "requires true or false."))
        }
        value = value == "true"
      } else {
        if (
          !type %in%
            c(
              "integer",
              "float",
              "point2",
              "vector2",
              "point3",
              "vector3",
              "normal3",
              "point",
              "vector",
              "normal",
              "rgb",
              "color",
              "spectrum",
              "blackbody"
            )
        ) {
          pbrt_error(command, paste("Unknown parameter type:", type))
        }
        value = suppressWarnings(as.numeric(value))
        if (any(!is.finite(value))) {
          pbrt_error(command, paste("Non-numeric or nonfinite value for", name))
        }
        if (type == "integer" && any(value != floor(value))) {
          pbrt_error(command, paste(name, "requires integer values."))
        }
        tuple = if (type %in% c("point2", "vector2")) {
          2L
        } else if (
          type %in%
            c(
              "point3",
              "vector3",
              "normal3",
              "point",
              "vector",
              "normal",
              "rgb",
              "color"
            )
        ) {
          3L
        } else {
          1L
        }
        if (length(value) %% tuple != 0) {
          pbrt_error(command, paste("Wrong tuple length for", name))
        }
      }
      parameter = list(type = type, value = value)
      if (
        name %in%
          names(command$params) &&
          !identical(command$params[[name]], parameter)
      ) {
        pbrt_error(command, paste("Duplicate parameter:", name))
      }
      command$params[[name]] = parameter
    }
    commands[[i]] = command
  }
  commands
}

#' @param filename PBRT input, optionally gzip-compressed.
#' @return Compact grammar tokens and typed arrays.
#' @keywords internal
#' @noRd
pbrt_file_tokens = function(filename) {
  if (!grepl("[.]gz$", filename, ignore.case = TRUE)) {
    return(pbrt_lex_cpp(filename))
  }
  temporary = tempfile(fileext = ".pbrt")
  on.exit(unlink(temporary), add = TRUE)
  input = gzfile(filename, "rb")
  on.exit(close(input), add = TRUE)
  output = file(temporary, "wb")
  tryCatch(
    {
      repeat {
        bytes = readBin(input, what = "raw", n = 16L * 1024L * 1024L)
        if (!length(bytes)) {
          break
        }
        writeBin(bytes, output)
      }
    },
    finally = close(output)
  )
  pbrt_lex_cpp(temporary, filename)
}

#' @param command Parsed directive.
#' @param state Current graphics state.
#' @param context Import resources and output accumulator.
#' @return Updated graphics state.
#' @keywords internal
#' @noRd
pbrt_execute = function(command, state, context) {
  name = command$name
  args = command$args
  if (context$ended) {
    pbrt_error(command, "Directive after WorldEnd.")
  }
  if (
    !context$world &&
      name %in%
        c(
          "Shape",
          "Material",
          "MakeNamedMaterial",
          "NamedMaterial",
          "Texture",
          "LightSource",
          "AreaLightSource",
          "ObjectBegin",
          "ObjectInstance"
        )
  ) {
    pbrt_error(command, paste(name, "must follow WorldBegin."))
  }
  if (name == "Include") {
    return(pbrt_execute_file(
      pbrt_asset(args[1], context, command),
      state,
      context,
      continuation = command$params
    ))
  }
  if (name == "Import") {
    if (!context$world) {
      pbrt_error(command, "Import must follow WorldBegin.")
    }
    pbrt_import_file(pbrt_asset(args[1], context, command), state, context)
    return(state)
  }
  if (
    name %in%
      c(
        "Identity",
        "Translate",
        "Scale",
        "Rotate",
        "LookAt",
        "Transform",
        "ConcatTransform"
      )
  ) {
    transform = pbrt_transform(command)
    fields = switch(
      state$active,
      All = c("transform", "end_transform"),
      StartTime = "transform",
      EndTime = "end_transform"
    )
    for (field in fields) {
      state[[field]] = if (name %in% c("Identity", "Transform")) {
        transform
      } else {
        state[[field]] %*% transform
      }
    }
    return(state)
  }
  if (name == "ActiveTransform") {
    if (!args[1] %in% c("All", "StartTime", "EndTime")) {
      pbrt_error(command, "Invalid ActiveTransform operand.")
    }
    state$active = args[1]
    return(state)
  }
  if (name == "CoordinateSystem") {
    context$coordinates[[args[1]]] = list(
      start = state$transform,
      end = state$end_transform
    )
  } else if (name == "CoordSysTransform") {
    if (is.null(context$coordinates[[args[1]]])) {
      pbrt_error(command, paste("Undefined coordinate system:", args[1]))
    }
    state$transform = context$coordinates[[args[1]]]$start
    state$end_transform = context$coordinates[[args[1]]]$end
  } else if (name %in% c("AttributeBegin", "TransformBegin", "ObjectBegin")) {
    if (name == "ObjectBegin") {
      if (!is.null(context$object)) {
        pbrt_error(command, "Nested ObjectBegin is not supported.")
      }
      if (exists(args[1], context$objects, inherits = FALSE)) {
        pbrt_error(command, paste("Duplicate object:", args[1]))
      }
      context$object = args[1]
      context$outer_rows = context$rows
      context$rows = pbrt_row_buffer(context$mesh_cache)
    }
    context$stack[[length(context$stack) + 1L]] = list(
      kind = name,
      state = state
    )
  } else if (name %in% c("AttributeEnd", "TransformEnd", "ObjectEnd")) {
    expected = switch(
      name,
      AttributeEnd = "AttributeBegin",
      TransformEnd = "TransformBegin",
      ObjectEnd = "ObjectBegin"
    )
    if (
      !length(context$stack) ||
        utils::tail(context$stack, 1)[[1]]$kind != expected
    ) {
      pbrt_error(command, paste("Unmatched", name))
    }
    saved = utils::tail(context$stack, 1)[[1]]$state
    context$stack = utils::head(context$stack, -1L)
    if (name == "TransformEnd") {
      state$transform = saved$transform
      state$end_transform = saved$end_transform
      state$active = saved$active
    } else {
      state = saved
    }
    if (name == "ObjectEnd") {
      definition = new.env(parent = emptyenv())
      definition$rows = pbrt_buffer_rows(context$rows)
      definition$resolving = FALSE
      definition$prototype = NULL
      assign(context$object, definition, context$objects)
      context$rows = context$outer_rows
      context$object = NULL
    }
  } else if (name == "WorldBegin") {
    if (context$world || length(context$stack)) {
      pbrt_error(
        command,
        "WorldBegin must occur once, outside attribute scopes."
      )
    }
    context$world = TRUE
    state$transform = diag(4)
    state$end_transform = diag(4)
    state$active = "All"
    context$coordinates$world = list(start = diag(4), end = diag(4))
  } else if (name == "WorldEnd") {
    if (!context$world || length(context$stack)) {
      pbrt_error(command, "Unmatched WorldEnd.")
    }
    context$ended = TRUE
  } else if (
    name %in%
      c(
        "Camera",
        "Film",
        "Sampler",
        "Integrator",
        "PixelFilter",
        "Accelerator",
        "TransformTimes"
      )
  ) {
    if (context$world) {
      pbrt_error(command, paste(name, "must precede WorldBegin."))
    }
    if (name == "TransformTimes") {
      command$args = suppressWarnings(as.numeric(args))
      if (any(!is.finite(command$args)) || command$args[2] <= command$args[1]) {
        pbrt_error(command, "TransformTimes requires finite start < end.")
      }
    }
    command$transform = state$transform
    command$end_transform = state$end_transform
    context$settings[[name]] = command
    if (name == "Camera") {
      context$coordinates$camera = list(
        start = solve(state$transform),
        end = solve(state$end_transform)
      )
    }
  } else if (name == "ColorSpace") {
    if (args[1] != "srgb") {
      pbrt_note(
        context,
        command,
        "Only linear sRGB values are supported; other ColorSpace values are not converted."
      )
    }
  } else if (name == "ReverseOrientation") {
    state$reverse = !state$reverse
  } else if (name == "Attribute") {
    if (!args[1] %in% c("shape", "light", "material", "medium", "texture")) {
      pbrt_error(command, "Unknown Attribute target.")
    }
    previous = state$defaults[[args[1]]]
    state$defaults[[args[1]]] = utils::modifyList(
      if (is.null(previous)) list() else previous,
      command$params
    )
  } else if (name %in% c("Material", "MakeNamedMaterial", "NamedMaterial")) {
    if (name == "NamedMaterial") {
      material = context$materials[[args[1]]]
      if (is.null(material)) {
        pbrt_note(
          context,
          command,
          paste(
            "Undefined material:",
            args[1],
            "(possibly a forward reference); replaced by neutral diffuse."
          )
        )
        material = list(
          material = diffuse(color = rep(.5, 3)),
          interface = FALSE
        )
      }
    } else {
      material = pbrt_material(command, state, context)
      if (name == "MakeNamedMaterial") {
        if (args[1] %in% names(context$materials)) {
          pbrt_error(command, paste("Duplicate named material:", args[1]))
        }
        context$materials[[args[1]]] = material
      }
    }
    if (name != "MakeNamedMaterial") {
      state$material = material$material
      state$interface = material$interface
    }
  } else if (name == "Texture") {
    if (args[1] %in% names(context$textures)) {
      pbrt_note(
        context,
        command,
        paste(
          "Duplicate named texture:",
          args[1],
          "is redefined for subsequent references; lexical texture scope is not preserved."
        )
      )
    }
    context$textures[[args[1]]] = pbrt_texture(command, state, context)
  } else if (name == "MakeNamedMedium") {
    if (args[1] %in% names(context$media)) {
      pbrt_error(command, paste("Duplicate named medium:", args[1]))
    }
    context$media[[args[1]]] = pbrt_medium(command, state, context)
  } else if (name == "MediumInterface") {
    for (medium in args[nzchar(args)]) {
      if (!medium %in% names(context$media)) {
        pbrt_error(command, paste("Undefined medium:", medium))
      }
    }
    if (nzchar(args[2])) {
      pbrt_note(
        context,
        command,
        "Exterior media are unsupported; only bounded interior media are converted."
      )
    }
    state$inside = args[1]
    state$outside = args[2]
  } else if (name == "AreaLightSource") {
    defaults = state$defaults$light
    command$params = utils::modifyList(
      if (is.null(defaults)) list() else defaults,
      command$params
    )
    state$area = if (args[1] %in% c("", "none")) NULL else command
  } else if (name == "LightSource") {
    if (!is.null(context$object)) {
      pbrt_error(
        command,
        "Non-area lights inside ObjectBegin cannot be instanced."
      )
    }
    pbrt_light(command, state, context)
  } else if (name %in% c("Shape", "ObjectInstance")) {
    if (!context$world) {
      pbrt_error(command, paste(name, "must follow WorldBegin."))
    }
    if (name == "ObjectInstance") {
      # PBRT resolves names after all imports finish, including references to
      # objects in sibling imports. Keep only placement state until then.
      reference = structure(
        list(
          command = command,
          transform = state$transform,
          end_transform = state$end_transform,
          times = context$settings$TransformTimes
        ),
        class = "pbrt_instance_reference"
      )
      pbrt_buffer_append(context$rows, reference)
      return(state)
    } else {
      rows = pbrt_shape(command, state, context)
    }
    if (!is.null(rows) && nrow(rows)) {
      pbrt_buffer_append(context$rows, rows)
    }
  } else {
    pbrt_note(context, command, paste("Unsupported directive omitted:", name))
  }
  state
}

#' @param filename Resolved imported PBRT file.
#' @param state Graphics state inherited at the Import directive.
#' @param context Parent import resources and output buffers.
#' @return Nothing; appends geometry and lights without changing parent state.
#' @keywords internal
#' @noRd
pbrt_import_file = function(filename, state, context) {
  # Import inherits graphics state without exporting subsequent changes.
  # Object definitions belong to the shared scene and resolve after parsing;
  # other translator symbol tables retain their existing lexical behavior.
  child = list2env(as.list(context, all.names = TRUE), parent = emptyenv())
  child$stack = list()
  on.exit(
    {
      # Retain provenance and ownership even on failure, so read_pbrt's cleanup
      # removes assets generated in nested imports as well as the parent file.
      context$source_files = unique(c(context$source_files, child$source_files))
      context$assets = unique(c(context$assets, child$assets))
    },
    add = TRUE
  )
  pbrt_execute_file(filename, state, child)
  if (length(child$stack) || !identical(child$object, context$object)) {
    stop("Unclosed PBRT scope in imported file: ", filename, call. = FALSE)
  }
  invisible(NULL)
}
#' @param command Parsed directive.
#' @param context Import resources.
#' @param defaults Default `list()`. Inherited typed parameter declarations.
#' @return Parameter reader recording consumed fields for unsupported diagnostics.
#' @keywords internal
#' @noRd
pbrt_parameters = function(command, context, defaults = list()) {
  p = new.env(parent = emptyenv())
  p$command = command
  p$context = context
  p$values = utils::modifyList(
    if (is.null(defaults)) list() else defaults,
    command$params
  )
  p$used = character()
  p
}

#' @param p Parameter reader.
#' @param name Parameter name.
#' @param default Default `NULL`. Value when absent.
#' @param size Default `NULL`. Required scalar/vector length, when applicable.
#' @return Parameter value, checked against a non-NULL default's type.
#' @keywords internal
#' @noRd
pbrt_get = function(p, name, default = NULL, size = NULL) {
  entry = p$values[[name]]
  if (is.null(entry)) {
    return(default)
  }
  p$used = c(p$used, name)
  value = entry$value
  if (is.numeric(default) && entry$type == "texture") {
    if (length(value) != 1L || is.null(p$context$textures[[value]])) {
      pbrt_error(p$command, paste("Undefined texture for", name))
    }
    texture = p$context$textures[[value]]
    if (!is.null(texture$graph)) {
      pbrt_note(
        p$context,
        p$command,
        paste(
          "Procedural texture replaced by its constant default for unsupported input:",
          name
        )
      )
      value = default
    } else if (is.null(texture$image)) {
      value = texture$value
    } else {
      pbrt_note(
        p$context,
        p$command,
        paste("Image-valued parameter replaced by its constant default:", name)
      )
      value = default
    }
  } else if (
    is.numeric(default) && entry$type %in% c("spectrum", "blackbody")
  ) {
    value = pbrt_spectrum(p, name, default)
  }
  if (!is.null(size) && length(value) != size) {
    pbrt_error(p$command, paste(name, "requires", size, "value(s)."))
  }
  if (
    !is.null(default) &&
      typeof(value) != typeof(default) &&
      !(is.numeric(default) && is.numeric(value))
  ) {
    pbrt_error(p$command, paste("Unexpected parameter type for", name))
  }
  value
}

#' @param buffer An import-local row buffer.
#' @param row One scene fragment or diagnostic data frame.
#' @return NULL after appending the row.
#' @keywords internal
#' @noRd
pbrt_buffer_append = function(buffer, row) {
  if (is.null(row)) {
    return(invisible(NULL))
  }
  if (inherits(row, "pbrt_instance_reference")) {
    buffer$new_instances = buffer$new_instances + 1L
  }
  # Environment-bound list updates copy their container. Bound that work to a
  # small chunk; completed chunks live in separate bindings and never grow.
  # A fresh shallow container also avoids [[<- recursively scanning a shared
  # instance graph to rule out a cycle. The child objects stay shared.
  buffer$pending = c(buffer$pending, list(row))
  if (length(buffer$pending) == 256L) {
    buffer$chunks = buffer$chunks + 1L
    buffer[[as.character(buffer$chunks)]] = pbrt_store_mesh_rows(
      pbrt_batch_curve_rows(buffer$pending),
      buffer$mesh_cache
    )
    buffer$pending = list()
  }
  invisible(NULL)
}

#' @param rows Bounded scene fragments.
#' @param cache Shared asset directory, or NULL for in-memory storage.
#' @return Scene fragments retaining lazy references to eligible inline meshes.
#' @keywords internal
#' @noRd
pbrt_store_mesh_rows = function(rows, cache) {
  rows = pbrt_batch_mesh_rows(rows, file_backed = !is.null(cache))
  if (is.null(cache)) {
    return(rows)
  }
  lapply(rows, function(row) {
    if (
      !inherits(row, 'ray_scene') ||
        nrow(row) != 1L ||
        !identical(row$shape, 'raymesh')
    ) {
      return(row)
    }
    info = row$shape_info[[1]]
    settings = info$shape_properties
    mat = row$material[[1]]
    # Keep transport boundaries, per-face materials and topology operations on
    # their existing path. This cache changes storage, not material semantics.
    if (
      !isTRUE(settings$override_material) ||
        isTRUE(settings$importance_sample_lights) ||
        !identical(settings$subdivision_levels, 1) ||
        !identical(settings$displacement_texture, '') ||
        isTRUE(settings$recalculate_normals) ||
        !is.null(info$medium) ||
        !is.null(mat$subsurface) ||
        isTRUE(mat$fog) ||
        isTRUE(mat$implicit_sample) ||
        (!is.null(mat$openpbr) && mat$openpbr$emission_luminance != 0) ||
        !(mat$type %in% c(1L, 4L, 11L, 12L))
    ) {
      return(row)
    }
    mesh = info$mesh_info[[1]]
    normals = vapply(
      mesh$shapes,
      function(x) all(x$has_vertex_normals),
      logical(1)
    )
    uv = vapply(mesh$shapes, function(x) all(x$has_vertex_tex), logical(1))
    if (!length(normals) || any(normals != normals[1]) || any(uv != uv[1])) {
      return(row)
    }
    for (i in seq_along(mesh$shapes)) {
      shape = mesh$shapes[[i]]
      count = nrow(mesh$vertices[[i]])
      if (
        any(shape$has_vertex_normals) != normals[i] ||
          any(shape$has_vertex_tex) != uv[i] ||
          (normals[i] &&
            (nrow(mesh$normals[[i]]) != count ||
              !identical(shape$indices, shape$norm_indices))) ||
          (uv[i] &&
            (nrow(mesh$texcoords[[i]]) != count ||
              !identical(shape$indices, shape$tex_indices)))
      ) {
        return(row)
      }
    }
    if (
      !dir.exists(cache$directory) &&
        !dir.create(cache$directory, recursive = TRUE)
    ) {
      stop('Unable to create PBRT mesh cache.', call. = FALSE)
    }
    filename = tempfile(
      'mesh-',
      normalizePath(cache$directory),
      fileext = '.ply'
    )
    pbrt_write_mesh_ply(mesh, filename, normals[1], uv[1])
    # Keep material, transforms, animation, orientation and row metadata intact.
    # Only the geometry payload and its loader settings change.
    row$shape = 'ply'
    info$mesh_info = list(NA)
    info$fileinfo = filename
    info$shape_properties = list(
      scale_ply = 1,
      basename = dirname(filename),
      subdivision_levels = 1,
      recalculate_normals = FALSE,
      calculate_consistent_normals = settings$calculate_consistent_normals
    )
    row$shape_info[[1]] = info
    row
  })
}

#' @param mesh Uniformly indexed, compatible mesh shapes.
#' @param filename Destination PLY file.
#' @param normals Whether every shape supplies normals.
#' @param uv Whether every shape supplies texture coordinates.
#' @return NULL after atomically writing a binary PLY asset.
#' @keywords internal
#' @noRd
pbrt_write_mesh_ply = function(mesh, filename, normals, uv) {
  vertex_counts = vapply(mesh$vertices, nrow, integer(1))
  face_counts = vapply(mesh$shapes, function(x) nrow(x$indices), integer(1))
  has_ptex = any(vapply(
    mesh$shapes,
    function(x) !is.null(x$ptex_face_indices),
    logical(1)
  ))
  if (sum(vertex_counts) > .Machine$integer.max) {
    stop('PLY mesh exceeds signed 32-bit index range.', call. = FALSE)
  }
  temporary = paste0(filename, '.partial')
  output = file(temporary, 'wb')
  on.exit({
    if (!is.null(output)) {
      close(output)
    }
    unlink(temporary)
  })
  header = c(
    'ply',
    'format binary_little_endian 1.0',
    paste('element vertex', sum(vertex_counts)),
    paste('property float', c('x', 'y', 'z')),
    if (normals) paste('property float', c('nx', 'ny', 'nz')),
    if (uv) paste('property float', c('u', 'v')),
    paste('element face', sum(face_counts)),
    'property list uchar int vertex_indices',
    if (has_ptex) 'property int face_indices',
    'end_header'
  )
  writeBin(charToRaw(paste0(paste(header, collapse = '\n'), '\n')), output)
  # Bound temporary interleaving buffers, including a single very large shape.
  # PLY's float attributes match the native loader's Float32 input conversion.
  for (i in seq_along(vertex_counts)) {
    for (begin in seq.int(1L, vertex_counts[i], by = 65536L)) {
      selected = seq.int(begin, min(vertex_counts[i], begin + 65535L))
      values = mesh$vertices[[i]][selected, , drop = FALSE]
      if (normals) {
        values = cbind(values, mesh$normals[[i]][selected, , drop = FALSE])
      }
      if (uv) {
        values = cbind(values, mesh$texcoords[[i]][selected, , drop = FALSE])
      }
      writeBin(as.double(t(values)), output, size = 4L, endian = 'little')
    }
  }
  offset = 0L
  for (i in seq_along(face_counts)) {
    for (begin in seq.int(1L, face_counts[i], by = 65536L)) {
      selected = seq.int(begin, min(face_counts[i], begin + 65535L))
      indices = mesh$shapes[[i]]$indices[selected, , drop = FALSE] + offset
      bytes = writeBin(
        as.integer(t(indices)),
        raw(),
        size = 4L,
        endian = 'little'
      )
      # Each face is one uchar count followed by three little-endian int32s.
      records = rbind(as.raw(3L), matrix(bytes, nrow = 12L))
      if (has_ptex) {
        face_indices = mesh$shapes[[i]]$ptex_face_indices
        values = if (is.null(face_indices)) {
          rep(0L, length(selected))
        } else {
          face_indices[selected]
        }
        face_bytes = writeBin(
          as.integer(values),
          raw(),
          size = 4L,
          endian = 'little'
        )
        records = rbind(records, matrix(face_bytes, nrow = 4L))
      }
      writeBin(as.vector(records), output)
    }
    offset = offset + vertex_counts[i]
  }
  close(output)
  output = NULL
  if (!file.rename(temporary, filename)) {
    stop('Unable to finish PBRT mesh asset.', call. = FALSE)
  }
  invisible(NULL)
}

#' @param rows A bounded list of imported scene rows or diagnostics.
#' @return Rows with adjacent compatible nonemissive curves packed numerically.
#' @keywords internal
#' @noRd
pbrt_batch_curve_rows = function(rows) {
  inputs = lapply(rows, pbrt_curve_batch_input)
  result = list()
  begin = 1L
  while (begin <= length(rows)) {
    input = inputs[[begin]]
    end = begin
    count = if (is.null(input)) 0L else ncol(input$values)
    if (!is.null(input)) {
      while (end < length(rows)) {
        next_input = inputs[[end + 1L]]
        if (
          is.null(next_input) ||
            !identical(input$common, next_input$common) ||
            count + ncol(next_input$values) > 4096L
        ) {
          break
        }
        count = count + ncol(next_input$values)
        end = end + 1L
      }
    }
    if (count >= 16L) {
      row = rows[[begin]][1L, ]
      values = do.call(cbind, lapply(inputs[begin:end], `[[`, 'values'))
      row$shape_info[[1]]$shape_properties = list(curve_data = values)
      result = c(result, list(row))
    } else {
      result = c(result, rows[begin:end])
    }
    begin = end + 1L
  }
  result
}

#' @param row An imported scene fragment, possibly a multi-span curve chain.
#' @return Common row settings and a 24-by-span descriptor matrix, or NULL when ineligible.
#' @keywords internal
#' @noRd
pbrt_curve_batch_input = function(row) {
  if (!inherits(row, 'ray_scene') || !nrow(row) || !all(row$shape == 'curve')) {
    return(NULL)
  }
  materials = vctrs::vec_data(row$material)
  # Diffuse/rough diffuse, hair, OpenPBR and translucent surfaces only. Keep
  # other transport families and all emitting/medium-bearing rows independent.
  eligible = vapply(
    materials,
    function(material) {
      material$type %in%
        c(1L, 4L, 9L, 11L, 12L) &&
        !material$fog &&
        !material$implicit_sample &&
        is.null(material$subsurface) &&
        (is.null(material$openpbr) || material$openpbr$emission_luminance == 0)
    },
    logical(1)
  )
  if (!all(eligible)) {
    return(NULL)
  }
  info = vctrs::vec_data(row$shape_info)
  if (any(vapply(info, function(x) !is.null(x$medium), logical(1)))) {
    return(NULL)
  }
  # Drop only per-span geometry before checking every remaining scene field.
  # In particular, chain rows may not silently inherit another row's transform,
  # animation, material, flip, texture or medium settings.
  properties = lapply(info, `[[`, 'shape_properties')
  for (i in seq_along(info)) {
    info[[i]]$shape_properties = list()
  }
  common = row
  common$shape_info = vctrs::vec_restore(info, row$shape_info)
  if (nrow(common) > 1L) {
    common = vctrs::vec_unique(common)
  }
  if (nrow(common) != 1L) {
    return(NULL)
  }
  fields = c(
    'p1',
    'p2',
    'p3',
    'p4',
    'width',
    'width_end',
    'u_min',
    'u_max',
    'split_depth',
    'curvetype',
    'normal',
    'normal_end'
  )
  if (
    !all(vapply(
      properties,
      function(x) {
        identical(names(x), fields) || identical(names(x), 'curve_data')
      },
      logical(1)
    ))
  ) {
    return(NULL)
  }
  # Column-major schema shared with native CURVE construction: control points
  # 1:12, widths 13:14, parameter bounds 15:16, depth/type 17:18, normals 19:24.
  values = do.call(
    cbind,
    lapply(properties, function(x) {
      if (!is.null(x$curve_data)) {
        return(x$curve_data)
      }
      matrix(as.numeric(unlist(x[fields], use.names = FALSE)), nrow = 24L)
    })
  )
  list(common = common, values = values)
}

#' @param rows A bounded list of scene rows or diagnostics.
#' @param file_backed Default `FALSE`. Also batch compatible OpenPBR/translucent surfaces.
#' @return Rows with adjacent compatible diffuse triangle meshes batched.
#' @keywords internal
#' @noRd
pbrt_batch_mesh_rows = function(rows, file_backed = FALSE) {
  if (length(rows) < 16L) {
    return(rows)
  }
  result = list()
  begin = 1L
  while (begin <= length(rows)) {
    row = rows[[begin]]
    end = begin
    compatible = NULL
    if (
      inherits(row, 'ray_scene') &&
        nrow(row) == 1L &&
        identical(row$shape, 'raymesh')
    ) {
      info = row$shape_info[[1]]
      settings = info$shape_properties
      mat = row$material[[1]]
      mesh = info$mesh_info[[1]]
      normal_mode = vapply(
        mesh$shapes,
        function(shape) {
          if (all(shape$has_vertex_normals)) {
            1L
          } else if (!any(shape$has_vertex_normals)) {
            0L
          } else {
            -1L
          }
        },
        integer(1)
      )
      # Keep transport boundaries and topology-changing operations separate.
      # Mixing complete and missing normals would disable the consistency table
      # on formerly smooth meshes, so each batch has one normal-coverage mode.
      if (
        mat$type %in%
          c(
            get_material_enum('diffuse'),
            get_material_enum('oren-nayar'),
            if (file_backed) {
              c(get_material_enum('openpbr'), get_material_enum('translucent'))
            }
          ) &&
          is.null(mat$subsurface) &&
          is.null(info$medium) &&
          !isTRUE(mat$fog) &&
          !isTRUE(mat$implicit_sample) &&
          (is.null(mat$openpbr) || mat$openpbr$emission_luminance == 0) &&
          isTRUE(settings$override_material) &&
          !isTRUE(settings$importance_sample_lights) &&
          identical(settings$subdivision_levels, 1) &&
          identical(settings$displacement_texture, '') &&
          !isTRUE(settings$recalculate_normals) &&
          length(normal_mode) &&
          all(normal_mode == normal_mode[1]) &&
          normal_mode[1] >= 0L
      ) {
        shape_count = length(mesh$shapes)
        compatible = row
        compatible$shape_info[[1]]$mesh_info = list(NULL)
        while (end < length(rows)) {
          next_row = rows[[end + 1L]]
          if (
            !inherits(next_row, 'ray_scene') ||
              nrow(next_row) != 1L ||
              !identical(next_row$shape, 'raymesh')
          ) {
            break
          }
          next_mesh = next_row$shape_info[[1]]$mesh_info[[1]]
          if (shape_count + length(next_mesh$shapes) > 256L) {
            break
          }
          next_row$shape_info[[1]]$mesh_info = list(NULL)
          if (!identical(compatible, next_row)) {
            break
          }
          same_normals = all(vapply(
            next_mesh$shapes,
            function(shape) {
              if (normal_mode[1] == 1L) {
                all(shape$has_vertex_normals)
              } else {
                !any(shape$has_vertex_normals)
              }
            },
            logical(1)
          ))
          if (!same_normals) {
            break
          }
          shape_count = shape_count + length(next_mesh$shapes)
          end = end + 1L
        }
      }
    }
    if (end - begin + 1L >= 16L) {
      meshes = lapply(rows[begin:end], function(x) {
        x$shape_info[[1]]$mesh_info[[1]]
      })
      combined = meshes[[1]]
      # Concatenate each container once; retain separate vertices and indices.
      # No welding, reindexing, per-shape material changes or repeated rbind.
      for (field in names(combined)) {
        combined[[field]] = vctrs::list_unchop(lapply(meshes, `[[`, field))
      }
      attr(combined, 'material_hashes') = unlist(
        lapply(meshes, attr, 'material_hashes'),
        use.names = FALSE
      )
      row$shape_info[[1]]$mesh_info[[1]] = combined
    } else {
      end = begin
    }
    result = c(result, list(row))
    begin = end + 1L
  }
  result
}

#' @param context Import resources.
#' @param command Source directive.
#' @param message Diagnostic text describing the actual fallback.
#' @return NULL, or an error in strict mode.
#' @keywords internal
#' @noRd
pbrt_note = function(context, command, message) {
  if (context$strict) {
    pbrt_error(
      command,
      paste0(message, " Use strict = FALSE to allow this conversion.")
    )
  }
  pbrt_buffer_append(
    context$diagnostics,
    data.frame(
      file = command$file,
      line = command$line,
      directive = command$name,
      message = message
    )
  )
  invisible(NULL)
}

#' @param p Parameter reader.
#' @return NULL, after reporting unconsumed fields.
#' @keywords internal
#' @noRd
pbrt_unused = function(p) {
  unused = setdiff(names(p$values), p$used)
  if (length(unused)) {
    pbrt_note(
      p$context,
      p$command,
      paste("Unsupported parameter(s) ignored:", paste(unused, collapse = ", "))
    )
  }
  invisible(NULL)
}

#' @param command Source directive.
#' @param message Error text.
#' @return Does not return.
#' @keywords internal
#' @noRd
pbrt_error = function(command, message) {
  stop(
    command$file,
    ":",
    command$line,
    ": ",
    command$name,
    ": ",
    message,
    call. = FALSE
  )
}

#' @param filename Asset name relative to the main scene, or an absolute path.
#' @param context Import resources.
#' @param command Source directive.
#' @return Existing canonical absolute path.
#' @keywords internal
#' @noRd
pbrt_asset = function(filename, context, command) {
  if (!is.character(filename) || length(filename) != 1 || !nzchar(filename)) {
    pbrt_error(command, "Expected one nonempty asset filename.")
  }
  filename = path.expand(filename)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", filename)) {
    filename = file.path(context$root, filename)
  }
  if (!file.exists(filename) || dir.exists(filename)) {
    pbrt_error(command, paste("Asset not found:", filename))
  }
  normalizePath(filename, mustWork = TRUE)
}

#' @param tokens PBRT tokens, possibly quoted.
#' @param command Source directive.
#' @return Strings decoded without evaluating R expressions.
#' @keywords internal
#' @noRd
pbrt_unquote = function(tokens, command) {
  quoted = startsWith(tokens, '"')
  tokens[quoted] = substring(tokens[quoted], 2L, nchar(tokens[quoted]) - 1L)
  escaped = which(quoted & grepl("\\", tokens, fixed = TRUE))
  for (i in escaped) {
    chars = strsplit(tokens[i], "", fixed = TRUE)[[1]]
    result = character(length(chars))
    at = 1L
    to = 0L
    while (at <= length(chars)) {
      ch = chars[at]
      if (ch == "\\") {
        at = at + 1L
        if (at > length(chars)) {
          pbrt_error(command, "Incomplete string escape.")
        }
        ch = switch(
          chars[at],
          n = "\n",
          r = "\r",
          t = "\t",
          b = "\b",
          f = "\f",
          '"' = '"',
          "\\" = "\\",
          NULL
        )
        if (is.null(ch)) pbrt_error(command, "Unsupported string escape.")
      }
      to = to + 1L
      result[to] = ch
      at = at + 1L
    }
    tokens[i] = paste0(result[seq_len(to)], collapse = "")
  }
  tokens
}

#' @param command Transform directive.
#' @return Invertible affine 4x4 transform using column vectors.
#' @keywords internal
#' @noRd
pbrt_transform = function(command) {
  value = suppressWarnings(as.numeric(command$args))
  if (any(!is.finite(value))) {
    pbrt_error(command, "Transform operands must be finite numbers.")
  }
  m = diag(4)
  if (command$name == "Translate") {
    m[1:3, 4] = value
  } else if (command$name == "Scale") {
    diag(m) = c(value, 1)
  } else if (command$name == "Rotate") {
    axis = value[2:4]
    if (sum(axis^2) == 0) {
      pbrt_error(command, "Rotation axis cannot be zero.")
    }
    axis = axis / sqrt(sum(axis^2))
    a = value[1] * pi / 180
    cross = matrix(
      c(0, axis[3], -axis[2], -axis[3], 0, axis[1], axis[2], -axis[1], 0),
      3
    )
    m[1:3, 1:3] = cos(a) *
      diag(3) +
      (1 - cos(a)) * tcrossprod(axis) +
      sin(a) * cross
  } else if (command$name == "LookAt") {
    eye = value[1:3]
    forward = value[4:6] - eye
    up = value[7:9]
    if (sum(forward^2) == 0 || sum(up^2) == 0) {
      pbrt_error(command, "Degenerate LookAt frame.")
    }
    forward = forward / sqrt(sum(forward^2))
    right = c(
      up[2] * forward[3] - up[3] * forward[2],
      up[3] * forward[1] - up[1] * forward[3],
      up[1] * forward[2] - up[2] * forward[1]
    )
    if (sum(right^2) < 1e-20 * sum(up^2)) {
      pbrt_error(command, "LookAt up is parallel to the view direction.")
    }
    right = right / sqrt(sum(right^2))
    up = c(
      forward[2] * right[3] - forward[3] * right[2],
      forward[3] * right[1] - forward[1] * right[3],
      forward[1] * right[2] - forward[2] * right[1]
    )
    m[1:3, ] = cbind(right, up, forward, eye)
    m = solve(m) # PBRT stores world-to-camera at Camera directives.
  } else if (command$name %in% c("Transform", "ConcatTransform")) {
    m = matrix(value, 4, 4) # PBRT's on-disk matrices are column-major.
  }
  if (
    any(!is.finite(m)) ||
      max(abs(m[4, ] - c(0, 0, 0, 1))) > 1e-10 ||
      rcond(m) < .Machine$double.eps
  ) {
    pbrt_error(command, "Transform must be a nonsingular affine matrix.")
  }
  m
}

#' @param command Material directive.
#' @param state Current graphics state.
#' @param context Import resources.
#' @return A ray material and an interface-only flag.
#' @keywords internal
#' @noRd
pbrt_material = function(command, state, context) {
  p = pbrt_parameters(command, context, state$defaults$material)
  type = if (command$name == "MakeNamedMaterial") {
    pbrt_get(p, "type", "diffuse", 1L)
  } else {
    command$args[1]
  }
  roughness = NULL
  color = NULL
  interface = type %in% c("", "none", "interface")
  material = diffuse(color = rep(.5, 3))
  bump_name = if ("displacement" %in% names(p$values)) {
    "displacement"
  } else {
    "bumpmap"
  }
  bump = pbrt_texture_parameter(p, bump_name, 0, context, allow_graph = TRUE)
  if (is.null(bump$image) && is.null(bump$graph) && any(bump$value != 0)) {
    bump$image = pbrt_write_asset(
      array(mean(bump$value), c(2, 2, 3)),
      ".exr",
      context
    )
  }
  # Texture evaluation has already decoded and scaled the height into linear
  # EXR. Preserve its precision and absolute values: PBRT's curved-surface
  # tangent construction also uses height * dndu/dndv, not only height slopes.
  bump_args = if (!is.null(bump$image)) {
    list(bump_texture = bump$image, bump_intensity = 1)
  } else {
    list()
  }
  if (
    type %in% c("diffuse", "matte", "coateddiffuse", "plastic", "substrate")
  ) {
    color = pbrt_texture_parameter(
      p,
      if (type %in% c("matte", "plastic", "substrate")) "Kd" else "reflectance",
      rep(.5, 3),
      context,
      allow_graph = TRUE
    )
    color_args = if (!is.null(color$image)) {
      list(
        image_texture = color$image,
        image_repeat = color$uv_repeat,
        image_offset = color$uv_offset
      )
    } else {
      list()
    }
    if (type %in% c("diffuse", "matte")) {
      sigma = pbrt_get(p, "sigma", 0, 1L)
      if (sigma != 0) {
        pbrt_note(
          context,
          command,
          "PBRT Oren-Nayar is converted to rayrender's energy-preserving rough diffuse model."
        )
      }
      material = do.call(
        diffuse,
        c(
          list(color = pbrt_rgb(color$value, p), sigma = sigma),
          color_args,
          bump_args
        )
      )
    } else {
      roughness = pbrt_roughness(p, if (type == "coateddiffuse") 0 else .1)
      eta = pbrt_get(p, "eta", 1.5, 1L)
      pbrt_note(
        context,
        command,
        "Layered material approximated with OpenPBR's dielectric coat over diffuse base."
      )
      material = do.call(
        openpbr,
        c(
          list(
            base_color = pbrt_rgb(color$value, p),
            specular_weight = 0,
            coat_weight = 1,
            coat_ior = eta,
            coat_roughness = mean(roughness)
          ),
          color_args,
          bump_args
        )
      )
    }
  } else if (type == "diffusetransmission") {
    reflection = pbrt_texture_parameter(
      p,
      "reflectance",
      rep(.25, 3),
      context,
      allow_graph = TRUE
    )
    transmission = pbrt_texture_parameter(
      p,
      "transmittance",
      rep(.25, 3),
      context,
      allow_graph = TRUE
    )
    scale = pbrt_get(p, "scale", 1, 1L)
    if (!is.finite(scale) || scale < 0) {
      pbrt_error(
        command,
        "Diffuse transmission scale must be finite and nonnegative."
      )
    }
    # Bake only scalar radiometric scaling, retaining each texture's independent
    # UV repeat/offset for native lookup. The two BSDF lobes remain independent.
    lobes = list(reflection, transmission)
    for (i in seq_along(lobes)) {
      if (!is.null(lobes[[i]]$image) && scale != 1) {
        image = rayimage::ray_read_image(lobes[[i]]$image)[,, 1:3, drop = FALSE]
        lobes[[i]]$image = pbrt_write_asset(
          pmin(pmax(image * scale, 0), 1),
          ".exr",
          context
        )
      }
      lobes[[i]]$value = pmin(
        pmax(rep(lobes[[i]]$value, length.out = 3L) * scale, 0),
        1
      )
    }
    reflection = lobes[[1]]
    transmission = lobes[[2]]
    if (!is.null(reflection$image)) {
      reflection$value = rep(0, 3)
    }
    if (!is.null(transmission$image)) {
      transmission$value = rep(0, 3)
    }
    total = reflection$value + transmission$value
    if (any(total > 1)) {
      pbrt_note(
        context,
        command,
        "Diffuse reflection plus transmission exceeds one; combined energy is normalized per channel."
      )
      reflection$value = reflection$value / pmax(total, 1)
      transmission$value = transmission$value / pmax(total, 1)
    }
    material = translucent(
      reflectance = reflection$value,
      transmittance = transmission$value,
      image_texture = if (is.null(reflection$image)) "" else reflection$image,
      transmission_texture = if (is.null(transmission$image)) {
        ""
      } else {
        transmission$image
      },
      image_repeat = reflection$uv_repeat,
      image_offset = reflection$uv_offset
    )
    material[[1]]$transmission_repeat = transmission$uv_repeat
    material[[1]]$transmission_offset = transmission$uv_offset
    if (!is.null(reflection$graph)) {
      material[[1]]$texture_graphs$color = texture_scale(
        reflection$graph,
        scale
      )
    }
    if (!is.null(transmission$graph)) {
      material[[1]]$texture_graphs$transmission = texture_scale(
        transmission$graph,
        scale
      )
    }
    if (!is.null(bump$image) || !is.null(bump$graph)) {
      pbrt_note(
        context,
        command,
        "Diffuse transmission uses geometric normals; its bump map is omitted."
      )
    }
  } else if (type == "disney") {
    material = pbrt_disney_material(p, context, bump_args)
  } else if (type %in% c("dielectric", "glass", "thindielectric")) {
    roughness = pbrt_roughness(p, 0)
    eta = pbrt_get(p, "eta", 1.5, 1L)
    if (eta <= 0) {
      pbrt_error(command, "Dielectric eta must be positive.")
    }
    if (type == "thindielectric") {
      pbrt_note(
        context,
        command,
        "Thin dielectric is approximated by OpenPBR's thin-walled transmission model."
      )
      material = do.call(
        openpbr,
        c(
          list(
            transmission_weight = 1,
            specular_ior = eta,
            specular_roughness = mean(roughness),
            geometry_thin_walled = TRUE
          ),
          bump_args
        )
      )
    } else if (all(roughness == 0)) {
      material = do.call(dielectric, c(list(refraction = eta), bump_args))
    } else {
      material = do.call(
        microfacet,
        c(
          list(transmission = TRUE, eta = eta, roughness = roughness),
          bump_args
        )
      )
    }
  } else if (type == "coatedconductor") {
    conductor = pbrt_roughness(p, 0, prefix = "conductor.")
    coat = pbrt_roughness(p, 0, prefix = "interface.")
    eta = pbrt_get(p, "interface.eta", 1.5, 1L)
    color = pbrt_texture_parameter(p, "conductor.reflectance", NULL, context)
    if (is.null(color$value)) {
      n = pbrt_spectrum(p, "conductor.eta", c(.2, .92, 1.1))
      k = pbrt_spectrum(p, "conductor.k", c(3.91, 2.45, 2.14))
      color$value = ((n - 1)^2 + k^2) / ((n + 1)^2 + k^2)
    }
    pbrt_note(
      context,
      command,
      "Coated conductor approximated with OpenPBR's dielectric coat over a metallic base; PBRT's layered multiple-scattering model differs."
    )
    material = do.call(
      openpbr,
      c(
        list(
          base_metalness = 1,
          base_color = pbrt_rgb(color$value, p),
          specular_roughness = mean(conductor),
          coat_weight = 1,
          coat_ior = eta,
          coat_roughness = mean(coat)
        ),
        if (!is.null(color$image)) {
          list(
            image_texture = color$image,
            image_repeat = color$uv_repeat,
            image_offset = color$uv_offset
          )
        } else {
          list()
        },
        bump_args
      )
    )
  } else if (type == "hair") {
    args = list(
      eta = pbrt_get(p, "eta", 1.55, 1L),
      beta_m = pbrt_get(p, "beta_m", .3, 1L),
      beta_n = pbrt_get(p, "beta_n", .3, 1L),
      alpha = pbrt_get(p, "alpha", 2, 1L)
    )
    if ("sigma_a" %in% names(p$values)) {
      args$sigma_a = pbrt_spectrum(p, "sigma_a", rep(0, 3))
    } else if ("color" %in% names(p$values)) {
      args$color = pbrt_rgb(pbrt_spectrum(p, "color", rep(.5, 3)), p)
    } else {
      args$pigment = pbrt_get(p, "eumelanin", 1.3, 1L)
      args$red_pigment = pbrt_get(p, "pheomelanin", 0, 1L)
    }
    material = do.call(hair, args)
  } else if (type %in% c("conductor", "metal", "mirror")) {
    roughness = pbrt_roughness(p, 0)
    if ("reflectance" %in% names(p$values) || type == "mirror") {
      color = pbrt_texture_parameter(
        p,
        if (type == "mirror") "Kr" else "reflectance",
        rep(.9, 3),
        context
      )
      reflectance = pbrt_rgb(color$value, p)
      eta = rep(1, 3)
      k = 2 * sqrt(pmin(reflectance, .9999) / pmax(1 - reflectance, .0001))
      if (!is.null(color$image)) {
        pbrt_note(
          context,
          command,
          "Textured conductor reflectance is reduced to its fallback constant value."
        )
      }
    } else {
      eta = pbrt_spectrum(p, "eta", c(.2, .92, 1.1))
      k = pbrt_spectrum(p, "k", c(3.91, 2.45, 2.14))
      if (!"eta" %in% names(p$values)) {
        pbrt_note(
          context,
          command,
          "Default copper spectrum approximated with RGB optical constants."
        )
      }
    }
    material = do.call(
      if (all(roughness == 0)) metal else microfacet,
      c(
        list(eta = eta, kappa = k),
        if (all(roughness == 0)) list() else list(roughness = roughness),
        bump_args
      )
    )
  } else if (type == "subsurface") {
    roughness = pbrt_roughness(p, 0)
    eta = pbrt_get(p, "eta", 1.33, 1L)
    sa = pbrt_spectrum(p, "sigma_a", rep(.0011, 3))
    ss = pbrt_spectrum(p, "sigma_s", rep(2.55, 3))
    scale = pbrt_get(p, "scale", 1, 1L)
    pbrt_note(
      context,
      command,
      "PBRT's diffusion BSSRDF is replaced by rayrender's volumetric random-walk subsurface model."
    )
    material = subsurface(
      sigma_a = sa * scale,
      sigma_s = ss * scale,
      g = pbrt_get(p, "g", 0, 1L),
      refraction = eta,
      roughness = mean(roughness)
    )
  } else if (!interface) {
    pbrt_note(
      context,
      command,
      paste("Unsupported material replaced by neutral diffuse:", type)
    )
  }
  if (!is.null(color$graph)) {
    material[[1]]$texture_graphs$color = color$graph
  }
  alpha_graphs = attr(roughness, "alpha_graphs")
  if (!is.null(alpha_graphs)) {
    if (type %in% c("dielectric", "glass", "conductor", "metal")) {
      material[[1]]$texture_graphs$roughness = alpha_graphs[[1]]
      material[[1]]$texture_graphs$roughness_v = alpha_graphs[[2]]
      material[[1]]$texture_graphs$roughness_is_alpha = TRUE
    } else if (type == "thindielectric") {
      # OpenPBR uses perceptual roughness; anisotropic thin-glass import keeps
      # the existing mean approximation and reports that limitation.
      material[[1]]$texture_graphs$roughness = texture_node(
        "power",
        "scalar",
        child = texture_mix(alpha_graphs[[1]], alpha_graphs[[2]]),
        exponent = .5
      )
    } else {
      pbrt_note(
        context,
        command,
        paste(
          "Procedural roughness is not yet supported for PBRT material:",
          type
        )
      )
    }
  }
  if (!is.null(bump$image)) {
    material[[1]]$texture_offsets$bump = bump$uv_offset
    material[[1]]$bump_repeat = bump$uv_repeat
  }
  if (!is.null(bump$graph) && type != "diffusetransmission") {
    material[[1]]$texture_graphs$bump = texture_scalar(bump$graph)
    material[[1]]$bump_intensity = 1
  }
  pbrt_unused(p)
  list(material = material, interface = interface)
}

#' @param command Texture declaration.
#' @param state Current graphics state.
#' @param context Import resources.
#' @return Constant value or image descriptor with UV repeat.
#' @keywords internal
#' @noRd
pbrt_texture = function(command, state, context) {
  p = pbrt_parameters(command, context, state$defaults$texture)
  type = command$args[3]
  result = list(
    value = if (command$args[2] == "float") .5 else rep(.5, 3),
    uv_repeat = c(1, 1),
    uv_offset = c(0, 0)
  )
  if (!command$args[2] %in% c("float", "spectrum", "color")) {
    pbrt_error(command, "Unknown texture value type.")
  }
  if (type == "ptex") {
    filename = pbrt_asset(pbrt_get(p, "filename", "", 1L), context, command)
    encoding = pbrt_get(p, "encoding", "gamma 2.2", 1L)
    scale = pbrt_get(p, "scale", 1, 1L)
    node = texture_ptex(
      filename,
      type = if (command$args[2] == "float") "scalar" else "color",
      encoding = encoding
    )
    node$pbrt_encoding = TRUE
    result$graph = texture_scale(node, scale)
  } else if (type == "directionmix") {
    default = if (command$args[2] == "float") 0 else rep(0, 3)
    a = pbrt_texture_parameter(p, "tex1", default, context, allow_graph = TRUE)
    b = pbrt_texture_parameter(
      p,
      "tex2",
      default + 1,
      context,
      allow_graph = TRUE
    )
    direction = as.vector(
      state$transform[1:3, 1:3] %*% pbrt_get(p, "dir", c(0, 1, 0), 3L)
    )
    result$graph = texture_direction_mix(
      pbrt_texture_graph(
        a,
        if (command$args[2] == "float") "scalar" else "color"
      ),
      pbrt_texture_graph(
        b,
        if (command$args[2] == "float") "scalar" else "color"
      ),
      direction = direction
    )
  } else if (type == "constant") {
    result = pbrt_texture_parameter(
      p,
      "value",
      if (command$args[2] == "float") 1 else rep(1, 3),
      context
    )
  } else if (type == "scale") {
    result = pbrt_texture_parameter(p, "tex", 1, context, allow_graph = TRUE)
    scale = pbrt_get(p, "scale", 1, 1L)
    if (!is.null(result$graph)) {
      result$graph = texture_scale(result$graph, scale)
    }
    result$value = result$value * scale
    if (!is.null(result$image)) {
      image = rayimage::ray_read_image(result$image, convert_to_array = TRUE)
      result$image = pbrt_write_asset(image * scale, ".exr", context)
    }
  } else if (type == "imagemap") {
    filename = pbrt_asset(pbrt_get(p, "filename", "", 1L), context, command)
    mapping = pbrt_get(p, "mapping", "uv", 1L)
    if (mapping != "uv") {
      pbrt_note(
        context,
        command,
        "Non-UV image mapping replaced by UV mapping."
      )
    }
    result$uv_repeat = c(
      pbrt_get(p, "uscale", 1, 1L),
      pbrt_get(p, "vscale", 1, 1L)
    )
    offset = c(
      pbrt_get(p, "udelta", 0, 1L),
      pbrt_get(p, "vdelta", 0, 1L)
    )
    if (any(!is.finite(offset))) {
      pbrt_error(command, "Image UV offsets must be finite.")
    }
    wrap = pbrt_get(p, "wrap", "repeat", 1L)
    if (wrap != "repeat") {
      pbrt_note(context, command, "Image wrap mode replaced by repeat.")
    }
    pbrt_get(p, "filter")
    pbrt_get(p, "maxanisotropy")
    encoding = pbrt_get(
      p,
      "encoding",
      if (tolower(tools::file_ext(filename)) %in% c("exr", "hdr", "pfm")) {
        "linear"
      } else {
        "sRGB"
      },
      1L
    )
    if (!encoding %in% c("sRGB", "linear")) {
      pbrt_note(
        context,
        command,
        "Unsupported image encoding replaced by sRGB."
      )
    }
    image = pbrt_read_image(
      filename,
      source_linear = encoding == "linear"
    )
    result$uv_offset = offset
    scale = pbrt_get(p, "scale", 1, 1L)
    if (any(!is.finite(image)) || !is.finite(scale)) {
      pbrt_error(
        command,
        "Image values and texture scale must be finite."
      )
    }
    # Decode once to linear EXR, avoiding rayrender's approximate LDR gamma
    # decoder and allowing both PBRT's explicit linear and sRGB encodings.
    # PBRT float imagemaps use a nontrivial alpha channel when present. RGB
    # float maps are averaged when consumed as an opacity mask below.
    if (
      command$args[2] == "float" && dim(image)[3] == 4L && any(image[,, 4] != 1)
    ) {
      for (channel in 1:3) {
        image[,, channel] = image[,, 4]
      }
    }
    result$image = pbrt_write_asset(
      image[,, 1:3, drop = FALSE] * scale,
      ".exr",
      context
    )
  } else if (type == "checkerboard") {
    dimension = pbrt_get(p, "dimension", 2, 1L)
    mapping = pbrt_get(p, "mapping", "uv", 1L)
    if (dimension != 2 || mapping != "uv") {
      pbrt_note(
        context,
        command,
        "Checkerboard is baked as a 2D UV checker texture."
      )
    }
    a = pbrt_texture_parameter(
      p,
      "tex1",
      rep(1, 3),
      context,
      allow_graph = TRUE
    )
    b = pbrt_texture_parameter(
      p,
      "tex2",
      rep(0, 3),
      context,
      allow_graph = TRUE
    )
    nested = !is.null(a$image) ||
      !is.null(b$image) ||
      !is.null(a$graph) ||
      !is.null(b$graph)
    uscale = pbrt_get(p, "uscale", 1, 1L)
    vscale = pbrt_get(p, "vscale", 1, 1L)
    udelta = pbrt_get(p, "udelta", 0, 1L)
    vdelta = pbrt_get(p, "vdelta", 0, 1L)
    pbrt_get(p, "aamode")
    if (nested) {
      output_type = if (command$args[2] == "float") "scalar" else "color"
      result$graph = texture_checker(
        pbrt_texture_graph(a, output_type),
        pbrt_texture_graph(b, output_type),
        coordinates = texture_coordinates(
          "uv",
          scale = c(uscale, vscale, 1),
          offset = c(udelta, vdelta, 0)
        )
      )
    }
    # Bake the unit UV square at a fixed 256 pixels. This bounds asset size;
    # fractional offsets and negative scales are included before sampling.
    uv = (seq_len(256) - .5) / 256
    parity = outer(
      floor(rev(uv) * vscale + vdelta),
      floor(uv * uscale + udelta),
      "+"
    ) %%
      2
    values = array(0, c(256, 256, 3))
    a = rep(a$value, length.out = 3)
    b = rep(b$value, length.out = 3)
    for (channel in 1:3) {
      values[,, channel] = ifelse(parity == 0, a[channel], b[channel])
    }
    result$image = pbrt_write_asset(values, ".exr", context)
  } else {
    pbrt_note(
      context,
      command,
      paste("Unsupported texture replaced by a constant:", type)
    )
  }
  if (any(result$uv_repeat <= 0)) {
    pbrt_note(
      context,
      command,
      "Nonpositive texture repeat replaced by positive repeat."
    )
    result$uv_repeat = pmax(abs(result$uv_repeat), 1e-8)
  }
  pbrt_unused(p)
  result
}

#' @param command Medium declaration.
#' @param state Current graphics state.
#' @param context Import resources.
#' @return A homogeneous medium or NULL for unsupported media.
#' @keywords internal
#' @noRd
pbrt_medium = function(command, state, context) {
  p = pbrt_parameters(command, context, state$defaults$medium)
  type = pbrt_get(p, "type", "homogeneous", 1L)
  if (!type %in% c("homogeneous", "nanovdb", "uniformgrid")) {
    pbrt_note(context, command, paste("Unsupported medium omitted:", type))
    return(list(medium = NULL))
  }
  args = list(
    sigma_a = pbrt_spectrum(p, "sigma_a", rep(1, 3)),
    sigma_s = pbrt_spectrum(p, "sigma_s", rep(1, 3)),
    density_scale = pbrt_get(p, "scale", 1, 1L),
    g = pbrt_get(p, "g", 0, 1L),
    emission = pbrt_spectrum(p, "Le", rep(0, 3)),
    emission_scale = pbrt_get(
      p,
      if ("LeScale" %in% names(p$values)) "LeScale" else "Lescale",
      1,
      if (type == "uniformgrid") NULL else 1L
    )
  )
  if (type != "homogeneous") {
    # Both renderers apply scale * (temperature - offset), in kelvin.
    args$temperature_scale = pbrt_get(p, "temperaturescale", 1, 1L)
    offset = pbrt_get(
      p,
      "temperatureoffset",
      pbrt_get(p, "temperaturecutoff", 0, 1L),
      1L
    )
    args$temperature_offset = offset
  }
  if (type == "nanovdb") {
    args$filename = pbrt_asset(
      pbrt_get(p, "filename", "", 1L),
      context,
      command
    )
    args$density_grid = pbrt_get(p, "gridname", "density", 1L)
    args$temperature_grid = pbrt_get(p, "temperaturename", "temperature", 1L)
    # PBRT's PixelSensor integrates CIE response without dividing by its Y
    # integral. BlackbodyRGB() divides by 106.856895, unlike PBRT's blackbody
    # volume emission (which bypasses its usual photometric light normalization).
    # Convert here so native rayrender temperature units remain unchanged.
    args$emission_scale = args$emission_scale * 106.856895
    result = do.call(nanovdb_medium, args)
    result$temperature_optional = TRUE
    pbrt_note(
      context,
      command,
      "NanoVDB temperature emission uses rayrender's RGB blackbody model; files must contain float grids, with density required and temperature optional. Uncompressed, ZIP and BLOSC files are supported through openvdbr."
    )
  } else if (type == "uniformgrid") {
    dims = vapply(
      c("nx", "ny", "nz"),
      function(name) pbrt_get(p, name, 1, 1L),
      numeric(1)
    )
    if (any(dims < 1 | dims %% 1 != 0)) {
      pbrt_error(command, "Grid dimensions must be positive integers.")
    }
    density = pbrt_get(p, "density", size = prod(dims))
    if (is.null(density)) {
      pbrt_error(command, "uniformgrid requires density samples.")
    }
    args$density = array(density, dims)
    if (length(args$emission_scale) != 1L) {
      if (length(args$emission_scale) != prod(dims)) {
        pbrt_error(command, "Lescale requires one value or one per grid cell.")
      }
      args$emission_scale = array(args$emission_scale, dims)
    }
    args$bounds = rbind(
      pbrt_get(p, "p0", c(0, 0, 0), 3L),
      pbrt_get(p, "p1", c(1, 1, 1), 3L)
    )
    temperature = pbrt_get(p, "temperature", size = prod(dims))
    if (!is.null(temperature)) {
      args$temperature = array(temperature, dims)
      args$emission_scale = args$emission_scale * 106.856895
    }
    result = do.call(grid_medium, args)
  } else {
    result = do.call(homogeneous_medium, args)
  }
  pbrt_unused(p)
  if (!isTRUE(all.equal(state$transform, state$end_transform))) {
    pbrt_note(
      context,
      command,
      "Medium's definition transform is frozen at StartTime."
    )
  }
  list(medium = result, transform = state$transform)
}

#' @param command Light directive.
#' @param state Current graphics state.
#' @param context Import resources and light accumulator.
#' @return NULL, after adding a light to the import.
#' @keywords internal
#' @noRd
pbrt_light = function(command, state, context) {
  if (!isTRUE(all.equal(state$transform, state$end_transform))) {
    pbrt_note(
      context,
      command,
      "Animated light transforms are frozen at StartTime."
    )
  }
  p = pbrt_parameters(command, context, state$defaults$light)
  type = command$args[1]
  if (type == "infinite") {
    filename = pbrt_get(p, "filename", "", 1L)
    radiance = pbrt_spectrum(p, "L", rep(1, 3)) * pbrt_get(p, "scale", 1, 1L)
    if (nzchar(filename)) {
      source = pbrt_asset(filename, context, command)
      image = pbrt_read_environment(source)
      if (
        length(dim(image)) != 3L || dim(image)[3] < 3L || any(!is.finite(image))
      ) {
        pbrt_error(command, "Environment image must contain finite RGB pixels.")
      }
      image = pbrt_environment_image(
        image[,, 1:3, drop = FALSE],
        state$transform
      )
      for (channel in 1:3) {
        image[,, channel] = image[,, channel] * radiance[channel]
      }
      light = infinite_light(pbrt_write_asset(image, ".exr", context))
      pbrt_buffer_append(context$lights, light)
    } else {
      image = array(rep(radiance, each = 8), c(2, 4, 3))
      light = infinite_light(pbrt_write_asset(image, ".exr", context))
      pbrt_buffer_append(context$lights, light)
    }
  } else if (type == "distant") {
    radiance = pbrt_spectrum(p, "L", rep(1, 3)) * pbrt_get(p, "scale", 1, 1L)
    from = pbrt_get(p, "from", c(0, 0, 0), 3L)
    to = pbrt_get(p, "to", c(0, 0, 1), 3L)
    direction = as.vector(state$transform[1:3, 1:3] %*% (from - to))
    if (sum(direction^2) == 0) {
      pbrt_error(command, "Distant light direction is zero.")
    }
    pbrt_note(
      context,
      command,
      "Delta distant light approximated by a 0.53-degree disk with matched integrated irradiance."
    )
    solid_angle = 2 * pi * (1 - cos(.53 * pi / 360))
    intensity = max(radiance)
    light = disk_light(
      color = if (intensity > 0) radiance / intensity else rep(0, 3),
      intensity = intensity / solid_angle,
      direction = direction
    )
    pbrt_buffer_append(context$lights, light)
  } else if (type %in% c("point", "spot")) {
    intensity = pbrt_spectrum(p, "I", rep(1, 3)) * pbrt_get(p, "scale", 1, 1L)
    from = pbrt_get(p, "from", c(0, 0, 0), 3L)
    from = as.vector(state$transform %*% c(from, 1))[1:3]
    light_args = list(
      color = if (max(intensity) > 0) intensity / max(intensity) else rep(0, 3),
      intensity = max(intensity),
      position = from
    )
    if (type == "spot") {
      to = as.vector(
        state$transform %*% c(pbrt_get(p, "to", c(0, 0, 1), 3L), 1)
      )[1:3]
      angle = pbrt_get(p, "coneangle", 30, 1L)
      delta = pbrt_get(p, "conedeltaangle", 5, 1L)
      light_args$direction = to - from
      light_args$cone_angle = angle
      light_args$falloff_angle = delta
    }
    pbrt_buffer_append(
      context$point_lights,
      do.call(if (type == "point") point_light else spot_light, light_args)
    )
  } else {
    pbrt_note(context, command, paste("Unsupported light omitted:", type))
  }
  pbrt_unused(p)
  invisible(NULL)
}

#' @param rows Scene rows to place.
#' @param state Start/end graphics transforms.
#' @param local Primitive coordinate conversion.
#' @param context Import settings.
#' @return Scene rows with static or animated placements.
#' @keywords internal
#' @noRd
pbrt_place = function(rows, state, local, context) {
  moving = !identical(state$transform, state$end_transform) &&
    !isTRUE(all.equal(state$transform, state$end_transform))
  times = if (is.null(context$settings$TransformTimes)) {
    c(0, 1)
  } else {
    context$settings$TransformTimes$args
  }
  for (i in seq_len(nrow(rows))) {
    rows$transforms[[i]]$group_transform = list(
      if (moving) local else state$transform %*% local
    )
    if (moving) {
      rows$animation_info[[i]]$start_transform_animation = list(state$transform)
      rows$animation_info[[i]]$end_transform_animation = list(
        state$end_transform
      )
      rows$animation_info[[i]]$start_time = times[1]
      rows$animation_info[[i]]$end_time = times[2]
    }
  }
  rows
}

#' @param command Shape directive.
#' @param state Current graphics state.
#' @param context Import resources.
#' @return Scene rows or NULL when a shape is explicitly omitted.
#' @keywords internal
#' @noRd
pbrt_shape = function(command, state, context) {
  p = pbrt_parameters(command, context, state$defaults$shape)
  type = command$args[1]
  material = state$material
  if (!is.null(state$area)) {
    area = pbrt_parameters(state$area, context)
    if (state$area$args[1] != "diffuse") {
      pbrt_note(
        context,
        state$area,
        "Area light converted to diffuse emission."
      )
    }
    radiance = pbrt_spectrum(area, "L", rep(1, 3)) *
      pbrt_get(area, "scale", 1, 1L)
    if (pbrt_get(area, "twosided", FALSE, 1L)) {
      pbrt_note(
        context,
        state$area,
        "Two-sided area emission reduced to front-side emission."
      )
    }
    intensity = max(radiance)
    if (
      !material[[1]]$type %in% c(1L, 4L) ||
        any(material[[1]]$properties[[1]][1:3] != 0)
    ) {
      pbrt_note(
        context,
        state$area,
        "Area emitter's reflecting surface is replaced by an emission-only material."
      )
    }
    material = light(
      color = if (intensity > 0) radiance / intensity else rep(0, 3),
      intensity = intensity
    )
    pbrt_unused(area)
  }
  if (!is.null(p$values$alpha)) {
    alpha = pbrt_texture_parameter(p, "alpha", 1, context)
    if (is.null(alpha$image)) {
      if (length(alpha$value) != 1L || !is.finite(alpha$value)) {
        pbrt_error(
          command,
          "Shape alpha must be a finite scalar or float texture."
        )
      }
      if (alpha$value <= 0) {
        if (alpha$value != 0 || is.null(state$area)) {
          return(NULL)
        }
        # PBRT treats constant zero alpha on an area emitter specially: the
        # surface disappears but still contributes sampled illumination.
        # rayrender's existing invisible light retains that illumination and
        # primary/specular invisibility, with different diffuse-ray visibility.
        material = light(
          color = if (intensity > 0) radiance / intensity else rep(0, 3),
          intensity = intensity,
          invisible = TRUE
        )
        alpha$value = 1
        pbrt_note(
          context,
          command,
          paste(
            "Constant-zero-alpha area emitter converted to an invisible light;",
            "it may still be intersected after a diffuse bounce."
          )
        )
      }
    }
    if (type == "curve") {
      if (!is.null(alpha$image) || alpha$value < 1) {
        pbrt_note(
          context,
          command,
          "Curve alpha masks are not supported; mask omitted."
        )
      }
    } else {
      material = pbrt_alpha_material(material, alpha, context)
    }
  }
  local = diag(4)
  flipped = state$reverse
  if (type == "sphere") {
    radius = pbrt_get(p, "radius", 1, 1L)
    zmin = pbrt_get(p, "zmin", -radius, 1L)
    zmax = pbrt_get(p, "zmax", radius, 1L)
    phimax = pbrt_get(p, "phimax", 360, 1L)
    if (radius <= 0 || zmin >= zmax || phimax <= 0 || phimax > 360) {
      pbrt_error(command, "Invalid sphere dimensions.")
    }
    if (zmin > -radius || zmax < radius || phimax < 360) {
      pbrt_note(
        context,
        command,
        "Clipped sphere omitted; rayrender's analytic sphere must be complete."
      )
      return(NULL)
    }
    rows = sphere(radius = radius, material = material, flipped = flipped)
  } else if (type == "disk") {
    radius = pbrt_get(p, "radius", 1, 1L)
    inner = pbrt_get(p, "innerradius", 0, 1L)
    height = pbrt_get(p, "height", 0, 1L)
    phimax = pbrt_get(p, "phimax", 360, 1L)
    if (radius <= 0 || inner < 0 || inner >= radius) {
      pbrt_error(command, "Invalid disk radii.")
    }
    if (phimax != 360) {
      pbrt_note(context, command, "Partial disk omitted.")
      return(NULL)
    }
    # rayrender's disk and cylinder are y-up; PBRT's primitives are z-up.
    local[1:3, 1:3] = matrix(c(1, 0, 0, 0, 0, 1, 0, -1, 0), 3)
    local[3, 4] = height
    rows = disk(
      radius = radius,
      inner_radius = inner,
      material = material,
      flipped = flipped
    )
  } else if (type == "cylinder") {
    radius = pbrt_get(p, "radius", 1, 1L)
    zmin = pbrt_get(p, "zmin", -1, 1L)
    zmax = pbrt_get(p, "zmax", 1, 1L)
    phimax = pbrt_get(p, "phimax", 360, 1L)
    if (radius <= 0 || zmax <= zmin || phimax <= 0 || phimax > 360) {
      pbrt_error(command, "Invalid cylinder dimensions.")
    }
    local[1:3, 1:3] = matrix(c(1, 0, 0, 0, 0, 1, 0, -1, 0), 3)
    local[3, 4] = (zmin + zmax) / 2
    rows = cylinder(
      radius = radius,
      length = zmax - zmin,
      # The y-up to z-up rotation reverses the azimuth sign. Reflect the
      # clipping interval, rather than reflecting the cylinder's geometry.
      phi_min = 360 - phimax,
      phi_max = 360,
      capped = FALSE,
      material = material,
      flipped = flipped
    )
  } else if (type == "plymesh") {
    file = pbrt_asset(pbrt_get(p, "filename", "", 1L), context, command)
    if (file.info(file)$size == 0) {
      pbrt_error(command, paste("Empty PLY asset:", file))
    }
    if ("displacement" %in% names(p$values)) {
      rows = pbrt_displaced_ply(file, p, state, material, flipped)
    } else {
      rows = ply_model(file, material = material, flipped = flipped)
    }
  } else if (type %in% c("trianglemesh", "bilinearmesh", "loopsubdiv")) {
    rows = pbrt_mesh(p, type, material, flipped)
  } else if (type == "curve") {
    points = pbrt_get(p, "P")
    degree = pbrt_get(p, "degree", 3, 1L)
    basis = pbrt_get(p, "basis", "bezier", 1L)
    if (
      degree != 3 ||
        !basis %in% c("bezier", "bspline") ||
        length(points) %% 3 ||
        length(points) < 12 ||
        (basis == "bezier" && (length(points) / 3 - 1) %% 3)
    ) {
      pbrt_note(
        context,
        command,
        "Only cubic Bezier and uniform B-spline curve chains are supported; curve omitted."
      )
      return(NULL)
    }
    points = matrix(points, ncol = 3, byrow = TRUE)
    width = pbrt_get(p, "width", 1, 1L)
    widths = c(
      pbrt_get(p, "width0", width, 1L),
      pbrt_get(p, "width1", width, 1L)
    )
    split_depth = pbrt_get(p, "splitdepth", 3L, 1L)
    curve_type = pbrt_get(p, "type", "flat", 1L)
    if (!curve_type %in% c("flat", "cylinder", "ribbon")) {
      pbrt_error(command, "Invalid curve type.")
    }
    normals = pbrt_get(p, "N")
    count = if (basis == "bezier") {
      (nrow(points) - 1L) / 3L
    } else {
      nrow(points) - 3L
    }
    if (curve_type == "ribbon" && length(normals) != 3L * (count + 1L)) {
      pbrt_error(
        command,
        "Ribbon curves require one normal at each chain endpoint."
      )
    }
    rows = vector("list", count)
    for (i in seq_len(count)) {
      if (basis == "bezier") {
        control = points[seq.int(3L * (i - 1L) + 1L, length.out = 4L), ]
      } else {
        # Exact uniform cubic B-spline conversion, as in PBRT's
        # CubicBSplineToBezier. Adjacent spans share position and tangent.
        control = points[seq.int(i, length.out = 4L), ]
        control = rbind(
          (control[1, ] + 4 * control[2, ] + control[3, ]) / 6,
          (2 * control[2, ] + control[3, ]) / 3,
          (control[2, ] + 2 * control[3, ]) / 3,
          (control[2, ] + 4 * control[3, ] + control[4, ]) / 6
        )
      }
      args = list(
        p1 = control,
        width = widths[1] + diff(widths) * (i - 1) / count,
        width_end = widths[1] + diff(widths) * i / count,
        type = curve_type,
        split_depth = split_depth,
        material = material,
        flipped = flipped
      )
      if (curve_type == "ribbon") {
        args$normal = normals[3L * (i - 1L) + 1:3]
        args$normal_end = normals[3L * i + 1:3]
      }
      rows[[i]] = do.call(bezier_curve, args)
    }
    # Most PBRT hair directives contain one segment. Its constructor already
    # returned a scene row; inferring and casting its type again is unnecessary.
    rows = if (count == 1L) {
      rows[[1]]
    } else {
      vctrs::list_unchop(rows, ptype = rows[[1]][FALSE, ])
    }
  } else {
    pbrt_note(context, command, paste("Unsupported shape omitted:", type))
    return(NULL)
  }
  pbrt_unused(p)
  if (
    type %in%
      c("sphere", "disk", "cylinder") &&
      (nzchar(material[[1]]$image) ||
        nzchar(material[[1]]$alphaimage) ||
        nzchar(material[[1]]$bump_texture) ||
        (!is.null(material[[1]]$transmission_texture) &&
          nzchar(material[[1]]$transmission_texture)))
  ) {
    pbrt_note(
      context,
      command,
      "Textured primitive uses rayrender's native UV convention, which differs from PBRT."
    )
  }
  rows = pbrt_place(rows, state, local, context)
  if (nzchar(state$inside)) {
    medium = context$media[[state$inside]]$medium
    if (!is.null(medium)) {
      if (
        !type %in%
          c(
            "sphere",
            "plymesh",
            "trianglemesh",
            "bilinearmesh",
            "loopsubdiv"
          ) ||
          state$reverse
      ) {
        pbrt_note(
          context,
          command,
          "Medium omitted: boundaries must be closed, outward-facing spheres or meshes."
        )
      } else {
        definition = context$media[[state$inside]]
        medium$medium_transform = solve(state$transform %*% local) %*%
          definition$transform
        if (
          !isTRUE(all.equal(state$transform, state$end_transform)) &&
            medium$type != "homogeneous"
        ) {
          pbrt_note(
            context,
            command,
            "Spatial medium moves with its animated boundary after matching its definition transform at StartTime."
          )
        }
        rows = set_medium(rows, medium, keep_surface = !state$interface)
      }
    }
  } else if (state$interface) {
    return(NULL)
  }
  rows
}

#' @param material Shape-local material copy.
#' @param alpha Constant/image opacity descriptor with its own UV mapping.
#' @param context Import resources and opacity image cache.
#' @return Material carrying an opacity scalar or a cached linear PNG mask.
#' @keywords internal
#' @noRd
pbrt_alpha_material = function(material, alpha, context) {
  if (is.null(alpha$image)) {
    material[[1]]$alpha_value = min(max(alpha$value, 0), 1)
    return(material)
  }
  # Native alpha textures read the PNG alpha channel as coverage, without a
  # color transfer function. Convert decoded linear values once per texture,
  # not once per shape/instance, and leave UV scaling/translation to lookup.
  filename = context$alpha_images[[alpha$image]]
  if (is.null(filename)) {
    image = rayimage::ray_read_image(alpha$image, convert_to_array = TRUE)
    coverage = pmin(pmax((image[,, 1] + image[,, 2] + image[,, 3]) / 3, 0), 1)
    mask = array(1, c(dim(image)[1:2], 4L))
    mask[,, 4] = coverage
    filename = pbrt_write_asset(mask, ".png", context)
    context$alpha_images[[alpha$image]] = filename
  }
  material[[1]]$alphaimage = filename
  material[[1]]$alpha_repeat = alpha$uv_repeat
  material[[1]]$texture_offsets$alpha = alpha$uv_offset
  material
}

#' @param p Parameter reader.
#' @param name Parameter name.
#' @param default Fallback numeric value.
#' @return Nonnegative numeric scalar/RGB value; unsupported spectra diagnosed.
#' @keywords internal
#' @noRd
pbrt_spectrum = function(p, name, default) {
  entry = p$values[[name]]
  if (is.null(entry)) {
    return(default)
  }
  value = pbrt_get(p, name)
  # Constant wavelength tables describe wavelength-independent coefficients.
  # A neutral RGB coefficient represents them directly, including values > 1.
  if (
    name %in%
      c("sigma_a", "sigma_s") &&
      entry$type == "spectrum" &&
      is.numeric(value) &&
      length(value) >= 4L &&
      length(value) %% 2L == 0L
  ) {
    wavelengths = value[seq.int(1L, length(value), 2L)]
    amplitudes = value[seq.int(2L, length(value), 2L)]
    if (
      all(is.finite(value)) &&
        all(diff(wavelengths) > 0) &&
        wavelengths[1] <= 360 &&
        utils::tail(wavelengths, 1) >= 830 &&
        all(amplitudes == amplitudes[1]) &&
        amplitudes[1] >= 0
    ) {
      return(rep(amplitudes[1], length.out = length(default)))
    }
  }
  if (entry$type %in% c("spectrum", "blackbody")) {
    pbrt_note(
      p$context,
      p$command,
      paste("Spectral parameter replaced by its RGB default:", name)
    )
    return(default)
  }
  if (!is.numeric(value) || !length(value) %in% c(1, 3) || any(value < 0)) {
    pbrt_error(
      p$command,
      paste(name, "must be a nonnegative scalar or RGB value.")
    )
  }
  if (length(default) == 3) {
    value = rep(value, length.out = 3)
  }
  value
}

#' @param p Parameter reader.
#' @param name Parameter name.
#' @param default Constant fallback value.
#' @param context Import resources.
#' @param allow_graph Default `FALSE`. Whether this consumer supports procedural graphs.
#' @return Constant/image descriptor, resolving named textures.
#' @keywords internal
#' @noRd
pbrt_texture_parameter = function(
  p,
  name,
  default,
  context,
  allow_graph = FALSE
) {
  if (!is.null(p$values[[name]]) && p$values[[name]]$type == "texture") {
    reference = pbrt_get(p, name, size = 1L)
    result = context$textures[[reference]]
    if (is.null(result)) {
      pbrt_error(p$command, paste("Undefined texture:", reference))
    }
    if (!is.null(result$graph) && !allow_graph) {
      pbrt_note(
        context,
        p$command,
        paste("Procedural texture is not supported for this input:", name)
      )
      return(list(value = default, uv_repeat = c(1, 1), uv_offset = c(0, 0)))
    }
    return(result)
  }
  if (!is.null(p$values[[name]]) && p$values[[name]]$type == "float") {
    return(list(
      value = pbrt_get(p, name, default, 1L),
      uv_repeat = c(1, 1),
      uv_offset = c(0, 0)
    ))
  }
  list(
    value = pbrt_spectrum(p, name, default),
    uv_repeat = c(1, 1),
    uv_offset = c(0, 0)
  )
}

#' @param descriptor Imported constant, image, or procedural texture.
#' @param type Expected scalar or color output.
#' @return Public texture graph retaining independent UV mapping.
#' @keywords internal
#' @noRd
pbrt_texture_graph = function(descriptor, type) {
  if (!is.null(descriptor$graph)) {
    node = descriptor$graph
  } else if (!is.null(descriptor$image)) {
    node = texture_image(
      descriptor$image,
      type = type,
      encoding = "linear",
      coordinates = texture_coordinates(
        "uv",
        scale = c(descriptor$uv_repeat, 1),
        offset = c(descriptor$uv_offset, 0)
      )
    )
  } else {
    node = texture_constant(descriptor$value)
  }
  if (type == "scalar") texture_scalar(node) else node
}

#' @param filename Environment image file path.
#' @return Linear sRGB radiance array, using EXR primaries when present.
#' @keywords internal
#' @noRd
pbrt_read_environment = function(filename) {
  if (tolower(tools::file_ext(filename)) != "exr") {
    return(pbrt_read_image(filename))
  }
  # Read the header together with the pixels. An untagged PBRT image defaults
  # to sRGB, whereas rayimage's general-purpose EXR fallback is ACEScg.
  source = libopenexr::read_exr(filename)
  image = array(0, c(source$height, source$width, 3L))
  channels = list(source$r, source$g, source$b)
  chromaticities = source$metadata$chromaticities
  conversion = diag(3)
  if (!is.null(chromaticities)) {
    xy = do.call(cbind, chromaticities[c("red", "green", "blue", "white")])
    if (
      !identical(dim(xy), c(2L, 4L)) || any(!is.finite(xy)) || any(xy[2, ] == 0)
    ) {
      stop("Invalid environment EXR chromaticities: ", filename, call. = FALSE)
    }
    xyz = rbind(xy[1, ], xy[2, ], 1 - colSums(xy))
    xyz = sweep(xyz, 2, xy[2, ], "/")
    primaries = xyz[, 1:3, drop = FALSE]
    rgb_to_xyz = sweep(primaries, 2, solve(primaries, xyz[, 4]), "*")
    # PBRT reconstructs spectra in the image's color space, then integrates
    # them into the film's XYZ response. It does not chromatically adapt the
    # image white point to the film white point. Match that linear conversion.
    conversion = rayimage::CS_SRGB$xyz_to_rgb %*% rgb_to_xyz
  }
  for (channel in 1:3) {
    image[,, channel] = pmax(
      0,
      conversion[channel, 1] *
        channels[[1]] +
        conversion[channel, 2] * channels[[2]] +
        conversion[channel, 3] * channels[[3]]
    )
  }
  image
}

#' @param filename Image file path.
#' @param source_linear Default `NULL`. Detect encoding from the format, or override it.
#' @return Linear RGB image array, including floating-point PFM assets.
#' @keywords internal
#' @noRd
pbrt_read_image = function(filename, source_linear = NULL) {
  if (tolower(tools::file_ext(filename)) != "pfm") {
    args = list(filename, convert_to_array = TRUE)
    if (!is.null(source_linear)) {
      args$source_linear = source_linear
    }
    return(do.call(rayimage::ray_read_image, args))
  }
  input = file(filename, "rb")
  on.exit(close(input))
  header = readLines(input, n = 3L, warn = FALSE)
  fields = strsplit(paste(header, collapse = " "), "[[:space:]]+")[[1]]
  if (length(fields) != 4L || !fields[1] %in% c("PF", "Pf")) {
    stop("Invalid PFM header: ", filename, call. = FALSE)
  }
  size_scale = suppressWarnings(as.numeric(fields[2:4]))
  size = size_scale[1:2]
  scale = size_scale[3]
  channels = if (fields[1] == "PF") 3L else 1L
  count = prod(size) * channels
  if (
    any(!is.finite(size_scale)) ||
      any(size < 1 | size != floor(size)) ||
      scale == 0 ||
      !is.finite(count) ||
      count * 4 > file.info(filename)$size
  ) {
    stop("Invalid PFM dimensions or scale: ", filename, call. = FALSE)
  }
  pixels = readBin(
    input,
    numeric(),
    n = count,
    size = 4L,
    endian = if (scale < 0) "little" else "big"
  )
  if (length(pixels) != count || any(!is.finite(pixels))) {
    stop("Truncated or nonfinite PFM pixels: ", filename, call. = FALSE)
  }
  # PFM stores interleaved channels, left to right, with the bottom row first.
  image = aperm(array(pixels * abs(scale), c(channels, size)), c(3, 2, 1))
  image = image[seq.int(size[2], 1L), , , drop = FALSE]
  if (channels == 1L) {
    image = image[,, rep(1L, 3L), drop = FALSE]
  }
  if (identical(source_linear, FALSE)) {
    image = ifelse(
      image <= 0.04045,
      image / 12.92,
      ((image + 0.055) / 1.055)^2.4
    )
  }
  image
}

#' @param image RGB array to write.
#' @param extension File extension.
#' @param context Import resources and generated asset list.
#' @return Absolute path to a newly generated image.
#' @keywords internal
#' @noRd
pbrt_write_asset = function(image, extension, context) {
  if (
    !dir.exists(context$asset_dir) &&
      !dir.create(context$asset_dir, recursive = TRUE)
  ) {
    stop("Unable to create PBRT asset directory.", call. = FALSE)
  }
  filename = tempfile(
    "pbrt-",
    tmpdir = normalizePath(context$asset_dir),
    fileext = extension
  )
  if (extension == ".png") {
    png::writePNG(image, filename)
  } else {
    # rayimage's EXR writer drops singleton channel dimensions before passing
    # them to libopenexr, which requires matrices. Duplicating a constant axis
    # preserves the texture while keeping both spatial dimensions at least two.
    dimensions = dim(image)
    if (any(dimensions[1:2] == 1L)) {
      rows = rep(seq_len(dimensions[1]), length.out = max(2L, dimensions[1]))
      columns = rep(seq_len(dimensions[2]), length.out = max(2L, dimensions[2]))
      image = image[rows, columns, , drop = FALSE]
    }
    # Generated RGB assets contain rayrender's linear sRGB values. Supply
    # their primaries explicitly rather than tagging plain arrays as ACEScg.
    image = rayimage::ray_read_image(
      image,
      source_linear = TRUE,
      assume_colorspace = rayimage::CS_SRGB,
      assume_white = "D65"
    )
    rayimage::ray_write_image(image, filename, write_linear = TRUE)
  }
  context$assets = c(context$assets, filename)
  filename
}

#' @param value Numeric scalar or RGB color.
#' @param p Parameter reader for source diagnostics.
#' @return Three nonnegative linear RGB values bounded by one.
#' @keywords internal
#' @noRd
pbrt_rgb = function(value, p) {
  if (
    !is.numeric(value) ||
      !length(value) %in% c(1, 3) ||
      any(!is.finite(value) | value < 0 | value > 1)
  ) {
    pbrt_error(
      p$command,
      "Reflectance must be a scalar or RGB vector between zero and one."
    )
  }
  rep(value, length.out = 3)
}

#' @param p Parameter reader.
#' @param default Default roughness in PBRT units.
#' @param prefix Default `""`. Layer prefix for coated conductor parameters.
#' @return Two legacy roughness values, optionally carrying native alpha graphs.
#' @keywords internal
#' @noRd
pbrt_roughness = function(p, default, prefix = "") {
  names = paste0(prefix, c("roughness", "uroughness", "vroughness"))
  textured = any(vapply(
    names,
    function(name) {
      entry = p$values[[name]]
      if (is.null(entry) || entry$type != "texture") {
        return(FALSE)
      }
      descriptor = p$context$textures[[entry$value]]
      !is.null(descriptor$graph) || !is.null(descriptor$image)
    },
    logical(1)
  ))
  if (textured) {
    base = pbrt_texture_parameter(
      p,
      names[1],
      default,
      p$context,
      allow_graph = TRUE
    )
    u = if (is.null(p$values[[names[2]]])) {
      base
    } else {
      pbrt_texture_parameter(
        p,
        names[2],
        default,
        p$context,
        allow_graph = TRUE
      )
    }
    v = if (is.null(p$values[[names[3]]])) {
      base
    } else {
      pbrt_texture_parameter(
        p,
        names[3],
        default,
        p$context,
        allow_graph = TRUE
      )
    }
    remap = pbrt_get(p, "remaproughness", TRUE, 1L)
    graphs = lapply(list(u, v), function(x) {
      node = texture_scalar(pbrt_texture_graph(x, "scalar"))
      texture_node(
        "power",
        "scalar",
        child = node,
        exponent = if (remap) 0.5 else 1
      )
    })
    # Nonzero placeholders select a rough BSDF. Native evaluation uses the
    # complete graph; these numbers never replace its spatial variation.
    return(structure(c(.3, .3), alpha_graphs = graphs))
  }
  roughness = pbrt_get(p, paste0(prefix, "roughness"), default, 1L)
  roughness = c(
    pbrt_get(p, paste0(prefix, "uroughness"), roughness, 1L),
    pbrt_get(p, paste0(prefix, "vroughness"), roughness, 1L)
  )
  remap = pbrt_get(p, "remaproughness", TRUE, 1L)
  if (any(roughness < 0)) {
    pbrt_error(p$command, "Roughness cannot be negative.")
  }
  # PBRT v4 RoughnessToAlpha(r) = sqrt(r); rayrender uses alpha = r^2.
  alpha = if (remap) sqrt(roughness) else roughness
  if (any(alpha > 1)) {
    pbrt_note(
      p$context,
      p$command,
      "GGX alpha above one clamped to rayrender's supported range."
    )
    alpha = pmin(alpha, 1)
  }
  sqrt(alpha)
}

#' @param p Disney material parameter reader.
#' @param context Import resources and conversion diagnostics.
#' @param bump_args Prepared bump-map arguments.
#' @return An OpenPBR material approximating PBRT's Disney BSDF.
#' @keywords internal
#' @noRd
pbrt_disney_material = function(p, context, bump_args) {
  # PBRT v3 materials/disney.cpp defines the legacy Disney parameters/defaults:
  # https://github.com/mmp/pbrt-v3/blob/master/src/materials/disney.cpp
  # Keep this separate from pbrt_roughness(): Disney already squares perceptual
  # roughness, unlike PBRT v4's generic RoughnessToAlpha conversion.
  color = pbrt_texture_parameter(p, "color", rep(.5, 3), context)
  base = pbrt_rgb(color$value, p)
  roughness = pbrt_texture_parameter(p, "roughness", .5, context)
  defaults = list(
    metallic = 0,
    eta = 1.5,
    speculartint = 0,
    anisotropic = 0,
    sheen = 0,
    sheentint = .5,
    clearcoat = 0,
    clearcoatgloss = 1,
    spectrans = 0,
    flatness = 0,
    difftrans = 1
  )
  values = lapply(names(defaults), function(name) {
    value = pbrt_get(p, name, defaults[[name]], 1L)
    if (
      !is.finite(value) ||
        value < 0 ||
        (name != "eta" && value > 1) ||
        (name == "eta" && value == 0)
    ) {
      pbrt_error(p$command, paste("Invalid Disney parameter:", name))
    }
    value
  })
  names(values) = names(defaults)
  thin = pbrt_get(p, "thin", FALSE, 1L)
  distance = pbrt_get(p, "scatterdistance", rep(0, 3))
  if (
    !is.numeric(distance) ||
      !length(distance) %in% c(1, 3) ||
      any(!is.finite(distance) | distance < 0)
  ) {
    pbrt_error(
      p$command,
      "Disney scatterdistance must be a nonnegative scalar or RGB distance."
    )
  }
  distance = rep(distance, length.out = 3)
  fitted_surface = !thin &&
    !any(distance > 0) &&
    values$spectrans == 0 &&
    values$eta >= 1
  r = if (is.null(roughness$image)) {
    roughness$value
  } else {
    rayimage::ray_read_image(roughness$image, convert_to_array = TRUE)[,, 1]
  }
  if (
    !is.numeric(r) ||
      any(!is.finite(r) | r < 0 | r > 1) ||
      (is.null(roughness$image) && length(r) != 1L)
  ) {
    pbrt_error(
      p$command,
      "Disney roughness must be a scalar or scalar image in [0, 1]."
    )
  }
  pbrt_note(
    context,
    p$command,
    "Disney material converted to OpenPBR; diffuse/retroreflection, metallic Fresnel, additive sheen/clearcoat, and subsurface transport are approximations."
  )

  # Match the two GGX alpha widths, including Disney's .001 alpha floor.
  # OpenPBR preserves RMS alpha rather than Disney's geometric mean, and uses
  # alpha_y/alpha_x = 1 - anisotropy rather than 1 - .9 * anisotropic.
  aspect = sqrt(1 - .9 * values$anisotropic)
  alpha_x = pmax(r^2 / aspect, .001)
  alpha_y = pmax(r^2 * aspect, .001)
  converted_roughness = ((alpha_x^2 + alpha_y^2) / 2)^.25
  anisotropy = if (is.null(roughness$image)) {
    1 - alpha_y / alpha_x
  } else {
    .9 * values$anisotropic
  }
  if (any(converted_roughness > 1)) {
    pbrt_note(
      context,
      p$command,
      "Anisotropic Disney roughness exceeds OpenPBR's range and is clamped to one."
    )
  }
  converted_roughness = pmin(converted_roughness, 1)

  # Disney normalizes tint by luminance. The fitted surface mapping transfers
  # sheen color gain into its weight; other color clipping is diagnosed.
  luminance = sum(base * c(.212671, .715160, .072169))
  tint = if (luminance > 0) base / luminance else rep(1, 3)
  # PBRT's DisneyFresnel mixes dielectric Fresnel with a metallic Schlick lobe.
  # Its tint affects only partial metalness: at m=0 that Schlick lobe has zero
  # weight, and at m=1 its F0 is the base color, independent of speculartint.
  specular_tint_weight = if (values$metallic > 0 && values$metallic < 1) {
    values$speculartint
  } else {
    0
  }
  specular_tint = (1 - specular_tint_weight) + specular_tint_weight * tint
  sheen_tint = (1 - values$sheentint) + values$sheentint * tint
  if (
    (specular_tint_weight > 0 && any(specular_tint > 1)) ||
      (!fitted_surface && values$sheen > 0 && any(sheen_tint > 1))
  ) {
    pbrt_note(
      context,
      p$command,
      "Luminance-normalized Disney tint exceeds one and is clamped to OpenPBR's color range."
    )
  }
  args = list(
    base_color = base,
    base_metalness = values$metallic,
    base_diffuse_roughness = if (is.null(roughness$image)) r else .5,
    specular_ior = values$eta,
    specular_color = pmin(specular_tint, 1),
    specular_roughness = if (is.null(roughness$image)) {
      converted_roughness
    } else {
      .5
    },
    specular_roughness_anisotropy = anisotropy,
    transmission_weight = values$spectrans,
    # PBRT tints each solid refraction by sqrt(color). Thin OpenPBR represents
    # both interfaces together, so use the complete through-sheet color there.
    transmission_color = if (thin) base else sqrt(base),
    transmission_depth = 0,
    # Disney's clearcoat BRDF has an additional factor of 1/4 relative to a
    # conventional microfacet lobe with the same D, F, and Smith masking.
    coat_weight = .25 * values$clearcoat,
    coat_ior = 1.5,
    # Disney uses GTR1 alpha = lerp(gloss, .1, .001); OpenPBR uses GGX r^2.
    coat_roughness = sqrt(
      .1 * (1 - values$clearcoatgloss) + .001 * values$clearcoatgloss
    ),
    fuzz_weight = values$sheen * (1 - values$metallic) * (1 - values$spectrans),
    fuzz_color = pmin(sheen_tint, 1),
    geometry_thin_walled = thin,
    subsurface_color = base
  )
  if (thin) {
    # At zero phase anisotropy OpenPBR's thin subsurface lobe splits energy
    # equally forward/backward. Weight difftrans reproduces Disney's dt/2 split
    # for the opaque dielectric base without allocating a volume.
    args$subsurface_weight = values$difftrans
    args$subsurface_scatter_anisotropy = 0
    if (values$flatness > 0) {
      pbrt_note(
        context,
        p$command,
        "Disney thin-surface flatness has no matching OpenPBR fake-subsurface lobe and is omitted."
      )
    }
  } else if (any(distance > 0)) {
    volume = pbrt_disney_subsurface(base, distance)
    args$subsurface_weight = 1
    args$subsurface_color = volume$color
    args$subsurface_radius = volume$radius
    args$subsurface_radius_scale = volume$radius_scale
    args$subsurface_scatter_anisotropy = 0
    pbrt_note(
      context,
      p$command,
      "Disney subsurface uses the Hyperion albedo/mean-free-path fit (internal-reflection IOR 1.4); surface eta is retained. Random walks may require a larger max_depth than PBRT diffusion."
    )
  }
  if (fitted_surface) {
    f0 = ((values$eta - 1) / (values$eta + 1))^2
    if (f0 * (1 + values$metallic^3) > .99) {
      pbrt_note(
        context,
        p$command,
        "Fitted Disney metallic Fresnel exceeds the supported range and is clamped to .99."
      )
    }
    if (.112 * args$fuzz_weight * max(1, sheen_tint) > 1) {
      pbrt_note(
        context,
        p$command,
        "Fitted Disney sheen weight exceeds one and is clamped to OpenPBR's range."
      )
    }
    args = pbrt_disney_surface(args, specular_tint, sheen_tint)
  }
  if (!is.null(color$image)) {
    args$image_texture = color$image
    args$image_repeat = color$uv_repeat
    args$image_offset = color$uv_offset
    if (
      specular_tint_weight > 0 ||
        values$sheen > 0 ||
        values$spectrans > 0 ||
        (!thin && any(distance > 0)) ||
        (thin && values$difftrans > 0)
    ) {
      pbrt_note(
        context,
        p$command,
        "Disney base-color image is retained; tint, transmission, and subsurface colors use its constant fallback because OpenPBR exposes no separate maps for those controls."
      )
    }
  }
  if (!is.null(roughness$image)) {
    # The native OpenPBR roughness slot reads raw bytes, so export linear PNG.
    args$roughness_texture = pbrt_write_asset(
      converted_roughness,
      ".png",
      context
    )
    if (is.null(color$image)) {
      args$image_repeat = roughness$uv_repeat
    } else if (any(roughness$uv_repeat != color$uv_repeat)) {
      pbrt_note(
        context,
        p$command,
        "Disney roughness-map UV repeat is replaced by the base-color repeat."
      )
    }
    pbrt_note(
      context,
      p$command,
      "Disney roughness image drives OpenPBR specular roughness; diffuse roughness remains constant and the minimum-alpha anisotropy clamp is approximate."
    )
  }
  material = do.call(openpbr, c(args, bump_args))
  if (!is.null(roughness$image)) {
    material[[1]]$texture_offsets$roughness = roughness$uv_offset
  }
  material
}

#' @param color Linear RGB Disney surface albedo in [0, 1].
#' @param distance Nonnegative RGB Disney scatterdistance, in world units.
#' @return OpenPBR subsurface color, radius, and RGB radius scale.
#' @keywords internal
#' @noRd
pbrt_disney_subsurface = function(color, distance) {
  # Burley et al., The Design and Evolution of Disney's Hyperion Renderer
  # (2018), section 4.4.2: paired fits for isotropic volume albedo and extinction
  # MFP, assuming internal reflections at IOR 1.4. Do not pair this color fit
  # with the different, index-matched 2016 distance fit.
  # https://www.janwalter.org/Publications/HyperionTog2018.pdf
  absorption_fraction = exp(
    color * (-11.43 + color * (15.38 - 13.91 * color))
  )
  distance_scale = 4.012 +
    color * (-15.21 + color * (32.34 + color * (-34.68 + 13.91 * color)))

  # PBRT's DisneyBSSRDF multiplies scatterdistance by .2 before using it as
  # the normalized diffusion profile distance d. Hyperion gives MFP = d * s(A).
  mean_free_path = .2 * distance * distance_scale
  radius = max(mean_free_path)

  # The native OpenPBR reference inverts van de Hulst to get volume albedo.
  # Re-encode the fitted alpha = 1 - absorption_fraction using its forward
  # expression (g=0), rather than passing the Disney surface color unchanged.
  s = sqrt(absorption_fraction)
  encoded_color = (1 - s) * (1 - .139 * s) / (1 + 1.17 * s)
  # Preserve a perfectly white, nonabsorbing input; the empirical polynomial
  # otherwise leaves a small absorption floor at A=1. Black already maps to 0.
  encoded_color[color == 1] = 1
  list(
    color = encoded_color,
    radius = radius,
    radius_scale = if (radius > 0) mean_free_path / radius else rep(0, 3)
  )
}

#' @param args Direct Disney-to-OpenPBR arguments for an opaque solid surface.
#' @param specular_tint Unclipped luminance-normalized Disney specular tint.
#' @param sheen_tint Unclipped luminance-normalized Disney sheen tint.
#' @return OpenPBR arguments with fitted surface lobe controls.
#' @keywords internal
#' @noRd
pbrt_disney_surface = function(args, specular_tint, sheen_tint) {
  # PBRT blends metallic twice: its normal-incidence specular response is
  # (1-m^2)*F0 + m^2*C, but diffuse remains (1-m)*C. With q=1-m+m^2,
  # base_weight=q and metalness=m^2/q preserve both color contributions;
  # F0'=(1+m^3)*F0 preserves the untinted dielectric contribution.
  m = args$base_metalness
  q = 1 - m + m^2
  args$base_weight = q
  args$base_metalness = m^2 / q
  f0 = ((args$specular_ior - 1) / (args$specular_ior + 1))^2
  f0 = min(.99, f0 * (1 + m^3))
  args$specular_ior = (1 + sqrt(f0)) / (1 - sqrt(f0))

  # Rounded fits to the native PBRT/OpenPBR directional response, validated
  # with separate materials, angles, and three render lighting conditions.
  # See dev/experiments/*/2026-09-30-disney-surface-fit for data and reproduction.
  # EON roughness does not reproduce Disney's diffuse/retroreflection terms.
  args$base_diffuse_roughness = 0
  args$specular_color = pmin(
    1,
    pmax(
      0,
      1 + .890 * m * (1 - m) * (specular_tint - 1)
    )
  )
  gain = max(1, sheen_tint)
  args$fuzz_color = sheen_tint / gain
  args$fuzz_weight = min(1, .112 * args$fuzz_weight * gain)
  args$fuzz_roughness = .557

  # GTR1 and GGX cannot match at every angle. Fit their width and amplitude
  # together, then compensate OpenPBR's layered base interface/darkening.
  alpha = args$coat_roughness^2
  args$coat_weight = .242 * args$coat_weight
  args$coat_roughness = .936 * alpha^.43
  args$coat_darkening = 0
  args$specular_ior = args$specular_ior * (1 + .75 * args$coat_weight)
  args
}

#' @param image Linear RGB PBRT v4 equal-area environment image.
#' @param transform Light-to-world affine transform.
#' @return A latitude-longitude RGB array in rayrender's environment convention.
#' @keywords internal
#' @noRd
pbrt_environment_image = function(image, transform) {
  # Invert Clarberg's equal-area sphere map analytically, using atan2 instead
  # of PBRT's SIMD polynomial approximation. See PBRT v4 util/math.cpp,
  # EqualAreaSphereToSquare(). Pixel rows here run top to bottom.
  height = dim(image)[1]
  width = dim(image)[2]
  # Extract channels once; copying an entire HDR plane for every output row
  # makes large environment conversion cubic in image height.
  channels = lapply(1:3, function(channel) image[,, channel])
  output = array(0, c(height, 2L * height, 3L))
  inverse = solve(transform[1:3, 1:3])
  phi = -2 * pi * (seq_len(dim(output)[2]) - .5) / dim(output)[2]
  for (row in seq_len(height)) {
    theta = pi * (row - .5) / height
    world = rbind(
      sin(theta) * sin(phi),
      rep(cos(theta), length(phi)),
      sin(theta) * cos(phi)
    )
    local = inverse %*% world
    local = sweep(local, 2, sqrt(colSums(local^2)), "/")
    radius = sqrt(pmax(0, 1 - abs(local[3, ])))
    angle = atan2(abs(local[2, ]), abs(local[1, ])) * 2 / pi
    v = angle * radius
    u = radius - v
    south = local[3, ] < 0
    old_u = u[south]
    u[south] = 1 - v[south]
    v[south] = 1 - old_u
    u = (.5 + .5 * ifelse(local[1, ] < 0, -u, u))
    v = (.5 + .5 * ifelse(local[2, ] < 0, -v, v))
    x = pmin(width - 1, pmax(0, u * width - .5))
    y = pmin(height - 1, pmax(0, v * height - .5))
    x0 = floor(x)
    x1 = pmin(x0 + 1, width - 1)
    y0 = floor(y)
    y1 = pmin(y0 + 1, height - 1)
    dx = x - x0
    dy = y - y0
    for (channel in 1:3) {
      pixels = channels[[channel]]
      output[row, , channel] = (1 - dy) *
        ((1 - dx) *
          pixels[cbind(y0 + 1, x0 + 1)] +
          dx * pixels[cbind(y0 + 1, x1 + 1)]) +
        dy *
          ((1 - dx) *
            pixels[cbind(y1 + 1, x0 + 1)] +
            dx * pixels[cbind(y1 + 1, x1 + 1)])
    }
  }
  output
}

#' @param filename Resolved PLY path.
#' @param p Shape parameter reader.
#' @param state Graphics transforms.
#' @param material Converted material.
#' @param flipped Whether to reverse orientation.
#' @return Displaced raymesh scene row.
#' @keywords internal
#' @noRd
pbrt_displaced_ply = function(filename, p, state, material, flipped) {
  displacement = pbrt_texture_parameter(p, "displacement", 0, p$context)
  edge = pbrt_get(p, "edgelength", 1, 1L)
  if (edge <= 0) {
    pbrt_error(p$command, "Displacement edgelength must be positive.")
  }
  mesh = rayvertex::ply_mesh(filename)
  for (i in seq_along(mesh$shapes)) {
    if (!nrow(mesh$normals[[i]])) {
      mesh = rayvertex::smooth_normals_mesh(mesh, id = i)
    }
    if (!nrow(mesh$texcoords[[i]])) {
      if (!is.null(displacement$image)) {
        pbrt_note(
          p$context,
          p$command,
          "PLY displacement has no UV coordinates; using the texture value at (0,0)."
        )
      }
      mesh$texcoords[[i]] = matrix(0, nrow(mesh$vertices[[i]]), 2)
      mesh$shapes[[i]]$tex_indices = mesh$shapes[[i]]$indices
      mesh$shapes[[i]]$has_vertex_tex[] = TRUE
    }
  }
  if (any(displacement$uv_repeat != 1)) {
    pbrt_note(
      p$context,
      p$command,
      "Mesh displacement UV repeat requires manual adjustment."
    )
  }
  # Refine the whole mesh with simple subdivision to preserve its silhouette.
  # PBRT refines individual edges adaptively, so topology is an approximation.
  maximum = 0
  for (i in seq_along(mesh$shapes)) {
    vertices = mesh$vertices[[i]]
    world = cbind(vertices, 1) %*% t(state$transform)
    indices = mesh$shapes[[i]]$indices + 1L
    for (pair in list(c(1, 2), c(2, 3), c(3, 1))) {
      maximum = max(
        maximum,
        sqrt(rowSums(
          (world[indices[, pair[1]], 1:3, drop = FALSE] -
            world[indices[, pair[2]], 1:3, drop = FALSE])^2
        ))
      )
    }
  }
  levels = max(0, ceiling(log2(max(maximum / edge, 1))))
  # rayvertex treats level one as disabled; two refinements meet that bound.
  if (levels == 1) {
    levels = 2
  }
  faces = sum(vapply(
    mesh$shapes,
    function(shape) nrow(shape$indices),
    integer(1)
  ))
  if (faces * 4^levels > 2e6) {
    pbrt_error(
      p$command,
      "Displacement refinement would exceed two million triangles; increase edgelength before importing."
    )
  }
  if (levels > 0) {
    mesh = rayvertex::subdivide_mesh(
      mesh,
      subdivision_levels = levels,
      simple = TRUE
    )
  }
  pbrt_note(
    p$context,
    p$command,
    "PLY displacement uses rayvertex UV displacement and uniform simple subdivision in place of PBRT's adaptive edge refinement."
  )
  image = if (is.null(displacement$image)) {
    array(displacement$value[1], c(2, 2, 3))
  } else {
    rayimage::ray_read_image(displacement$image, convert_to_array = TRUE)
  }
  mesh = rayvertex::displace_mesh(mesh, image, verbose = FALSE)
  raymesh_model(
    mesh,
    material = material,
    flipped = flipped,
    importance_sample_lights = material[[1]]$type == 5L
  )
}

#' @param p Parameter reader.
#' @param type PBRT mesh type.
#' @param material Rayrender material.
#' @param flipped Whether to reverse surface orientation.
#' @return One raymesh scene row with zero-based topology.
#' @keywords internal
#' @noRd
pbrt_mesh = function(p, type, material, flipped) {
  points = pbrt_get(p, "P")
  if (!is.numeric(points) || !length(points) || length(points) %% 3) {
    pbrt_error(p$command, "Mesh P must contain xyz triples.")
  }
  vertices = matrix(points, ncol = 3, byrow = TRUE)
  arity = if (type == "bilinearmesh") 4L else 3L
  indices = pbrt_get(
    p,
    "indices",
    if (nrow(vertices) == arity) seq_len(arity) - 1L else NULL
  )
  if (
    !is.numeric(indices) ||
      !length(indices) ||
      length(indices) %% arity ||
      any(indices != floor(indices) | indices < 0 | indices >= nrow(vertices))
  ) {
    pbrt_error(p$command, "Invalid zero-based mesh indices.")
  }
  # Validated indices are bounded by an R matrix dimension, so they fit in an
  # integer. This also avoids repeated native coercions of each index channel.
  indices = matrix(as.integer(indices), ncol = arity, byrow = TRUE)
  face_indices = pbrt_get(p, "faceIndices")
  if (!is.null(face_indices)) {
    if (
      !is.numeric(face_indices) ||
        length(face_indices) != nrow(indices) ||
        any(
          !is.finite(face_indices) |
            face_indices < 0 |
            face_indices > .Machine$integer.max |
            face_indices != floor(face_indices)
        )
    ) {
      pbrt_error(
        p$command,
        "faceIndices must contain one nonnegative integer per source face."
      )
    }
    face_indices = as.integer(face_indices)
  }
  if (type == "bilinearmesh") {
    # PBRT patch order is p00,p10,p01,p11, not a perimeter quad order.
    a = vertices[indices[, 2] + 1, , drop = FALSE] -
      vertices[indices[, 1] + 1, , drop = FALSE]
    b = vertices[indices[, 3] + 1, , drop = FALSE] -
      vertices[indices[, 1] + 1, , drop = FALSE]
    d = vertices[indices[, 4] + 1, , drop = FALSE] -
      vertices[indices[, 1] + 1, , drop = FALSE]
    normal = cbind(
      a[, 2] * b[, 3] - a[, 3] * b[, 2],
      a[, 3] * b[, 1] - a[, 1] * b[, 3],
      a[, 1] * b[, 2] - a[, 2] * b[, 1]
    )
    if (
      any(
        abs(rowSums(normal * d)) >
          1e-8 * pmax(sqrt(rowSums(normal^2) * rowSums(d^2)), 1e-30)
      )
    ) {
      pbrt_note(
        p$context,
        p$command,
        "Nonplanar bilinear patches approximated by two triangles each."
      )
    }
    indices = rbind(
      indices[, c(1, 2, 4), drop = FALSE],
      indices[, c(1, 4, 3), drop = FALSE]
    )
    if (!is.null(face_indices)) face_indices = rep(face_indices, 2L)
  }
  normals = pbrt_get(p, "N")
  if (!is.null(normals)) {
    if (!is.numeric(normals) || length(normals) != length(points)) {
      pbrt_error(p$command, "Mesh N must contain one normal per vertex.")
    }
    normals = matrix(normals, ncol = 3, byrow = TRUE)
    if (any(rowSums(normals^2) == 0)) {
      pbrt_error(p$command, "Mesh normals must be nonzero.")
    }
  }
  uv = pbrt_get(p, "uv")
  if (is.null(uv)) {
    uv = pbrt_get(p, "st")
  }
  if (!is.null(uv)) {
    if (!is.numeric(uv) || length(uv) != 2L * nrow(vertices)) {
      pbrt_error(p$command, "Mesh UV must contain one pair per vertex.")
    }
    uv = matrix(uv, ncol = 2, byrow = TRUE)
  }
  levels = 1
  if (type == "loopsubdiv") {
    pbrt_note(
      p$context,
      p$command,
      "Loop subdivision limit surface approximated with rayrender's mesh subdivision."
    )
    levels = pbrt_get(p, "levels", 3, 1L)
    if (levels < 0 || levels != floor(levels)) {
      pbrt_error(p$command, "Subdivision levels must be a nonnegative integer.")
    }
  }
  mesh = pbrt_uniform_mesh(
    vertices = vertices,
    indices = indices,
    normals = normals,
    uv = uv,
    context = p$context
  )
  if (!is.null(face_indices)) {
    mesh$shapes[[1]]$ptex_face_indices = face_indices
  }
  if (type == "bilinearmesh" && is.null(uv)) {
    # Two triangles retain one default patch chart, rather than each receiving
    # the default triangle chart. UV indices need not equal vertex indices.
    patches = nrow(indices) %/% 2L
    mesh$texcoords[[1]] = rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
    mesh$shapes[[1]]$tex_indices = rbind(
      matrix(rep(c(0L, 1L, 3L), patches), ncol = 3L, byrow = TRUE),
      matrix(rep(c(0L, 3L, 2L), patches), ncol = 3L, byrow = TRUE)
    )
    mesh$shapes[[1]]$has_vertex_tex = rep(TRUE, nrow(indices))
  }
  if (!is.null(face_indices) && type == "loopsubdiv" && levels > 0) {
    pbrt_error(
      p$command,
      "Subdivision of Ptex-addressed meshes is not supported."
    )
  }
  if (type == "loopsubdiv" && levels <= 1) {
    # rayrender's lazy subdivision uses level one to mean 'disabled'. Use the
    # mesh API here so an explicitly requested single refinement is retained.
    mesh = rayvertex::subdivide_mesh(mesh, subdivision_levels = levels)
    levels = 1
  }
  raymesh_model(
    mesh,
    material = material,
    flipped = flipped,
    subdivision_levels = levels,
    # Coordinates, indices, normals and UVs were validated above; construct_mesh
    # supplies the remaining structure. Revalidating every tiny PBRT patch is
    # a substantial part of import time in scenes with many separate meshes.
    validate_mesh = FALSE,
    importance_sample_lights = material[[1]]$type == 5L
  )
}

#' @param vertices Validated PBRT vertex matrix.
#' @param indices Validated, zero-based triangle indices.
#' @param normals One normal per vertex, or NULL.
#' @param uv One UV pair per vertex, or NULL.
#' @param context Import resources and shared mesh template.
#' @return A rayvertex mesh with uniform attribute coverage.
#' @keywords internal
#' @noRd
pbrt_uniform_mesh = function(vertices, indices, normals, uv, context) {
  if (is.null(context$mesh_template)) {
    # Derive containers and material metadata from the public constructor.
    # PBRT's validation above guarantees uniform attribute coverage, so avoid
    # rechecking every face and rehashing the same material for every patch.
    context$mesh_template = rayvertex::construct_mesh(
      vertices = matrix(numeric(), 0L, 3L),
      indices = matrix(integer(), 0L, 3L),
      material = context$mesh_material
    )
  }
  mesh = context$mesh_template
  mesh$vertices[[1]] = vertices
  if (!is.null(normals)) {
    mesh$normals[[1]] = normals
  }
  if (!is.null(uv)) {
    mesh$texcoords[[1]] = uv
  }
  faces = nrow(indices)
  missing_indices = if (is.null(normals) || is.null(uv)) {
    matrix(-1, faces, 3L)
  } else {
    NULL
  }
  shape = mesh$shapes[[1]]
  shape$indices = indices
  shape$norm_indices = if (is.null(normals)) missing_indices else indices
  shape$tex_indices = if (is.null(uv)) missing_indices else indices
  shape$material_ids = rep(0, faces)
  shape$has_vertex_tex = rep(!is.null(uv), faces)
  shape$has_vertex_normals = rep(!is.null(normals), faces)
  mesh$shapes[[1]] = shape
  mesh
}
