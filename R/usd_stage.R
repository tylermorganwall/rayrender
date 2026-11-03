#' @param scene Scene with asset paths already remapped.
#' @param assets Complete asset dependency paths.
#' @param root Stage directory, used to read model geometry.
#' @return Nested text fragments and deferred numeric data for native streaming.
#' @keywords internal
#' @noRd
usd_stage_document = function(scene, assets, root) {
  writer = new.env(parent = emptyenv())
  writer$root = root
  writer$fps = if (is.null(attr(scene, "ray_render_settings")$fps)) {
    30
  } else {
    attr(scene, "ray_render_settings")$fps
  }
  writer$last_time = 0
  writer$warnings = character()
  writer$materials = list()
  writer$material_paths = new.env(hash = TRUE, parent = emptyenv())
  objects = lapply(seq_len(nrow(scene)), function(i) {
    row = lapply(scene, function(column) {
      if (is.list(column)) column[[i]] else column[i]
    })
    usd_object_text(
      row,
      sprintf("Object%08d", i - 1L),
      paste0("/World/Objects/", sprintf("Object%08d", i - 1L)),
      writer
    )
  })
  cameras = usd_camera_text(attributes(scene), writer)
  lights = usd_light_text(attributes(scene), writer)
  metadata = list(
    "customData = { dictionary rayrender = {\nint version = 1\ndictionary attributes = ",
    usd_metadata(attributes(scene), TRUE),
    "\ndictionary columns = ",
    usd_metadata(as.list(scene), TRUE, templates = TRUE),
    "\n} }"
  )
  world = usd_prim(
    "Xform",
    "World",
    c(
      usd_property("asset[]", "rayrender:assets", assets, custom = TRUE),
      usd_prim("Scope", "Objects", objects),
      usd_prim("Scope", "Materials", writer$materials),
      cameras,
      usd_prim("Scope", "Lights", lights)
    ),
    metadata
  )
  if (length(writer$warnings)) {
    warning(
      paste(sort(unique(writer$warnings)), collapse = "\n"),
      call. = FALSE
    )
  }
  list(
    '#usda 1.0\n(\ndefaultPrim = "World"\nupAxis = "Y"\nmetersPerUnit = 1\nstartTimeCode = 0\nendTimeCode = ',
    as.character(writer$last_time),
    "\nframesPerSecond = ",
    as.character(writer$fps),
    "\ntimeCodesPerSecond = ",
    as.character(writer$fps),
    "\n)\n",
    world
  )
}

#' @param row One object's descriptors.
#' @param name Prim name.
#' @param path Absolute prim path.
#' @param writer Shared authoring state.
#' @param archive Default `TRUE`. Attach the exact renderer record.
#' @return An object Xform and its geometry/material bindings.
#' @keywords internal
#' @noRd
usd_object_text = function(row, name, path, writer, archive = TRUE) {
  metadata = if (archive) {
    list(
      "customData = { dictionary rayrender = { dictionary record = ",
      usd_metadata(row, TRUE),
      " } }"
    )
  } else {
    ""
  }
  transform = usd_object_matrix(row)
  animation = row$animation_info
  transforms = usd_property("matrix4d", "xformOp:transform", transform)
  start = animation$start_transform_animation[[1]]
  if (is.matrix(start) && identical(dim(start), c(4L, 4L))) {
    end = animation$end_transform_animation[[1]]
    times = c(animation$start_time, animation$end_time) * writer$fps
    transforms = c(
      transforms,
      paste0(
        "matrix4d xformOp:transform.timeSamples = {\n",
        paste(
          c(
            paste0(
              times[1],
              ": ",
              usd_literal_native(transform %*% t(start), "matrix4d")
            ),
            paste0(
              times[2],
              ": ",
              usd_literal_native(transform %*% t(end), "matrix4d")
            )
          ),
          collapse = ",\n"
        ),
        "\n}"
      )
    )
    writer$last_time = max(writer$last_time, times)
  }
  if (!is.null(animation$track)) {
    writer$warnings = c(
      writer$warnings,
      "Authored animation tracks are retained in rayrender metadata; interchange samples only endpoint transforms."
    )
  }
  material = usd_material_text(row$material, writer)
  geometry = usd_geometry_text(row, paste0(path, "/Geometry"), material, writer)
  volume = NULL
  medium = row$shape_info$medium
  if (is.list(medium) && identical(medium$type, "nanovdb")) {
    is_vdb = tolower(tools::file_ext(medium$filename)) == "vdb"
    field = usd_prim(
      if (is_vdb) "OpenVDBAsset" else "RayrenderNanoVDBAsset",
      "Density",
      c(
        usd_property("asset", "filePath", medium$filename, custom = !is_vdb),
        usd_property(
          "token",
          "fieldName",
          medium$density_grid,
          custom = !is_vdb
        )
      )
    )
    volume = usd_prim(
      "Volume",
      "Volume",
      c(paste0("rel field:density = <", path, "/Volume/Density>"), field)
    )
    if (!is_vdb) {
      writer$warnings = c(
        writer$warnings,
        "NanoVDB assets use the rayrender volume extension; other USD renderers may require OpenVDB conversion."
      )
    }
  }
  usd_prim(
    "Xform",
    name,
    c(
      transforms,
      usd_property(
        "token[]",
        "xformOpOrder",
        "xformOp:transform",
        uniform = TRUE
      ),
      geometry,
      volume
    ),
    metadata
  )
}

#' @param attrs Scene attributes.
#' @param writer Shared authoring state.
#' @return Camera and render-settings prims.
#' @keywords internal
#' @noRd
usd_camera_text = function(attrs, writer) {
  cameras = attrs$ray_cameras
  settings = attrs$ray_render_settings
  width = if (is.null(settings$width)) 400 else settings$width
  height = if (is.null(settings$height)) 400 else settings$height
  paths = character()
  result = lapply(seq_along(cameras), function(i) {
    cam = cameras[[i]]
    name = paste0("Camera", i - 1L)
    paths[names(cameras)[i]] <<- paste0("/World/Cameras/", name)
    n = length(cam$motion$x)
    properties = list()
    for (frame in seq_len(n)) {
      sample = lapply(cam$motion, `[`, frame)
      from = unlist(sample[c("x", "y", "z")], use.names = FALSE)
      at = unlist(sample[c("dx", "dy", "dz")], use.names = FALSE)
      up = unlist(sample[c("upx", "upy", "upz")], use.names = FALSE)
      # USD cameras look down local -Z; matrix rows are right, up and backward.
      backward = from - at
      backward = backward / sqrt(sum(backward^2))
      right = usd_cross(up, backward)
      right = right / sqrt(sum(right^2))
      transform = diag(4)
      transform[1:3, 1:3] = rbind(right, usd_cross(backward, right), backward)
      transform[4, 1:3] = from
      fov = sample$fov
      vertical = if (fov == 0) sample$orthoy * 10 else 100 * tan(fov * pi / 360)
      horizontal = if (fov == 0) {
        sample$orthox * 10
      } else {
        vertical * width / height
      }
      values = list(
        "matrix4d xformOp:transform" = usd_literal_native(
          transform,
          "matrix4d"
        ),
        "token projection" = usd_literal_native(
          if (fov == 0) "orthographic" else "perspective",
          "token"
        ),
        "float focalLength" = "50",
        "float focusDistance" = usd_literal_native(sample$focal, "float"),
        "float fStop" = usd_literal_native(
          if (sample$aperture == 0) 0 else 5 / sample$aperture,
          "float"
        )
      )
      if (fov >= 0 && fov < 180) {
        values[["float verticalAperture"]] = usd_literal_native(
          vertical,
          "float"
        )
        values[["float horizontalAperture"]] = usd_literal_native(
          horizontal,
          "float"
        )
      } else {
        writer$warnings = c(
          writer$warnings,
          "Panoramic/physical camera projection uses rayrender camera metadata."
        )
      }
      for (key in names(values)) {
        properties[[key]] = c(properties[[key]], values[[key]])
      }
    }
    body = lapply(names(properties), function(key) {
      values = properties[[key]]
      if (n == 1) {
        return(list(key, " = ", values))
      }
      list(
        key,
        ".timeSamples = {",
        lapply(seq_along(values), function(i) {
          list(if (i > 1L) ", ", as.character(i - 1L), ": ", values[[i]])
        }),
        "}"
      )
    })
    writer$last_time = max(writer$last_time, n - 1)
    usd_prim(
      "Camera",
      name,
      c(
        body,
        usd_property(
          "token[]",
          "xformOpOrder",
          "xformOp:transform",
          uniform = TRUE
        )
      ),
      paste0(
        "customData = { dictionary rayrender = { string name = ",
        usd_literal_native(names(cameras)[i], "string"),
        " } }"
      )
    )
  })
  active = attrs$active_camera
  render = usd_prim(
    "RenderSettings",
    "RenderSettings",
    c(
      usd_property("int2", "resolution", c(width, height)),
      if (length(active) && active %in% names(paths)) {
        paste0("rel camera = <", paths[active], ">")
      }
    )
  )
  c(usd_prim("Scope", "Cameras", result), render)
}

#' @param attrs Scene attributes.
#' @param writer Shared authoring state.
#' @return Light prims.
#' @keywords internal
#' @noRd
usd_light_text = function(attrs, writer) {
  result = list()
  for (i in seq_along(attrs$ray_infinite_lights)) {
    light = attrs$ray_infinite_lights[[i]]
    type = "Xform"
    body = character()
    if (light$type == "image") {
      type = "DomeLight"
      body = c(
        usd_property("asset", "inputs:texture:file", light$filename),
        usd_property(
          "float",
          "inputs:intensity",
          if (is.null(light$intensity)) 1 else light$intensity
        ),
        usd_property(
          "float",
          "xformOp:rotateY",
          if (is.null(light$rotation)) 0 else light$rotation
        ),
        usd_property(
          "token[]",
          "xformOpOrder",
          "xformOp:rotateY",
          uniform = TRUE
        )
      )
    } else if (
      light$type %in% c("disk", "uniform_disk") && !is.null(light$direction)
    ) {
      type = "DistantLight"
      body = c(
        usd_property(
          "float",
          "inputs:angle",
          if (is.null(light$angular_diameter)) 0.53 else light$angular_diameter
        ),
        usd_property(
          "color3f",
          "inputs:color",
          if (is.null(light$color)) c(1, 1, 1) else light$color
        ),
        usd_property(
          "float",
          "inputs:intensity",
          if (is.null(light$intensity)) 1 else light$intensity
        ),
        usd_property(
          "matrix4d",
          "xformOp:transform",
          usd_direction_matrix(-light$direction)
        ),
        usd_property(
          "token[]",
          "xformOpOrder",
          "xformOp:transform",
          uniform = TRUE
        )
      )
    } else {
      writer$warnings = c(
        writer$warnings,
        "Procedural/celestial lighting uses rayrender metadata; sky datasets remain external."
      )
    }
    result[[length(result) + 1L]] = usd_prim(
      type,
      paste0("Light", i - 1L),
      body
    )
  }
  for (i in seq_along(attrs$ray_lights)) {
    light = attrs$ray_lights[[i]]
    transform = if (light$type == "spot") {
      usd_direction_matrix(light$direction)
    } else {
      diag(4)
    }
    transform[4, 1:3] = light$position
    body = c(
      usd_property("bool", "inputs:treatAsPoint", TRUE),
      usd_property("float", "inputs:radius", 0),
      usd_property("color3f", "inputs:color", light$color),
      usd_property("float", "inputs:intensity", light$intensity),
      usd_property("matrix4d", "xformOp:transform", transform),
      usd_property(
        "token[]",
        "xformOpOrder",
        "xformOp:transform",
        uniform = TRUE
      )
    )
    if (light$type == "spot") {
      body = c(
        body,
        usd_property("float", "inputs:shaping:cone:angle", light$cone_angle),
        usd_property(
          "float",
          "inputs:shaping:cone:softness",
          min(1, light$falloff_angle / max(light$cone_angle, 1e-9))
        )
      )
    }
    result[[length(result) + 1L]] = usd_prim(
      "SphereLight",
      paste0("Point", i - 1L),
      body,
      if (light$type == "spot") 'prepend apiSchemas = ["ShapingAPI"]' else ""
    )
  }
  result
}

#' @param type USD prim type.
#' @param name Prim name.
#' @param body Property declarations and child prims.
#' @param metadata Default `""`. Prim metadata declarations.
#' @return USDA prim definition.
#' @keywords internal
#' @noRd
usd_prim = function(type, name, body, metadata = "") {
  # Retain nested fragments: neither a prim nor its parent becomes one R string.
  list(list(
    "def ",
    type,
    " ",
    usd_literal_native(name, "string"),
    if (length(metadata) && !identical(metadata, "")) {
      list(" (\n", metadata, "\n)")
    },
    " {\n",
    lapply(body, function(fragment) list(fragment, "\n")),
    "\n}\n"
  ))
}

#' @param type USD value type.
#' @param name Property name.
#' @param value Property value.
#' @param custom Default `FALSE`. Mark a renderer-specific property.
#' @param uniform Default `FALSE`. Mark a uniform property.
#' @return USDA property declaration.
#' @keywords internal
#' @noRd
usd_property = function(type, name, value, custom = FALSE, uniform = FALSE) {
  list(list(
    paste0(
      if (custom) "custom ",
      if (uniform) "uniform ",
      type,
      " ",
      name,
      " = "
    ),
    structure(
      list(value = value, type = type),
      class = "rayrender_usd_literal"
    )
  ))
}

#' @param value Renderer metadata, retained by reference until writing.
#' @param fields Default `FALSE`. Encode named fields as a USD dictionary.
#' @param templates Default `FALSE`. Retain column attributes without their data.
#' @return A deferred metadata fragment for the native stream writer.
#' @keywords internal
#' @noRd
usd_metadata = function(value, fields = FALSE, templates = FALSE) {
  # The outer list keeps c() from flattening the typed node into ordinary fields.
  list(structure(
    list(value = value, fields = fields, templates = templates),
    class = "rayrender_usd_metadata"
  ))
}

#' @param row Object row.
#' @return USD row-vector transformation matrix.
#' @keywords internal
#' @noRd
usd_object_matrix = function(row) {
  transform = row$transforms
  result = diag(c(transform$scale[[1]], 1)) %*%
    usd_rotation_matrix(transform$angle[[1]], transform$order_rotation[[1]])
  translation = diag(4)
  translation[4, 1:3] = unlist(row[c("x", "y", "z")], use.names = FALSE)
  result = result %*% translation
  group = transform$group_transform[[1]]
  if (is.matrix(group) && identical(dim(group), c(4L, 4L))) {
    result = result %*% t(group)
  }
  result
}

#' @param descriptor Material descriptor.
#' @param writer Shared authoring state.
#' @return Absolute path of a deduplicated preview material.
#' @keywords internal
#' @noRd
usd_material_text = function(descriptor, writer) {
  key = digest::digest(descriptor, algo = "sha256")
  if (exists(key, writer$material_paths, inherits = FALSE)) {
    return(writer$material_paths[[key]])
  }
  name = paste0("m", key)
  path = paste0("/World/Materials/", name)
  props = descriptor$properties[[1]]
  color = if (length(props) >= 3L) props[1:3] else c(1, 1, 1)
  preview = c(
    usd_property("token", "info:id", "UsdPreviewSurface", uniform = TRUE),
    usd_property("color3f", "inputs:diffuseColor", color),
    usd_property("float", "inputs:roughness", 0.5),
    "token outputs:surface"
  )
  kind = descriptor$type
  if (kind %in% c("metal", 2)) {
    preview = c(preview, usd_property("float", "inputs:metallic", 1))
  }
  if (kind %in% c("dielectric", 3, "mf-t", 10)) {
    writer$warnings = c(
      writer$warnings,
      "Dielectric transmission uses rayrender metadata; USD Preview Surface is only an approximation."
    )
  }
  if (kind %in% c("light", 5)) {
    preview = c(preview, usd_property("color3f", "inputs:emissiveColor", color))
  }
  texture = NULL
  image = descriptor$image
  if (
    is.character(image) && length(image) == 1L && !is.na(image) && nzchar(image)
  ) {
    preview = c(
      preview,
      paste0(
        "color3f inputs:diffuseColor.connect = <",
        path,
        "/Texture.outputs:rgb>"
      )
    )
    uv = usd_prim(
      "Shader",
      "UV",
      c(
        usd_property(
          "token",
          "info:id",
          "UsdPrimvarReader_float2",
          uniform = TRUE
        ),
        usd_property("token", "inputs:varname", "st"),
        "float2 outputs:result"
      )
    )
    texture = c(
      uv,
      usd_prim(
        "Shader",
        "Texture",
        c(
          usd_property("token", "info:id", "UsdUVTexture", uniform = TRUE),
          usd_property("asset", "inputs:file", image),
          paste0("float2 inputs:st.connect = <", path, "/UV.outputs:result>"),
          "float3 outputs:rgb"
        )
      )
    )
  } else if (is.array(image)) {
    writer$warnings = c(
      writer$warnings,
      "Inline image arrays are retained in typed rayrender data; file-backed textures provide the USD preview image."
    )
  }
  if (
    any(!vapply(descriptor$texture_graphs, is.null, TRUE)) ||
      isTRUE(descriptor$noise != 0)
  ) {
    writer$warnings = c(
      writer$warnings,
      "Procedural materials retain their exact rayrender descriptors; interchange uses preview materials."
    )
  }
  writer$materials[[length(writer$materials) + 1L]] = usd_prim(
    "Material",
    name,
    c(
      paste0(
        "token outputs:surface.connect = <",
        path,
        "/Preview.outputs:surface>"
      ),
      usd_prim("Shader", "Preview", preview),
      texture
    ),
    list(
      "customData = { dictionary rayrender = { dictionary material = ",
      usd_metadata(descriptor),
      " } }"
    )
  )
  writer$material_paths[[key]] = path
  path
}

#' @param row Object descriptor.
#' @param path Geometry prim path.
#' @param material Material path.
#' @param writer Shared authoring state.
#' @return Geometry prim text.
#' @keywords internal
#' @noRd
usd_geometry_text = function(row, path, material, writer) {
  shapes = c(
    "sphere",
    "xy_rect",
    "xz_rect",
    "yz_rect",
    "box",
    "obj",
    "disk",
    "cylinder",
    "ellipsoid",
    "curve",
    "csg_object",
    "ply",
    "mesh3d",
    "raymesh",
    "instance"
  )
  kind = if (is.numeric(row$shape)) shapes[row$shape] else row$shape
  info = row$shape_info
  props = info$shape_properties
  type = "Xform"
  body = c(paste0("rel material:binding = <", material, ">"))
  if (isTRUE(info$flipped)) {
    body = c(
      body,
      usd_property("token", "orientation", "leftHanded", uniform = TRUE)
    )
  }
  if (kind == "sphere") {
    type = "Sphere"
    body = c(body, usd_property("double", "radius", props$radius))
  } else if (kind == "box") {
    type = "Cube"
    body = c(
      body,
      usd_property("double", "size", 1),
      usd_property("float3", "xformOp:scale", props$boxinfo),
      usd_property("token[]", "xformOpOrder", "xformOp:scale", uniform = TRUE)
    )
  } else if (
    kind == "cylinder" &&
      (if (is.null(props$phi_max)) (2 * pi) else props$phi_max) >= 2 * pi
  ) {
    type = "Cylinder"
    body = c(
      body,
      usd_property("double", "radius", props$radius),
      usd_property("double", "height", props$length),
      usd_property("token", "axis", "Y", uniform = TRUE)
    )
  } else if (kind %in% c("xy_rect", "xz_rect", "yz_rect")) {
    a = props$rectinfo[1] / 2
    b = props$rectinfo[2] / 2
    points = rbind(c(-a, -b, 0), c(a, -b, 0), c(a, b, 0), c(-a, b, 0))
    if (kind == "xz_rect") {
      points = points[, c(1, 3, 2)]
    }
    if (kind == "yz_rect") {
      points = points[, c(3, 1, 2)]
    }
    return(usd_mesh_text("Geometry", points, matrix(0:3, 1), body = body))
  } else if (kind == "mesh3d") {
    mesh = info$mesh_info[[1]]
    return(usd_mesh_text(
      "Geometry",
      mesh$vertices,
      mesh$indices,
      mesh$normals,
      mesh$texcoords,
      body = body
    ))
  } else if (kind == "raymesh") {
    mesh = info$mesh_info[[1]]
    for (i in seq_along(mesh$shapes)) {
      body = c(
        body,
        usd_mesh_text(
          paste0("Mesh", i - 1L),
          mesh$vertices[[i]],
          mesh$shapes[[i]]$indices,
          mesh$normals[[i]],
          mesh$texcoords[[i]],
          mesh$shapes[[i]]$norm_indices,
          mesh$shapes[[i]]$tex_indices
        )
      )
    }
    writer$warnings = c(
      writer$warnings,
      "Imported per-face raymesh materials are retained for rayrender; interchange uses the object preview material."
    )
  } else if (kind == "obj") {
    source = info$fileinfo
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", source)) {
      source = file.path(writer$root, source)
    }
    mesh = usd_obj_mesh(
      source,
      if (is.null(props$scale_obj)) 1 else props$scale_obj
    )
    writer$warnings = c(
      writer$warnings,
      "OBJ material libraries remain available for rayrender; interchange uses the object preview material."
    )
    return(usd_mesh_text(
      "Geometry",
      mesh$points,
      mesh$faces,
      mesh$normals,
      mesh$uv,
      mesh$normal_indices,
      mesh$uv_indices,
      body
    ))
  } else if (kind == "instance") {
    type = "PointInstancer"
    original = props$original_scene[[1]]
    prototype_path = paste0(path, "/Prototypes/Prototype")
    children = lapply(seq_len(nrow(original)), function(i) {
      record = lapply(original, function(column) {
        if (is.list(column)) column[[i]] else column[i]
      })
      name = sprintf("Object%08d", i - 1L)
      usd_object_text(
        record,
        name,
        paste0(prototype_path, "/", name),
        writer,
        archive = FALSE
      )
    })
    n = length(props$x_values)
    angles = cbind(props$angle_x, props$angle_y, props$angle_z)
    rotations = usd_instance_quaternions_native(
      angles,
      as.integer(row$transforms$order_rotation[[1]])
    )
    body = c(
      body,
      usd_prim("Scope", "Prototypes", usd_prim("Xform", "Prototype", children)),
      paste0("rel prototypes = <", prototype_path, ">"),
      usd_property("int[]", "protoIndices", rep.int(0L, n)),
      usd_property(
        "point3f[]",
        "positions",
        cbind(props$x_values, props$y_values, props$z_values)
      ),
      usd_property(
        "float3[]",
        "scales",
        cbind(props$scale_x, props$scale_y, props$scale_z)
      ),
      usd_property("quath[]", "orientations", rotations)
    )
  } else {
    if (
      is.character(info$fileinfo) &&
        !is.na(info$fileinfo) &&
        nzchar(info$fileinfo)
    ) {
      body = c(
        body,
        usd_property(
          "asset",
          "rayrender:sourceAsset",
          info$fileinfo,
          custom = TRUE
        )
      )
    }
    writer$warnings = c(
      writer$warnings,
      paste0(
        kind,
        " geometry has a rayrender-only representation in this export."
      )
    )
  }
  usd_prim(
    type,
    "Geometry",
    body,
    'prepend apiSchemas = ["MaterialBindingAPI"]'
  )
}

#' @param a First vector.
#' @param b Second vector.
#' @return Their three-dimensional cross product.
#' @keywords internal
#' @noRd
usd_cross = function(a, b) {
  c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
}

#' @param direction Direction corresponding to local negative Z.
#' @return Row-vector rotation matrix.
#' @keywords internal
#' @noRd
usd_direction_matrix = function(direction) {
  backward = -direction / sqrt(sum(direction^2))
  up = if (abs(backward[2]) > 0.99) c(1, 0, 0) else c(0, 1, 0)
  right = usd_cross(up, backward)
  right = right / sqrt(sum(right^2))
  result = diag(4)
  result[1:3, 1:3] = rbind(right, usd_cross(backward, right), backward)
  result
}

#' @param angles Euler angles in degrees.
#' @param order Rotation axis order.
#' @return Row-vector rotation matrix.
#' @keywords internal
#' @noRd
usd_rotation_matrix = function(angles, order) {
  result = diag(4)
  for (axis in order) {
    theta = angles[axis] * pi / 180
    cosine = cos(theta)
    sine = sin(theta)
    rotation = diag(4)
    pair = switch(axis, c(2L, 3L), c(3L, 1L), c(1L, 2L))
    rotation[pair, pair] = matrix(
      c(cosine, -sine, sine, cosine),
      2,
      byrow = TRUE
    )
    result = result %*% t(rotation)
  }
  result
}

#' @param name Prim name.
#' @param points Vertex matrix.
#' @param faces Zero-based face indices, as a matrix or list.
#' @param normals Default `NULL`. Vertex normals.
#' @param uv Default `NULL`. Texture coordinates.
#' @param normal_indices Default `NULL`. Face-varying normal indices.
#' @param uv_indices Default `NULL`. Face-varying UV indices.
#' @param body Default `character()`. Additional properties.
#' @return Mesh prim text.
#' @keywords internal
#' @noRd
usd_mesh_text = function(
  name,
  points,
  faces,
  normals = NULL,
  uv = NULL,
  normal_indices = NULL,
  uv_indices = NULL,
  body = character()
) {
  counts = if (is.matrix(faces)) {
    rep.int(ncol(faces), nrow(faces))
  } else {
    lengths(faces)
  }
  indices = if (is.matrix(faces)) {
    as.vector(t(faces))
  } else {
    unlist(faces, use.names = FALSE)
  }
  body = c(
    body,
    usd_property("point3f[]", "points", points[, 1:3, drop = FALSE]),
    usd_property("int[]", "faceVertexCounts", counts),
    usd_property("int[]", "faceVertexIndices", indices),
    usd_property("token", "subdivisionScheme", "none", uniform = TRUE)
  )
  if (is.matrix(normals) && length(normals)) {
    varying = length(normal_indices) > 0
    if (varying) {
      ni = if (is.matrix(normal_indices)) {
        as.vector(t(normal_indices))
      } else {
        unlist(normal_indices, use.names = FALSE)
      }
      normals = normals[ni + 1L, , drop = FALSE]
    }
    body = c(
      body,
      list(list(
        usd_property("normal3f[]", "normals", normals[, 1:3, drop = FALSE]),
        ' (interpolation = "',
        if (varying) "faceVarying" else "vertex",
        '")'
      ))
    )
  }
  if (is.matrix(uv) && length(uv)) {
    body = c(
      body,
      list(list(
        usd_property("texCoord2f[]", "primvars:st", uv[, 1:2, drop = FALSE]),
        ' (interpolation = "',
        if (length(uv_indices)) "faceVarying" else "vertex",
        '")'
      ))
    )
    if (length(uv_indices)) {
      body = c(
        body,
        usd_property(
          "int[]",
          "primvars:st:indices",
          if (is.matrix(uv_indices)) {
            as.vector(t(uv_indices))
          } else {
            unlist(uv_indices, use.names = FALSE)
          }
        )
      )
    }
  }
  usd_prim("Mesh", name, body, 'prepend apiSchemas = ["MaterialBindingAPI"]')
}

#' @param filename OBJ file.
#' @param scale Vertex scale.
#' @return Mesh vertices, faces, normals, UVs and their indices.
#' @keywords internal
#' @noRd
usd_obj_mesh = function(filename, scale) {
  points = normals = uv = faces = ni = ui = list()
  for (line in readLines(filename, warn = FALSE)) {
    words = strsplit(trimws(sub("#.*$", "", line)), "[[:space:]]+")[[1]]
    if (!length(words) || !nzchar(words[1])) {
      next
    }
    kind = words[1]
    if (kind == "v") {
      points[[length(points) + 1L]] = as.numeric(words[2:4]) * scale
    }
    if (kind == "vn") {
      normals[[length(normals) + 1L]] = as.numeric(words[2:4])
    }
    if (kind == "vt") {
      uv[[length(uv) + 1L]] = as.numeric(words[2:3])
    }
    if (kind == "f") {
      indices = strsplit(words[-1], "/", fixed = TRUE)
      resolved = lapply(seq_len(3), function(column) {
        values = vapply(
          indices,
          function(index) {
            if (length(index) >= column && nzchar(index[column])) {
              as.integer(index[column])
            } else {
              NA_integer_
            }
          },
          1L
        )
        size = c(length(points), length(uv), length(normals))[column]
        ifelse(values > 0, values - 1L, size + values)
      })
      if (
        anyNA(resolved[[1]]) ||
          any(resolved[[1]] < 0 | resolved[[1]] >= length(points))
      ) {
        stop("Invalid OBJ face in ", filename)
      }
      faces[[length(faces) + 1L]] = resolved[[1]]
      ui[[length(ui) + 1L]] = resolved[[2]]
      ni[[length(ni) + 1L]] = resolved[[3]]
    }
  }
  if (!length(points) || !length(faces)) {
    stop("OBJ contains no mesh: ", filename)
  }
  has_normals = length(normals) && !anyNA(unlist(ni))
  has_uv = length(uv) && !anyNA(unlist(ui))
  list(
    points = do.call(rbind, points),
    faces = faces,
    normals = if (has_normals) do.call(rbind, normals),
    normal_indices = if (has_normals) ni,
    uv = if (has_uv) do.call(rbind, uv),
    uv_indices = if (has_uv) ui
  )
}
