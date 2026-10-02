#'@title Process a scene
#'
#' @param scene Scene to process.
#' @param process_material_ids Default `TRUE`. Renumber shared material IDs.
#' @param subsurface_prepared Default `FALSE`. Whether this call's scene already has resolved subsurface interiors.
#' @return Processed scene and native shape/material metadata.
#'@keywords internal
process_scene = function(
  scene,
  process_material_ids = TRUE,
  subsurface_prepared = FALSE
) {
  if (!subsurface_prepared) {
    scene = prepare_subsurface(scene)
  }
  # Keep list-column updates local through all texture and material passes.
  materials = vctrs::vec_data(scene$material)
  shape_info = vctrs::vec_data(scene$shape_info)
  shapevec = unlist(lapply(
    tolower(scene$shape),
    switch,
    "sphere" = 1,
    "xy_rect" = 2,
    "xz_rect" = 3,
    "yz_rect" = 4,
    "box" = 5,
    "obj" = 6,
    "disk" = 7,
    "cylinder" = 8,
    "ellipsoid" = 9,
    "curve" = 10,
    "csg_object" = 11,
    "ply" = 12,
    "mesh3d" = 13,
    "raymesh" = 14,
    "instance" = 15,
    stop(sprintf("Shape not recognized"))
  ))
  scene$shape = shapevec
  # Convert shapes to enums
  typevec = rep(0, nrow(scene))

  for (i in seq_len(nrow(scene))) {
    typevec[i] = materials[[i]]$type
  }

  #alpha texture handler -- need to do this before images to override alpha if present
  # Allocate paths only for images we actually write. Millions of untextured
  # primitives otherwise create five unused filename strings apiece.
  for (i in seq_len(nrow(scene))) {
    alpha_input = materials[[i]]$alphaimage
    alpha_tex_bool = is.array(alpha_input)
    alpha_is_filename = is.character(alpha_input) && nchar(alpha_input) > 0
    if (alpha_tex_bool || alpha_is_filename) {
      alpha_file = tempfile("alphatemp", fileext = ".png")
    }
    if (alpha_tex_bool) {
      if (length(dim(alpha_input)) == 2) {
        png::writePNG(fliplr(t(alpha_input)), alpha_file)
      } else if (dim(alpha_input)[3] == 4) {
        alpha_input[,, 1] = alpha_input[,, 4]
        alpha_input[,, 2] = alpha_input[,, 4]
        alpha_input[,, 3] = alpha_input[,, 4]
        png::writePNG(
          fliplr(aperm(alpha_input[,, 1:3], c(2, 1, 3))),
          alpha_file
        )
      } else if (dim(alpha_input)[3] == 3) {
        png::writePNG(
          fliplr(aperm(alpha_input, c(2, 1, 3))),
          alpha_file
        )
      } else {
        stop(
          "alpha texture dims: c(",
          paste(dim(alpha_input), collapse = ", "),
          ") not valid for texture."
        )
      }
      materials[[i]]$alphaimage = alpha_file
    } else if (alpha_is_filename) {
      if (any(!file.exists(path.expand(alpha_input)))) {
        stop(paste0(
          "Cannot find the following texture file:\n",
          paste(alpha_input, collapse = "\n")
        ))
      }
      temp_array = png::readPNG(alpha_input)
      if (dim(temp_array)[3] == 4 && any(temp_array[,, 4] != 1)) {
        temp_array[,, 1] = temp_array[,, 4]
        temp_array[,, 2] = temp_array[,, 4]
        temp_array[,, 3] = temp_array[,, 4]
      }
      if (
        dim(temp_array)[3] == 4 &&
          all(temp_array[,, 4] == 1) &&
          any(temp_array[,, 1:3] != 1)
      ) {
        temp_array[,, 4] = temp_array[,, 1]
      }
      materials[[i]]$alphaimage = alpha_file
      png::writePNG(temp_array, alpha_file)
    }
  }

  #texture handler
  for (i in seq_len(nrow(scene))) {
    image_input = materials[[i]]$image
    image_tex_bool = is.array(image_input)
    image_is_filename = is.character(image_input) && nchar(image_input) > 0
    if (image_tex_bool) {
      image_file = tempfile("imagetemp", fileext = ".png")
      if (dim(image_input)[3] == 4) {
        png::writePNG(
          fliplr(aperm(image_input[,, 1:4], c(2, 1, 3))),
          image_file
        )
        #Handle PNG with alpha
        if (!alpha_tex_bool && any(image_input[,, 4] != 1)) {
          materials[[i]]$alphaimage = image_file
        }
      } else if (dim(image_input)[3] == 3) {
        png::writePNG(
          fliplr(aperm(image_input, c(2, 1, 3))),
          image_file
        )
      }
      materials[[i]]$image = image_file
    } else if (image_is_filename) {
      image_input = path.expand(image_input)
      if (!file.exists(image_input)) {
        stop(paste0("Cannot find the following texture file:\n", image_input))
      }
      # Automatic alpha extraction is a PNG operation. JPEG and floating-point
      # HDR/EXR textures are loaded directly by the renderer's image decoder.
      tmp_image = if (tolower(tools::file_ext(image_input)) == "png") {
        png::readPNG(image_input)
      } else {
        NULL
      }
      if (length(dim(tmp_image)) == 3 && dim(tmp_image)[3] == 4) {
        if (any(tmp_image[,, 4] != 1)) {
          tmp_image[,, 1] = tmp_image[,, 4]
          tmp_image[,, 2] = tmp_image[,, 4]
          tmp_image[,, 3] = tmp_image[,, 4]
          tmp_image[,, 4] = tmp_image[,, 4]
          alpha_file = materials[[i]]$alphaimage
          if (!nzchar(alpha_file)) {
            alpha_file = tempfile("alphatemp", fileext = ".png")
          }
          png::writePNG(tmp_image[,, 1:4], alpha_file)
          materials[[i]]$alphaimage = alpha_file
        }
      }
      materials[[i]]$image = image_input
    }
  }

  #displacement texture handler
  for (i in seq_len(nrow(scene))) {
    if (!scene$shape[[i]] %in% c(6, 13, 14)) {
      #obj, mesh3d, raymesh
      next
    }
    if (scene$shape[[i]] == 13) {
      image_input = shape_info[[i]]$mesh_info[[1]]$displacement_texture
    } else {
      image_input = shape_info[[i]]$shape_properties$displacement_texture
    }
    image_tex_bool = is.array(image_input)
    image_is_filename = is.character(image_input) && nchar(image_input) > 0
    if (image_tex_bool) {
      displacement_file = tempfile("disp_image_temp", fileext = ".png")
      if (dim(image_input)[3] == 4) {
        png::writePNG(
          fliplr(aperm(image_input[,, 1:3], c(2, 1, 3))),
          displacement_file
        )
      } else if (dim(image_input)[3] == 3) {
        png::writePNG(
          fliplr(aperm(image_input, c(2, 1, 3))),
          displacement_file
        )
      }
      shape_info[[
        i
      ]]$shape_properties$displacement_texture = displacement_file
    } else if (image_is_filename) {
      if (any(!file.exists(path.expand(image_input)))) {
        stop(paste0(
          "Cannot find the following displacement texture file:\n",
          paste(image_input, collapse = "\n")
        ))
      }
      shape_info[[
        i
      ]]$shape_properties$displacement_texture = path.expand(image_input)
    }
  }

  #bump texture handler
  for (i in seq_len(nrow(scene))) {
    bump_input = materials[[i]]$bump_texture
    bump_tex_bool = is.array(bump_input)
    bump_is_filename = is.character(bump_input) && nchar(bump_input) > 0
    if (bump_tex_bool) {
      bump_file = tempfile("bumptemp", fileext = ".png")
      bump_dims = dim(bump_input)
      if (length(bump_dims) == 2) {
        temp_array = array(0, dim = c(bump_dims, 3))
        temp_array[,, 1] = bump_input
        temp_array[,, 2] = bump_input
        temp_array[,, 3] = bump_input
        bump_dims = c(bump_dims, 3)
      } else {
        temp_array = bump_input
      }
      if (bump_dims[3] == 4) {
        png::writePNG(
          fliplr(aperm(temp_array[,, 1:3], c(2, 1, 3))),
          bump_file
        )
      } else if (bump_dims[3] == 3) {
        png::writePNG(
          fliplr(aperm(temp_array, c(2, 1, 3))),
          bump_file
        )
      }
      materials[[i]]$bump_texture = bump_file
    } else if (bump_is_filename) {
      if (any(!file.exists(path.expand(bump_input)))) {
        stop(paste0(
          "Cannot find the following texture file:\n",
          paste(bump_input, collapse = "\n")
        ))
      }
      materials[[i]]$bump_texture = path.expand(bump_input)
    }
  }

  #roughness texture handler
  for (i in seq_len(nrow(scene))) {
    roughness_input = materials[[i]]$roughness_texture
    rough_tex_bool = is.array(roughness_input)
    roughness_is_filename = is.character(roughness_input) &&
      nchar(roughness_input) > 0
    if (rough_tex_bool) {
      roughness_file = tempfile("roughtemp", fileext = ".png")
      if (length(dim(roughness_input)) == 2) {
        png::writePNG(fliplr(t(roughness_input)), roughness_file)
      } else if (dim(roughness_input)[3] == 3) {
        png::writePNG(
          fliplr(aperm(roughness_input, c(2, 1, 3))),
          roughness_file
        )
      } else {
        stop(
          "alpha texture dims: c(",
          paste(dim(roughness_input), collapse = ", "),
          ") not valid for texture."
        )
      }
      materials[[i]]$roughness_texture = roughness_file
    } else if (roughness_is_filename) {
      if (any(!file.exists(path.expand(roughness_input)))) {
        stop(paste0(
          "Cannot find the following texture file:\n",
          paste(roughness_input, collapse = "\n")
        ))
      }
      materials[[i]]$roughness_texture = path.expand(roughness_input)
    }
  }

  for (i in seq_len(nrow(scene))) {
    fileinfovec = shape_info[[i]]$fileinfo
    if (!is.na(fileinfovec)) {
      if (
        any(
          !file.exists(shape_info[[i]]$fileinfo) &
            nchar(shape_info[[i]]$fileinfo) > 0
        )
      ) {
        stop(paste0(
          "Cannot find the following obj/ply file:\n",
          fileinfovec,
          collapse = "\n"
        ))
      }
    }
  }

  #Spotlight handler
  if (any(typevec == 8)) {
    if (any(shapevec[typevec == 8] > 4)) {
      stop("spotlights are only supported for spheres and rects")
    }
    for (i in seq_len(nrow(scene))) {
      if (typevec[i] == 8) {
        materials[[i]]$properties[[1]][4:6] = materials[[
          i
        ]]$properties[[1]][4:6] -
          c(scene$x[i], scene$y[i], scene$z[i])
      }
    }
  }

  # Only process material IDs right before rendering
  if (process_material_ids) {
    #Material ID handler; these must show up in increasing order.  Note, this will
    #cause problems if `match` is ever changed to return doubles when matching in
    #long vectors as has happened with `which` recently.
    material_id = unlist(lapply(shape_info, \(x) x$material_id))
    is_na_mat = is.na(material_id)
    material_id_increasing = as.integer(
      match(material_id, unique(material_id)) - 1L
    )
    for (i in seq_len(nrow(scene))) {
      if (!is_na_mat[i]) {
        shape_info[[i]]$material_id = material_id_increasing[i]
      }
    }
  }

  #Detect any importance sampling
  any_light = FALSE
  for (i in seq_len(nrow(scene))) {
    any_light = any_light ||
      (materials[[i]]$type %in% c(5, 8)) ||
      (identical(materials[[i]]$type, 11L) &&
        materials[[i]]$openpbr$emission_luminance > 0)
    if (scene$shape[[i]] == 15) {
      #instance
      any_light = any_light || shape_info[[i]]$shape_properties$any_light
    }
  }

  scene$material = vctrs::vec_restore(materials, scene$material)
  scene$shape_info = vctrs::vec_restore(shape_info, scene$shape_info)

  scene_info = list()
  scene_info$scene = scene
  scene_info$shape = shapevec
  scene_info$typevec = typevec
  scene_info$any_light = any_light

  return(scene_info)
}
