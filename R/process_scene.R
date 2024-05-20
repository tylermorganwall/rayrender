#'@title Process a scene
#'
#'@keywords internal
#'#internal
process_scene = function(scene, process_material_ids = TRUE) {
  shapevec = unlist(lapply(tolower(scene$shape),switch,
                           "sphere" = 1,"xy_rect" = 2, "xz_rect" = 3,"yz_rect" = 4,"box" = 5, 
                           "obj" = 6, "disk" = 7, "cylinder" = 8, "ellipsoid" = 9,
                           "curve" = 10, "csg_object" = 11, "ply" = 12,
                           "mesh3d" = 13, "raymesh" = 14, "instance" = 15, stop(sprintf("Shape not recognized"))))
  scene$shape = shapevec
  # Convert shapes to enums
  typevec = rep(0,nrow(scene))
  
  for(i in seq_len(nrow(scene))) {
    typevec[i] = switch(scene$material[[i]]$type,
                        "diffuse" = 1L,"metal" = 2L,"dielectric" = 3L,
                        "oren-nayar" = 4L, "light" = 5L, "mf" = 6L,
                        "glossy" = 7L, "spotlight" = 8L, "hair" = 9L, "mf-t" = 10L,
                        stop(sprintf("Material type `%s` not found",scene$material[[i]]$type)))
    scene$material[[i]]$type = typevec[i]
  }
  
  #alpha texture handler -- need to do this before images to override alpha if present
  alpha_temp_file_names = tempfile(sprintf("alphatemp%i",seq_len(nrow(scene))),fileext = ".png")
  for(i in seq_len(nrow(scene))) {
    alpha_input = scene$material[[i]]$alphaimage[[1]]
    alpha_tex_bool = is.array(alpha_input)
    alpha_is_filename = is.character(alpha_input) && !is.na(alpha_input)
    if(alpha_tex_bool) {
      if(length(dim(alpha_input)) == 2) {
        png::writePNG(fliplr(t(alpha_input)), alpha_temp_file_names[i])
      } else if(dim(alpha_input)[3] == 4) {
        alpha_input[,,1] = alpha_input[,,4]
        alpha_input[,,2] = alpha_input[,,4]
        alpha_input[,,3] = alpha_input[,,4]
        png::writePNG(fliplr(aperm(alpha_input[,,1:3],c(2,1,3))), alpha_temp_file_names[i])
      } else if(dim(alpha_input)[3] == 3) {
        png::writePNG(fliplr(aperm(alpha_input,c(2,1,3))), alpha_temp_file_names[i])
      } else {
        stop("alpha texture dims: c(", paste(dim(alpha_input),collapse=", "), ") not valid for texture.")
      }
      scene$material[[i]]$alphaimage = alpha_temp_file_names[i]
    } else if(alpha_is_filename) {
      if(any(!file.exists(path.expand(alpha_input)) & nchar(alpha_input) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(alpha_input, collapse="\n")))
      }
      temp_array = png::readPNG(alpha_input) 
      if(dim(temp_array)[3] == 4 && any(temp_array[,,4] != 1)) {
        temp_array[,,1] = temp_array[,,4]
        temp_array[,,2] = temp_array[,,4]
        temp_array[,,3] = temp_array[,,4]
      }
      scene$material[[i]]$alphaimage = alpha_temp_file_names[i]
      png::writePNG(temp_array,alpha_temp_file_names[i])
    } else {
      scene$material[[i]]$alphaimage = ""
    }
  }
  
  #texture handler
  temp_file_names = tempfile(sprintf("imagetemp%i",seq_len(nrow(scene))), fileext = ".png")
  for(i in seq_len(nrow(scene))) {
    image_input = scene$material[[i]]$image[[1]]
    image_tex_bool = is.array(image_input)
    image_is_filename = is.character(image_input) && !is.na(image_input)
    if(image_tex_bool) {
      if(dim(image_input)[3] == 4) {
        png::writePNG(fliplr(aperm(image_input[,,1:3],c(2,1,3))),temp_file_names[i])
        #Handle PNG with alpha
        if(!alpha_tex_bool && any(image_input[,,4] != 1)) {
          image_input[,,1] = image_input[,,4]
          image_input[,,2] = image_input[,,4]
          image_input[,,3] = image_input[,,4]
          png::writePNG(fliplr(aperm(image_input[,,1:3],c(2,1,3))), alpha_temp_file_names[i])
          scene$material[[i]]$alphaimage = alpha_temp_file_names[i]
        }
      } else if(dim(image_input)[3] == 3){
        png::writePNG(fliplr(aperm(image_input,c(2,1,3))),temp_file_names[i])
      }
      scene$material[[i]]$image = temp_file_names[i]
    } else if(image_is_filename) {
      if(any(!file.exists(path.expand(image_input)) & nchar(image_input) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(image_input, collapse="\n")))
      }
      temp_file_names[i] = path.expand(image_input)
      scene$material[[i]]$image = temp_file_names[i]
    } else {
      scene$material[[i]]$image = ""
    }
  }
  
  #displacement texture handler
  disp_temp_file_names = tempfile(sprintf("disp_image_temp%i",seq_len(nrow(scene))), fileext = ".png")
  for(i in seq_len(nrow(scene))) {
    if(!scene$shape[[i]] %in% c(6,13,14)) { #obj, mesh3d, raymesh
      next
    }
    image_input = scene$shape_info[[i]]$shape_properties$displacement_texture[[1]]
    image_tex_bool = is.array(image_input)
    image_is_filename = is.character(image_input) && !is.na(image_input)
    if(image_tex_bool) {
      if(dim(image_input)[3] == 4) {
        png::writePNG(fliplr(aperm(image_input[,,1:3],c(2,1,3))),disp_temp_file_names[i])
      } else if(dim(image_input)[3] == 3){
        png::writePNG(fliplr(aperm(image_input,c(2,1,3))),disp_temp_file_names[i])
      }
      scene$shape_info[[i]]$shape_properties$displacement_texture = disp_temp_file_names[i]
    } else if(image_is_filename) {
      if(any(!file.exists(path.expand(image_input)) & nchar(image_input) > 0)) {
        stop(paste0("Cannot find the following displacement texture file:\n",
                    paste(image_input, collapse="\n")))
      }
      disp_temp_file_names[i] = path.expand(image_input)
      scene$shape_info[[i]]$shape_properties$displacement_texture = disp_temp_file_names[i]
    } else {
      scene$shape_info[[i]]$shape_properties$displacement_texture = ""
    }
  }

  #bump texture handler
  bump_temp_file_names = tempfile(sprintf("bumptemp%i",seq_len(nrow(scene))),fileext = ".png")
  for(i in seq_len(nrow(scene))) {
    bump_input = scene$material[[i]]$bump_texture[[1]]
    bump_tex_bool = is.array(bump_input)
    bump_is_filename = is.character(bump_input)  && !is.na(bump_input)
    if(bump_tex_bool) {
      bump_dims = dim(bump_input)
      if(length(bump_dims) == 2) {
        temp_array = array(0, dim = c(bump_dims,3))
        temp_array[,,1] = bump_input
        temp_array[,,2] = bump_input
        temp_array[,,3] = bump_input
        bump_dims = c(bump_dims,3)
      } else {
        temp_array = bump_input
      }
      if(bump_dims[3] == 4) {
        png::writePNG(fliplr(aperm(temp_array[,,1:3],c(2,1,3))),bump_temp_file_names[i])
      } else if(bump_dims[3] == 3){
        png::writePNG(fliplr(aperm(temp_array,c(2,1,3))),bump_temp_file_names[i])
      }
      scene$material[[i]]$bump_texture = bump_temp_file_names[i]
    } else if(bump_is_filename) {
      if(any(!file.exists(path.expand(bump_input)) & nchar(bump_input) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(bump_input, collapse="\n")))
      }
      bump_temp_file_names[i] = path.expand(bump_input)
      scene$material[[i]]$bump_texture = bump_temp_file_names[i]
    } else {
      scene$material[[i]]$bump_texture = ""
    }
  }
  
  #roughness texture handler
  rough_temp_file_names = tempfile(sprintf("roughtemp%i",seq_len(nrow(scene))),fileext = ".png")
  for(i in seq_len(nrow(scene))) {
    roughness_input = scene$material[[i]]$roughness_texture[[1]]
    rough_tex_bool = is.array(roughness_input)
    roughness_is_filename = is.character(roughness_input)  && !is.na(roughness_input)
    if(rough_tex_bool) {
      if(length(dim(roughness_input)) == 2) {
        png::writePNG(fliplr(t(roughness_input)), rough_temp_file_names[i])
      } else if(dim(roughness_input)[3] == 3) {
        png::writePNG(fliplr(aperm(roughness_input,c(2,1,3))), rough_temp_file_names[i])
      } else {
        stop("alpha texture dims: c(", paste(dim(roughness_input),collapse=", "), ") not valid for texture.")
      }
      scene$material[[i]]$roughness_texture = rough_temp_file_names[i]
    } else if(roughness_is_filename) {
      if(any(!file.exists(path.expand(roughness_input)) & nchar(roughness_input) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(roughness_input, collapse="\n")))
      }
      rough_temp_file_names[i] = path.expand(roughness_input)
      scene$material[[i]]$roughness_texture = rough_temp_file_names[i]
    } else{
      scene$material[[i]]$roughness_texture = ""
    }
  }
  
  for(i in seq_len(nrow(scene))) {
    fileinfovec = scene$shape_info[[i]]$fileinfo
    if(!is.na(fileinfovec)) {
      if(any(!file.exists(scene$shape_info[[i]]$fileinfo) & nchar(scene$shape_info[[i]]$fileinfo) > 0)) {
        stop(paste0("Cannot find the following obj/ply file:\n",
                    fileinfovec, collapse="\n"))
      }
    } 
  }

  #Spotlight handler
  if(any(typevec == 8)) {
    if(any(shapevec[typevec == 8] > 4)) {
      stop("spotlights are only supported for spheres and rects")
    }
    for(i in seq_len(nrow(scene))) {
      if(typevec[i] == 8) {
        scene$material[[i]]$properties[[1]][4:6] = scene$material[[i]]$properties[[1]][4:6] - 
          c(scene$x[i],
            scene$y[i],
            scene$z[i]) 
      }
    }
  }
  
  # Only process material IDs right before rendering
  if(process_material_ids) {
    #Material ID handler; these must show up in increasing order.  Note, this will
    #cause problems if `match` is ever changed to return doubles when matching in
    #long vectors as has happened with `which` recently.
    material_id = unlist(lapply(scene$shape_info, \(x) x$material_id))
    is_na_mat = is.na(material_id)
    material_id_increasing = as.integer(match(material_id, unique(material_id)) - 1L)
    for(i in seq_len(nrow(scene))) {
      if(!is_na_mat[i]) {
        scene$shape_info[[i]]$material_id = material_id_increasing[i]
      } 
    }
  }
  
  #Detect any importance sampling
  any_light = FALSE
  for(i in seq_len(nrow(scene))) {
    any_light = any_light || (scene$material[[i]]$type %in% c(5, 7)) #light and spotlight
    if(scene$shape[[i]] == 15) { #instance
      any_light = any_light || scene$shape_info[[i]]$shape_properties$any_light 
    }
  }
  
  
  scene_info = list()
  scene_info$scene = scene
  scene_info$shape = shapevec
  scene_info$typevec = typevec
  scene_info$any_light = any_light
  
  return(scene_info)
}
