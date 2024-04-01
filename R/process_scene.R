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
  alpha_array_list = vector(mode="list", length(nrow(scene)))
  alpha_tex_bool = vector(mode = "logical", length(nrow(scene)))
  alpha_filename_bool = vector(mode = "logical", length(nrow(scene)))
  alpha_temp_file_names = tempfile(sprintf("alphatemp%i",seq_len(nrow(scene))),fileext = ".png")
  for(i in seq_len(nrow(scene))) {
    alpha_array_list[i] = scene$material[[i]]$alphaimage
    alpha_tex_bool[i] = is.array(alpha_array_list[[i]])
    alpha_filename_bool[i] = is.character(alpha_array_list[[i]]) && !is.na(alpha_array_list[[i]])
    if(alpha_tex_bool[i]) {
      if(length(dim(alpha_array_list[[i]][[1]])) == 2) {
        png::writePNG(fliplr(t(alpha_array_list[[i]][[1]])), alpha_temp_file_names[i])
      } else if(dim(alpha_array_list[[i]][[1]])[3] == 4) {
        alpha_array_list[[i]][[1]][,,1] = alpha_array_list[[i]][[1]][,,4]
        alpha_array_list[[i]][[1]][,,2] = alpha_array_list[[i]][[1]][,,4]
        alpha_array_list[[i]][[1]][,,3] = alpha_array_list[[i]][[1]][,,4]
        png::writePNG(fliplr(aperm(alpha_array_list[[i]][[1]][,,1:3],c(2,1,3))), alpha_temp_file_names[i])
      } else if(dim(alpha_array_list[[i]][[1]])[3] == 3) {
        png::writePNG(fliplr(aperm(alpha_array_list[[i]][[1]],c(2,1,3))), alpha_temp_file_names[i])
      } else {
        stop("alpha texture dims: c(", paste(dim(alpha_array_list[[i]][[1]]),collapse=", "), ") not valid for texture.")
      }
      scene$material[[i]]$alphaimage = TRUE
    } else if(alpha_filename_bool[i]) {
      if(any(!file.exists(path.expand(alpha_array_list[[i]][[1]])) & nchar(alpha_array_list[[i]][[1]]) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(alpha_array_list[[i]][[1]], collapse="\n")))
      }
      temp_array = png::readPNG(alpha_array_list[[i]][[1]]) 
      if(dim(temp_array)[3] == 4 && any(temp_array[,,4] != 1)) {
        temp_array[,,1] = temp_array[,,4]
        temp_array[,,2] = temp_array[,,4]
        temp_array[,,3] = temp_array[,,4]
      }
      scene$material[[i]]$alphaimage = TRUE
      png::writePNG(temp_array,alpha_temp_file_names[i])
    } else {
      scene$material[[i]]$alphaimage = FALSE
    }
  }
  alpha_tex_bool = alpha_tex_bool | alpha_filename_bool
  alpha_list = list()
  alpha_list$alpha_temp_file_names = alpha_temp_file_names
  alpha_list$alpha_tex_bool = alpha_tex_bool
  
  #texture handler
  image_array_list = vector(mode="list", length(nrow(scene))) #scene$image
  image_tex_bool = vector(mode = "logical", length(nrow(scene)))
  image_filename_bool = vector(mode = "logical", length(nrow(scene)))
  temp_file_names = tempfile(sprintf("imagetemp%i",seq_len(nrow(scene))),fileext = ".png")
  for(i in seq(nrow(scene))) {
    image_array_list[i] = scene$material[[i]]$image
    image_tex_bool[i] = is.array(image_array_list[[i]])
    image_filename_bool[i] = is.character(image_array_list[[i]]) && !is.na(image_array_list[[i]])
    if(image_tex_bool[i]) {
      if(dim(image_array_list[[i]])[3] == 4) {
        png::writePNG(fliplr(aperm(image_array_list[[i]][,,1:3],c(2,1,3))),temp_file_names[i])
        #Handle PNG with alpha
        if(!alpha_tex_bool[i] && any(image_array_list[[i]][,,4] != 1)) {
          image_array_list[[i]][,,1] = image_array_list[[i]][,,4]
          image_array_list[[i]][,,2] = image_array_list[[i]][,,4]
          image_array_list[[i]][,,3] = image_array_list[[i]][,,4]
          png::writePNG(fliplr(aperm(image_array_list[[i]][,,1:3],c(2,1,3))), alpha_temp_file_names[i])
          alpha_list$alpha_tex_bool[i] = TRUE
          scene$material[[i]]$alphaimage = TRUE
        }
      } else if(dim(image_array_list[[i]])[3] == 3){
        png::writePNG(fliplr(aperm(image_array_list[[i]],c(2,1,3))),temp_file_names[i])
      }
      scene$material[[i]]$image = TRUE
    } else if(image_filename_bool[i]) {
      if(any(!file.exists(path.expand(image_array_list[[i]])) & nchar(image_array_list[[i]]) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(image_array_list[[i]], collapse="\n")))
      }
      temp_file_names[i] = path.expand(image_array_list[[i]])
      scene$material[[i]]$image = TRUE
    } else {
      scene$material[[i]]$image = FALSE
    }
  }
  image_list = list()
  image_tex_bool = image_tex_bool | image_filename_bool
  image_list$image_tex_bool = image_tex_bool
  image_list$image_temp_file_names = temp_file_names

  #bump texture handler
  bump_array_list = vector(mode="list", length(nrow(scene)))
  bump_tex_bool = vector(mode = "logical", length(nrow(scene)))
  bump_filename_bool = vector(mode = "logical", length(nrow(scene)))
  bump_temp_file_names = tempfile(sprintf("bumptemp%i",seq_len(nrow(scene))),fileext = ".png")
  for(i in seq(nrow(scene))) {
    bump_array_list[i] = scene$material[[i]]$bump_texture
    bump_tex_bool[i] = is.array(bump_array_list[[i]])
    bump_filename_bool[i] = is.character(bump_array_list[[i]])  && !is.na(bump_array_list[[i]])
    if(bump_tex_bool[i]) {
      bump_dims = dim(bump_array_list[[i]][[1]])
      if(length(bump_dims) == 2) {
        temp_array = array(0, dim = c(bump_dims,3))
        temp_array[,,1] = bump_array_list[[i]][[1]]
        temp_array[,,2] = bump_array_list[[i]][[1]]
        temp_array[,,3] = bump_array_list[[i]][[1]]
        bump_dims = c(bump_dims,3)
      } else {
        temp_array = bump_array_list[[i]][[1]]
      }
      if(bump_dims[3] == 4) {
        png::writePNG(fliplr(aperm(temp_array[,,1:3],c(2,1,3))),bump_temp_file_names[i])
      } else if(bump_dims[3] == 3){
        png::writePNG(fliplr(aperm(temp_array,c(2,1,3))),bump_temp_file_names[i])
      }
      scene$material[[i]]$bump_texture = TRUE
    } else if(bump_filename_bool[i]) {
      if(any(!file.exists(path.expand(bump_array_list[[i]][[1]])) & nchar(bump_array_list[[i]][[1]]) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(bump_array_list[[i]][[1]], collapse="\n")))
      }
      scene$material[[i]]$bump_texture = TRUE
      bump_temp_file_names[i] = path.expand(bump_array_list[[i]][[1]])
    } else {
      scene$material[[i]]$bump_texture = FALSE
    }
  }
  bump_tex_bool = bump_tex_bool | bump_filename_bool
  bump_list = list()
  bump_list$bump_temp_file_names = bump_temp_file_names
  bump_list$bump_tex_bool = bump_tex_bool
  
  #roughness texture handler
  roughness_array_list = vector(mode="list", length(nrow(scene)))
  rough_tex_bool = vector(mode = "logical", length(nrow(scene)))
  rough_filename_bool = vector(mode = "logical", length(nrow(scene)))
  rough_temp_file_names = tempfile(sprintf("roughtemp%i",seq_len(nrow(scene))),fileext = ".png")
  for(i in seq(nrow(scene))) {
    roughness_array_list[i] = scene$material[[i]]$roughness_texture
    rough_tex_bool[i] = is.array(roughness_array_list[[i]])
    rough_filename_bool[i] = is.character(roughness_array_list[[i]])  && !is.na(roughness_array_list[[i]])
    if(rough_tex_bool[i]) {
      if(length(dim(roughness_array_list[[i]][[1]])) == 2) {
        png::writePNG(fliplr(t(roughness_array_list[[i]][[1]])), rough_temp_file_names[i])
      } else if(dim(roughness_array_list[[i]][[1]])[3] == 3) {
        png::writePNG(fliplr(aperm(roughness_array_list[[i]][[1]],c(2,1,3))), rough_temp_file_names[i])
      } else {
        stop("alpha texture dims: c(", paste(dim(roughness_array_list[[i]][[1]]),collapse=", "), ") not valid for texture.")
      }
      scene$material[[i]]$roughness_texture = TRUE
    } else if(rough_filename_bool[i]) {
      if(any(!file.exists(path.expand(roughness_array_list[[i]][[1]])) & nchar(roughness_array_list[[i]][[1]]) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(roughness_array_list[[i]][[1]], collapse="\n")))
      }
      rough_temp_file_names[i] = path.expand(roughness_array_list[[i]][[1]])
      scene$material[[i]]$roughness_texture = TRUE
    } else{
      scene$material[[i]]$roughness_texture = FALSE
    }
  }
  rough_tex_bool = rough_tex_bool | rough_filename_bool
  roughness_list = list()
  roughness_list$rough_temp_file_names = rough_temp_file_names
  roughness_list$rough_tex_bool = rough_tex_bool
  
  for(i in seq(nrow(scene))) {
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
    
    # material_id = list()
    # for(i in seq_len(nrow(scene))) {
    #   #Handle instances
    #   if(typevec[i] == 15) {
    #     material_id[[i]] = unlist(lapply(scene$shape_info, \(x) x$material_id))
    #   }
    # }
    # material_id = unlist(material_id)
    # is_na_mat = is.na(material_id)
    # material_id_increasing = as.integer(match(material_id, unique(material_id)) - 1L)
    # for(i in seq_len(nrow(scene))) {
    #   if(typevec[i] == 15) {
    #     instance_scene = scene$shape_info[[i]]$shape_properties
    #     for(i in seq_len(length(scene))) {
    # 
    #     }
    # 
    #   } else {
    #     if(!is_na_mat[cntr]) {
    #       scene$shape_info[[i]]$material_id = material_id_increasing[cntr]
    #     }
    #     cntr = cntr + 1
    #   }
    # }
  }
  
  
  scene_info = list()
  scene_info$scene = scene
  scene_info$shape = shapevec
  scene_info$typevec = typevec
  scene_info$image_list = image_list
  scene_info$bump_list = bump_list
  scene_info$alpha_list = alpha_list
  scene_info$roughness_list = roughness_list
  
  return(scene_info)
}
