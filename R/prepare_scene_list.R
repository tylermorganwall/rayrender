#'@title Prepare the scene list
#'
#'@keywords internal
#'@examples
#'#internal
prepare_scene_list = function(scene, width = 400, height = 400, fov = 20, 
                              lookfrom = c(0,1,10), lookat = c(0,0,0), camera_up = c(0,1,0), 
                              samples = 100,  camera_description_file = NA, 
                              camera_scale = 1, iso = 100, film_size = 22,
                              min_variance = 0.00005, min_adaptive_size = 8,
                              sample_method = "sobol", 
                              max_depth = NA, roulette_active_depth = 100,
                              ambient_light = FALSE, 
                              aperture = 0.1, clamp_value = Inf,
                              filename = NULL, backgroundhigh = "#80b4ff",backgroundlow = "#ffffff",
                              shutteropen = 0.0, shutterclose = 1.0, focal_distance=NULL, ortho_dimensions = c(1,1),
                              tonemap ="gamma", bloom = TRUE, parallel=TRUE, bvh_type = "sah",
                              environment_light = NULL, rotate_env = 0, intensity_env = 1,
                              debug_channel = "none", return_raw_array = FALSE,
                              progress = interactive(), verbose = FALSE, sample_dist = Inf,
                              keep_colors = FALSE) {
  if(verbose) {
    currenttime = proc.time()
    cat("Building Scene: ")
  }
  if(debug_channel == "none" && is.na(max_depth)) {
    max_depth = 50
  }
  if(debug_channel != "none" && is.na(max_depth)) {
    max_depth = 1
  }
  iso = iso/100
  film_size = film_size/1000
  lookvec = (lookat - lookfrom)
  i1 = c(2,3,1)
  i2 = c(3,1,2)
  if(all(lookvec[i1]*camera_up[i2] - lookvec[i2]*camera_up[i1] == 0)) {
    stop("camera_up value c(", paste(camera_up, collapse=","), ") is aligned exactly with camera vector (lookat - lookfrom). Choose a different value for camera_up.")
  }
  backgroundhigh = convert_color(backgroundhigh)
  backgroundlow = convert_color(backgroundlow)
  position_list = list()
  position_list$xvec = scene$x 
  position_list$yvec = scene$y
  position_list$zvec = scene$z
  rvec = scene$radius
  shapevec = unlist(lapply(tolower(scene$shape),switch,
                           "sphere" = 1,"xy_rect" = 2, "xz_rect" = 3,"yz_rect" = 4,"box" = 5, "triangle" = 6, 
                           "obj" = 7, "objcolor" = 8, "disk" = 9, "cylinder" = 10, "ellipsoid" = 11,
                           "objvertexcolor" = 12, "cone" = 13, "curve" = 14, "csg_object" = 15, "ply" = 16,
                           "mesh3d" = 17))
  typevec = unlist(lapply(tolower(scene$type),switch,
                          "diffuse" = 1,"metal" = 2,"dielectric" = 3, 
                          "oren-nayar" = 4, "light" = 5, "microfacet" = 6, 
                          "glossy" = 7, "spotlight" = 8, "hair" = 9, "microfacet_transmission" = 10))
  sigmavec = unlist(scene$sigma)
  
  if(!tonemap %in% c("gamma","reinhold","uncharted", "hbd", "raw")) {
    stop("tonemap value ", tonemap, " not recognized")
  }
  toneval = switch(tonemap, "gamma" = 1,"reinhold" = 2,"uncharted" = 3,"hbd" = 4, "raw" = 5)
  proplist = scene$properties
  
  checkeredlist = scene$checkercolor
  checkeredbool = purrr::map_lgl(checkeredlist,.f = ~all(!is.na(.x)))
  
  #glossy 
  glossyinfo = scene$glossyinfo
  
  #gradient handler
  gradient_info = list()
  gradient_info$gradient_colors = scene$gradient_color
  gradient_info$isgradient = purrr::map_lgl(gradient_info$gradient_colors,.f = ~all(!is.na(.x)))
  gradient_info$gradient_trans = scene$gradient_transpose
  gradient_info$is_world_gradient = scene$world_gradient
  gradient_info$gradient_control_points = scene$gradient_point_info
  gradient_info$type = unlist(lapply(tolower(scene$gradient_type),switch,
                                     "hsv" = TRUE, "rgb" = FALSE, FALSE))
  #noise handler
  noisebool = purrr::map_lgl(scene$noise, .f = ~.x > 0)
  noisevec = scene$noise
  noisephasevec = scene$noisephase * pi/180
  noiseintvec = scene$noiseintensity
  noisecolorlist = scene$noisecolor
  
  #rotation handler
  rot_angle_list = scene$angle
  
  #fog handler
  fog_bool = scene$fog
  fog_vec = scene$fogdensity
  
  #flip handler
  flip_vec = scene$flipped
  
  #light handler
  light_prop_vec =  scene$lightintensity
  
  if(!any(typevec == 5) && !any(typevec == 8) && is.null(ambient_light) && is.null(environment_light)) {
    ambient_light = TRUE
  } else {
    if(is.null(ambient_light)) {
      ambient_light = FALSE
    }
  }
  
  #alpha texture handler
  alpha_array_list = scene$alphaimage
  alpha_tex_bool = purrr::map_lgl(alpha_array_list,.f = ~is.array(.x[[1]]))
  alpha_filename_bool = purrr::map_lgl(alpha_array_list,.f = ~is.character(.x[[1]]))
  alpha_temp_file_names = purrr::map_chr(alpha_tex_bool, .f = (function(.x) tempfile(fileext = ".png")))
  for(i in 1:length(alpha_array_list)) {
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
    }
    if(alpha_filename_bool[i]) {
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
      png::writePNG(temp_array,alpha_temp_file_names[i])
    }
  }
  alpha_tex_bool = alpha_tex_bool | alpha_filename_bool
  alphalist = list()
  alphalist$alpha_temp_file_names = alpha_temp_file_names
  alphalist$alpha_tex_bool = alpha_tex_bool
  
  #texture handler
  image_array_list = scene$image
  image_tex_bool = purrr::map_lgl(image_array_list,.f = ~is.array(.x))
  image_filename_bool = purrr::map_lgl(image_array_list,.f = ~is.character(.x))
  temp_file_names = purrr::map_chr(image_tex_bool,.f = ~ifelse(.x, tempfile(fileext = ".png"),""))
  for(i in 1:length(image_array_list)) {
    if(image_tex_bool[i]) {
      if(dim(image_array_list[[i]])[3] == 4) {
        png::writePNG(fliplr(aperm(image_array_list[[i]][,,1:3],c(2,1,3))),temp_file_names[i])
        #Handle PNG with alpha
        if(!alpha_tex_bool[i] && any(image_array_list[[i]][,,4] != 1)) {
          image_array_list[[i]][,,1] = image_array_list[[i]][,,4]
          image_array_list[[i]][,,2] = image_array_list[[i]][,,4]
          image_array_list[[i]][,,3] = image_array_list[[i]][,,4]
          png::writePNG(fliplr(aperm(image_array_list[[i]][,,1:3],c(2,1,3))), alpha_temp_file_names[i])
          alphalist$alpha_tex_bool[i] = TRUE
        }
      } else if(dim(image_array_list[[i]])[3] == 3){
        png::writePNG(fliplr(aperm(image_array_list[[i]],c(2,1,3))),temp_file_names[i])
      }
    }
    if(image_filename_bool[i]) {
      if(any(!file.exists(path.expand(image_array_list[[i]])) & nchar(image_array_list[[i]]) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(image_array_list[[i]], collapse="\n")))
      }
      temp_file_names[i] = path.expand(image_array_list[[i]])
    }
  }
  image_tex_bool = image_tex_bool | image_filename_bool
  image_repeat = scene$image_repeat
  
  #bump texture handler
  bump_array_list = scene$bump_texture
  bump_tex_bool = purrr::map_lgl(bump_array_list,.f = ~is.array(.x[[1]]))
  bump_filename_bool = purrr::map_lgl(bump_array_list,.f = ~is.character(.x[[1]]))
  bump_temp_file_names = purrr::map_chr(bump_tex_bool,.f = ~ifelse(.x, tempfile(fileext = ".png"),""))
  for(i in 1:length(bump_array_list)) {
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
    }
    if(bump_filename_bool[i]) {
      if(any(!file.exists(path.expand(bump_array_list[[i]][[1]])) & nchar(bump_array_list[[i]][[1]]) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(bump_array_list[[i]][[1]], collapse="\n")))
      }
      bump_temp_file_names[i] = path.expand(bump_array_list[[i]][[1]])
    }
  }
  bump_tex_bool = bump_tex_bool | bump_filename_bool
  bump_intensity = scene$bump_intensity
  alphalist$bump_temp_file_names = bump_temp_file_names
  alphalist$bump_tex_bool = bump_tex_bool
  alphalist$bump_intensity = bump_intensity
  
  #roughness texture handler
  roughness_array_list = scene$roughness_texture
  rough_tex_bool = purrr::map_lgl(roughness_array_list,.f = ~is.array(.x[[1]]))
  rough_filename_bool = purrr::map_lgl(roughness_array_list,.f = ~is.character(.x[[1]]))
  rough_temp_file_names = purrr::map_chr(rough_tex_bool, .f = (function(.x) tempfile(fileext = ".png")))
  for(i in 1:length(roughness_array_list)) {
    if(rough_tex_bool[i]) {
      tempgloss = glossyinfo[[i]]
      if(length(dim(roughness_array_list[[i]][[1]])) == 2) {
        png::writePNG(fliplr(t(roughness_array_list[[i]][[1]])), rough_temp_file_names[i])
      } else if(dim(roughness_array_list[[i]][[1]])[3] == 3) {
        png::writePNG(fliplr(aperm(roughness_array_list[[i]][[1]],c(2,1,3))), rough_temp_file_names[i])
      } else {
        stop("alpha texture dims: c(", paste(dim(roughness_array_list[[i]][[1]]),collapse=", "), ") not valid for texture.")
      }
    }
    if(rough_filename_bool[i]) {
      if(any(!file.exists(path.expand(roughness_array_list[[i]][[1]])) & nchar(roughness_array_list[[i]][[1]]) > 0)) {
        stop(paste0("Cannot find the following texture file:\n",
                    paste(roughness_array_list[[i]][[1]], collapse="\n")))
      }
      rough_temp_file_names[i] = path.expand(roughness_array_list[[i]][[1]])
    }
  }
  rough_tex_bool = rough_tex_bool | rough_filename_bool
  roughness_list = list()
  roughness_list$rough_temp_file_names = rough_temp_file_names
  roughness_list$rough_tex_bool = rough_tex_bool
  
  
  #implicit sampling handler
  implicit_vec = scene$implicit_sample
  
  #order rotation handler
  order_rotation_list = scene$order_rotation
  
  #group handler
  group_bool = purrr::map_lgl(scene$group_transform,.f = ~all(!is.na(.x)))
  group_transform = scene$group_transform
  
  
  #animation handler
  animation_bool = purrr::map_lgl(scene$start_transform_animation,.f = ~all(!is.na(.x))) & 
    purrr::map_lgl(scene$end_transform_animation,.f = ~all(!is.na(.x)))
  start_transform_animation = scene$start_transform_animation
  end_transform_animation = scene$end_transform_animation
  animation_start_time = scene$start_time
  animation_end_time = scene$end_time
  
  #triangle normal handler
  tri_normal_bools = purrr::map2_lgl(shapevec,proplist,.f = ~.x == 6 && all(!is.na(.y)))
  tri_color_vert = scene$tricolorinfo
  is_tri_color = purrr::map_lgl(tri_color_vert,.f = ~all(!is.na(.x)))
  
  #obj handler
  fileinfovec = scene$fileinfo
  fileinfovec[is.na(fileinfovec)] = ""
  objfilenamevec = purrr::map_chr(fileinfovec, path.expand)
  if(any(!file.exists(objfilenamevec) & nchar(objfilenamevec) > 0)) {
    stop(paste0("Cannot find the following obj/ply files:\n",
                paste(objfilenamevec[!file.exists(objfilenamevec) & nchar(objfilenamevec) > 0], 
                      collapse="\n")
    ))
  }
  base_dir = function(x) {
    dirname_processed = dirname(x)
    if(dirname_processed == ".") {
      return("")
    } else {
      return(dirname_processed)
    }
  }
  objbasedirvec = purrr::map_chr(objfilenamevec, base_dir)
  
  #bg image handler
  if(!is.null(environment_light)) {
    hasbackground = TRUE
    backgroundstring = path.expand(environment_light)
    if(!file.exists(environment_light)) {
      hasbackground = FALSE
      warning("file '", environment_light, "' cannot be found, not using background image.")
    }
    if(dir.exists(environment_light)) {
      stop("environment_light argument '", environment_light, "' is a directory, not a file.")
    }
  } else {
    hasbackground = FALSE
    backgroundstring = ""
  }
  
  #scale handler
  scale_factor = scene$scale_factor
  
  if(length(lookfrom) != 3) {
    stop("lookfrom must be length-3 numeric vector")
  }
  if(length(lookat) != 3) {
    stop("lookat must be length-3 numeric vector")
  }
  if(is.null(focal_distance)) {
    focal_distance = sqrt(sum((lookfrom-lookat)^2))
  }
  if(!is.null(options("cores")[[1]])) {
    numbercores = options("cores")[[1]]
  } else {
    numbercores = parallel::detectCores()
  }
  if(!parallel) {
    numbercores = 1
  }
  if(!is.numeric(debug_channel)) {
    debug_channel = unlist(lapply(tolower(debug_channel),switch,
                                  "none" = 0,"depth" = 1,"normals" = 2, "uv" = 3, "bvh" = 4,
                                  "variance" = 5, "normal" = 2, "dpdu" = 6, "dpdv" = 7, "color" = 8, 
                                  "position" = 10, "direction" = 11, "time" = 12, "shape" = 13,
                                  "pdf" = 14, "error" = 15, "bounces" = 16, "camera" = 17,
                                  "ao" = 18, 0))
    light_direction = c(0,1,0)
  } else {
    light_direction = debug_channel
    debug_channel = 9
  }
  if(debug_channel == 4) {
    message("rayrender must be compiled with option DEBUGBVH for this debug option to work")
  }
  
  if(fov == 0) {
    if(length(ortho_dimensions) != 2) {
      stop("ortho_dimensions must be length-2 numeric vector")
    }
  }
  if(verbose) {
    buildingtime = proc.time() - currenttime
    cat(sprintf("%0.3f seconds \n",buildingtime[3]))
  }
  sample_method = unlist(lapply(tolower(sample_method),switch,
                                "random" = 0,"stratified" = 1, "sobol" = 2,"sobol_blue" = 3, 0))
  
  camera_info = list()
  strat_dim = c()
  if(length(samples) == 2) {
    strat_dim = samples
    samples = samples[1]*samples[2]
  } else {
    strat_dim = rep(min(floor(sqrt(samples)),8),2)
  }
  
  fov = ifelse(fov < 0, 0, fov);
  if(!is.na(camera_description_file)) {
    camera_description_file = switch(camera_description_file, 
                                     "50mm" = system.file("extdata","dgauss.50mm.txt",
                                                          package="rayrender"),
                                     "wide" = system.file("extdata","wide.22mm.txt",
                                                          package="rayrender"),
                                     "fisheye" = system.file("extdata","fisheye.10mm.txt",
                                                             package="rayrender"),
                                     "telephoto" = system.file("extdata","telephoto.250mm.txt",
                                                               package="rayrender"),
                                     camera_description_file)
    
    if(file.exists(camera_description_file)) {
      real_camera_info = as.matrix(utils::read.delim(camera_description_file, header=FALSE, comment.char="#"))
      fov = -1
    } else {
      warning("Camera description file `", camera_description_file, "` not found. Ignoring.")
      camera_description_file = NA
      real_camera_info = matrix(nrow=0,ncol=4)
    }
  } else {
    real_camera_info = matrix(nrow=0,ncol=4)
  }
  
  camera_info$nx = width
  camera_info$ny = height
  camera_info$ns = samples
  camera_info$fov = fov
  camera_info$lookfrom = lookfrom
  camera_info$lookat = lookat
  camera_info$aperture = aperture
  camera_info$camera_up = camera_up
  camera_info$shutteropen = shutteropen
  camera_info$shutterclose = shutterclose
  camera_info$ortho_dimensions = ortho_dimensions
  camera_info$focal_distance = focal_distance
  camera_info$max_depth = max_depth
  camera_info$roulette_active_depth = roulette_active_depth
  camera_info$sample_method = sample_method
  camera_info$stratified_dim = strat_dim
  camera_info$light_direction = light_direction
  camera_info$bvh = switch(bvh_type,"sah" = 1, "equal" = 2, 1)
  camera_info$real_camera_info = real_camera_info
  camera_info$film_size = film_size
  camera_info$camera_scale = camera_scale
  camera_info$sample_dist = sample_dist
  camera_info$keep_colors = keep_colors
  camera_info$iso = iso
  
  animation_info = list()
  animation_info$animation_bool            = animation_bool            
  animation_info$start_transform_animation = start_transform_animation 
  animation_info$end_transform_animation   = end_transform_animation   
  animation_info$animation_start_time      = animation_start_time      
  animation_info$animation_end_time        = animation_end_time        
  
  if(max_depth <= 0) {
    stop("max_depth must be greater than zero")
  }
  if(roulette_active_depth <= 0) {
    stop("roulette_active_depth must be greater than zero")
  }
  
  #Spotlight handler
  if(any(typevec == 8)) {
    if(any(shapevec[typevec == 8] > 4)) {
      stop("spotlights are only supported for spheres and rects")
    }
    for(i in 1:length(proplist)) {
      if(typevec[i] == 8) {
        proplist[[i]][4:6] = proplist[[i]][4:6] - c(position_list$xvec[i],position_list$yvec[i],position_list$zvec[i]) 
      }
    }
  }
  
  
  #Material ID handler; these must show up in increasing order.  Note, this will
  #cause problems if `match` is every changed to return doubles when matching in
  #long vectors as has happened with `which` recently.
  material_id = scene$material_id
  material_id = as.integer(match(material_id, unique(material_id)) - 1L)
  material_id_bool = !is.na(scene$material_id)
  
  if(min_adaptive_size < 1) {
    warning("min_adaptive_size cannot be less than one: setting to one")
    min_adaptive_size = 1
  }
  
  if(min_variance < 0) {
    stop("min_variance cannot be less than zero")
  }
  
  #CSG handler
  csg_list = scene$csg_object
  csg_info = list()
  csg_info$csg = csg_list
  
  #mesh3d handler
  mesh_list = scene$mesh_info
  
  
  scene_info = list()
  scene_info$ambient_light = ambient_light
  scene_info$type = typevec
  scene_info$shape = shapevec
  scene_info$radius = rvec
  scene_info$position_list = position_list
  scene_info$properties = proplist
  scene_info$n = length(typevec)
  scene_info$bghigh = backgroundhigh
  scene_info$bglow = backgroundlow
  scene_info$ischeckered = checkeredbool
  scene_info$checkercolors = checkeredlist
  scene_info$gradient_info = gradient_info
  scene_info$noise=noisevec
  scene_info$isnoise=noisebool
  scene_info$noisephase=noisephasevec
  scene_info$noiseintensity=noiseintvec
  scene_info$noisecolorlist = noisecolorlist
  scene_info$angle = rot_angle_list
  scene_info$isimage = image_tex_bool
  scene_info$filelocation = temp_file_names
  scene_info$alphalist = alphalist
  scene_info$lightintensity = light_prop_vec
  scene_info$isflipped = flip_vec
  scene_info$isvolume=fog_bool
  scene_info$voldensity = fog_vec
  scene_info$implicit_sample = implicit_vec
  scene_info$order_rotation_list = order_rotation_list
  scene_info$clampval = clamp_value
  scene_info$isgrouped = group_bool  
  scene_info$group_transform= group_transform
  scene_info$tri_normal_bools = tri_normal_bools
  scene_info$is_tri_color = is_tri_color
  scene_info$tri_color_vert= tri_color_vert
  scene_info$fileinfo = objfilenamevec
  scene_info$filebasedir = objbasedirvec
  scene_info$progress_bar = progress
  scene_info$numbercores = numbercores
  scene_info$hasbackground = hasbackground
  scene_info$background = backgroundstring
  scene_info$scale_list = scale_factor
  scene_info$sigmavec = sigmavec
  scene_info$rotate_env = rotate_env
  scene_info$intensity_env = intensity_env
  scene_info$verbose = verbose
  scene_info$debug_channel = debug_channel
  scene_info$shared_id_mat=material_id
  scene_info$is_shared_mat=material_id_bool
  scene_info$min_variance = min_variance
  scene_info$min_adaptive_size = min_adaptive_size
  scene_info$glossyinfo = glossyinfo
  scene_info$image_repeat = image_repeat
  scene_info$csg_info = csg_info
  scene_info$mesh_list=mesh_list
  scene_info$roughness_list = roughness_list
  scene_info$animation_info = animation_info
  
  all_info = list()
  all_info$scene_info = scene_info
  all_info$camera_info = camera_info
  
  return(all_info)
}