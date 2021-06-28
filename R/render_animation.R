#' Render Animation
#' 
#' Takes the scene description and renders an image, either to the device or to a filename. 
#'
#' @param scene Tibble of object locations and properties. 
#' @param camera_motion Data frame of camera motion vectors, calculated with `generate_camera_motion()`.
#' @param start_frame Default `1`. Frame to start the animation. 
#' @param width Default `400`. Width of the render, in pixels.
#' @param height Default `400`. Height of the render, in pixels.
#' @param samples Default `100`. The maximum number of samples for each pixel. If this is a length-2
#' vector and the `sample_method` is `stratified`, this will control the number of strata in each dimension.
#' The total number of samples in this case will be the product of the two numbers.
#' @param min_variance Default `0.00005`. Minimum acceptable variance for a block of pixels for the 
#' adaptive sampler. Smaller numbers give higher quality images, at the expense of longer rendering times.
#' If this is set to zero, the adaptive sampler will be turned off and the renderer
#' will use the maximum number of samples everywhere.
#' @param min_adaptive_size Default `8`. Width of the minimum block size in the adaptive sampler.
#' @param sample_method Default `sobol`. The type of sampling method used to generate
#' random numbers. The other options are `random` (worst quality but simple), 
#' `stratified` (only implemented for completion), 
#' and `sobol_blue` (best option for sample counts below 256). 
#' @param max_depth Default `50`. Maximum number of bounces a ray can make in a scene.
#' @param roulette_active_depth Default `10`. Number of ray bounces until a ray can stop bouncing via
#' Russian roulette.
#' @param ambient_light Default `FALSE`, unless there are no emitting objects in the scene. 
#' If `TRUE`, the background will be a gradient varying from `backgroundhigh` directly up (+y) to 
#' `backgroundlow` directly down (-y).
#' @param clamp_value Default `Inf`. If a bright light or a reflective material is in the scene, occasionally
#' there will be bright spots that will not go away even with a large number of samples. These 
#' can be removed (at the cost of slightly darkening the image) by setting this to a small number greater than 1. 
#' @param filename Default `NULL`. If present, the renderer will write to the filename instead
#' of the current device.
#' @param backgroundhigh Default `#80b4ff`. The "high" color in the background gradient. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param backgroundlow Default `#ffffff`. The "low" color in the background gradient. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param shutteropen Default `0`. Time at which the shutter is open. Only affects moving objects.
#' @param shutterclose Default `1`. Time at which the shutter is open. Only affects moving objects.
#' @param focal_distance Default `NULL`, automatically set to the `lookfrom-lookat` distance unless
#' otherwise specified.
#' @param ortho_dimensions Default `c(1,1)`. Width and height of the orthographic camera. Will only be used if `fov = 0`. 
#' @param tonemap Default `gamma`. Choose the tone mapping function,
#' Default `gamma` solely adjusts for gamma and clamps values greater than 1 to 1. 
#' `reinhold` scales values by their individual color channels `color/(1+color)` and then performs the 
#' gamma adjustment. `uncharted` uses the mapping developed for Uncharted 2 by John Hable. `hbd` uses an
#' optimized formula by Jim Hejl and Richard Burgess-Dawson. Note: If set to anything other than `gamma`,
#' objects with material `light()` may not be anti-aliased. If `raw`, the raw array of HDR values will be
#' returned, rather than an image or a plot.
#' @param bloom Default `TRUE`. Set to `FALSE` to get the raw, pathtraced image. Otherwise,
#' this performs a convolution of the HDR image of the scene with a sharp, long-tailed
#' exponential kernel, which does not visibly affect dimly pixels, but does result in emitters light
#' slightly bleeding into adjacent pixels. This provides an antialiasing effect for lights, even when
#' tonemapping the image. Pass in a matrix to specify the convolution kernel manually, or a positive number
#' to control the intensity of the bloom (higher number = more bloom).
#' @param environment_light Default `NULL`. An image to be used for the background for rays that escape
#' the scene. Supports both HDR (`.hdr`) and low-dynamic range (`.png`, `.jpg`) images.
#' @param rotate_env Default `0`. The number of degrees to rotate the environment map around the scene.
#' @param intensity_env Default `1`. The amount to increase the intensity of the environment lighting. Useful
#' if using a LDR (JPEG or PNG) image as an environment map.
#' @param debug_channel Default `none`. If `depth`, function will return a depth map of rays into the scene 
#' instead of an image. If `normals`, function will return an image of scene normals, mapped from 0 to 1.
#' If `uv`, function will return an image of the uv coords. If `variance`, function will return an image 
#' showing the number of samples needed to take for each block to converge. If `dpdu` or `dpdv`, function will return
#' an image showing the differential `u` and `u` coordinates. If `color`, function will return the raw albedo
#' values (with white for `metal` and `dielectric` materials). If `preview`, an image rendered with `render_preview()` 
#' will be returned.
#' @param return_raw_array Default `FALSE`. If `TRUE`, function will return raw array with RGB intensity
#' information.
#' @param parallel Default `FALSE`. If `TRUE`, it will use all available cores to render the image
#'  (or the number specified in `options("cores")` if that option is not `NULL`).
#' @param bvh_type Default `"sah"`, "surface area heuristic". Method of building the bounding volume
#' hierarchy structure used when rendering. Other option is "equal", which splits tree into groups
#' of equal size.
#' @param progress Default `TRUE` if interactive session, `FALSE` otherwise. 
#' @param preview_light_direction Default `c(0,-1,0)`. Vector specifying the orientation for the global light using for phong shading.
#' @param preview_exponent Default `6`. Phong exponent.  
#' @param verbose Default `FALSE`. Prints information and timing information about scene
#' construction and raytracing progress.
#' @export
#' @importFrom  grDevices col2rgb
#' @return Raytraced plot to current device, or an image saved to a file. 
#'
#' @examples
#' #Create and animate flying through a scene on a simulated roller coaster
#' \donttest{
#' set.seed(3)
#' elliplist = list()
#' ellip_colors = rainbow(8)
#' for(i in 1:1200) {
#'   elliplist[[i]] = ellipsoid(x=10*runif(1)-5,y=10*runif(1)-5,z=10*runif(1)-5,
#'                              angle = 360*runif(3), a=0.1,b=0.05,c=0.1,
#'                              material=glossy(color=sample(ellip_colors,1)))
#' }
#' ellip_scene = do.call(rbind, elliplist)
#' 
#' camera_pos = list(c(0,1,15),c(5,-5,5),c(-5,5,-5),c(0,1,-15))
#' 
#' #Plot the camera path and render from above using the path object:
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(0,20,0), width=800,height=800,samples=32,
#'                camera_up = c(0,0,1),
#'                fov=80)
#'             
#' #Side view     
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(20,0,0),width=800,height=800,samples=32,
#'                  fov=80)
#'  
#' #View from the start        
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(0,1.5,16),width=800,height=800,samples=32,
#'                  fov=80)
#'                  
#' #Generate Camera movement, setting the lookat position to be same as camera position, but offset
#' #slightly in front. We'll render 12 frames, but you'd likely want more in a real animation.
#' 
#' camera_motion =  generate_camera_motion(positions = camera_pos, lookats = camera_pos, 
#'                                         offset_lookat = 1, fovs=80, frames=12) 
#'                                         
#' #This returns a data frame of individual camera positions, interpolated by cubic bezier curves.
#' camera_motion
#' 
#' #Pass NA filename to plot to the device. We'll keep the path and offset it slightly to see
#' #where we're going. This results in a "roller coaster" effect.
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(obj_model(r_obj(),x=10,y=-10,scale_obj=3, angle=c(0,-45,0),
#'                        material=dielectric(attenuation=c(1,1,0.3)))) %>% 
#'   add_object(pig(x=-7,y=10,z=-5,scale=1,angle=c(0,-45,80),emotion="angry")) %>% 
#'   add_object(pig(x=0,y=-0.25,z=-15,scale=1,angle=c(0,225,-20), order_rotation=c(3,2,1),
#'                  emotion="angry", spider=TRUE)) %>% 
#'   add_object(path(camera_pos, y=-0.2,material=diffuse(color="red"))) %>% 
#'   render_animation(filename = NA, camera_motion = camera_motion, samples=100,
#'                    sample_method="sobol_blue",  
#'                    clamp_value=10, width=400, height=400)
#' 
#' }
render_animation = function(scene, camera_motion, start_frame = 1,
                            width = 400, height = 400, 
                            samples = 100, min_variance = 0.00005, min_adaptive_size = 8,
                            sample_method = "sobol",
                            max_depth = 50, roulette_active_depth = 10,
                            ambient_light = FALSE, 
                            clamp_value = Inf,
                            filename = "rayimage", backgroundhigh = "#80b4ff",backgroundlow = "#ffffff",
                            shutteropen = 0.0, shutterclose = 1.0, focal_distance=NULL, ortho_dimensions = c(1,1),
                            tonemap ="gamma", bloom = TRUE, parallel=TRUE, bvh_type = "sah",
                            environment_light = NULL, rotate_env = 0, intensity_env = 1,
                            debug_channel = "none", return_raw_array = FALSE,
                            progress = interactive(), verbose = FALSE,
                            preview_light_direction = c(0,-1,0), preview_exponent = 6) { 
  if(verbose) {
    currenttime = proc.time()
    cat("Building Scene: ")
  }
  if(is.na(filename)) {
    filename = ""
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
                          "glossy" = 7, "spotlight" = 8, "hair" = 9))
  sigmavec = unlist(scene$sigma)
  
  if(!tonemap %in% c("gamma","reinhold","uncharted", "hbd", "raw")) {
    stop("tonemap value ", tonemap, " not recognized")
  }
  toneval = switch(tonemap, "gamma" = 1,"reinhold" = 2,"uncharted" = 3,"hbd" = 4, "raw" = 5)
  movingvec = purrr::map_lgl(scene$velocity,.f = ~any(.x != 0))
  proplist = scene$properties
  vel_list = scene$velocity
  
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
  
  if(!any(typevec == 5) && !any(typevec == 8) && missing(ambient_light) && missing(environment_light)) {
    ambient_light = TRUE
  }
  
  #texture handler
  image_array_list = scene$image
  image_tex_bool = purrr::map_lgl(image_array_list,.f = ~is.array(.x))
  image_filename_bool = purrr::map_lgl(image_array_list,.f = ~is.character(.x))
  temp_file_names = purrr::map_chr(image_tex_bool,.f = ~ifelse(.x, tempfile(fileext = ".png"),""))
  for(i in 1:length(image_array_list)) {
    if(image_tex_bool[i]) {
      if(dim(image_array_list[[i]])[3] == 4) {
        png::writePNG(fliplr(aperm(image_array_list[[i]][,,1:3],c(2,1,3))),temp_file_names[i])
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
  
  #movement handler
  if(shutteropen == shutterclose) {
    movingvec = rep(FALSE,length(movingvec))
  }
  
  #implicit sampling handler
  implicit_vec = scene$implicit_sample
  
  #order rotation handler
  order_rotation_list = scene$order_rotation
  
  #group handler
  group_bool = purrr::map_lgl(scene$pivot_point,.f = ~all(!is.na(.x)))
  group_pivot = scene$pivot_point 
  group_angle = scene$group_angle 
  group_order_rotation = scene$group_order_rotation 
  group_translate = scene$group_translate 
  group_scale = scene$group_scale
  
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
  objbasedirvec = purrr::map_chr(objfilenamevec, dirname)
  
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
  
  if(!is.null(options("cores")[[1]])) {
    numbercores = options("cores")[[1]]
  } else {
    numbercores = parallel::detectCores()
  }
  if(!parallel) {
    numbercores = 1
  }
  if(debug_channel != "preview") {
    debug_channel = unlist(lapply(tolower(debug_channel),switch,
                                  "none" = 0,"depth" = 1,"normals" = 2, "uv" = 3, "bvh" = 4,
                                  "variance" = 5, "normal" = 2, "dpdu" = 6, "dpdv" = 7, "color" = 8, 
                                  0))
    light_direction = c(0,1,0,0)
  } else {
    light_direction = c(-preview_light_direction, preview_exponent)
    debug_channel = 9
  }
  if(debug_channel == 4) {
    message("rayrender must be compiled with option DEBUGBVH for this debug option to work")
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
  
  camera_info$nx = width
  camera_info$ny = height
  camera_info$ns = samples
  camera_info$shutteropen = shutteropen
  camera_info$shutterclose = shutterclose
  camera_info$max_depth = max_depth
  camera_info$roulette_active_depth = roulette_active_depth
  camera_info$sample_method = sample_method
  camera_info$stratified_dim = strat_dim
  camera_info$light_direction = light_direction
  camera_info$bvh = switch(bvh_type,"sah" = 1, "equal" = 2, 1)
  
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
  scene_info$velocity = vel_list
  scene_info$moving = movingvec
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
  scene_info$group_pivot=group_pivot
  scene_info$group_translate = group_translate
  scene_info$group_angle = group_angle
  scene_info$group_order_rotation = group_order_rotation
  scene_info$group_scale = group_scale
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
  
  #Camera Movement Info
  if(filename != "") {
    filename_str = paste0(filename,str_pad(1:nrow(camera_motion),width=ceiling(log10(nrow(camera_motion))),pad="0"),".png")
  } else {
    filename_str = rep("", length(1:nrow(camera_motion)))
  }
  
  #Pathrace Scene
  rgb_mat = render_animation_rcpp(camera_info = camera_info, scene_info = scene_info, camera_movement = camera_motion,
                              start_frame = start_frame - 1, filenames = filename_str, post_process_frame  = post_process_frame,
                              toneval=toneval, bloom = bloom) 
}
