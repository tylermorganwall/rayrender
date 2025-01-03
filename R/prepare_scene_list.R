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
                              keep_colors = FALSE, integrator_type = "nee", denoise = TRUE,
                              print_debug_info = FALSE) {
  #Process images, convert shapes and materials to enums, extract positions, and 
  scene_info = process_scene(scene)
  if(!is.numeric(debug_channel)) {
    debug_channel = unlist(lapply(tolower(debug_channel),switch,
                                  "none" = 0,"depth" = 1,"normals" = 2, "uv" = 3, "bvh" = 4,
                                  "variance" = 5, "normal" = 2, "dpdu" = 6, "dpdv" = 7, "color" = 8, 
                                  "position" = 10, "direction" = 11, "time" = 12, "shape" = 13,
                                  "pdf" = 14, "error" = 15, "bounces" = 16, "camera" = 17,
                                  "ao" = 18, "material" = 19, 0))
    if(debug_channel != 0) {
      denoise = FALSE
    }
    light_direction = c(0,1,0)
  } else {
    denoise = FALSE
    light_direction = debug_channel
    debug_channel = 9
  }
  integrator_type = switch(integrator_type, "nee" = 1L, "rtiow" = 2L, "basic" = 3L, 
                    stop(integrator_type, " not recognized as valid `integrator_type`"))
  if(debug_channel == 4) {
    message("rayrender must be compiled with option DEBUGBVH for this debug option to work")
  }
  if(debug_channel == 0 && is.na(max_depth)) {
    max_depth = 50
  }
  if(debug_channel == 16 && is.na(max_depth)) {
    max_depth = 10000
  }
  if(debug_channel != 0 && is.na(max_depth)) {
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
  
  if(!tonemap %in% c("gamma","reinhold","uncharted", "hbd", "raw")) {
    stop("tonemap value ", tonemap, " not recognized")
  }
  
  if(!scene_info$any_light && 
     is.null(ambient_light) && 
     is.null(environment_light)) {
    ambient_light = TRUE
  } else {
    if(is.null(ambient_light)) {
      ambient_light = FALSE
    }
  }
  
  #Background image handler
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
  if(length(lookfrom) != 3) {
    stop("lookfrom must be length-3 numeric vector")
  }
  if(length(lookat) != 3) {
    stop("lookat must be length-3 numeric vector")
  }
  if(is.null(focal_distance)) {
    focal_distance = sqrt(sum((lookfrom-lookat)^2))
  }
  numbercores = getOption("cores", default = getOption("Ncpus", default = parallel::detectCores()))
  if(!parallel) {
    numbercores = 1
  }
  
  if(fov == 0) {
    if(length(ortho_dimensions) != 2) {
      stop("ortho_dimensions must be length-2 numeric vector")
    }
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
  
  fov = ifelse(fov < 0, 0, fov)
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
  
  if(max_depth <= 0) {
    stop("max_depth must be greater than zero")
  }
  if(roulette_active_depth <= 0) {
    stop("roulette_active_depth must be greater than zero")
  }
  
  if(min_adaptive_size < 1) {
    warning("min_adaptive_size cannot be less than one: setting to one")
    min_adaptive_size = 1
  }
  
  if(min_variance < 0) {
    stop("min_variance cannot be less than zero")
  } else if (min_variance == 0) {
    min_variance = 0L
  }
  
  render_info = list()
  render_info$ambient_light = ambient_light
  render_info$bghigh = backgroundhigh
  render_info$bglow = backgroundlow
  render_info$clampval = clamp_value
  render_info$progress_bar = progress
  render_info$numbercores = numbercores
  render_info$hasbackground = hasbackground
  render_info$background = backgroundstring
  render_info$rotate_env = rotate_env
  render_info$intensity_env = intensity_env
  render_info$verbose = verbose
  render_info$debug_channel = debug_channel
  render_info$min_variance = min_variance
  render_info$min_adaptive_size = min_adaptive_size
  render_info$integrator_type = integrator_type
  render_info$denoise = denoise
  render_info$print_debug_info = print_debug_info

  all_info = list()
  all_info$scene_info = scene_info
  all_info$camera_info = camera_info
  all_info$render_info = render_info
  all_info$scene = scene_info$scene
  
  return(all_info)
}
