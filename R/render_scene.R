#' Render Scene
#' 
#' Takes the scene description and renders an image, either to the device or to a filename. 
#'
#' @param scene Tibble of object locations and properties. 
#' @param width Default `400`. Width of the render, in pixels.
#' @param height Default `400`. Height of the render, in pixels.
#' @param fov Default `20`. Field of view, in degrees. If this is zero, the camera will use an orthographic projection. The size of the plane
#' used to create the orthographic projection is given in argument `ortho_dimensions`.
#' @param samples Default `100`. The maximum number of samples for each pixel. If this is a length-2
#' vector and the `sample_method` is `stratified`, this will control the number of strata in each dimension.
#' The total number of samples in this case will be the product of the two numbers.
#' @param min_variance Default `0.00005`. Minimum acceptable variance for a block of pixels for the 
#' adaptive sampler. Smaller numbers give higher quality images, at the expense of longer rendering times.
#' If this is set to zero, the adaptive sampler will be turned off and the renderer
#' will use the maximum number of samples everywhere.
#' @param min_adaptive_size Default `8`. Width of the minimum block size in the adaptive sampler.
#' @param sample_method Default `random`. The type of sampling method used to generate
#' random numbers. The other option is `stratified`, which can improve the render quality (at the cost
#' of increased time allocating the random samples).
#' @param max_depth Default `50`. Maximum number of bounces a ray can make in a scene.
#' @param roulette_active_depth Default `10`. Number of ray bounces until a ray can stop bouncing via
#' Russian roulette.
#' @param ambient_light Default `FALSE`, unless there are no emitting objects in the scene. 
#' If `TRUE`, the background will be a gradient varying from `backgroundhigh` directly up (+y) to 
#' `backgroundlow` directly down (-y).
#' @param lookfrom Default `c(0,1,10)`. Location of the camera.
#' @param lookat Default `c(0,0,0)`. Location where the camera is pointed.
#' @param camera_up Default `c(0,1,0)`. Vector indicating the "up" position of the camera.
#' @param aperture Default `0.1`. Aperture of the camera. Smaller numbers will increase depth of field, causing
#' less blurring in areas not in focus.
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
#' showing the number of samples needed to take for each block to converge (when the 
#' @param parallel Default `FALSE`. If `TRUE`, it will use all available cores to render the image
#'  (or the number specified in `options("cores")` if that option is not `NULL`).
#' @param progress Default `TRUE` if interactive session, `FALSE` otherwise. 
#' @param verbose Default `FALSE`. Prints information and timing information about scene
#' construction and raytracing progress.
#' @export
#' @importFrom  grDevices col2rgb
#' @return Raytraced plot to current device, or an image saved to a file. 
#'
#' @examples
#' #Generate a large checkered sphere as the ground
#' \donttest{
#' scene = generate_ground(depth=-0.5, material = diffuse(color="white", checkercolor="darkgreen"))
#' render_scene(scene,parallel=TRUE,samples=500)
#' }
#' 
#' #Add a sphere to the center
#' \donttest{
#' scene = scene %>%
#'   add_object(sphere(x=0,y=0,z=0,radius=0.5,material = diffuse(color=c(1,0,1))))
#' render_scene(scene,fov=20,parallel=TRUE,samples=500)
#' }
#' 
#' #Add a marbled cube 
#' \donttest{
#' scene = scene %>%
#'   add_object(cube(x=1.1,y=0,z=0,material = diffuse(noise=3)))
#' render_scene(scene,fov=20,parallel=TRUE,samples=500)
#' }
#' 
#' #Add a metallic gold sphere, using stratified sampling for a higher quality render
#' \donttest{
#' scene = scene %>%
#'   add_object(sphere(x=-1.1,y=0,z=0,radius=0.5,material = metal(color="gold",fuzz=0.1)))
#' render_scene(scene,fov=20,parallel=TRUE,samples=500, sample_method = "stratified")
#' }
#' 
#' #Lower the number of samples to render more quickly (here, we also use only one core).
#' \donttest{
#' render_scene(scene, samples=4)
#' }
#' 
#' #Add a floating R plot using the iris dataset as a png onto a floating 2D rectangle
#' 
#' \donttest{
#' tempfileplot = tempfile()
#' png(filename=tempfileplot,height=400,width=800)
#' plot(iris$Petal.Length,iris$Sepal.Width,col=iris$Species,pch=18,cex=4)
#' dev.off()
#' 
#' image_array = aperm(png::readPNG(tempfileplot),c(2,1,3))
#' scene = scene %>%
#'   add_object(xy_rect(x=0,y=1.1,z=0,xwidth=2,angle = c(0,180,0),
#'                      material = diffuse(image = image_array)))
#' render_scene(scene,fov=20,parallel=TRUE,samples=500)
#' }
#' 
#' #Move the camera
#' \donttest{
#' render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15,parallel=TRUE)
#' }
#' 
#' #Change the background gradient to a night time ambiance
#' \donttest{
#' render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15,
#'              backgroundhigh = "#282375", backgroundlow = "#7e77ea", parallel=TRUE,
#'              samples=500)
#' }
#'                  
#'#Increase the aperture to blur objects that are further from the focal plane.
#' \donttest{
#' render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15,
#'              aperture = 0.5,parallel=TRUE,samples=500)
#' }
#'                  
#'#Spin the camera around the scene, decreasing the number of samples to render faster. To make 
#'#an animation, specify the a filename in `render_scene` for each frame and use the `av` package
#'#or ffmpeg to combine them all into a movie.
#'
#'t=1:30 
#'xpos = 10 * sin(t*12*pi/180+pi/2)
#'zpos = 10 * cos(t*12*pi/180+pi/2)
#'\donttest{
#'#Save old par() settings
#'old.par = par(no.readonly = TRUE)
#'on.exit(par(old.par))
#'par(mfrow=c(5,6))
#'for(i in 1:30) {
#'  render_scene(scene, samples=16, 
#'    lookfrom = c(xpos[i],1.5,zpos[i]),lookat = c(0,0.5,0), parallel=TRUE)
#'}
#'}
render_scene = function(scene, width = 400, height = 400, fov = 20, 
                        samples = 100, min_variance = 0.00005, min_adaptive_size = 8,
                        sample_method = "random",
                        max_depth = 50, roulette_active_depth = 10,
                        ambient_light = FALSE, 
                        lookfrom = c(0,1,10), lookat = c(0,0,0), camera_up = c(0,1,0), 
                        aperture = 0.1, clamp_value = Inf,
                        filename = NULL, backgroundhigh = "#80b4ff",backgroundlow = "#ffffff",
                        shutteropen = 0.0, shutterclose = 1.0, focal_distance=NULL, ortho_dimensions = c(1,1),
                        tonemap ="gamma", bloom = TRUE, parallel=TRUE, 
                        environment_light = NULL, rotate_env = 0, intensity_env = 1,
                        debug_channel = "none",
                        progress = interactive(), verbose = FALSE) { 
  if(verbose) {
    currenttime = proc.time()
    cat("Building Scene: ")
  }
  #Check if Cornell Box scene and set camera if user did not:
  if(!is.null(attr(scene,"cornell"))) {
    corn_message = "Setting default values for Cornell box: "
    missing_corn = FALSE
    if(missing(lookfrom)) {
      lookfrom = c(278, 278, -800)
      corn_message = paste0(corn_message, "lookfrom `c(278,278,-800)` ")
      missing_corn = TRUE
    }
    if(missing(lookat)) {
      lookat = c(278, 278, 0)
      corn_message = paste0(corn_message, "lookat `c(278,278,0)` ")
      missing_corn = TRUE
    }
    if(missing(fov)) {
      fov=40
      corn_message = paste0(corn_message, "fov `40` ")
      missing_corn = TRUE
    }
    if(fov == 0 && missing(ortho_dimensions)) {
      ortho_dimensions = c(580,580)
      corn_message = paste0(corn_message, "ortho_dimensions `c(580, 580)` ")
      missing_corn = TRUE
    }
    corn_message = paste0(corn_message,".")
    if(missing_corn) {
      message(corn_message)
    }
  }
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
                          "objvertexcolor" = 12, "cone" = 13))
  typevec = unlist(lapply(tolower(scene$type),switch,
                          "diffuse" = 1,"metal" = 2,"dielectric" = 3, "oren-nayar" = 4, "light" = 5, "microfacet" = 6, "glossy" = 7))
  sigmavec = unlist(scene$sigma)
  
  assertthat::assert_that(tonemap %in% c("gamma","reinhold","uncharted", "hbd", "raw"))
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
  
  if(!any(typevec == 5) && missing(ambient_light) && missing(environment_light)) {
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
    stop(paste0("Cannot find the following .obj files:\n",
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
  
  assertthat::assert_that(all(c(length(position_list$xvec),length(position_list$yvec),length(position_list$zvec),length(rvec),length(typevec),length(proplist)) == length(position_list$xvec)))
  assertthat::assert_that(all(!is.null(typevec)))
  assertthat::assert_that(length(lookfrom) == 3)
  assertthat::assert_that(length(lookat) == 3)
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

  debug_channel = unlist(lapply(tolower(debug_channel),switch,
                          "none" = 0,"depth" = 1,"normals" = 2, "uv" = 3, "bvh" = 4,
                          "variance" = 5))
  if(debug_channel == 4) {
    message("rayrender must be compiled with option DEBUGBVH for this debug option to work")
  }
  
  if(fov == 0) {
    assertthat::assert_that(length(ortho_dimensions) == 2)
  }
  if(verbose) {
    buildingtime = proc.time() - currenttime
    cat(sprintf("%0.3f seconds \n",buildingtime[3]))
  }
  sample_method = unlist(lapply(tolower(sample_method),switch,
                                "random" = 0,"stratified" = 1, 0))
  
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
  
  assertthat::assert_that(max_depth > 0)
  assertthat::assert_that(roulette_active_depth > 0)
  
  
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
  
  rgb_mat = render_scene_rcpp(camera_info = camera_info, ambient_light = ambient_light,
                             type = typevec, shape = shapevec, radius = rvec, 
                             position_list = position_list,
                             properties = proplist, velocity = vel_list, moving = movingvec,
                             n = length(typevec), 
                             bghigh = backgroundhigh, bglow = backgroundlow,
                             ischeckered = checkeredbool, checkercolors = checkeredlist,
                             gradient_info = gradient_info, 
                             noise=noisevec,isnoise=noisebool,noisephase=noisephasevec, 
                             noiseintensity=noiseintvec, noisecolorlist = noisecolorlist,
                             angle = rot_angle_list, isimage = image_tex_bool, filelocation = temp_file_names,
                             alphalist = alphalist,
                             lightintensity = light_prop_vec,isflipped = flip_vec,
                             isvolume=fog_bool, voldensity = fog_vec, 
                             implicit_sample = implicit_vec, order_rotation_list = order_rotation_list, clampval = clamp_value,
                             isgrouped = group_bool, group_pivot=group_pivot, group_translate = group_translate,
                             group_angle = group_angle, group_order_rotation = group_order_rotation, group_scale = group_scale,
                             tri_normal_bools = tri_normal_bools, is_tri_color = is_tri_color, tri_color_vert= tri_color_vert,
                             fileinfo = objfilenamevec, filebasedir = objbasedirvec,
                             progress_bar = progress, numbercores = numbercores, 
                             hasbackground = hasbackground, background = backgroundstring, scale_list = scale_factor,
                             sigmavec = sigmavec, rotate_env = rotate_env, intensity_env = intensity_env,
                             verbose = verbose, debug_channel = debug_channel,
                             shared_id_mat=material_id, is_shared_mat=material_id_bool, 
                             min_variance = min_variance, min_adaptive_size = min_adaptive_size, 
                             glossyinfo = glossyinfo, image_repeat = image_repeat) 
  full_array = array(0,c(ncol(rgb_mat$r),nrow(rgb_mat$r),3))
  full_array[,,1] = flipud(t(rgb_mat$r))
  full_array[,,2] = flipud(t(rgb_mat$g))
  full_array[,,3] = flipud(t(rgb_mat$b))
  if(debug_channel == 1) {
    returnmat = full_array[,,1]
    returnmat[is.infinite(returnmat)] = NA
    if(is.null(filename)) {
      plot_map(full_array)
      return(invisible(full_array))
    } else {
      save_png(full_array,filename)
      return(invisible(full_array))
    }
  } else if (debug_channel %in% c(2,3,4,5)) {
    if(is.null(filename)) {
      plot_map(full_array)
      return(invisible(full_array))
    } else {
      save_png(full_array,filename)
      return(invisible(full_array))
    }
  } 
  if(!is.matrix(bloom)) {
    if(is.numeric(bloom) && length(bloom) == 1) {
      kernel = rayimage::generate_2d_exponential(0.1,11,3*1/bloom)
      full_array = rayimage::render_convolution(image = full_array, kernel = kernel, min_value = 1, preview=FALSE)
    } else {
      if(bloom) {
        kernel = rayimage::generate_2d_exponential(0.1,11,3)
        full_array = rayimage::render_convolution(image = full_array, kernel = kernel, min_value = 1, preview=FALSE)
      }
    }
  } else {
    kernel = bloom
    if(ncol(kernel) %% 2 == 0) {
      newkernel = matrix(0, ncol = ncol(kernel) + 1, nrow = nrow(kernel))
      newkernel[,1:ncol(kernel)] = kernel
      kernel = newkernel
    }
    if(nrow(kernel) %% 2 == 0) {
      newkernel = matrix(0, ncol = ncol(kernel), nrow = nrow(kernel) + 1)
      newkernel[1:nrow(kernel),] = kernel
      kernel = newkernel
    }
    full_array = rayimage::render_convolution(image = full_array, kernel = kernel,  min_value = 1, preview=FALSE)
  }
  tonemapped_channels = tonemap_image(height,width,full_array[,,1],full_array[,,2],full_array[,,3],toneval)
  full_array = array(0,c(nrow(tonemapped_channels$r),ncol(tonemapped_channels$r),3))
  full_array[,,1] = tonemapped_channels$r
  full_array[,,2] = tonemapped_channels$g
  full_array[,,3] = tonemapped_channels$b
  if(toneval == 5) {
    return(full_array)
  }

  array_from_mat = array(full_array,dim=c(nrow(full_array),ncol(full_array),3))
  if(any(is.na(array_from_mat ))) {
    array_from_mat[is.na(array_from_mat)] = 0
  }
  if(any(array_from_mat > 1 | array_from_mat < 0,na.rm = TRUE)) {
    array_from_mat[array_from_mat > 1] = 1
    array_from_mat[array_from_mat < 0] = 0
  }
  if(is.null(filename)) {
    plot_map(array_from_mat)
  } else {
    save_png(array_from_mat,filename)
  }
  return(invisible(array_from_mat))
}
