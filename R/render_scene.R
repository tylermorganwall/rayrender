#' Render Scene
#' 
#' Takes the scene description and renders an image, either to the device or to a filename.
#'
#' @param scene Tibble of object locations and properties. 
#' @param width Default `400`. Width of the render, in pixels.
#' @param height Default `400`. Height of the render, in pixels.
#' @param fov Default `20`. Field of view, in degrees. If this is zero, the camera will use an orthographic projection. The size of the plane
#' used to create the orthographic projection is given in argument `ortho_dimensions`.
#' @param samples Default `100`. Number of samples for each pixel.
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
#' objects with material `light()` may not be anti-aliased.
#' @param environment_light Default `NULL`. An image to be used for the background for rays that escape
#' the scene. Supports both HDR (`.hdr`) and low-dynamic range (`.png`, `.jpg`) images.
#' @param rotate_env Default `0`. The number of degrees to rotate the environment map around the scene.
#' @param parallel Default `FALSE`. If `TRUE`, it will use all available cores to render the image
#'  (or the number specified in `options("cores")` if that option is not `NULL`).
#' @param progress Default `TRUE` if interactive session, `FALSE` otherwise. 
#' @param verbose Default `FALSE`. Prints information and timing information about scene
#' construction and raytracing progress.
#' @param debug Default `NULL`. If `bvh`, will return an image indicated the number of BVH lookups.
#' @export
#' @importFrom  grDevices col2rgb
#' @return Raytraced plot to current device, or an image saved to a file.
#'
#' @examples
#' #Generate a large checkered sphere as the ground
#' scene = generate_ground(depth=-0.5, material = diffuse(color="white", checkercolor="darkgreen"))
#' \donttest{
#' render_scene(scene,parallel=TRUE,samples=500)
#' }
#' 
#' #Add a sphere to the center
#' scene = scene %>%
#'   add_object(sphere(x=0,y=0,z=0,radius=0.5,material = diffuse(color=c(1,0,1))))
#' \donttest{
#' render_scene(scene,fov=20,parallel=TRUE,samples=500)
#' }
#' 
#' #Add a marbled cube 
#' scene = scene %>%
#'   add_object(cube(x=1.1,y=0,z=0,material = diffuse(noise=3)))
#' \donttest{
#' render_scene(scene,fov=20,parallel=TRUE,samples=500)
#' }
#' 
#' #Add a metallic gold sphere
#' scene = scene %>%
#'   add_object(sphere(x=-1.1,y=0,z=0,radius=0.5,material = metal(color="gold",fuzz=0.1)))
#' \donttest{
#' render_scene(scene,fov=20,parallel=TRUE,samples=500)
#' }
#' 
#' #Lower the number of samples to render more quickly (here, we also use only one core).
#' render_scene(scene, samples=4)
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
#'  render_scene(scene, samples=5,
#'    lookfrom = c(xpos[i],1.5,zpos[i]),lookat = c(0,0.5,0), parallel=TRUE)
#'}
#'}
render_scene = function(scene, width = 400, height = 400, fov = 20, samples = 100, ambient_light = FALSE,
                        lookfrom = c(0,1,10), lookat = c(0,0,0), camera_up = c(0,1,0), 
                        aperture = 0.1, clamp_value = Inf,
                        filename = NULL, backgroundhigh = "#80b4ff",backgroundlow = "#ffffff",
                        shutteropen = 0.0, shutterclose = 1.0, focal_distance=NULL, ortho_dimensions = c(1,1),
                        tonemap ="gamma", parallel=TRUE,
                        environment_light = NULL, rotate_env = 0,
                        progress = interactive(), verbose = FALSE, debug = NULL) { 
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
  xvec = scene$x 
  yvec = scene$y
  zvec = scene$z
  rvec = scene$radius
  shapevec = unlist(lapply(tolower(scene$shape),switch,
                          "sphere" = 1,"xy_rect" = 2, "xz_rect" = 3,"yz_rect" = 4,"box" = 5, "triangle" = 6, 
                          "obj" = 7, "objcolor" = 8, "disk" = 9, "cylinder" = 10, "ellipsoid" = 11))
  typevec = unlist(lapply(tolower(scene$type),switch,
                          "diffuse" = 1,"metal" = 2,"dielectric" = 3, "oren-nayar" = 4, "light" = 5))
  sigmavec = unlist(scene$sigma)
  assertthat::assert_that(tonemap %in% c("gamma","reinhold","uncharted", "hbd"))
  toneval = switch(tonemap, "gamma" = 1,"reinhold" = 2,"uncharted" = 3,"hbd" = 4)
  movingvec = purrr::map_lgl(scene$velocity,.f = ~any(.x != 0))
  proplist = scene$properties
  vel_list = scene$velocity
  
  checkeredlist = scene$checkercolor
  checkeredbool = purrr::map_lgl(checkeredlist,.f = ~all(!is.na(.x)))
  
  #gradient handler
  gradient_info = list()
  gradient_info$gradient_colors = scene$gradient_color
  gradient_info$isgradient = purrr::map_lgl(gradient_info$gradient_colors,.f = ~all(!is.na(.x)))
  gradient_info$gradient_trans = scene$gradient_transpose
  
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
  temp_file_names = purrr::map_chr(image_tex_bool,.f = ~ifelse(.x, tempfile(fileext = ".png"),""))
  for(i in 1:length(image_array_list)) {
    if(image_tex_bool[i]) {
      if(dim(image_array_list[[i]])[3] == 4) {
        png::writePNG(fliplr(aperm(image_array_list[[i]][,,1:3],c(2,1,3))),temp_file_names[i])
      } else if(dim(image_array_list[[i]])[3] == 3){
        png::writePNG(fliplr(aperm(image_array_list[[i]],c(2,1,3))),temp_file_names[i])
      }
    }
  }
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
  } else {
    hasbackground = FALSE
    backgroundstring = ""
  }
  
  #scale handler
  scale_factor = scene$scale_factor
  
  assertthat::assert_that(all(c(length(xvec),length(yvec),length(zvec),length(rvec),length(typevec),length(proplist)) == length(xvec)))
  assertthat::assert_that(all(!is.null(typevec)))
  for(i in 1:length(xvec)) {
    if(typevec[i] == 1) {
      if(shapevec[i] == 1) {
        # assertthat::assert_that(length(proplist[[i]]) == 3)
      }
    } else if (typevec[i] == 2) {
      # assertthat::assert_that(length(proplist[[i]]) == 4)
    } else if (typevec[i] == 3) {
      # assertthat::assert_that(length(proplist[[i]]) == 7)
    } 
  }
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
  if(is.null(debug)) {
    debugval = 0
  } else if(debug == "bvh") {
    debugval = 1
  } else {
    debugval = 0
  }
  if(fov == 0) {
    assertthat::assert_that(length(ortho_dimensions) == 2)
  }
  rgb_mat = render_scene_rcpp(nx = width, ny = height, ns = samples, fov = fov, ambient_light = ambient_light,
                             lookfromvec = lookfrom, lookatvec = lookat, aperture=aperture,camera_up = camera_up,
                             type = typevec, shape = shapevec, radius = rvec, 
                             x = xvec, y = yvec, z = zvec,
                             properties = proplist, velocity = vel_list, moving = movingvec,
                             n = length(typevec), 
                             bghigh = backgroundhigh, bglow = backgroundlow,
                             shutteropen = shutteropen, shutterclose = shutterclose,
                             ischeckered = checkeredbool, checkercolors = checkeredlist,
                             gradient_info = gradient_info,
                             noise=noisevec,isnoise=noisebool,noisephase=noisephasevec, 
                             noiseintensity=noiseintvec, noisecolorlist = noisecolorlist,
                             angle = rot_angle_list, isimage = image_tex_bool, filelocation = temp_file_names,
                             lightintensity = light_prop_vec,isflipped = flip_vec,
                             focus_distance=focal_distance,
                             isvolume=fog_bool, voldensity = fog_vec , parallel=parallel,
                             implicit_sample = implicit_vec, order_rotation_list = order_rotation_list, clampval = clamp_value,
                             isgrouped = group_bool, group_pivot=group_pivot, group_translate = group_translate,
                             group_angle = group_angle, group_order_rotation = group_order_rotation, group_scale = group_scale,
                             tri_normal_bools = tri_normal_bools, is_tri_color = is_tri_color, tri_color_vert= tri_color_vert,
                             fileinfo = objfilenamevec, filebasedir = objbasedirvec, toneval = toneval,
                             progress_bar = progress, numbercores = numbercores, debugval = debugval,
                             hasbackground = hasbackground, background = backgroundstring, scale_list = scale_factor,
                             ortho_dimensions = ortho_dimensions, sigmavec = sigmavec, rotate_env = rotate_env,
                             verbose = verbose) 
  full_array = array(0,c(ncol(rgb_mat$r),nrow(rgb_mat$r),3))
  full_array[,,1] = t(rgb_mat$r)
  full_array[,,2] = t(rgb_mat$g)
  full_array[,,3] = t(rgb_mat$b)
  array_from_mat = array(full_array,dim=c(nrow(full_array),ncol(full_array),3))
  if(any(is.na(array_from_mat ))) {
    array_from_mat[is.na(array_from_mat)] = 0
  }
  if(any(array_from_mat > 1 | array_from_mat < 0,na.rm = TRUE)) {
    array_from_mat[array_from_mat > 1] = 1
    array_from_mat[array_from_mat < 0] = 0
  }
  if(is.null(filename)) {
    plot_map(flipud(array_from_mat))
  } else {
    save_png(flipud(array_from_mat),filename)
  }
}
