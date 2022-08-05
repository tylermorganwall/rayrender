#' Render Ambient Occlusion
#' 
#' Takes the scene description and renders an image using ambient occlusion, 
#' either to the device or to a filename. 
#'
#' @param scene Tibble of object locations and properties. 
#' @param width Default `400`. Width of the render, in pixels.
#' @param height Default `400`. Height of the render, in pixels.
#' @param fov Default `20`. Field of view, in degrees. If this is `0`, the camera will use an orthographic projection. The size of the plane
#' used to create the orthographic projection is given in argument `ortho_dimensions`. From `0` to `180`, this uses a perspective
#' projections. If this value is `360`, a 360 degree environment image will be rendered. 
#' @param sample_dist Default `10`. Ambient occlusion sampling distance.
#' @param keep_colors Default `FALSE`. Whether to keep the diffuse material colors.
#' @param background_color Default `"white"`. Background color.
#' @param samples Default `100`. The maximum number of samples for each pixel. If this is a length-2
#' vector and the `sample_method` is `stratified`, this will control the number of strata in each dimension.
#' The total number of samples in this case will be the product of the two numbers.
#' @param camera_description_file Default `NA`. Filename of a camera description file for rendering with
#' a realistic camera. Several camera files are built-in: `"50mm"`,`"wide"`,`"fisheye"`, and `"telephoto"`.
#' @param camera_scale Default `1`. Amount to scale the camera up or down in size. Use this rather than scaling a 
#' scene.
#' @param iso Default `100`. Camera exposure.
#' @param film_size Default `22`, in `mm` (scene units in `m`. Size of the film if using a realistic camera, otherwise
#' ignored.
#' @param min_variance Default `0.00005`. Minimum acceptable variance for a block of pixels for the 
#' adaptive sampler. Smaller numbers give higher quality images, at the expense of longer rendering times.
#' If this is set to zero, the adaptive sampler will be turned off and the renderer
#' will use the maximum number of samples everywhere.
#' @param min_adaptive_size Default `8`. Width of the minimum block size in the adaptive sampler.
#' @param sample_method Default `sobol`. The type of sampling method used to generate
#' random numbers. The other options are `random` (worst quality but fastest), 
#' `stratified` (only implemented for completion), `sobol_blue` (best option for sample counts below 256), 
#' and `sobol` (slowest but best quality, better than `sobol_blue` for sample counts greater than 256).
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
#' @param shutteropen Default `0`. Time at which the shutter is open. Only affects moving objects.
#' @param shutterclose Default `1`. Time at which the shutter is open. Only affects moving objects.
#' @param focal_distance Default `NULL`, automatically set to the `lookfrom-lookat` distance unless
#' otherwise specified.
#' @param ortho_dimensions Default `c(1,1)`. Width and height of the orthographic camera. Will only be used if `fov = 0`. 
#' @param parallel Default `FALSE`. If `TRUE`, it will use all available cores to render the image
#'  (or the number specified in `options("cores")` if that option is not `NULL`).
#' @param bvh_type Default `"sah"`, "surface area heuristic". Method of building the bounding volume
#' hierarchy structure used when rendering. Other option is "equal", which splits tree into groups
#' of equal size.
#' @param progress Default `TRUE` if interactive session, `FALSE` otherwise. 
#' @param verbose Default `FALSE`. Prints information and timing information about scene
#' construction and raytracing progress.
#' @export
#' @importFrom  grDevices col2rgb
#' @return Raytraced plot to current device, or an image saved to a file. Invisibly returns the
#' array (containing either debug data or the RGB)
#'
#' @examples
#' #Generate and render a regular scene and an ambient occlusion version of that scene
#' if(rayrender:::run_documentation()) {
#'angles = seq(0,360,by=36)
#'xx = rev(c(rep(c(1,0.5),5),1) * sinpi(angles/180))
#'yy = rev(c(rep(c(1,0.5),5),1) * cospi(angles/180))
#'star_polygon = data.frame(x=xx,y=yy)
#'hollow_star = rbind(star_polygon,0.8*star_polygon)
#'
#'generate_ground(material = diffuse(color="grey20", checkercolor = "grey50",sigma=90)) %>%
#'  add_object(sphere(material=metal())) %>%
#'  add_object(obj_model(y=-1,x=-1.8,r_obj(), angle=c(0,135,0),material = diffuse(sigma=90))) %>%
#'  add_object(pig(x=1.8,y=-1.2,scale=0.5,angle=c(0,90,0),diffuse_sigma = 90)) %>%
#'  add_object(extruded_polygon(hollow_star,top=-0.5,bottom=-1, z=-2,
#'                              hole = nrow(star_polygon),
#'                              material=diffuse(color="red",sigma=90))) %>%
#'  render_scene(parallel = TRUE,width=800,height=800,
#'               fov=70,clamp_value=10,samples=128, aperture=0.1,
#'               lookfrom=c(-0.9,1.2,-4.5),lookat=c(0,-1,0))
#'}
#' if(rayrender:::run_documentation()) {
#'
#'#Render the scene with ambient occlusion
#'generate_ground(material = diffuse(color="grey20", checkercolor = "grey50",sigma=90)) %>%
#'  add_object(sphere(material=metal())) %>%
#'  add_object(obj_model(y=-1,x=-1.8,r_obj(), angle=c(0,135,0),material = diffuse(sigma=90))) %>%
#'  add_object(pig(x=1.8,y=-1.2,scale=0.5,angle=c(0,90,0),diffuse_sigma = 90)) %>%
#'  add_object(extruded_polygon(hollow_star,top=-0.5,bottom=-1, z=-2,
#'                              hole = nrow(star_polygon),
#'                              material=diffuse(color="red",sigma=90))) %>%
#'  render_ao(parallel = TRUE,width=800,height=800, sample_dist=10,
#'            fov=70,samples=128, aperture=0.1,
#'            lookfrom=c(-0.9,1.2,-4.5),lookat=c(0,-1,0))
#'  }
#' if(rayrender:::run_documentation()) {
#'#Decrease the ray occlusion search distance
#'generate_ground(material = diffuse(color="grey20", checkercolor = "grey50",sigma=90)) %>%
#'  add_object(sphere(material=metal())) %>%
#'  add_object(obj_model(y=-1,x=-1.8,r_obj(), angle=c(0,135,0),material = diffuse(sigma=90))) %>%
#'  add_object(pig(x=1.8,y=-1.2,scale=0.5,angle=c(0,90,0),diffuse_sigma = 90)) %>%
#'  add_object(extruded_polygon(hollow_star,top=-0.5,bottom=-1, z=-2,
#'                              hole = nrow(star_polygon),
#'                              material=diffuse(color="red",sigma=90))) %>%
#'  render_ao(parallel = TRUE,width=800,height=800, sample_dist=1,
#'            fov=70,samples=128, aperture=0.1,
#'            lookfrom=c(-0.9,1.2,-4.5),lookat=c(0,-1,0))
#' }
#' if(rayrender:::run_documentation()) {
#'#Turn on colors
#'generate_ground(material = diffuse(color="grey20", checkercolor = "grey50",sigma=90)) %>%
#'  add_object(sphere(material=metal())) %>%
#'  add_object(obj_model(y=-1,x=-1.8,r_obj(), angle=c(0,135,0),material = diffuse(sigma=90))) %>%
#'  add_object(pig(x=1.8,y=-1.2,scale=0.5,angle=c(0,90,0),diffuse_sigma = 90)) %>%
#'  add_object(extruded_polygon(hollow_star,top=-0.5,bottom=-1, z=-2,
#'                              hole = nrow(star_polygon),
#'                              material=diffuse(color="red",sigma=90))) %>%
#'  render_ao(parallel = TRUE,width=800,height=800, sample_dist=1,
#'            fov=70,samples=128, aperture=0.1, keep_colors = TRUE,
#'            lookfrom=c(-0.9,1.2,-4.5),lookat=c(0,-1,0))
#'
#' }
render_ao = function(scene, width = 400, height = 400, fov = 20, 
                     sample_dist = 10, keep_colors = FALSE,
                     samples = 100,  camera_description_file = NA, 
                     camera_scale = 1, iso = 100, film_size = 22,
                     min_variance = 0.0000, min_adaptive_size = 8,
                     sample_method = "sobol", background_color = "white",
                     lookfrom = c(0,1,10), lookat = c(0,0,0), camera_up = c(0,1,0), 
                     aperture = 0.1, clamp_value = Inf,
                     filename = NULL, 
                     shutteropen = 0.0, shutterclose = 1.0, focal_distance=NULL, ortho_dimensions = c(1,1),
                     parallel=TRUE, bvh_type = "sah", 
                     progress = interactive(), verbose = FALSE) { 
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
      lookat = c(278, 278, 555/2)
      corn_message = paste0(corn_message, "lookat `c(278,278,555/2)` ")
      missing_corn = TRUE
    }
    if(missing(fov) && is.na(camera_description_file)) {
      fov=40
      corn_message = paste0(corn_message, "fov `40` ")
      missing_corn = TRUE
    }
    if(fov == 0 && missing(ortho_dimensions) && is.na(camera_description_file)) {
      ortho_dimensions = c(580,580)
      corn_message = paste0(corn_message, "ortho_dimensions `c(580, 580)` ")
      missing_corn = TRUE
    }
    corn_message = paste0(corn_message,".")
    if(missing_corn) {
      message(corn_message)
    }
  }
  debug_channel = "ao"
  tonemap = "gamma"
  bloom = FALSE
  
  background_color = convert_color(background_color)
  scene_list = prepare_scene_list(scene = scene, width = width, height = height, fov = fov, 
                                  samples = samples,  camera_description_file = camera_description_file, 
                                  camera_scale = camera_scale, iso = iso, film_size = film_size,
                                  min_variance = min_variance, min_adaptive_size = min_adaptive_size,
                                  sample_method = sample_method, 
                                  max_depth = 1000, roulette_active_depth = 1000,
                                  ambient_light = FALSE, 
                                  lookfrom = lookfrom, lookat = lookat, camera_up = camera_up, 
                                  aperture = aperture, clamp_value = clamp_value,
                                  filename = filename, backgroundhigh = background_color, backgroundlow = background_color,
                                  shutteropen = shutteropen, shutterclose = shutterclose, 
                                  focal_distance = focal_distance, ortho_dimensions = ortho_dimensions,
                                  tonemap = "gamma", bloom = FALSE, parallel=parallel, bvh_type = bvh_type,
                                  environment_light = NULL, rotate_env = 0, 
                                  intensity_env = 1,
                                  debug_channel = debug_channel, return_raw_array = FALSE,
                                  progress = progress, verbose = verbose, sample_dist = sample_dist,
                                  keep_colors = keep_colors)
  
  camera_info = scene_list$camera_info
  scene_info = scene_list$scene_info
  
  camera_info$preview = FALSE
  camera_info$interactive = FALSE
  
  #Pathrace Scene
  rgb_mat = render_scene_rcpp(camera_info = camera_info, scene_info = scene_info) 
  
  return_array = post_process_scene(rgb_mat, iso, tonemap, debug_channel, filename, FALSE, bloom)
  return(invisible(return_array))
}
