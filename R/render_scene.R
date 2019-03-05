#' Title
#'
#' @param scene Tibble of sphere locations and properties. Generate each row using one of the
#' material functions: [lambertian()], [metal()], and [dielectric()].
#' @param width Default `200`. Width of the render, in pixels.
#' @param height Default `200`. Height of the render, in pixels.
#' @param fov Default `20`. Field of view, in degrees.
#' @param samples Default `100`. Number of samples for each pixel.
#' @param lookfrom Default `c(10,1,0)`. Location of the camera.
#' @param lookat Default `c(0,0,0)`. Location where the camera is pointed.
#' @param aperture Default `0.1`. Aperture of the camera. Higher numbers will increase depth of field.
#' @param filename Default `NULL`. If present, the renderer will write to the filename instead
#' of the current device.
#' @param backgroundhigh Default `#80b4ff`. The "high" color in the background gradient. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param backgroundlow Default `#ffffff`. The "low" color in the background gradient. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param shutteropen Default `0`. Time at which the shutter is open. Only affects moving objects.
#' @param shutterclose Default `1`. Time at which the shutter is open. Only affects moving objects.
#' @export
#' @importFrom  grDevices col2rgb
#' @return Raytraced plot to current device, or an image saved to a file.
#'
#' @examples
#' #Start with the ground
#' scene = generate_ground(depth=-0.5)
#' render_scene(scene)
#' 
#' #Add a sphere to the center
#' scene = add_sphere(scene, lambertian(x=0,y=0,z=0,radius=0.5,color=c(1,0,1)))
#' render_scene(scene)
#' 
#' #Add a metal ball (using hexcode color representation)
#' scene = add_sphere(scene, metal(x=0,y=0,z=1,radius=0.5,color="#ffffff",fuzz=0))
#' render_scene(scene)
#' 
#' #Add a brushed metal ball 
#' scene = add_sphere(scene, metal(x=0,y=1,z=0,radius=0.5,color=c(0.3,0.6,1),fuzz=0.25))
#' render_scene(scene)
#' 
#' #Add a dielectric (glass) ball
#' scene = add_sphere(scene, dielectric(x=0,y=0,z=-1,radius=0.5,refraction=1.6))
#' render_scene(scene)
#' 
#' #Add a grid of glass balls in front
#' glass_array_list = list() 
#' 
#' yloc = seq(-0.33,+0.33,length.out=3)
#' zloc = seq(-0.33,0.33,length.out=3)
#' locations = expand.grid(y=yloc,z=zloc)
#' for(i in 1:9) {
#'   glass_array_list[[i]] = dielectric(x=1,y=locations$y[i],z=locations$z[i], radius=0.15)
#' }
#' glass_array = do.call(rbind,glass_array_list)
#' 
#' scene = add_sphere(scene, glass_array)
#' render_scene(scene)
#' 
#' #Move the camera
#' render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15)
#' 
#' #Change the background gradient to a night time ambience
#' render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15,
#'                  backgroundhigh = "#282375", backgroundlow = "#7e77ea")
#'                  
#'#Increase the aperture to give more depth of field.
#' render_scene(scene,lookfrom = c(7,1.5,10),lookat = c(0,0.5,0),fov=15,
#'                  aperture = 1)
render_scene = function(scene, width = 400, height = 400, fov = 20, samples = 100,
                        lookfrom = c(10,1,0), lookat = c(0,0,0), aperture = 0.1,
                        filename = NULL, backgroundhigh = "#80b4ff",backgroundlow = "#ffffff",
                        shutteropen = 0.0, shutterclose = 1.0) { 
  backgroundhigh = convert_color(backgroundhigh)
  backgroundlow = convert_color(backgroundlow)
  xvec = scene$x 
  yvec = scene$y
  zvec = scene$z
  rvec = scene$radius
  typevec = unlist(lapply(tolower(scene$type),switch,"lambertian" = 1,"metal" = 2,"dielectric" = 3))
  movingvec = purrr::map_lgl(scene$velocity,.f = ~any(.x != 0))
  proplist = scene$properties
  vel_list = scene$velocity
  checkeredlist = scene$checkercolor
  checkeredbool = purrr::map_lgl(checkeredlist,.f = ~all(!is.na(.x)))
  noisebool = purrr::map_lgl(scene$noise, .f = ~.x > 0)
  noisevec = scene$noise
  if(shutteropen == shutterclose) {
    movingvec = rep(FALSE,length(movingvec))
  }
  assertthat::assert_that(all(c(length(xvec),length(yvec),length(zvec),length(rvec),length(typevec),length(proplist)) == length(xvec)))
  assertthat::assert_that(all(!is.null(typevec)))
  for(i in 1:length(xvec)) {
    if(typevec[i] == 1) {
      assertthat::assert_that(length(proplist[[i]]) == 3)
    } else if (typevec[i] == 2) {
      assertthat::assert_that(length(proplist[[i]]) == 4)
    } else if (typevec[i] == 3) {
      assertthat::assert_that(length(proplist[[i]]) == 1)
    }
  }
  rgb_mat = generate_initial(nx = width, ny = height, ns = samples, fov = fov,
                             lookfromvec = lookfrom, lookatvec = lookat, aperture=aperture,
                             type = typevec, radius = rvec,
                             x = xvec, y = yvec, z = zvec,
                             properties = proplist, velocity = vel_list, moving = movingvec,
                             n = length(typevec), 
                             bghigh = backgroundhigh, bglow = backgroundlow,
                             shutteropen = shutteropen, shutterclose = shutterclose,
                             ischeckered = checkeredbool, checkercolors = checkeredlist,
                             noise=noisevec,isnoise=noisebool) 
  full_array = array(0,c(ncol(rgb_mat$r),nrow(rgb_mat$r),3))
  full_array[,,1] = t(rgb_mat$r)
  full_array[,,2] = t(rgb_mat$g)
  full_array[,,3] = t(rgb_mat$b)
  array_from_mat = array(full_array,dim=c(nrow(full_array),ncol(full_array),3))
  if(is.null(filename)) {
    rayshader::plot_map(flipud(array_from_mat))
  } else {
    rayshader::save_png(flipud(array_from_mat),filename)
  }
}