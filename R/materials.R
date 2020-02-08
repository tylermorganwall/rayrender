#' Diffuse Material
#'
#' @param color Default `white`. The color of the surface. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param checkercolor Default `NA`. If not `NA`, determines the secondary color of the checkered surface. 
#' Can be either a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param checkerperiod Default `3`. The period of the checker pattern. Increasing this value makes the checker 
#' pattern bigger, and decreasing it makes it smaller
#' @param noise Default `0`. If not `0`, covers the surface in a turbulent marble pattern. This value will determine
#' the amount of turbulence in the texture.
#' @param noisephase Default `0`. The phase of the noise. The noise will repeat at `360`.
#' @param noiseintensity Default `10`. Intensity of the noise.
#' @param noisecolor Default `#000000`. The secondary color of the noise pattern.
#' Can be either a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param gradient_color Default `NA`. If not `NA`, creates a secondary color for a linear gradient 
#' between the this color and color specified in `color`. Direction is determined by `gradient_transpose`.
#' @param gradient_transpose Default `FALSE`. If `TRUE`, this will use the `v` coordinate texture instead
#' of the `u` coordinate texture to map the gradient.
#' @param image_array A 3-layer RGB array to be used as the texture on the surface of the object.
#' @param fog Default `FALSE`. If `TRUE`, the object will be a volumetric scatterer.
#' @param fogdensity Default `0.01`. The density of the fog. Higher values will produce more opaque objects.
#' @param sigma Default `NULL`. A number between 0 and Infinity specifying the roughness of the surface using the Oren-Nayar microfacet model.
#' Higher numbers indicate a roughed surface, where sigma is the standard deviation of the microfacet orientation angle. When 0, this reverts
#' to the default lambertian behavior.
#' @param importance_sample Default `FALSE`. If `TRUE`, the object will be sampled explicitly during 
#' the rendering process. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#'
#' @return Single row of a tibble describing the diffuse material.
#' @export
#' @importFrom  grDevices col2rgb
#'
#' @examples
#' #Generate the cornell box and add a single white sphere to the center
#' scene = generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,material=diffuse()))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Add a checkered rectangular cube below             
#' scene = scene %>%
#'   add_object(cube(x=555/2,y=555/8,z=555/2,xwidth=555/2,ywidth=555/4,zwidth=555/2,
#'   material = diffuse(checkercolor="purple",checkerperiod=20)))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#'   
#' #Add a marbled sphere           
#' scene = scene %>%
#'   add_object(sphere(x=555/2+555/4,y=555/2,z=555/2,radius=555/8,
#'   material = diffuse(noise=1/20)))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Add an orange volumetric (fog) cube           
#' scene = scene %>%
#'   add_object(cube(x=555/2-555/4,y=555/2,z=555/2,xwidth=555/4,ywidth=555/4,zwidth=555/4,
#'   material = diffuse(fog=TRUE, fogdensity=0.05,color="orange")))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #' #Add an line segment with a color gradient        
#' scene = scene %>%
#'   add_object(segment(start = c(555,450,450),end=c(0,450,450),radius = 50, 
#'                      material = diffuse(color="#1f7326", gradient_color = "#a60d0d")))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
diffuse = function(color = "#ffffff", 
                   checkercolor = NA, checkerperiod = 3,
                   noise = 0, noisephase = 0, noiseintensity = 10, noisecolor = "#000000",
                   gradient_color = NA, gradient_transpose = FALSE,
                   image_array = NA, fog = FALSE, fogdensity = 0.01, 
                   sigma = NULL, importance_sample = FALSE) {
  if(all(!is.na(checkercolor))) {
    checkercolor = convert_color(checkercolor)
  } else {
    checkercolor = NA
  }
  if(all(!is.na(gradient_color))) {
    gradient_color = convert_color(gradient_color)
  } else {
    gradient_color = NA
  }
  
  info = convert_color(color)
  noisecolor = convert_color(noisecolor)
  if(!is.array(image_array) && !is.na(image_array)) {
    image = NA
    warning("Image not in recognized format (array or matrix), ignoring")
  }
  type = "diffuse"
  if(!is.null(sigma) && is.numeric(sigma)) {
    if(sigma < 0) {
      warning("sigma must be greater than 0 (input: ", sigma, ")--ignoring and using lambertian model")
    } else {
      if(sigma == 0) {
        type = "diffuse"
      } else {
        type = "oren-nayar"
        sigma = sigma*pi/180
      }
    }
  } else {
    sigma = 0
  }
  assertthat::assert_that(checkerperiod != 0)
  tibble::tibble(type = type, 
                 properties = list(info), checkercolor=list(c(checkercolor,checkerperiod)), 
                 gradient_color = list(gradient_color), gradient_transpose = gradient_transpose,
                 noise=noise, noisephase = noisephase, noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                 image = list(image_array), lightintensity = NA,
                 fog=fog, fogdensity=fogdensity,implicit_sample = importance_sample, sigma = sigma)
}

#' Metallic Material
#'
#' @param color Default `white`. The color of the sphere. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param fuzz  Default `0`. The roughness of the metallic surface. Maximum `1`.
#' @param importance_sample Default `FALSE`. If `TRUE`, the object will be sampled explicitly during 
#' the rendering process. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the metallic material.
#' @export
#'
#' @examples
#' #Generate the cornell box with a single metal sphere in the center
#' scene = generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,material=metal()))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' #Add a rotated shiny metal cube     
#' scene = scene %>%
#'   add_object(cube(x=380,y=150/2,z=200,xwidth=150,ywidth=150,zwidth=150,
#'   material = metal(color="#8B4513"),angle=c(0,45,0)))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' #Add a brushed metal cube (setting the fuzz variable)           
#' scene = scene %>%
#'   add_object(cube(x=150,y=150/2,z=300,xwidth=150,ywidth=150,zwidth=150,
#'   material = metal(color="#FAFAD2",fuzz=0.1),angle=c(0,-30,0)))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
metal = function(color = "#ffffff", fuzz = 0,  importance_sample = FALSE) {
  color = convert_color(color)
  tibble::tibble(type = "metal", 
                 properties = list(c(color,fuzz)), 
                 checkercolor=list(NA), 
                 gradient_color = list(NA), gradient_transpose = FALSE,
                 noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 lightinfo = list(NA),
                 image = list(NA), lightintensity = NA,fog=FALSE,fogdensity=0.01,
                 implicit_sample = importance_sample, sigma = 0)
}

#' Dielectric (glass) Material
#'
#' @param color Default `white`. The color of the surface. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param refraction Default `1.5`. The index of refraction.
#' @param importance_sample Default `FALSE`. If `TRUE`, the object will be sampled explicitly during 
#' the rendering process. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#'
#' @return Single row of a tibble describing the dielectric material.
#' @export
#'
#' @examples
#' #Generate a checkered ground
#' scene = generate_ground(depth=-0.5,
#'                         material=diffuse(color="white", checkercolor="grey30",checkerperiod=2))
#' \donttest{
#' render_scene(scene,parallel=TRUE)
#' }
#' 
#' #Add a glass sphere
#' \donttest{
#' scene %>%
#'   add_object(sphere(x=-0.5,radius=0.5,material=dielectric())) %>%
#'   render_scene(parallel=TRUE,samples=400)
#' }
#' 
#' #Add a rotated colored glass cube
#' \donttest{
#' scene %>%
#'   add_object(sphere(x=-0.5,radius=0.5,material=dielectric())) %>%
#'   add_object(cube(x=0.5,xwidth=0.5,material=dielectric(color="darkgreen"),angle=c(0,-45,0))) %>%
#'   render_scene(parallel=TRUE,samples=40)
#' }
#' 
#' #Add an area light behind and at an angle and turn off the ambient lighting
#' \donttest{
#' scene %>%
#'   add_object(sphere(x=-0.5,radius=0.5,material=dielectric())) %>%
#'   add_object(cube(x=0.5,xwidth=0.5,material=dielectric(color="darkgreen"),angle=c(0,-45,0))) %>%
#'   add_object(yz_rect(z=-3,y=1,x=0,zwidth=3,ywidth=1.5,
#'                      material=light(intensity=15),
#'                      angle=c(0,-90,45), order_rotation = c(3,2,1))) %>%
#'   render_scene(parallel=TRUE,aperture=0, ambient_light=FALSE,samples=1000)
#' }
dielectric = function(color="white", refraction = 1.5, importance_sample = FALSE) {
  color = convert_color(color)
  tibble::tibble(type = "dielectric", 
                 properties = list(c(color,refraction)), 
                 checkercolor=list(NA), 
                 gradient_color = list(NA), gradient_transpose = FALSE,
                 noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 image = list(NA), lightintensity = NA, 
                 fog=FALSE, fogdensity=NA, implicit_sample = importance_sample, sigma = 0)
}

#' Light Material
#'
#' @param color Default `white`. The color of the light Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param intensity Default `NA`. If a positive value, this will turn this object into a light emitting the value specified
#' in `color` (ignoring other properties). Higher values will produce a brighter light.
#' @param importance_sample Default `TRUE`. Keeping this on for lights improves the convergence of the rendering 
#' algorithm, in most cases. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#'
#' @return Single row of a tibble describing the diffuse material.
#' @export
#' @importFrom  grDevices col2rgb
#'
#' @examples
#' #Generate the cornell box without a light and add a single white sphere to the center
#' scene = generate_cornell(light=FALSE) %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,material=light()))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #All gather around the orb
#' scene = generate_ground(material = diffuse(checkercolor="grey50")) %>%
#'   add_object(sphere(radius=0.5,material=light(intensity=5,color="red"))) %>%
#'   add_object(obj_model(r_obj(), z=-3,x=-1.5,y=-1, angle=c(0,45,0))) %>%
#'   add_object(pig(scale=0.3, x=1.5,z=-2,y=-1.5,angle=c(0,-135,0)))
#' \donttest{
#' render_scene(scene, samples=500, parallel=TRUE, clamp_value=10)
#' }
light = function(color = "#ffffff", intensity = 10, importance_sample = TRUE) {
  info = convert_color(color)
  tibble::tibble(type = "light", 
                 properties = list(info), checkercolor=list(NA), 
                 gradient_color = list(NA), gradient_transpose = FALSE,
                 noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 image = list(NA), lightintensity = intensity,
                 fog=FALSE, fogdensity=0.01, implicit_sample = importance_sample, sigma = 0)
}


#' Lambertian Material (deprecated)
#'
#' @param ... Arguments to pass to diffuse() function.
#'
#' @return Single row of a tibble describing the diffuse material.
#' @export
#' @importFrom  grDevices col2rgb
#'
#' @examples
#' #Deprecated lambertian material. Will display a warning.
#' \donttest{
#' scene = generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,material=lambertian()))
#'   render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=10,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
lambertian = function(...) {
  warning("lambertian() deprecated--use diffuse() instead.")
  diffuse(...)
}

