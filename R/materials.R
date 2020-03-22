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
#' @param image_texture Default `NA`. A 3-layer RGB array or filename to be used as the texture on the surface of the object.
#' @param alpha_texture Default `NA`. A matrix or filename (specifying a greyscale image) to be used to specify the transparency.
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
                   image_texture = NA, alpha_texture = NA,
                   fog = FALSE, fogdensity = 0.01, 
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
  if(!is.array(image_texture) && !is.na(image_texture) && !is.character(image_texture)) {
    image_texture = NA
    warning("Texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(alpha_texture) && !is.na(alpha_texture) && !is.character(alpha_texture)) {
    alpha_texture = NA
    warning("Alpha texture not in recognized format (array, matrix, or filename), ignoring.")
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
                 image = list(image_texture), alphaimage = list(alpha_texture), lightintensity = NA,
                 fog=fog, fogdensity=fogdensity,implicit_sample = importance_sample, sigma = sigma)
}

#' Metallic Material
#'
#' @param color Default `white`. The color of the sphere. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param fuzz  Default `0`. The roughness of the metallic surface. Maximum `1`.
#' @param alpha_texture Default `NA`. A matrix or filename (specifying a greyscale image) to be used to specify the transparency.
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
metal = function(color = "#ffffff", fuzz = 0,  alpha_texture = NA, 
                 importance_sample = FALSE) {
  color = convert_color(color)
  if(!is.array(alpha_texture) && !is.na(alpha_texture) && !is.character(alpha_texture)) {
    alpha_texture = NA
    warning("Alpha texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  tibble::tibble(type = "metal", 
                 properties = list(c(color,fuzz)), 
                 checkercolor=list(NA), 
                 gradient_color = list(NA), gradient_transpose = FALSE,
                 noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 lightinfo = list(NA), 
                 image = list(NA), alphaimage = list(alpha_texture), 
                 lightintensity = NA,fog=FALSE,fogdensity=0.01,
                 implicit_sample = importance_sample, sigma = 0)
}

#' Dielectric (glass) Material
#' 
#'
#' @param color Default `white`. The color of the surface. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param refraction Default `1.5`. The index of refraction.
#' @param attenuation Default `c(0,0,0)`. The Beer-Lambert color-channel specific exponential attenuation 
#' through the material. Higher numbers will result in less of that color making it through the material.
#' Note: This assumes the object has a closed surface. 
#' @param priority Default `0`. When two dielectric materials overlap, the one with the lower priority value
#' is used for intersection. NOTE: If the camera is placed inside a dielectric object, its priority value
#' will not be taken into account when determining hits to other objects also inside the object.
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
#' scene = generate_ground(depth=-0.5, material = diffuse(checkercolor="grey30",checkerperiod=2))
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
#'   render_scene(parallel=TRUE,samples=400)
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
#' 
#' #Color glass using Beer-Lambert attenuation, which attenuates light on a per-channel
#' #basis as it travels through the material. This effect is what gives some types of glass
#' #a green glow at the edges. We will get this effect by setting a lower attenuation value 
#' #for the `green` (second) channel in the dielectric `attenuation` argument.
#' \donttest{
#' generate_ground(depth=-0.5,material=diffuse(checkercolor="grey30",checkerperiod=2)) %>%
#'   add_object(sphere(z=-5,x=-0.5,y=1,material=light(intensity=10))) %>%
#'   add_object(cube(y=0.3,ywidth=0.1,xwidth=2,zwidth=2,
#'                   material=dielectric(attenuation=c(1.2,0.2,1.2)),angle=c(45,110,0))) %>%
#'   render_scene(parallel=TRUE, samples = 1000)
#' }
#' 
#' #If you have overlapping dielectrics, the `priority` value can help disambiguate what 
#' #object wins. Here, I place a bubble inside a cube by setting a lower priority value and
#' #making the inner sphere have a index of refraction of 1. I also place spheres at the corners.
#' \donttest{
#' generate_ground(depth=-0.51,material=diffuse(checkercolor="grey30",checkerperiod=2)) %>%
#'   add_object(cube(material = dielectric(priority=2, attenuation = c(10,3,10)))) %>%
#'   add_object(sphere(radius=0.49,material = dielectric(priority=1, refraction=1))) %>%
#'   add_object(sphere(radius=0.25,x=0.5,z=-0.5,y=0.5, 
#'                     material = dielectric(priority=0,attenuation = c(10,3,10) ))) %>%
#'   add_object(sphere(radius=0.25,x=-0.5,z=0.5,y=0.5,
#'                     material = dielectric(priority=0,attenuation = c(10,3,10)))) %>%
#'   render_scene(parallel=TRUE, samples = 400,lookfrom=c(5,1,5))
#' }
#' 
#' # We can also use this as a basic Constructive Solid Geometry interface by setting 
#' # the index of refraction equal to empty space, 1. This will subtract out those regions.
#' # Here I make a concave lens by subtracting two spheres from a cube.
#' \donttest{
#' generate_ground(depth=-0.51,material=diffuse(checkercolor="grey30",checkerperiod=2,sigma=90)) %>%
#'   add_object(cube(material = dielectric(attenuation = c(6,6,2),priority=1))) %>%
#'   add_object(sphere(radius=1,x=1.01,
#'                     material = dielectric(priority=0,refraction=1))) %>%
#'   add_object(sphere(radius=1,x=-1.01,
#'                     material = dielectric(priority=0,refraction=1))) %>%
#'   add_object(sphere(y=10,x=3,material=light(intensit=150))) %>%
#'   render_scene(parallel=TRUE, samples = 400,lookfrom=c(5,3,5))
#' }
dielectric = function(color="white", refraction = 1.5,  attenuation = c(0,0,0), 
                      priority = 0, importance_sample = FALSE) {
  color = convert_color(color)
  tibble::tibble(type = "dielectric", 
                 properties = list(c(color, refraction, attenuation, priority)), 
                 checkercolor=list(NA), 
                 gradient_color = list(NA), gradient_transpose = FALSE,
                 noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 image = list(NA), alphaimage = list(NA), lightintensity = NA, 
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
                 image = list(NA), alphaimage = list(NA), lightintensity = intensity,
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

