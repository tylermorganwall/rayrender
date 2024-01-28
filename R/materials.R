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
#' @param gradient_point_start Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `color`.
#' @param gradient_point_end Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `gradient_color`.
#' @param gradient_type Default `hsv`. Colorspace to calculate the gradient. Alternative `rgb`.
#' @param image_texture Default `NA`. A 3-layer RGB array or filename to be used as the texture on the surface of the object.
#' @param image_repeat Default `1`. Number of times to repeat the image across the surface.
#' `u` and `v` repeat amount can be set independently if user passes in a length-2 vector.
#' @param alpha_texture Default `NA`. A matrix or filename (specifying a greyscale image) to 
#' be used to specify the transparency.
#' @param bump_texture Default `NA`. A matrix, array, or filename (specifying a greyscale image) to 
#' be used to specify a bump map for the surface.
#' @param bump_intensity Default `1`. Intensity of the bump map. High values may lead to unphysical results.
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
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Add a checkered rectangular cube below             
#' scene = scene %>%
#'   add_object(cube(x=555/2,y=555/8,z=555/2,xwidth=555/2,ywidth=555/4,zwidth=555/2,
#'   material = diffuse(checkercolor="purple",checkerperiod=20)))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#'   
#' #Add a marbled sphere           
#' scene = scene %>%
#'   add_object(sphere(x=555/2+555/4,y=555/2,z=555/2,radius=555/8,
#'   material = diffuse(noise=1/20)))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Add an orange volumetric (fog) cube           
#' scene = scene %>%
#'   add_object(cube(x=555/2-555/4,y=555/2,z=555/2,xwidth=555/4,ywidth=555/4,zwidth=555/4,
#'   material = diffuse(fog=TRUE, fogdensity=0.05,color="orange")))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #' #Add an line segment with a color gradient        
#' scene = scene %>%
#'   add_object(segment(start = c(555,450,450),end=c(0,450,450),radius = 50, 
#'                      material = diffuse(color="#1f7326", gradient_color = "#a60d0d")))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
diffuse = function(color = "#ffffff", 
                   checkercolor = NA, checkerperiod = 3,
                   noise = 0, noisephase = 0, noiseintensity = 10, noisecolor = "#000000",
                   gradient_color = NA, gradient_transpose = FALSE, 
                   gradient_point_start = NA, gradient_point_end = NA, gradient_type = "hsv",
                   image_texture = NA_character_, image_repeat = 1, alpha_texture = NA,
                   bump_texture = NA, bump_intensity = 1,
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
  if(!any(is.na(gradient_point_start)) && !any(is.na(gradient_point_end)) && !any(is.na(gradient_color))) {
    stopifnot(length(gradient_point_start) == 3)
    stopifnot(length(gradient_point_end) == 3)
    stopifnot(is.numeric(gradient_point_start))
    stopifnot(is.numeric(gradient_point_end))
    gradient_point_info = c(gradient_point_start,gradient_point_end)
    is_world_gradient = TRUE
  } else {
    is_world_gradient = FALSE
    gradient_point_info = NA
  }
  
  info = convert_color(color)
  noisecolor = convert_color(noisecolor)
  if(!is.array(image_texture) && !is.na(image_texture) && !is.character(image_texture)) {
    image_texture = NA_character_
    warning("Texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(alpha_texture) && !is.na(alpha_texture) && !is.character(alpha_texture)) {
    alpha_texture = NA
    warning("Alpha texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(bump_texture) && !is.na(bump_texture) && !is.character(bump_texture)) {
    bump_texture = NA
    warning("Bump texture not in recognized format (array, matrix, or filename), ignoring.")
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
  if(length(image_repeat) == 1) {
    image_repeat = c(image_repeat,image_repeat)
  }
  stopifnot(checkerperiod != 0)
  ray_material(list(type = type, 
               properties = list(info), checkercolor=list(c(checkercolor,checkerperiod)), 
               gradient_color = list(gradient_color), gradient_transpose = gradient_transpose,
               world_gradient = is_world_gradient,gradient_point_info = list(gradient_point_info),
               gradient_type = gradient_type,
               noise=noise, noisephase = noisephase * pi / 180, noiseintensity = noiseintensity, noisecolor = list(noisecolor),
               image = list(image_texture), image_repeat = list(image_repeat), 
               alphaimage = list(alpha_texture), lightintensity = NA_real_,
               fog=fog, fogdensity=fogdensity,implicit_sample = importance_sample, 
               sigma = sigma, glossyinfo = list(NA), bump_texture = list(bump_texture),
               bump_intensity = bump_intensity, roughness_texture = list(NA_character_)))
}

#' Metallic Material
#'
#' @param color Default `white`. The color of the sphere. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param eta Default `0`. Wavelength dependent refractivity of the material (red, green, and blue channels).
#' If single number, will be repeated across all three channels.
#' @param kappa Default `0`. Wavelength dependent absorption of the material (red, green, and blue channels).
#' If single number, will be repeated across all three channels.
#' @param fuzz  Default `0`. Deprecated--Use the microfacet material instead, as it is designed for rough metals. 
#' The roughness of the metallic surface. Maximum `1`.
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
#' @param gradient_point_start Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `color`.
#' @param gradient_point_end Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `gradient_color`.
#' @param gradient_type Default `hsv`. Colorspace to calculate the gradient. Alternative `rgb`.
#' @param image_texture Default `NA`. A 3-layer RGB array or filename to be used as the texture on the surface of the object.
#' @param image_repeat Default `1`. Number of times to repeat the image across the surface.
#' `u` and `v` repeat amount can be set independently if user passes in a length-2 vector.
#' @param alpha_texture Default `NA`. A matrix or filename (specifying a greyscale image) to be used to specify the transparency.
#' @param bump_texture Default `NA`. A matrix, array, or filename (specifying a greyscale image) to 
#' be used to specify a bump map for the surface.
#' @param bump_intensity Default `1`. Intensity of the bump map. High values may lead to unphysical results.
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
#' # Generate the cornell box with a single chrome sphere in the center. For other metals,
#' # See the website refractiveindex.info for eta and k data, use wavelengths 5
#' # 80nm (R), 530nm (G), and 430nm (B).
#' scene = generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,
#'   material=metal(eta=c(3.2176,3.1029,2.1839), k = c(3.3018,3.33,3.0339))))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=50,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' #Add an aluminum rotated shiny metal block     
#' scene = scene %>%
#'   add_object(cube(x=380,y=150/2,z=200,xwidth=150,ywidth=150,zwidth=150,
#'   material = metal(eta = c(1.07,0.8946,0.523), k = c(6.7144,6.188,4.95)),angle=c(0,45,0)))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' #Add a copper metal cube      
#' scene = scene %>%
#'   add_object(cube(x=150,y=150/2,z=300,xwidth=150,ywidth=150,zwidth=150,
#'                   material = metal(eta = c(0.497,0.8231,1.338), 
#'                                    k = c(2.898,2.476,2.298)),
#'                   angle=c(0,-30,0)))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Finally, let's add a lead pipe
#' scene2 = scene %>%
#'   add_object(cylinder(x=450,y=200,z=400,length=400,radius=30,
#'                   material = metal(eta = c(1.44,1.78,1.9), 
#'                                    k = c(3.18,3.36,3.43)),
#'                   angle=c(0,-30,0)))
#' if(run_documentation()) {
#' render_scene(scene2, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
metal = function(color = "#ffffff", 
                 eta = 0, kappa = 0, fuzz = 0,  
                 checkercolor = NA, checkerperiod = 3,
                 noise = 0, noisephase = 0, noiseintensity = 10, noisecolor = "#000000",
                 gradient_color = NA, gradient_transpose = FALSE,
                 gradient_point_start = NA, gradient_point_end = NA, gradient_type = "hsv",
                 image_texture = NA_character_, image_repeat = 1, alpha_texture = NA,
                 bump_texture = NA, bump_intensity = 1,
                 importance_sample = FALSE) {
  color = convert_color(color)
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
  if(!any(is.na(gradient_point_start)) && !any(is.na(gradient_point_end)) && !any(is.na(gradient_color))) {
    stopifnot(length(gradient_point_start) == 3)
    stopifnot(length(gradient_point_end) == 3)
    stopifnot(is.numeric(gradient_point_start))
    stopifnot(is.numeric(gradient_point_end))
    gradient_point_info = c(gradient_point_start,gradient_point_end)
    is_world_gradient = TRUE
  } else {
    is_world_gradient = FALSE
    gradient_point_info = NA
  }
  noisecolor = convert_color(noisecolor)
  if(!is.array(image_texture) && !is.na(image_texture) && !is.character(image_texture)) {
    image_texture = NA_character_
    warning("Texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(alpha_texture) && !is.na(alpha_texture) && !is.character(alpha_texture)) {
    alpha_texture = NA
    warning("Alpha texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(bump_texture) && !is.na(bump_texture) && !is.character(bump_texture)) {
    bump_texture = NA
    warning("Bump texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(length(eta) == 1) {
    eta = c(eta,eta,eta)
  }
  if(length(kappa) == 1) {
    kappa = c(kappa,kappa,kappa)
  }
  if(length(eta) > 3 || length(eta) == 2) {
    stop("eta must be either single number or 3-component vector")
  }
  if(length(kappa) > 3 || length(kappa) == 2) {
    stop("kappa must be either single number or 3-component vector")
  }
  if(length(image_repeat) == 1) {
    image_repeat = c(image_repeat,image_repeat)
  }
  glossyinfo = list(c(1, 0, 0, eta, kappa));
  ray_material(list(type = "metal", 
                 properties = list(c(color,fuzz)), 
                 checkercolor=list(c(checkercolor,checkerperiod)), 
                 gradient_color = list(gradient_color), gradient_transpose = gradient_transpose,
                 world_gradient = is_world_gradient,gradient_point_info = list(gradient_point_info),
                 gradient_type = gradient_type,
                 noise=noise, noisephase = noisephase * pi / 180, 
                 noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                 lightinfo = list(NA), 
                 image = list(image_texture), image_repeat = list(image_repeat),
                 alphaimage = list(alpha_texture), 
                 lightintensity = NA_real_,fog=FALSE,fogdensity=0.01,
                 implicit_sample = importance_sample, 
                 sigma = 0, glossyinfo = glossyinfo, bump_texture = list(bump_texture),
                 bump_intensity = bump_intensity, roughness_texture = list(NA_character_)))
}

#' Dielectric (glass) Material
#' 
#'
#' @param color Default `white`. The color of the surface. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param refraction Default `1.5`. The index of refraction.
#' @param attenuation Default `c(0,0,0)`. The Beer-Lambert color-channel specific exponential attenuation 
#' through the material. Higher numbers will result in less of that color making it through the material. If a character string
#' is provided (either as a named R color or a hex string), this will be converted to a length-3 vector equal to one minus the RGB
#' color vector, which should approximate the color being passed.
#' Note: This assumes the object has a closed surface. 
#' @param attenuation_intensity Default `1`. Changes the attenuation by a multiplicative factor. Values lower than one will make the 
#' dielectric more transparent, while values greater than one will make the glass more opaque.
#' @param priority Default `0`. When two dielectric materials overlap, the one with the lower priority value
#' is used for intersection. NOTE: If the camera is placed inside a dielectric object, its priority value
#' will not be taken into account when determining hits to other objects also inside the object.
#' @param importance_sample Default `FALSE`. If `TRUE`, the object will be sampled explicitly during 
#' the rendering process. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#' @param bump_texture Default `NA`. A matrix, array, or filename (specifying a greyscale image) to 
#' be used to specify a bump map for the surface.
#' @param bump_intensity Default `1`. Intensity of the bump map. High values may lead to unphysical results.
#'
#' @return Single row of a tibble describing the dielectric material.
#' @export
#'
#' @examples
#' #Generate a checkered ground
#' scene = generate_ground(depth=-0.5, material = diffuse(checkercolor="grey30",checkerperiod=2))
#' if(run_documentation()) {
#' render_scene(scene,parallel=TRUE)
#' }
#' 
#' #Add a glass sphere
#' if(run_documentation()) {
#' scene %>%
#'   add_object(sphere(x=-0.5,radius=0.5,material=dielectric())) %>%
#'   render_scene(parallel=TRUE,samples=128)
#' }
#' 
#' #Add a rotated colored glass cube
#' if(run_documentation()) {
#' scene %>%
#'   add_object(sphere(x=-0.5,radius=0.5,material=dielectric())) %>%
#'   add_object(cube(x=0.5,xwidth=0.5,material=dielectric(color="darkgreen"),angle=c(0,-45,0))) %>%
#'   render_scene(parallel=TRUE,samples=128)
#' }
#' 
#' #Add an area light behind and at an angle and turn off the ambient lighting
#' if(run_documentation()) {
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
#' if(run_documentation()) {
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
#' if(run_documentation()) {
#' generate_ground(depth=-0.51,material=diffuse(checkercolor="grey30",checkerperiod=2)) %>%
#'   add_object(cube(material = dielectric(priority=2, attenuation = c(10,3,10)))) %>%
#'   add_object(sphere(radius=0.49,material = dielectric(priority=1, refraction=1))) %>%
#'   add_object(sphere(radius=0.25,x=0.5,z=-0.5,y=0.5, 
#'                     material = dielectric(priority=0,attenuation = c(10,3,10) ))) %>%
#'   add_object(sphere(radius=0.25,x=-0.5,z=0.5,y=0.5,
#'                     material = dielectric(priority=0,attenuation = c(10,3,10)))) %>%
#'   render_scene(parallel=TRUE, samples = 128,lookfrom=c(5,1,5)) 
#' }
#' 
#' # We can also use this as a basic Constructive Solid Geometry interface by setting 
#' # the index of refraction equal to empty space, 1. This will subtract out those regions.
#' # Here I make a concave lens by subtracting two spheres from a cube.
#' if(run_documentation()) {
#' generate_ground(depth=-0.51,material=diffuse(checkercolor="grey30",checkerperiod=2,sigma=90)) %>%
#'   add_object(cube(material = dielectric(attenuation = c(3,3,1),priority=1))) %>%
#'   add_object(sphere(radius=1,x=1.01,
#'                     material = dielectric(priority=0,refraction=1))) %>%
#'   add_object(sphere(radius=1,x=-1.01, 
#'                     material = dielectric(priority=0,refraction=1))) %>%
#'   add_object(sphere(y=10,x=3,material=light(intensit=150))) %>%
#'   render_scene(parallel=TRUE, samples = 128,lookfrom=c(5,3,5))
#' }
dielectric = function(color="white", refraction = 1.5,  
                      attenuation = c(0,0,0), attenuation_intensity = 1,
                      priority = 0, importance_sample = FALSE,
                      bump_texture = NA, bump_intensity = 1) {
  color = convert_color(color)
  stopifnot(attenuation_intensity >= 0 && 
            is.numeric(attenuation_intensity) && 
            !is.na(attenuation_intensity) &&
            !is.infinite(attenuation_intensity))
  if(is.character(attenuation)) {
    attenuation = 1 - convert_color(attenuation) 
  } 
  attenuation = attenuation * attenuation_intensity
  
  if(!is.array(bump_texture) && !is.na(bump_texture) && !is.character(bump_texture)) {
    bump_texture = NA
    warning("Bump texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  ray_material(list(type = "dielectric", 
                      properties = list(c(color, refraction, attenuation, priority)), 
                      checkercolor=list(NA), 
                      gradient_color = list(NA), gradient_transpose = FALSE,
                      world_gradient = FALSE,gradient_point_info = list(NA),
                      gradient_type = NA,
                      noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                      image = list(NA_character_), image_repeat = list(c(1,1)), 
                      alphaimage = list(NA_character_), lightintensity = NA_real_, 
                      fog=FALSE, fogdensity=NA_real_, implicit_sample = importance_sample, 
                      sigma = 0, glossyinfo = list(NA), bump_texture = list(bump_texture),
                      bump_intensity = bump_intensity, roughness_texture = list(NA_character_)))
}

#' Microfacet Material
#'
#' @param color Default `white`. The color of the surface. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param roughness Default `0.0001`. Roughness of the surface, between `0` (smooth) and `1` (diffuse). 
#' Can be either a single number, or two numbers indicating an anisotropic distribution of normals. `0` is a smooth surface, while
#' `1` is extremely rough. This can be used to create a wide-variety of materials (e.g. `0-0.01` is specular 
#' metal, `0.02`-`0.1` is brushed metal, `0.1`-`0.3` is a rough metallic surface , `0.3`-`0.5` is diffuse, and 
#' above that is a rough satin-like material). 
#' Two numbers will specify the x and y roughness separately (e.g. `roughness = c(0.01, 0.001)` gives an 
#' etched metal effect). If `0`, this defaults to the `metal()` material for faster evaluation.
#' @param transmission Default `FALSE`. If `TRUE`, this material will be a rough dielectric instead of a rough metallic surface.
#' @param eta Default `0`. Wavelength dependent refractivity of the material (red, green, and blue channels).
#' If single number, will be repeated across all three channels. If `transmission = TRUE`, this is a single value representing the
#' index of refraction of the material.
#' @param kappa Default `0`. Wavelength dependent absorption of the material (red, green, and blue channels).
#' If single number, will be repeated across all three channels. If `transmission = TRUE`, this length-3 vector specifies 
#' the attenuation of the dielectric (analogous to the dielectric `attenuation` argument).
#' @param microfacet Default `tbr`.  Type of microfacet distribution. Alternative option `beckmann`.
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
#' @param gradient_point_start Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `color`.
#' @param gradient_point_end Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `gradient_color`.
#' @param gradient_type Default `hsv`. Colorspace to calculate the gradient. Alternative `rgb`.
#' @param image_texture Default `NA`. A 3-layer RGB array or filename to be used as the texture on the surface of the object.
#' @param image_repeat Default `1`. Number of times to repeat the image across the surface.
#' `u` and `v` repeat amount can be set independently if user passes in a length-2 vector.
#' @param alpha_texture Default `NA`. A matrix or filename (specifying a greyscale image) to be used to specify the transparency.
#' @param bump_texture Default `NA`. A matrix, array, or filename (specifying a greyscale image) to 
#' be used to specify a bump map for the surface.
#' @param bump_intensity Default `1`. Intensity of the bump map. High values may lead to unphysical results.
#' @param roughness_texture Default `NA`. A matrix, array, or filename (specifying a greyscale image) to 
#' be used to specify a roughness map for the surface.
#' @param roughness_range Default ` c(0.0001, 0.2)`. This is a length-2 vector that specifies the range of roughness values 
#' that the `roughness_texture` can take. 
#' @param roughness_flip Default `FALSE`. Setting this to `TRUE` flips the roughness values specified in the `roughness_texture`
#' so high values are now low values and vice versa.
#' @param importance_sample Default `FALSE`. If `TRUE`, the object will be sampled explicitly during 
#' the rendering process. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#'
#' @return Single row of a tibble describing the microfacet material.
#' @export
#'
#' @examples
#' # Generate a golden egg, using eta and kappa taken from physical measurements
#' # See the website refractiveindex.info for eta and k data, use 
#' # wavelengths 580nm (R), 530nm (G), and 430nm (B).
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(ellipsoid(x=555/2,555/2,y=150, a=100,b=150,c=100,
#'              material=microfacet(roughness=0.1,
#'                                  eta=c(0.216,0.42833,1.3184), kappa=c(3.239,2.4599,1.8661)))) %>% 
#'  render_scene(lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, parallel=TRUE,clamp_value=10)
#'  }
#' if(run_documentation()) {
#' #Make the roughness anisotropic (either horizontal or vertical), adding an extra light in front
#' #to show off the different microfacet orientations
#' generate_cornell() %>%
#'   add_object(sphere(x=555/2,z=50,y=75,radius=20,material=light())) %>% 
#'   add_object(ellipsoid(x=555-150,555/2,y=150, a=100,b=150,c=100,
#'              material=microfacet(roughness=c(0.3,0.1),
#'                                  eta=c(0.216,0.42833,1.3184), kappa=c(3.239,2.4599,1.8661)))) %>% 
#'  add_object(ellipsoid(x=150,555/2,y=150, a=100,b=150,c=100,
#'              material=microfacet(roughness=c(0.1,0.3),
#'                                  eta=c(0.216,0.42833,1.3184), kappa=c(3.239,2.4599,1.8661)))) %>%  
#'  render_scene(lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40,  parallel=TRUE,clamp_value=10)
#'}
#' if(run_documentation()) {
#' #Render a rough silver R with a smaller golden egg in front
#' generate_cornell() %>%
#'   add_object(obj_model(r_obj(),x=555/2,z=350,y=0, scale_obj = 200, angle=c(0,200,0),
#'              material=microfacet(roughness=0.2,
#'                                  eta=c(1.1583,0.9302,0.5996), kappa=c(6.9650,6.396,5.332)))) %>% 
#'  add_object(ellipsoid(x=200,z=200,y=80, a=50,b=80,c=50,
#'              material=microfacet(roughness=0.1,
#'                                  eta=c(0.216,0.42833,1.3184), kappa=c(3.239,2.4599,1.8661)))) %>% 
#'  render_scene(lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, parallel=TRUE,clamp_value=10)
#'  }
#' if(run_documentation()) {
#' #Increase the roughness
#' generate_cornell() %>%
#'   add_object(obj_model(r_obj(),x=555/2,z=350,y=0, scale_obj = 200, angle=c(0,200,0),
#'              material=microfacet(roughness=0.5,
#'                                  eta=c(1.1583,0.9302,0.5996), kappa=c(6.9650,6.396,5.332)))) %>% 
#'  add_object(ellipsoid(x=200,z=200,y=80, a=50,b=80,c=50,
#'              material=microfacet(roughness=0.3,
#'                                  eta=c(0.216,0.42833,1.3184), kappa=c(3.239,2.4599,1.8661)))) %>% 
#'  render_scene(lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, parallel=TRUE,clamp_value=10)
#'  }
#' if(run_documentation()) {
#'  #Use transmission for a rough dielectric
#' generate_cornell() %>%
#'   add_object(obj_model(r_obj(),x=555/2,z=350,y=0, scale_obj = 200, angle=c(0,200,0),
#'              material=microfacet(roughness=0.3, transmission=T, eta=1.6))) %>% 
#'  add_object(ellipsoid(x=200,z=200,y=80, a=50,b=80,c=50,
#'              material=microfacet(roughness=0.3, transmission=T, eta=1.6))) %>% 
#'  render_scene(lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, parallel=TRUE,clamp_value=10, min_variance=1e-6)
#' }
microfacet = function(color="white", roughness = 0.0001, transmission = FALSE,
                      eta = 0, kappa = 0, microfacet = "tbr", 
                      checkercolor = NA, checkerperiod = 3,
                      noise = 0, noisephase = 0, noiseintensity = 10, noisecolor = "#000000",
                      gradient_color = NA, gradient_transpose = FALSE,
                      gradient_point_start = NA_real_, gradient_point_end = NA_real_, gradient_type = "hsv",
                      image_texture = NA_character_, image_repeat = 1, alpha_texture = NA_character_,
                      bump_texture = NA_character_, bump_intensity = 1, roughness_texture = NA_character_,
                      roughness_range = c(0.0001, 0.2), roughness_flip = FALSE,
                      importance_sample = FALSE) {
  microtype = switch(microfacet, "tbr" = 1,"beckmann" = 2, 1)
  roughness[roughness <= 0] = 0
  roughness[roughness > 1] = 1
  if(length(roughness) == 1) {
    alphax = roughness^2
    alphay = roughness^2
  } else {
    alphax = roughness[1]^2
    alphay = roughness[2]^2
  }
  if(length(roughness_range) != 2) {
    stop("length of roughness_range must be 2")
  }
  if(roughness_range[1] > roughness_range[2]) {
    roughness_range = rev(roughness_range)
  }
  if(roughness_range[1] == roughness_range[2] && !is.na(roughness_texture)) {
    roughness_texture = NA
    roughness = roughness_range[1]
  }
  if(length(eta) == 1) {
    eta = c(eta,eta,eta)
  }
  if(length(kappa) == 1) {
    kappa = c(kappa,kappa,kappa)
  }
  if(length(eta) > 3 || length(eta) == 2) {
    stop("eta must be either single number or 3-component vector")
  }
  if(length(kappa) > 3 || length(kappa) == 2) {
    stop("kappa must be either single number or 3-component vector")
  }
  color = convert_color(color)
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
  if(!any(is.na(gradient_point_start)) && !any(is.na(gradient_point_end)) && !any(is.na(gradient_color))) {
    stopifnot(length(gradient_point_start) == 3)
    stopifnot(length(gradient_point_end) == 3)
    stopifnot(is.numeric(gradient_point_start))
    stopifnot(is.numeric(gradient_point_end))
    gradient_point_info = c(gradient_point_start,gradient_point_end)
    is_world_gradient = TRUE
  } else {
    is_world_gradient = FALSE
    gradient_point_info = NA
  }
  noisecolor = convert_color(noisecolor)
  if(!is.array(image_texture) && !is.na(image_texture) && !is.character(image_texture)) {
    image_texture = NA_character_
    warning("Texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(alpha_texture) && !is.na(alpha_texture) && !is.character(alpha_texture)) {
    alpha_texture = NA
    warning("Alpha texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(bump_texture) && !is.na(bump_texture) && !is.character(bump_texture)) {
    bump_texture = NA
    warning("Bump texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(roughness_texture) && !is.na(roughness_texture) && !is.character(roughness_texture)) {
    roughness_texture = NA
    warning("Roughness texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(length(image_repeat) == 1) {
    image_repeat = c(image_repeat,image_repeat)
  }
  roughness_flip = ifelse(roughness_flip,1,0)
  glossyinfo = list(c(microtype, alphax, alphay, eta, kappa, roughness_range,roughness_flip))
  if(alphax == 0 && alphay == 0 ) {
    if(!transmission) {
      ray_material(list(type = "metal", 
                          properties = list(c(color, 0)), 
                          gradient_color = list(gradient_color), gradient_transpose = FALSE,
                          world_gradient = is_world_gradient, gradient_point_info = list(gradient_point_info),
                          gradient_type = gradient_type,
                          checkercolor=list(c(checkercolor,checkerperiod)), 
                          noise=noise, noisephase = noisephase * pi / 180, 
                        noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                          image = list(image_texture), image_repeat = list(image_repeat), 
                          alphaimage = list(alpha_texture), lightintensity = NA_real_, 
                          fog=FALSE, fogdensity=NA_real_, implicit_sample = importance_sample, 
                          sigma = 0, glossyinfo = glossyinfo, bump_texture = list(bump_texture),
                          bump_intensity = bump_intensity, roughness_texture = list(NA_character_)))
    } else {
      ray_material(list(type = "dielectric", 
                          properties = list(c(color, eta[1], c(0,0,0), 0)), 
                          gradient_color = list(gradient_color), gradient_transpose = FALSE,
                          world_gradient = is_world_gradient, gradient_point_info = list(gradient_point_info),
                          gradient_type = gradient_type,
                          checkercolor=list(c(checkercolor,checkerperiod)), 
                          noise=noise, noisephase = noisephase * pi / 180, noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                          image = list(image_texture), image_repeat = list(image_repeat), 
                          alphaimage = list(alpha_texture), lightintensity = NA_real_, 
                          fog=FALSE, fogdensity=NA_real_, implicit_sample = importance_sample, 
                          sigma = 0, glossyinfo = glossyinfo, bump_texture = list(bump_texture),
                          bump_intensity = bump_intensity, roughness_texture = list(NA_character_)))
    }
  } else {
    if(transmission) {
      typeval = "mf-t"
    } else {
      typeval = "mf"
    }
    ray_material(list(type = typeval, 
                 properties = list(c(color)), 
                 gradient_color = list(gradient_color), gradient_transpose = FALSE,
                 world_gradient = is_world_gradient,gradient_point_info = list(gradient_point_info),
                 gradient_type = gradient_type,
                 checkercolor=list(c(checkercolor,checkerperiod)), 
                 noise=noise, noisephase = noisephase * pi / 180, noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                 image = list(image_texture), image_repeat = list(image_repeat), 
                 alphaimage = list(alpha_texture), lightintensity = NA_real_, 
                 fog=FALSE, fogdensity=NA_real_, implicit_sample = importance_sample, 
                 sigma = 0, glossyinfo = glossyinfo, bump_texture = list(bump_texture),
                 bump_intensity = bump_intensity, roughness_texture = list(roughness_texture)))
  }
}

#' Light Material
#'
#' @param color Default `white`. The color of the light Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param intensity Default `10`. If a positive value, this will turn this object into a light emitting the value specified
#' in `color` (ignoring other properties). Higher values will produce a brighter light.
#' @param importance_sample Default `TRUE`. Keeping this on for lights improves the convergence of the rendering 
#' algorithm, in most cases. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#' @param spotlight_focus Default `NA`, no spotlight. Otherwise, a length-3 numeric vector specifying
#' the x/y/z coordinates that the spotlight should be focused on. Only works for spheres and rectangles.
#' @param spotlight_width Default `30`. Angular width of the spotlight.
#' @param spotlight_start_falloff Default `15`. Angle at which the light starts fading in intensity.
#' @param invisible Default `FALSE`. If `TRUE`, the light itself will be invisible.
#' @param image_texture Default `NA`. A 3-layer RGB array or filename to be used as the texture on the surface of the object.
#' @param image_repeat Default `1`. Number of times to repeat the image across the surface.
#' `u` and `v` repeat amount can be set independently if user passes in a length-2 vector.
#' @param gradient_color Default `NA`. If not `NA`, creates a secondary color for a linear gradient 
#' between the this color and color specified in `color`. Direction is determined by `gradient_transpose`.
#' @param gradient_transpose Default `FALSE`. If `TRUE`, this will use the `v` coordinate texture instead
#' of the `u` coordinate texture to map the gradient.
#' @param gradient_point_start Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `color`.
#' @param gradient_point_end Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `gradient_color`.
#' @param gradient_type Default `hsv`. Colorspace to calculate the gradient. Alternative `rgb`.
#' 
#' @return Single row of a tibble describing the light material.
#' @export
#' @importFrom  grDevices col2rgb
#'
#' @examples
#' #Generate the cornell box without a light and add a single white sphere to the center
#' scene = generate_cornell(light=FALSE) %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,material=light()))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Remove the light for direct camera rays, but keep the lighting
#' scene = generate_cornell(light=FALSE) %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,
#'              material=light(intensity=15,invisible=TRUE)))
#' if(run_documentation()) {
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=128,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #All gather around the orb
#' scene = generate_ground(material = diffuse(checkercolor="grey50")) %>%
#'   add_object(sphere(radius=0.5,material=light(intensity=5,color="red"))) %>%
#'   add_object(obj_model(r_obj(), z=-3,x=-1.5,y=-1, angle=c(0,45,0))) %>%
#'   add_object(pig(scale=0.3, x=1.5,z=-2,y=-1.5,angle=c(0,-135,0)))
#' if(run_documentation()) {
#' render_scene(scene, samples=128, parallel=TRUE, clamp_value=10)
#' }
light = function(color = "#ffffff", intensity = 10, importance_sample = TRUE, 
                 spotlight_focus = NA, spotlight_width = 30, spotlight_start_falloff = 15,
                 invisible = FALSE, image_texture = NA_character_, image_repeat = 1, 
                 gradient_color = NA, gradient_transpose = FALSE,
                 gradient_point_start = NA, gradient_point_end = NA, gradient_type = "hsv") {
  info = convert_color(color)
  if(!is.array(image_texture) && !is.na(image_texture) && !is.character(image_texture)) {
    image_texture = NA_character_
    warning("Texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(length(image_repeat) == 1) {
    image_repeat = c(image_repeat,image_repeat)
  }
  if(all(!is.na(gradient_color))) {
    gradient_color = convert_color(gradient_color)
  } else {
    gradient_color = NA_real_
  }
  if(!any(is.na(gradient_point_start)) && !any(is.na(gradient_point_end)) && !any(is.na(gradient_color))) {
    stopifnot(length(gradient_point_start) == 3)
    stopifnot(length(gradient_point_end) == 3)
    stopifnot(is.numeric(gradient_point_start))
    stopifnot(is.numeric(gradient_point_end))
    gradient_point_info = c(gradient_point_start,gradient_point_end)
    is_world_gradient = TRUE
  } else {
    is_world_gradient = FALSE
    gradient_point_info = NA_real_
  }
  if(invisible) {
    invisible = 1
  } else {
    invisible = 0
  }
  if(all(!is.na(spotlight_focus))) {
    stopifnot(length(spotlight_focus) == 3)
    spotlight_width = min(c(spotlight_width,180))
    spotlight_start_falloff = min(c(spotlight_start_falloff,90))
    info = c(info, spotlight_focus, cospi(spotlight_width/180), cospi(spotlight_start_falloff/180), invisible)
    ray_material(list(type = "spotlight", 
                        properties = list(info), checkercolor=list(NA_real_), 
                        gradient_color = list(NA_real_), gradient_transpose = FALSE,
                        world_gradient = FALSE, gradient_point_info = list(NA_real_),
                        gradient_type = NA,
                        noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                        image = list(image_texture), image_repeat = list(image_repeat),
                        alphaimage = list(NA_character_), lightintensity = intensity,
                        fog=FALSE, fogdensity = 0.01, implicit_sample = importance_sample, 
                        sigma = 0, glossyinfo = list(NA_real_), bump_texture = list(NA_character_),
                        bump_intensity = 1, roughness_texture = list(NA_character_)))
  } else {
    info = c(info, invisible)
    ray_material(list(type = "light", 
                        properties = list(info), checkercolor=list(NA_real_), 
                        gradient_color = list(gradient_color), gradient_transpose = gradient_transpose,
                        world_gradient = is_world_gradient, gradient_point_info = list(gradient_point_info),
                        gradient_type = gradient_type,
                        noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                        image = list(image_texture), image_repeat = list(image_repeat),
                        alphaimage = list(NA_character_), lightintensity = intensity,
                        fog=FALSE, fogdensity = 0.01, implicit_sample = importance_sample, 
                        sigma = 0, glossyinfo = list(NA_real_), bump_texture = list(NA_character_),
                        bump_intensity = 1, roughness_texture = list(NA_character_)))
  }
}

#' Glossy Material
#'
#' @param color Default `white`. The color of the surface. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param gloss Default `0.8`. Gloss of the surface, between `1` (completely glossy) and `0` (rough glossy). 
#' Can be either a single number, or two numbers indicating an anisotropic distribution of normals (as in `microfacet()`).
#' @param reflectance Default `0.03`. The reflectivity of the surface. `1` is a full mirror, `0` is diffuse with a glossy highlight.
#' @param microfacet Default `tbr`.  Type of microfacet distribution. Alternative option `beckmann`.
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
#' @param gradient_point_start Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `color`.
#' @param gradient_point_end Default `NA`. If not `NA`, this changes the behavior from mapping texture coordinates to 
#' mapping to world space coordinates. This should be a length-3 vector specifying the x,y, and z points where the gradient
#' begins with value `gradient_color`.
#' @param gradient_type Default `hsv`. Colorspace to calculate the gradient. Alternative `rgb`.
#' @param image_texture Default `NA`. A 3-layer RGB array or filename to be used as the texture on the surface of the object.
#' @param image_repeat Default `1`. Number of times to repeat the image across the surface.
#' `u` and `v` repeat amount can be set independently if user passes in a length-2 vector.
#' @param alpha_texture Default `NA`. A matrix or filename (specifying a greyscale image) to be used to specify the transparency.
#' @param bump_texture Default `NA`. A matrix, array, or filename (specifying a greyscale image) to 
#' be used to specify a bump map for the surface.
#' @param bump_intensity Default `1`. Intensity of the bump map. High values may lead to unphysical results.
#' @param roughness_texture Default `NA`. A matrix, array, or filename (specifying a greyscale image) to 
#' be used to specify a roughness map for the surface.
#' @param roughness_range Default ` c(0.0001, 0.2)`. This is a length-2 vector that specifies the range of roughness values 
#' that the `roughness_texture` can take. 
#' @param roughness_flip Default `FALSE`. Setting this to `TRUE` flips the roughness values specified in the `roughness_texture`
#' so high values are now low values and vice versa.
#' @param importance_sample Default `FALSE`. If `TRUE`, the object will be sampled explicitly during 
#' the rendering process. If the object is particularly important in contributing to the light paths
#' in the image (e.g. light sources, refracting glass ball with caustics, metal objects concentrating light),
#' this will help with the convergence of the image.
#'
#' @return Single row of a tibble describing the glossy material.
#' @export
#'
#' @examples
#' if(run_documentation()) {
#' #Generate a glossy sphere
#' generate_ground(material=diffuse(sigma=90)) %>%
#'   add_object(sphere(y=0.2,material=glossy(color="#2b6eff"))) %>% 
#'   add_object(sphere(y=2.8,material=light())) %>%
#'   render_scene(parallel=TRUE,clamp_value=10,samples=128,sample_method="sobol_blue")
#'  }
#' if(run_documentation()) {
#' #Change the color of the underlying diffuse layer
#' generate_ground(material=diffuse(sigma=90)) %>%
#'   add_object(sphere(y=0.2,x=-2.1,material=glossy(color="#fc3d03"))) %>% 
#'   add_object(sphere(y=0.2,material=glossy(color="#2b6eff"))) %>% 
#'   add_object(sphere(y=0.2,x=2.1,material=glossy(color="#2fed4f"))) %>% 
#'   add_object(sphere(y=8,z=-5,radius=3,material=light(intensity=20))) %>%
#'   render_scene(parallel=TRUE,clamp_value=10,samples=128,fov=40,sample_method="sobol_blue")
#'  }
#' if(run_documentation()) {
#' #Change the amount of gloss 
#' generate_ground(material=diffuse(sigma=90)) %>%
#'   add_object(sphere(y=0.2,x=-2.1,material=glossy(gloss=1,color="#fc3d03"))) %>% 
#'   add_object(sphere(y=0.2,material=glossy(gloss=0.5,color="#2b6eff"))) %>% 
#'   add_object(sphere(y=0.2,x=2.1,material=glossy(gloss=0,color="#2fed4f"))) %>% 
#'   add_object(sphere(y=8,z=-5,radius=3,material=light(intensity=20))) %>%
#'   render_scene(parallel=TRUE,clamp_value=10,samples=128,fov=40,sample_method="sobol_blue")
#'  }
#' if(run_documentation()) {
#' #Add gloss to a pattern 
#' generate_ground(material=diffuse(sigma=90)) %>%
#'   add_object(sphere(y=0.2,x=-2.1,material=glossy(noise=2,noisecolor="black"))) %>% 
#'   add_object(sphere(y=0.2,material=glossy(color="#ff365a",checkercolor="#2b6eff"))) %>% 
#'   add_object(sphere(y=0.2,x=2.1,material=glossy(color="blue",gradient_color="#2fed4f"))) %>% 
#'   add_object(sphere(y=8,z=-5,radius=3,material=light(intensity=20))) %>%
#'   render_scene(parallel=TRUE,clamp_value=10,samples=128,fov=40,sample_method="sobol_blue")
#'  }
#' if(run_documentation()) {
#' #Add an R and a fill light (this may look familiar)
#' generate_ground(material=diffuse()) %>%
#'   add_object(sphere(y=0.2,material=glossy(color="#2b6eff",reflectance=0.05))) %>% 
#'   add_object(obj_model(r_obj(),z=1,y=-0.05,scale=0.45,material=diffuse())) %>%
#'   add_object(sphere(y=6,z=1,radius=4,material=light(intensity=3))) %>%
#'   add_object(sphere(z=15,material=light(intensity=50))) %>%
#'   render_scene(parallel=TRUE,clamp_value=10,samples=128,sample_method="sobol_blue")
#' }
glossy = function(color="white", gloss = 1, reflectance = 0.05, microfacet = "tbr", 
                  checkercolor = NA, checkerperiod = 3,
                  noise = 0, noisephase = 0, noiseintensity = 10, noisecolor = "#000000",
                  gradient_color = NA, gradient_transpose = FALSE,
                  gradient_point_start = NA_real_, gradient_point_end = NA_real_, gradient_type = "hsv",
                  image_texture = NA_character_, image_repeat = 1, alpha_texture = NA_character_, 
                  bump_texture = NA_character_, roughness_texture = NA_character_, bump_intensity = 1,
                  roughness_range = c(0.0001, 0.2), roughness_flip = FALSE,
                  importance_sample = FALSE) {
  microtype = switch(microfacet, "tbr" = 1,"beckmann" = 2, 1)
  gloss[gloss <= 0] = 0
  gloss[gloss > 1] = 1
  gloss = 1 - gloss
  gloss = gloss/2
  if(length(gloss) == 1) {
    alphax = gloss^2
    alphay = gloss^2
  } else {
    alphax = gloss[1]^2
    alphay = gloss[2]^2
  }
  if(length(roughness_range) != 2) {
    stop("length of roughness_range must be 2")
  }
  if(roughness_range[1] > roughness_range[2]) {
    roughness_range = rev(roughness_range)
  }
  if(roughness_range[1] == roughness_range[2] && !is.na(roughness_texture)) {
    roughness_texture = NA_character_
    gloss = 1-roughness_range[1]
  }
  color = convert_color(color)
  reflectance = rep(reflectance,3)
  if(all(!is.na(checkercolor))) {
    checkercolor = convert_color(checkercolor)
  } else {
    checkercolor = NA_real_
  }
  if(all(!is.na(gradient_color))) {
    gradient_color = convert_color(gradient_color)
  } else {
    gradient_color = NA_real_
  }
  if(!any(is.na(gradient_point_start)) && !any(is.na(gradient_point_end)) && !any(is.na(gradient_color))) {
    stopifnot(length(gradient_point_start) == 3)
    stopifnot(length(gradient_point_end) == 3)
    stopifnot(is.numeric(gradient_point_start))
    stopifnot(is.numeric(gradient_point_end))
    gradient_point_info = c(gradient_point_start,gradient_point_end)
    is_world_gradient = TRUE
  } else {
    is_world_gradient = FALSE
    gradient_point_info = NA_real_
  }
  noisecolor = convert_color(noisecolor)
  if(!is.array(image_texture) && !is.na(image_texture) && !is.character(image_texture)) {
    image_texture = NA_character_
    warning("Texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(alpha_texture) && !is.na(alpha_texture) && !is.character(alpha_texture)) {
    alpha_texture = NA_character_
    warning("Alpha texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(bump_texture) && !is.na(bump_texture) && !is.character(bump_texture)) {
    bump_texture = NA_character_
    warning("Bump texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(!is.array(roughness_texture) && !is.na(roughness_texture) && !is.character(roughness_texture)) {
    roughness_texture = NA_character_
    warning("Roughness texture not in recognized format (array, matrix, or filename), ignoring.")
  }
  if(length(image_repeat) == 1) {
    image_repeat = c(image_repeat,image_repeat)
  }
  roughness_flip = ifelse(roughness_flip,1,0)
  
  glossyinfo = list(c(microtype, alphax, alphay, reflectance, c(1,1,1), roughness_range, roughness_flip));
  ray_material(list(type = "glossy", 
                      properties = list(c(color)), 
                      gradient_color = list(gradient_color), gradient_transpose = FALSE,
                      world_gradient = is_world_gradient, gradient_point_info = list(gradient_point_info),
                      gradient_type = gradient_type,
                      checkercolor=list(c(checkercolor,checkerperiod)), 
                      noise=noise, noisephase = noisephase * pi / 180, noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                      image = list(image_texture), image_repeat = list(image_repeat), 
                      alphaimage = list(alpha_texture), lightintensity = NA_real_, 
                      fog=FALSE, fogdensity = NA_real_, implicit_sample = importance_sample, 
                      sigma = 0, glossyinfo = glossyinfo, bump_texture = list(bump_texture),
                      bump_intensity = bump_intensity, roughness_texture = list(roughness_texture)))
}


SigmaAFromConcentration = function(ce, cp) {
  eumelaninSigmaA = c(0.419, 0.697, 1.37);
  pheomelaninSigmaA = c(0.187, 0.4, 1.05);
  sigma_a = ce * eumelaninSigmaA + cp * pheomelaninSigmaA
  return(sigma_a);
}

SigmaAFromReflectance = function(c, beta_n) {
  sigma_a = (log(c) / (5.969 - 0.215 * beta_n + 2.532 * (beta_n)^2 -
             10.73 * (beta_n)^3 + 5.574 * (beta_n)^4 +
             0.245 * (beta_n)^5))^2;
  return(sigma_a);
}

#' Hair Material
#'
#' @param pigment Default `1.3`. Concentration of the eumelanin pigment in the hair. Blonde hair has concentrations around 0.3, brown around 1.3, and black around 8.
#' @param red_pigment Default `0`.Concentration of the pheomelanin pigment in the hair. Pheomelanin makes red hair red.
#' @param color Default `NA`. Approximate color. Overrides `pigment`/`redness` arguments.
#' @param sigma_a Default `NA`. Attenuation. Overrides `color` and `pigment`/`redness` arguments.
#' @param eta Default `1.55`. Index of refraction of the hair medium.
#' @param beta_m Default `0.3`. Longitudinal roughness of the hair. Should be between 0 and 1. This roughness controls the size and shape of the hair highlight.
#' @param beta_n Default `0.3`. Azimuthal roughness of the hair. Should be between 0 and 1.
#' @param alpha Default `2`. Angle of scales on the hair surface, in degrees.
#'
#' @return Single row of a tibble describing the hair material.
#' @export
#' @importFrom  grDevices col2rgb
#'
#' @examples
#' #Create a hairball
#' if(run_documentation()) {
#' #Generate rendom points on a sphere
#' lengthval = 0.5
#' theta = acos(2*runif(10000)-1.0);
#' phi = 2*pi*(runif(10000))
#' bezier_list = list()
#' 
#' #Grow the hairs
#' for(i in 1:length(phi)) {
#'   pointval = c(sin(theta[i]) * sin(phi[i]),
#'                cos(theta[i]),
#'                sin(theta[i]) * cos(phi[i]))
#'   bezier_list[[i]] = bezier_curve(width=0.01, width_end=0.008,
#'                                   p1 = pointval,
#'                                   p2 = (1+(lengthval*0.33))*pointval, 
#'                                   p3 = (1+(lengthval*0.66))*pointval,
#'                                   p4 = (1+(lengthval)) * pointval,
#'                                   material=hair(pigment = 0.3, red_pigment = 1.3,
#'                                                 beta_m = 0.3, beta_n= 0.3),
#'                                   type="flat")
#' }
#' hairball = dplyr::bind_rows(bezier_list)
#' 
#' generate_ground(depth=-2,material=diffuse(color="grey20")) %>%
#'   add_object(sphere()) %>%
#'   add_object(hairball) %>%
#'   add_object(sphere(y=20,z=20,radius=5,material=light(color="white",intensity = 100))) %>%
#'   render_scene(samples=64, lookfrom=c(0,3,10),clamp_value = 10,
#'                fov=20)
#' }
#' if(run_documentation()) {         
#'                
#' #Specify the color directly and increase hair roughness
#' for(i in 1:length(phi)) {
#'   pointval = c(sin(theta[i]) * sin(phi[i]),
#'                cos(theta[i]),
#'                sin(theta[i]) * cos(phi[i]))
#'   bezier_list[[i]] = bezier_curve(width=0.01, width_end=0.008,
#'                                   p1 = pointval,
#'                                   p2 = (1+(lengthval*0.33))*pointval, 
#'                                   p3 = (1+(lengthval*0.66))*pointval,
#'                                   p4 = (1+(lengthval)) * pointval,
#'                                   material=hair(color="purple",
#'                                                 beta_m = 0.5, beta_n= 0.5),
#'                                   type="flat")
#' }
#' hairball = dplyr::bind_rows(bezier_list)
#' generate_ground(depth=-2,material=diffuse(color="grey20")) %>%
#'   add_object(sphere()) %>%
#'   add_object(hairball) %>%
#'   add_object(sphere(y=20,z=20,radius=5,material=light(color="white",intensity = 100))) %>%
#'   render_scene(samples=64, lookfrom=c(0,3,10),clamp_value = 10,
#'                fov=20)
#' }
hair = function(pigment = 1.3, red_pigment = 0, color = NA, sigma_a = NA, 
                eta = 1.55, beta_m = 0.3, beta_n = 0.3, alpha = 2) {
  if(!is.na(sigma_a)) {
    sigma_a = sigma_a
  } else if(!is.na(color)) {
    sigma_a = SigmaAFromReflectance(convert_color(color), beta_n)
  } else {
    stopifnot(pigment >= 0)
    stopifnot(red_pigment >= 0)
    sigma_a = SigmaAFromConcentration(pigment, red_pigment)
  }
  
  info = c(sigma_a, eta, beta_m, beta_n, alpha)
  ray_material(list(type = "hair", 
                      properties = list(info), checkercolor=list(NA_real_), 
                      gradient_color = list(NA_real_), gradient_transpose = NA,
                      world_gradient = FALSE, gradient_point_info = list(NA_real_),
                      gradient_type = NA,
                      noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                      image = list(NA_character_), image_repeat = list(NA_real_),
                      alphaimage = list(NA_character_), lightintensity = NA_real_,
                      fog=FALSE, fogdensity = NA_real_, implicit_sample = FALSE, 
                      sigma = NA_real_, glossyinfo = list(NA_real_), bump_texture = list(NA_character_),
                      bump_intensity = NA_real_, roughness_texture = list(NA_character_)))
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
#' if(run_documentation()) {
#' scene = generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,material=lambertian()))
#'   render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=10,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
lambertian = function(...) {
  warning("lambertian() deprecated--use diffuse() instead.")
  diffuse(...)
}

