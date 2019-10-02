#' Lambertian (diffuse) Material
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
#' @param image_array A 3-layer RGB array to be used as the texture on the surface of the object.
#' @param lightintensity Default `NA`. If a positive value, this will turn this object into a light emitting the value specified
#' in `color` (ignoring other properties). Higher values will produce a brighter light.
#' @param fog Default `FALSE`. If `TRUE`, the object will be a volumetric scatterer.
#' @param fogdensity Default `0.01`. The density of the fog. Higher values will produce more opaque objects.
#' @param implicit_sample Default `FALSE`, unless the object is a light. If `TRUE`, the object will
#' be sampled as part of the scattering probability density function.
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#' @importFrom  grDevices col2rgb
#'
#' @examples
#' #Generate the cornell box and add a single white sphere to the center
#' scene = generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=555/8,material=lambertian()))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Add a checkered rectangular cube below             
#' scene = scene %>%
#'   add_object(cube(x=555/2,y=555/8,z=555/2,xwidth=555/2,ywidth=555/4,zwidth=555/2,
#'   material = lambertian(checkercolor="purple",checkerperiod=20)))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#'   
#' #Add a marbled sphere           
#' scene = scene %>%
#'   add_object(sphere(x=555/2+555/4,y=555/2,z=555/2,radius=555/8,
#'   material = lambertian(noise=1/20)))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
#' 
#' #Add an orange volumetric (fog) cube           
#' scene = scene %>%
#'   add_object(cube(x=555/2-555/4,y=555/2,z=555/2,xwidth=555/4,ywidth=555/4,zwidth=555/4,
#'   material = lambertian(fog=TRUE, fogdensity=0.05,color="orange")))
#' \donttest{
#' render_scene(scene, lookfrom=c(278,278,-800),lookat = c(278,278,0), samples=500,
#'              aperture=0, fov=40, ambient_light=FALSE, parallel=TRUE)
#' }
lambertian = function(color = "#ffffff", checkercolor = NA, checkerperiod = 3,
                      noise = 0, noisephase = 0, noiseintensity = 10, noisecolor = "#000000",
                      image_array = NA, 
                      lightintensity = NA, fog = FALSE, fogdensity = 0.01, implicit_sample = FALSE) {
  if(all(!is.na(checkercolor))) {
    checkercolor = convert_color(checkercolor)
  } else {
    checkercolor = NA
  }
  color = convert_color(color)
  noisecolor = convert_color(noisecolor)
  if(!is.array(image_array) && !is.na(image_array)) {
    image = NA
    warning("Image not in recognized format (array or matrix), ignoring")
  }
  assertthat::assert_that(checkerperiod != 0)
  tibble::tibble(type = "lambertian", 
                 properties = list(color), checkercolor=list(c(checkercolor,checkerperiod)), 
                 noise=noise, noisephase = noisephase, noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                 image = list(image_array), lightintensity = lightintensity,
                 fog=fog, fogdensity=fogdensity,implicit_sample = implicit_sample)
}

#' Metallic Material
#'
#' @param color Default `white`. The color of the sphere. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param fuzz  Default `0`. The roughness of the metallic surface. Maximum `1`.
#' @param implicit_sample Default `FALSE`. If `TRUE`, the object will
#' be sampled as part of the scattering probability density function.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the sphere in the scene.
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
metal = function(color = "#ffffff", fuzz = 0,  implicit_sample = FALSE) {
  color = convert_color(color)
  tibble::tibble(type = "metal", 
                 properties = list(c(color,fuzz)), 
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 islight = FALSE, lightinfo = list(NA),
                 image = list(NA), lightintensity = NA,fog=FALSE,fogdensity=0.01,
                 implicit_sample = implicit_sample)
}

#' Dielectric (glass) Material
#'
#' @param color Default `white`. The color of the surface. Can be either
#' a hexadecimal code, R color string, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param refraction Default `1.5`. The index of refraction.
#' @param implicit_sample Default `TRUE`. If `FALSE`, the object will not 
#' be sampled as part of the scattering probability density function.
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' #Generate a checkered ground
#' scene = generate_ground(depth=-0.5,
#'                         material=lambertian(color="white", checkercolor="grey30",checkerperiod=2))
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
#'                      material=lambertian(lightintensity=15),
#'                      angle=c(0,-90,45), order_rotation = c(3,2,1))) %>%
#'   render_scene(parallel=TRUE,aperture=0, ambient_light=FALSE,samples=1000)
#' }
dielectric = function(color="white", refraction = 1.5, implicit_sample = FALSE) {
  color = convert_color(color)
  tibble::tibble(type = "dielectric", 
                 properties = list(c(color,refraction)), 
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 image = list(NA), lightintensity = NA, 
                 fog=FALSE, fogdensity=NA, implicit_sample = implicit_sample)
}
