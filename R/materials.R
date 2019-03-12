#' Lambertian
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param color Default `#ffffff`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param velocity Default `c(0,0,0)`. Velocity of the sphere, used for motion blur.
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#' @importFrom  grDevices col2rgb
#'
#' @examples
#' #create a single green sphere for the ground
#' scene = lambertian(x=0,y=-1000.5,z=0,radius=1000, color = c(0.8,1,0))
#' 
#' #add a single red sphere, using a hexcode
#' scene = rbind(scene, lambertian(x=0,y=0.5,z=0,radius = 1, color = "#00ff00"))
lambertian = function(color = "#ffffff", checkercolor = NA, checkerperiod = 3,
                      noise = 0, noisephase = 0, noiseintensity = 10, noisecolor = "#000000",
                      image_array = NA, 
                      lightintensity = NA, fog = FALSE, fogdensity = 0.01) {
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
  tibble::tibble(type = "lambertian", 
                 properties = list(color), checkercolor=list(c(checkercolor,checkerperiod)), 
                 noise=noise, noisephase = noisephase, noiseintensity = noiseintensity, noisecolor = list(noisecolor),
                 image = list(image_array), lightintensity = lightintensity,
                 fog=fog, fogcolor = list(color), fogdensity=fogdensity)
}

#' Metal
#'
#' @param color Default `#ffffff`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param fuzz  Default `0`. The roughness of the metallic surface. Maximum `1`.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' #create a single red metallic sphere
#' metal(x=0,y=1,z=0,radius=1, color = "#ff0000")
#' 
#' #create a single red brushed metallic sphere, using an RGB vector
#' metal(x=0,y=1,z=0,radius=1, color = c(1,0,0),fuzz=0)
metal = function(color = "#ffffff", fuzz = 0) {
  color = convert_color(color)
  tibble::tibble(type = "metal", 
                 properties = list(c(color,fuzz)), 
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 islight = FALSE, lightinfo = list(NA),
                 image = list(NA), fogcolor = list(NA),lightintensity = NA,fog=FALSE,fogdensity=0.01)
}

#' Dielectric (glass)
#'
#' @param refraction Default `1.5`.  Index of refraction.
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
dielectric = function(color="#ffffff", refraction = 1.5, fog = FALSE, fogcolor = "#ffffff", fogdensity = 0.01) {
  color = convert_color(color)
  fogcolor = convert_color(fogcolor)
  tibble::tibble(type = "dielectric", 
                 properties = list(c(color,refraction)), 
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0, noisecolor = list(c(0,0,0)),
                 image = list(NA), lightintensity = NA, 
                 fog=fog, fogcolor = list(fogcolor), fogdensity=fogdensity)
}

#' Fog (
