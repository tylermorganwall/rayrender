#' Lambertian Sphere
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
lambertian = function(x=0, y=0, z=0, radius=1, color = "#ffffff", checkercolor = NA, velocity = c(0,0,0),
                      noise = 0, noisephase = 0, noiseintensity = 10, angle = 0, image_array = NA) {
  if(all(!is.na(checkercolor))) {
    checkerlist = list(convert_color(checkercolor))
  } else {
    checkerlist = list(NA)
  }
  color = convert_color(color)
  if(!is.array(image_array) && !is.na(image_array)) {
    image = NA
    warning("Image not in recognized format (array or matrix), ignoring")
  }
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = "lambertian", properties = list(color), velocity = list(velocity), checkercolor=checkerlist, 
                 noise=noise, noisephase = noisephase, noiseintensity = noiseintensity, angle=angle, image = list(image_array))
}

#' Metal Sphere
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param color Default `#ffffff`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' @param fuzz  Default `0`. The roughness of the metallic surface. Maximum `1`.
#' @param velocity Default `c(0,0,0)`. Velocity of the sphere, used for motion blur.
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
metal = function(x=0, y=0, z=0, radius=1, color = "#ffffff", fuzz = 0,velocity = c(0,0,0)) {
  color = convert_color(color)
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = "metal", properties = list(c(color,fuzz)), velocity = list(velocity), 
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0,angle=0,image = list(NA))
}

#' Dielelectric Sphere
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param refraction Default `1.5`.  Index of refraction.
#' @param velocity Default `c(0,0,0)`. Velocity of the sphere, used for motion blur.
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
dielectric = function(x=0, y=0, z=0, radius=1, refraction = 1.5,velocity = c(0,0,0)) {
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = "dielectric", properties = list(refraction), velocity = list(velocity),
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0,angle=0,image = list(NA))
}