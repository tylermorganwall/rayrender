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
                      noise = 0, noisephase = 0, noiseintensity = 10, 
                      angle = 0, image_array = NA,
                      lightintensity = NA) {
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
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = "lambertian", 
                 properties = list(color), velocity = list(velocity), checkercolor=checkerlist, 
                 noise=noise, noisephase = noisephase, noiseintensity = noiseintensity, 
                 angle=angle, image = list(image_array), lightintensity = lightintensity,flipped=FALSE,
                 fog=FALSE,fogdensity=FALSE)
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
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = "metal", 
                 properties = list(c(color,fuzz)), velocity = list(velocity), 
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0,
                 islight = FALSE, lightinfo = list(NA),
                 angle=0,image = list(NA), lightintensity = NA,flipped=FALSE,fog=FALSE,fogdensity=FALSE)
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
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = "dielectric", 
                 properties = list(refraction), velocity = list(velocity),
                 checkercolor=list(NA), noise=0, noisephase = 0, noiseintensity = 0,
                 angle=0,image = list(NA),lightintensity = NA,flipped=FALSE,fog=FALSE,fogdensity=FALSE)
}

#' xy rectangle
#'
#' @param x0 Default `0`. x-coordinate of the center of the rectangle.
#' @param x1 Default `0`. y-coordinate of the center of the rectangle.
#' @param y0 Default `0`. x-coordinate of the center of the rectangle.
#' @param y1 Default `0`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param color Default `#ffffff`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' 
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
xy_rect = function(x0=0, x1=1, y0=0, y1=0, z=0, color = "#ffffff",lightintensity = NA,flipped=FALSE) {
  rectinfo = c(convert_color(color),x0,x1,y0,y1,z)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, type = "xy", 
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=list(NA), noise=NA, noisephase = 0, noiseintensity = NA,
                 angle=NA,image = list(NA),lightintensity = lightintensity,flipped=flipped,fog=FALSE,fogdensity=FALSE)
}

#' yz rectangle
#'
#' @param x0 Default `0`. x-coordinate of the center of the rectangle.
#' @param x1 Default `0`. y-coordinate of the center of the rectangle.
#' @param y0 Default `0`. x-coordinate of the center of the rectangle.
#' @param y1 Default `0`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param color Default `#ffffff`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' 
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
yz_rect = function(x0=0, x1=1, y0=0, y1=0, z=0, color = "#ffffff",lightintensity = NA,flipped=FALSE) {
  rectinfo = c(convert_color(color),x0,x1,y0,y1,z)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, type = "yz", 
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=list(NA), noise=NA, noisephase = 0, noiseintensity = NA,
                 angle=NA,image = list(NA),lightintensity = lightintensity,flipped=flipped,fog=FALSE,fogdensity=FALSE)
}

#' xz rectangle
#'
#' @param x0 Default `0`. x-coordinate of the center of the rectangle.
#' @param x1 Default `0`. y-coordinate of the center of the rectangle.
#' @param y0 Default `0`. x-coordinate of the center of the rectangle.
#' @param y1 Default `0`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param color Default `#ffffff`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' 
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
xz_rect = function(x0=0, x1=1, y0=0, y1=0, z=0, color = "#ffffff",lightintensity = NA,flipped=FALSE) {
  rectinfo = c(convert_color(color),x0,x1,y0,y1,z)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, type = "xz", 
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=list(NA), noise=NA, noisephase = 0, noiseintensity = NA,
                 angle=NA,image = list(NA),lightintensity = lightintensity,flipped=flipped,fog=FALSE,fogdensity=FALSE)
}

#' Box
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param xwidth Default `1`. y-coordinate of the center of the rectangle.
#' @param y Default `0`. x-coordinate of the center of the rectangle.
#' @param ywidth Default `1`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param zwidth Default `1`. z-coordinate of the center of the rectangle.
#' @param color Default `#ffffff`. The color of the sphere. Can be either
#' a hexadecimal code, or a numeric rgb vector listing three intensities between `0` and `1`.
#' 
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
box = function(x, x1, y, y1, z, z1, color = "#ffffff",lightintensity = NA,
                   noise = 0, noisephase = 0, noiseintensity = 10,
                   angle = 0, image_array = NA,fog=FALSE,fogdensity=0.01) {
  rectinfo = c(convert_color(color),x1,y1,z1)
  tibble::tibble(x=x,y=y,z=z,radius=NA, type = "box", 
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=list(NA), noise=noise, noisephase = noisephase, noiseintensity = noiseintensity,
                 angle=angle,image = list(image_array),lightintensity = lightintensity,flipped=FALSE,fog=fog,fogdensity=fogdensity)
}
