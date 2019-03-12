#' Sphere Object
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param material Default `[lambertian()]`.The material, called from one of the material functions [lambertian()], [metal()], or [dielectric()].
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
sphere = function(x=0, y=0, z=0, radius=1, material=lambertian(), angle = 0, velocity = c(0,0,0)) {
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = material$type, shape="sphere",
                 properties = material$properties, velocity = list(velocity), 
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=angle,image = material$image,lightintensity = material$lightintensity,
                 flipped=NA,fog=material$fog,fogcolor = material$fogcolor, fogdensity=material$fogdensity)
}

#' Cube Object
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param material Default `[lambertian()]`.The material, called from one of the material functions [lambertian()], [metal()], or [dielectric()].
#' @param angle Default `0`. Angle of rotation around the y-axis.
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
cube = function(x=0, x1=1, y=0, y1=1, z=0, z1=1, material=lambertian(), angle = 0, velocity = c(0,0,0)) {
  boxinfo = c(unlist(material$properties),x1,y1,z1)
  tibble::tibble(x=x,y=y,z=z,radius=NA, type = material$type, shape="box",
                 properties = list(boxinfo), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=angle,image = material$image,lightintensity = material$lightintensity,
                 flipped=NA,fog=material$fog,fogcolor = material$fogcolor, fogdensity=material$fogdensity)
}

#' Rectangular XY Plane Object 
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param xwidth Default `1`. y-coordinate of the center of the rectangle.
#' @param y Default `0`. x-coordinate of the center of the rectangle.
#' @param ywidth Default `1`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param material Default `[lambertian()]`.The material, called from one of the material functions [lambertian()], [metal()], or [dielectric()].
#' @param angle Default `0`. Angle of rotation around the y-axis.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
xy_rect = function(x=0, xwidth=1, y=0, ywidth=1, z=0, material = lambertian(), angle = 0, flipped=FALSE) {
  rectinfo = c(unlist(material$properties),x,xwidth,y,ywidth,z)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, type = material$type, shape="xy_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor=material$noisecolor,
                 angle=angle,image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog,fogcolor = material$fogcolor, fogdensity=material$fogdensity)
}

#' Rectangular YZ Plane Object
#'
#' @param y Default `0`. x-coordinate of the center of the rectangle.
#' @param ywidth Default `1`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. x-coordinate of the center of the rectangle.
#' @param zwidth Default `1`. y-coordinate of the center of the rectangle.
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param material Default [lambertian()].The material, called from one of the material functions [lambertian()], [metal()], or [dielectric()].
#' @param angle Default `0`. Angle of rotation around the y-axis.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
yz_rect = function(y=0, ywidth=1, z=0, zwidth=1, x=0, material = lambertian(), angle = 0, flipped=FALSE) {
  rectinfo = c(unlist(material$properties),y,ywidth,z,zwidth,x)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, type = material$type, shape="yz_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=angle,image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog,fogcolor = material$fogcolor, fogdensity=material$fogdensity)
}

#' Rectangular XZ Plane Object
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param xwidth Default `1`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. x-coordinate of the center of the rectangle.
#' @param zwidth Default `1`. y-coordinate of the center of the rectangle.
#' @param y Default `0`. z-coordinate of the center of the rectangle.
#' @param material Default [lambertian()].The material, called from one of the material functions [lambertian()], [metal()], or [dielectric()].
#' @param angle Default `0`. Angle of rotation around the y-axis.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' dielectric(x=0,y=1,z=0,radius=1, refraction = 2)
xz_rect = function(x=0, xwidth=1, z=0, zwidth=1, y=0, material = lambertian(), angle = 0, flipped=FALSE) {
  rectinfo = c(unlist(material$properties),x,xwidth,z,zwidth,y)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, 
                 type = material$type, shape="xz_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=angle,image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog,fogcolor = material$fogcolor, fogdensity=material$fogdensity)
}
