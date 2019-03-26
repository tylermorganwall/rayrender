#' Sphere Object
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param material Default  \code{\link{lambertian}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param velocity Default `c(0,0,0)`. Velocity of the sphere, used for motion blur.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
#' generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=100)) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
#' 
#' #Generate a GOLD sphere in the cornell box (it's GOLD!)
#' generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=100,z=555/2,radius=100,material=metal(color="gold",fuzz=0.2))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
#'   
#' #Add motion blur and show the sphere moving
#' generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=100,z=555/2,radius=100,
#'              material=metal(color="gold",fuzz=0.2),velocity=c(50,0,0))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
sphere = function(x=0, y=0, z=0, radius=1, material=lambertian(), 
                  angle = c(0,0,0), order_rotation = c(1,2,3), velocity = c(0,0,0), flipped=FALSE) {
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = material$type, shape="sphere",
                 properties = material$properties, velocity = list(velocity), 
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog,fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation))
}

#' Cube Object
#'
#' @param x Default `0`. x-coordinate of the center of the cube
#' @param y Default `0`. y-coordinate of the center of the cube
#' @param z Default `0`. z-coordinate of the center of the cube
#' @param width Default `1`. Cube width.
#' @param xwidth Default `1`. x-width of the cube. Overrides `width` argument for x-axis.
#' @param ywidth Default `1`. y-width of the cube. Overrides `width` argument for y-axis.
#' @param zwidth Default `1`. z-width of the cube. Overrides `width` argument for z-axis.
#' @param material Default  \code{\link{lambertian}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `0`. Angle of rotation around the y-axis.
#' @param velocity Default `c(0,0,0)`. Velocity of the cube, used for motion blur.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the cube in the scene.
#' @export
#'
#' @examples
#' #Generate a cube in the cornell box.
#' generate_cornell() %>%
#'   add_object(cube(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,zwidth=200,angle=c(0,30,0))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
#' 
#' #Generate a GOLD cube in the cornell box (it's GOLD!)
#' generate_cornell() %>%
#'   add_object(cube(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,zwidth=200,angle=c(0,30,0),
#'                   material = metal(color="gold", fuzz=0.2))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
cube = function(x=0, y=0, z=0, width=1, xwidth=1, ywidth=1, zwidth=1, 
                material=lambertian(), angle = c(0,0,0), order_rotation = c(1,2,3),velocity = c(0,0,0),
                flipped = FALSE) {
  xwidth = ifelse(missing(xwidth),width,xwidth)
  ywidth = ifelse(missing(ywidth),width,ywidth)
  zwidth = ifelse(missing(zwidth),width,zwidth)
  boxinfo = c(unlist(material$properties),xwidth,ywidth,zwidth)
  tibble::tibble(x=x,y=y,z=z,radius=NA, type = material$type, shape="box",
                 properties = list(boxinfo), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog, fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation))
}

#' Rectangular XY Plane Object 
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param y Default `0`. x-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param xwidth Default `1`. x-width of the rectangle.
#' @param ywidth Default `1`. y-width of the rectangle.
#' @param material Default  \code{\link{lambertian}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `0`. Angle of rotation around the y-axis.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the XY plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' generate_cornell() %>%
#'   add_object(xy_rect(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,
#'              material = lambertian(color="purple"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
#' 
#' #Generate a GOLD plane in the cornell box (it's GOLD!)
#' generate_cornell() %>%
#'   add_object(xy_rect(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,angle=c(0,30,0),
#'              material = metal(color="gold"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
xy_rect = function(x=0, y=0, z=0, xwidth=1, ywidth=1,  
                   material = lambertian(), angle = c(0,0,0), order_rotation = c(1,2,3), flipped=FALSE) {
  rectinfo = c(unlist(material$properties),x,xwidth,y,ywidth,z)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, type = material$type, shape="xy_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog,fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation))
}

#' Rectangular YZ Plane Object
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param y Default `0`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param ywidth Default `1`. y-width of the rectangle.
#' @param zwidth Default `1`. z-width of the rectangle.
#' @param material Default  \code{\link{lambertian}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `0`. Angle of rotation around the y-axis.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the YZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' generate_cornell() %>%
#'   add_object(yz_rect(x=100,y=100,z=555/2,ywidth=200,zwidth=200,
#'              material = lambertian(color="purple"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
#' 
#' #Generate a GOLD plane in the cornell box (it's GOLD!)
#' generate_cornell() %>%
#'   add_object(yz_rect(x=100,y=100,z=555/2,ywidth=200,zwidth=200, angle=c(0,30,0),
#'              material = metal(color="gold"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
yz_rect = function(x=0, y=0, z=0, ywidth=1, zwidth=1, material = lambertian(), 
                   angle = c(0,0,0), order_rotation = c(1,2,3), flipped=FALSE) {
  rectinfo = c(unlist(material$properties),y,ywidth,z,zwidth,x)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, type = material$type, shape="yz_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog, fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation))
}

#' Rectangular XZ Plane Object
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param y Default `0`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param xwidth Default `1`. x-width of the rectangle.
#' @param zwidth Default `1`. z-width of the rectangle.
#' @param material Default  \code{\link{lambertian}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `0`. Angle of rotation around the y-axis.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the XZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' generate_cornell() %>%
#'   add_object(xz_rect(x=555/2,y=100,z=555/2,xwidth=200,zwidth=200,
#'              material = lambertian(color="purple"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
#' 
#' #Generate a GOLD plane in the cornell box (it's GOLD!)
#' generate_cornell() %>%
#'   add_object(xz_rect(x=555/2,y=100,z=555/2,xwidth=200,zwidth=200,angle=c(0,30,0),
#'              material = metal(color="gold"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE)
xz_rect = function(x=0, xwidth=1, z=0, zwidth=1, y=0, material = lambertian(), 
                   angle = c(0,0,0), order_rotation = c(1,2,3), flipped=FALSE) {
  rectinfo = c(unlist(material$properties),x,xwidth,z,zwidth,y)
  tibble::tibble(x=NA,y=NA,z=NA,radius=NA, 
                 type = material$type, shape="xz_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog, fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation))
}
