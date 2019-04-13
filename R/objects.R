#' Sphere Object
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param material Default  \code{\link{lambertian}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param velocity Default `c(0,0,0)`. Velocity of the sphere, used for motion blur.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=555/2,z=555/2,radius=100)) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel=TRUE, clamp_value=5)
#' }
#' 
#' #Generate a GOLD sphere in the cornell box (it's GOLD!)
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=100,z=555/2,radius=100,material=metal(color="gold",fuzz=0.2))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel=TRUE, clamp_value=5)
#' }
#'   
#' #Add motion blur and show the sphere moving
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(sphere(x=555/2,y=100,z=555/2,radius=100,
#'              material=metal(color="gold",fuzz=0.2),velocity=c(50,0,0))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel=TRUE, clamp_value=5)
#' }
sphere = function(x=0, y=0, z=0, radius=1, material=lambertian(), 
                  angle = c(0,0,0), order_rotation = c(1,2,3), velocity = c(0,0,0), flipped=FALSE) {
  tibble::tibble(x=x,y=y,z=z,radius=radius, type = material$type, shape="sphere",
                 properties = material$properties, velocity = list(velocity), 
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog,fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA))
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
#' @param angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param velocity Default `c(0,0,0)`. Velocity of the cube, used for motion blur.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the cube in the scene.
#' @export
#'
#' @examples
#' #Generate a cube in the cornell box.
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(cube(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,zwidth=200,angle=c(0,30,0))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel=TRUE, clamp_value=5)
#' }
#' #Generate a GOLD cube in the cornell box (it's GOLD!)
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(cube(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,zwidth=200,angle=c(0,30,0),
#'                   material = metal(color="gold", fuzz=0.2))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel=TRUE, clamp_value=5)
#' }
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
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA))
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
#' @param angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the XY plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(xy_rect(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,
#'              material = lambertian(color="purple"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel = TRUE, clamp_value=5)
#' }
#' 
#' #Generate a GOLD plane in the cornell box (it's GOLD!)
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(xy_rect(x=555/2,y=100,z=555/2,xwidth=200,ywidth=200,angle=c(0,30,0),
#'              material = metal(color="gold"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel = TRUE, clamp_value=5)
#' }
xy_rect = function(x=0, y=0, z=0, xwidth=1, ywidth=1,  
                   material = lambertian(), angle = c(0,0,0), order_rotation = c(1,2,3), flipped=FALSE) {
  rectinfo = c(unlist(material$properties),x,xwidth,y,ywidth,z)
  tibble::tibble(x=x,y=y,z=z,radius=NA, type = material$type, shape="xy_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog,fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA))
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
#' @param angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the YZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(yz_rect(x=100,y=100,z=555/2,ywidth=200,zwidth=200,
#'              material = lambertian(color="purple"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel = TRUE, clamp_value=5)
#' }
#' #Generate a GOLD plane in the cornell box (it's GOLD!)
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(yz_rect(x=100,y=100,z=555/2,ywidth=200,zwidth=200, angle=c(0,30,0),
#'              material = metal(color="gold"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel = TRUE, clamp_value=5)
#' }
yz_rect = function(x=0, y=0, z=0, ywidth=1, zwidth=1, material = lambertian(), 
                   angle = c(0,0,0), order_rotation = c(1,2,3), flipped=FALSE) {
  rectinfo = c(unlist(material$properties),y,ywidth,z,zwidth,x)
  tibble::tibble(x=x,y=y,z=z,radius=NA, type = material$type, shape="yz_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog, fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA))
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
#' @param angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the XZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(xz_rect(x=555/2,y=100,z=555/2,xwidth=200,zwidth=200,
#'              material = lambertian(color="purple"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel = TRUE, clamp_value=5)
#' }
#' 
#' #Generate a GOLD plane in the cornell box (it's GOLD!)
#' \dontrun{
#' generate_cornell() %>%
#'   add_object(xz_rect(x=555/2,y=100,z=555/2,xwidth=200,zwidth=200,angle=c(0,30,0),
#'              material = metal(color="gold"))) %>%
#'   render_scene(lookfrom = c(278,278,-800) ,lookat = c(278,278,0), fov = 40, ambient_light=FALSE,
#'                samples=500, parallel = TRUE, clamp_value=5)
#' }
xz_rect = function(x=0, xwidth=1, z=0, zwidth=1, y=0, material = lambertian(), 
                   angle = c(0,0,0), order_rotation = c(1,2,3), flipped=FALSE) {
  rectinfo = c(unlist(material$properties),x,xwidth,z,zwidth,y)
  tibble::tibble(x=x,y=y,z=z,radius=NA, 
                 type = material$type, shape="xz_rect",
                 properties = list(rectinfo), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog, fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA))
}

#' Triangle Object
#'
#' @param v1 Default `c(1,0,0)`. Length-3 vector indicating the x, y, and z coordinate of the first triangle vertex.
#' @param v2 Default `c(0,1,0)`. Length-3 vector indicating the x, y, and z coordinate of the second triangle vertex.
#' @param v3 Default `c(-1,0,0)`. Length-3 vector indicating the x, y, and z coordinate of the third triangle vertex.
#' @param n1 Default `NA`. Length-3 vector indicating the normal vector associated with the first triangle vertex.
#' @param n2 Default `NA`. Length-3 vector indicating the normal vector associated with the second triangle vertex.
#' @param n3 Default `NA`. Length-3 vector indicating the normal vector associated with the third triangle vertex.
#' @param color1 Default `NA`. Length-3 vector or string indicating the color associated with the first triangle vertex. 
#' If NA but other vertices specified, color inherets from material.
#' @param color2 Default `NA`. Length-3 vector or string indicating the color associated with the second triangle vertex.
#' If NA but other vertices specified, color inherets from material.
#' @param color3 Default `NA`. Length-3 vector or string indicating the color associated with the third triangle vertex.
#' If NA but other vertices specified, color inherets from material.
#' @param material Default  \code{\link{lambertian}}.The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' 
#' @return Single row of a tibble describing the XZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
triangle = function(v1=c(1,0,0), v2=c(0,1,0),v3=c(-1,0,0), 
                    n1 = rep(NA,3), n2 = rep(NA,3), n3 = rep(NA,3),
                    color1 = rep(NA,3), color2 = rep(NA,3), color3 = rep(NA,3),
                    material = lambertian(), 
                    angle = c(0,0,0), order_rotation = c(1,2,3), flipped=FALSE) {
  info = c(unlist(material$properties),v1,v2,v3, n1, n2, n3)
  if(all(!is.na(color1))) {
    color1 = convert_color(color1)
  }
  if(all(!is.na(color2))) {
    color2 = convert_color(color2)
  }
  if(all(!is.na(color3))) {
    color3 = convert_color(color3)
  }
  if(any(is.na(color1)) && any(!is.na(c(color2,color3)))) {
    color1 = info[1:3]
  }
  if(any(is.na(color2)) && any(!is.na(c(color1,color3)))) {
    color2 = info[1:3]
  }
  if(any(is.na(color3)) && any(!is.na(c(color1,color2)))) {
    color3 = info[1:3]
  }
  colorvec = c(color1,color2,color3)
  tibble::tibble(x=0,y=0,z=0,radius=NA, 
                 type = material$type, shape="triangle",
                 properties = list(info), velocity = list(c(0,0,0)),
                 checkercolor=material$checkercolor, 
                 noise=material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity,noisecolor=material$noisecolor,
                 angle=list(angle),image = material$image,lightintensity = material$lightintensity,
                 flipped=flipped,fog=material$fog, fogdensity=material$fogdensity,
                 implicit_sample=material$implicit_sample,order_rotation=list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(colorvec))
}
