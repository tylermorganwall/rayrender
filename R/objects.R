#' Sphere Object
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @param material Default  \code{\link{diffuse}}. The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param velocity Default `c(0, 0, 0)`. Velocity of the sphere, used for motion blur.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(sphere(x = 555/2, y = 555/2, z = 555/2, radius = 100)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a gold sphere in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(sphere(x = 555/2, y = 100, z = 555/2, radius = 100, 
#'                     material = metal(color = "gold", fuzz = 0.2))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#'   
#' #Add motion blur and show the sphere moving
#' \donttest{
#' generate_cornell() %>%
#'   add_object(sphere(x = 555/2, y = 100, z = 555/2, radius = 100,
#'              material = metal(color = "gold", fuzz = 0.2), velocity = c(50, 0, 0))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
sphere = function(x = 0, y = 0, z = 0, radius = 1, material = diffuse(), 
                  angle = c(0, 0, 0), order_rotation = c(1, 2, 3), velocity = c(0, 0, 0), 
                  flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  tibble::tibble(x = x, y = y, z = z, radius = radius, type = material$type, shape = "sphere",
                 properties = material$properties, velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample, sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
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
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param velocity Default `c(0, 0, 0)`. Velocity of the cube.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the cube in the scene.
#' @export
#'
#' @examples
#' #Generate a cube in the cornell box.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 100, z = 555/2, 
#'                   xwidth = 200, ywidth = 200, zwidth = 200, angle = c(0, 30, 0))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' #Generate a gold cube in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 100, z = 555/2, 
#'                   xwidth = 200, ywidth = 200, zwidth = 200, angle = c(0, 30, 0),
#'                   material = metal(color = "gold", fuzz = 0.2))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a rotateddielectric box in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 200, z = 555/2, 
#'                   xwidth = 200, ywidth = 100, zwidth = 200, angle = c(30, 30, 30),
#'                   material = dielectric())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
cube = function(x = 0, y = 0, z = 0, width = 1, xwidth = 1, ywidth = 1, zwidth = 1, 
                material = diffuse(), angle = c(0, 0, 0), order_rotation = c(1, 2, 3), velocity = c(0, 0, 0),
                flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  xwidth = ifelse(missing(xwidth), width, xwidth)
  ywidth = ifelse(missing(ywidth), width, ywidth)
  zwidth = ifelse(missing(zwidth), width, zwidth)
  boxinfo = c(unlist(material$properties), xwidth, ywidth, zwidth)
  tibble::tibble(x = x, y = y, z = z, radius = NA, type = material$type, shape = "box",
                 properties = list(boxinfo), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' Rectangular XY Plane Object 
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param y Default `0`. x-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param xwidth Default `1`. x-width of the rectangle.
#' @param ywidth Default `1`. y-width of the rectangle.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Single row of a tibble describing the XY plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(xy_rect(x = 555/2, y = 100, z = 555/2, xwidth = 200, ywidth = 200,
#'              material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800), lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a gold plane in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(xy_rect(x = 555/2, y = 100, z = 555/2, 
#'                      xwidth = 200, ywidth = 200, angle = c(0, 30, 0),
#'                      material = metal(color = "gold"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
xy_rect = function(x = 0, y = 0, z = 0, xwidth = 1, ywidth = 1,  
                   material = diffuse(), angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  rectinfo = c(unlist(material$properties), x, xwidth, y, ywidth, z)
  tibble::tibble(x = x, y = y, z = z, radius = NA, type = material$type, shape = "xy_rect",
                 properties = list(rectinfo), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' Rectangular YZ Plane Object
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param y Default `0`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param ywidth Default `1`. y-width of the rectangle.
#' @param zwidth Default `1`. z-width of the rectangle.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Single row of a tibble describing the YZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(yz_rect(x = 100, y = 100, z = 555/2, ywidth = 200, zwidth = 200,
#'                      material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' #Generate a gold plane in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(yz_rect(x = 100, y = 100, z = 555/2, 
#'                      ywidth = 200, zwidth = 200, angle = c(0, 30, 0),
#'                      material = metal(color = "gold"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
yz_rect = function(x = 0, y = 0, z = 0, ywidth = 1, zwidth = 1, material = diffuse(), 
                   angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  rectinfo = c(unlist(material$properties), y, ywidth, z, zwidth, x)
  tibble::tibble(x = x, y = y, z = z, radius = NA, type = material$type, shape = "yz_rect",
                 properties = list(rectinfo), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' Rectangular XZ Plane Object
#'
#' @param x Default `0`. x-coordinate of the center of the rectangle.
#' @param y Default `0`. y-coordinate of the center of the rectangle.
#' @param z Default `0`. z-coordinate of the center of the rectangle.
#' @param xwidth Default `1`. x-width of the rectangle.
#' @param zwidth Default `1`. z-width of the rectangle.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Single row of a tibble describing the XZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a purple rectangle in the cornell box.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(xz_rect(x = 555/2, y = 100, z = 555/2, xwidth = 200, zwidth = 200,
#'              material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a gold plane in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(xz_rect(x = 555/2, y = 100, z = 555/2, 
#'              xwidth = 200, zwidth = 200, angle = c(0, 30, 0),
#'              material = metal(color = "gold"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
xz_rect = function(x = 0, xwidth = 1, z = 0, zwidth = 1, y = 0, material = diffuse(), 
                   angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  rectinfo = c(unlist(material$properties), x, xwidth, z, zwidth, y)
  tibble::tibble(x = x, y = y, z = z, radius = NA, 
                 type = material$type, shape = "xz_rect",
                 properties = list(rectinfo), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' Triangle Object
#'
#' @param v1 Default `c(1, 0, 0)`. Length-3 vector indicating the x, y, and z coordinate of the first triangle vertex.
#' @param v2 Default `c(0, 1, 0)`. Length-3 vector indicating the x, y, and z coordinate of the second triangle vertex.
#' @param v3 Default `c(-1, 0, 0)`. Length-3 vector indicating the x, y, and z coordinate of the third triangle vertex.
#' @param n1 Default `NA`. Length-3 vector indicating the normal vector associated with the first triangle vertex.
#' @param n2 Default `NA`. Length-3 vector indicating the normal vector associated with the second triangle vertex.
#' @param n3 Default `NA`. Length-3 vector indicating the normal vector associated with the third triangle vertex.
#' @param color1 Default `NA`. Length-3 vector or string indicating the color associated with the first triangle vertex. 
#' If NA but other vertices specified, color inherits from material.
#' @param color2 Default `NA`. Length-3 vector or string indicating the color associated with the second triangle vertex.
#' If NA but other vertices specified, color inherits from material.
#' @param color3 Default `NA`. Length-3 vector or string indicating the color associated with the third triangle vertex.
#' If NA but other vertices specified, color inherits from material.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Single row of a tibble describing the XZ plane in the scene.
#' @export
#'
#' @examples
#' #Generate a triangle in the Cornell box.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(triangle(v1 = c(100, 100, 100), v2 = c(555/2, 455, 455), v3 = c(455, 100, 100),
#'                       material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' #Pass individual colors to each vertex:
#' \donttest{
#' generate_cornell() %>%
#'   add_object(triangle(v1 = c(100, 100, 100), v2 = c(555/2, 455, 455), v3 = c(455, 100, 100),
#'                       color1 = "green", color2 = "yellow", color3 = "red")) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
triangle = function(v1 = c(1, 0, 0), v2 = c(0, 1, 0), v3 = c(-1, 0, 0), 
                    n1 = rep(NA, 3), n2 = rep(NA, 3), n3 = rep(NA, 3),
                    color1 = rep(NA, 3), color2 = rep(NA, 3), color3 = rep(NA, 3),
                    material = diffuse(), 
                    angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                    flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  info = c(unlist(material$properties), v1, v2, v3, n1, n2, n3)
  if(all(!is.na(color1))) {
    color1 = convert_color(color1)
  }
  if(all(!is.na(color2))) {
    color2 = convert_color(color2)
  }
  if(all(!is.na(color3))) {
    color3 = convert_color(color3)
  }
  if(any(is.na(color1)) && any(!is.na(c(color2, color3)))) {
    color1 = info[1:3]
  }
  if(any(is.na(color2)) && any(!is.na(c(color1, color3)))) {
    color2 = info[1:3]
  }
  if(any(is.na(color3)) && any(!is.na(c(color1, color2)))) {
    color3 = info[1:3]
  }
  colorvec = c(color1, color2, color3)
  tibble::tibble(x = 0, y = 0, z = 0, radius = NA, 
                 type = material$type, shape = "triangle",
                 properties = list(info), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(colorvec), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' Disk Object
#'
#' @param x Default `0`. x-coordinate of the center of the disk
#' @param y Default `0`. y-coordinate of the center of the disk
#' @param z Default `0`. z-coordinate of the center of the disk
#' @param radius Default `1`. Radius of the disk.
#' @param inner_radius Default `0`. Inner radius of the disk.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param velocity Default `c(0, 0, 0)`. Velocity of the disk.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the disk in the scene.
#' @export
#'
#' @examples
#' #Generate a disk in the cornell box.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(disk(x = 555/2, y = 50, z = 555/2, radius = 150, 
#'                   material = diffuse(color = "orange"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' #Rotate the disk.
#' \donttest{
#' generate_cornell() %>%
#'   add_object(disk(x = 555/2, y = 555/2, z = 555/2, radius = 150, angle = c(45, 0, 0), 
#'                   material = diffuse(color = "orange"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) , lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' #Pass a value for the inner radius.
#' \donttest{
#' generate_cornell() %>% 
#'   add_object(disk(x = 555/2, y = 555/2, z = 555/2, 
#'                   radius = 150, inner_radius = 75, angle = c(45, 0, 0), 
#'                   material = diffuse(color = "orange"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
disk = function(x = 0, y = 0, z = 0, radius = 1, inner_radius = 0, material = diffuse(), 
                angle = c(0, 0, 0), order_rotation = c(1, 2, 3), velocity = c(0, 0, 0), 
                flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  info = c(unlist(material$properties), inner_radius)
  tibble::tibble(x = x, y = y, z = z, radius = radius, type = material$type, shape = "disk",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' `obj` File Object
#' 
#' Load an obj file via a filepath. Currently only supports the diffuse texture with the `texture` argument. 
#' Note: light importance sampling currently not supported for this shape.
#'
#' @param filename Filename and path to the `obj` file. Can also be a `txt` file, if it's in the correct `obj` internally.
#' @param x Default `0`. x-coordinate to offset the model.
#' @param y Default `0`. y-coordinate to offset the model.
#' @param z Default `0`. z-coordinate to offset the model.
#' @param scale_obj Default `1`. Amount to scale the model. Use this to scale the object up or down on all axes, as it is
#' more robust to numerical precision errors than the generic scale option.
#' @param texture Default `FALSE`. Whether to load the obj file texture.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}. 
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Single row of a tibble describing the obj model in the scene.
#' @export
#'
#' @examples
#' #Load the included example R object file, by calling the r_obj() function. This
#' #returns the local file path to the `r.txt` obj file. The file extension is "txt" 
#' #due to package constraints, but the file contents are identical and it does not 
#' #affect the function.
#' 
#' \donttest{
#' generate_ground(material = diffuse(checkercolor = "grey50")) %>%
#'   add_object(obj_model(y = -0.8, filename = r_obj(),
#'                        material = metal(color = "gold", fuzz = 0.025))) %>%
#'   add_object(obj_model(x = 1.8, y = -0.8, filename = r_obj(), 
#'                        material = diffuse(color = "lightblue"))) %>%
#'   add_object(obj_model(x = -1.8, y = -0.8, filename = r_obj() , 
#'                        material = dielectric(color = "pink"))) %>%
#'   add_object(sphere(z = 20, x = 20, y = 20, radius = 10,
#'                     material = light(intensity = 20))) %>%
#'   render_scene(parallel = TRUE, samples = 400, 
#'                tonemap = "reinhold", aperture = 0.05, fov = 32, lookfrom = c(0, 2, 10))
#' }
#' 
#' #Use scale_obj to make objects bigger--this is more robust than the generic scale argument.
#' \donttest{
#' generate_ground(material = diffuse(checkercolor = "grey50")) %>%
#'   add_object(obj_model(y = -0.8, filename = r_obj(), scale_obj = 2,
#'                        material = diffuse(noise = TRUE, noiseintensity = 10,noisephase=45))) %>%
#'   add_object(sphere(z = 20, x = 20, y = 20, radius = 10,
#'                     material = light(intensity = 10))) %>%
#'   render_scene(parallel = TRUE, samples = 400, ambient = TRUE, 
#'                backgroundhigh="blue", backgroundlow="red",
#'                aperture = 0.05, fov = 32, lookfrom = c(0, 2, 10),
#'                lookat = c(0,1,0))
#' }
obj_model = function(filename, x = 0, y = 0, z = 0, scale_obj = 1, texture = FALSE,
                    material = diffuse(), 
                    angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                    flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  info = c(unlist(material$properties), scale_obj)
  if(texture) {
    shape = "objcolor"
  } else {
    shape = "obj"
  }
  tibble::tibble(x = x, y = y, z = z, radius = NA, 
                 type = material$type, shape = shape,
                 properties = list(info), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = filename, scale_factor = list(scale), group_scale = list(NA))
}

#' Cylinder Object
#'
#' @param x Default `0`. x-coordinate of the center of the cylinder
#' @param y Default `0`. y-coordinate of the center of the cylinder
#' @param z Default `0`. z-coordinate of the center of the cylinder
#' @param radius Default `1`. Radius of the cylinder.
#' @param length Default `1`. Length of the cylinder.
#' @param phi_min Default `0`. Minimum angle around the segment.
#' @param phi_max Default `360`. Maximum angle around the segment.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param velocity Default `c(0, 0, 0)`. Velocity of the cylinder.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the cylinder in the scene.
#' @export
#'
#' @examples
#' #Generate a cylinder in the cornell box. Add a cap to both ends.
#' 
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
#'                       length = 300, radius = 100, material = metal())) %>%
#'   add_object(disk(x = 555/2, y = 400, z = 555/2, 
#'                   radius = 100, material = metal())) %>%
#'   add_object(disk(x = 555/2, y = 100, z = 555/2, 
#'                   radius = 100, material = metal(), flipped = TRUE)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' #Rotate the cylinder
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
#'                       length = 300, radius = 100, angle = c(0, 0, 45),
#'                       material = diffuse())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' # Only render a subtended arc of the cylinder,
#' \donttest{
#' generate_cornell(lightintensity=3) %>%
#'   add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
#'                       length = 300, radius = 100, angle = c(45, 0, 0), phi_min = 0, phi_max = 180,
#'                       material = diffuse())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
cylinder = function(x = 0, y = 0, z = 0, radius = 1, length = 1,
                    phi_min = 0, phi_max = 360, material = diffuse(), 
                    angle = c(0, 0, 0), order_rotation = c(1, 2, 3), velocity = c(0, 0, 0), 
                    flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  assertthat::assert_that(phi_max > phi_min)
  info = c(unlist(material$properties), length, phi_min * pi / 180, phi_max * pi / 180)
  tibble::tibble(x = x, y = y, z = z, radius = radius, type = material$type, shape = "cylinder",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' Segment Object
#' 
#' Similar to the cylinder object, but specified by start and end points.
#'
#' @param start Default `c(0, -1, 0)`. Start point of the cylinder segment, specifing `x`, `y`, `z`.
#' @param end Default `c(0, 1, 0)`. End point of the cylinder segment, specifing `x`, `y`, `z`.
#' @param radius Default `1`. Radius of the segment.
#' @param phi_min Default `0`. Minimum angle around the segment.
#' @param phi_max Default `360`. Maximum angle around the segment.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param velocity Default `c(0, 0, 0)`. Velocity of the segment.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly. Notes: this will change the stated start/end position of the segment. 
#' Emissive objects may not currently function correctly when scaled.
#' 
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the segment in the scene.
#' @export
#'
#' @examples
#' #Generate a segment in the cornell box. 
#' \donttest{
#' generate_cornell() %>%
#'   add_object(segment(start = c(100, 100, 100), end = c(455, 455, 455), radius = 50)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#'
#' # Draw a line graph representing a normal distribution, but with metal:
#' xvals = seq(-3, 3, length.out = 30)
#' yvals = dnorm(xvals)
#' 
#' scene_list = list()
#' for(i in 1:(length(xvals) - 1)) {
#'   scene_list[[i]] = segment(start = c(555/2 + xvals[i] * 80, yvals[i] * 800, 555/2),
#'                             end = c(555/2 + xvals[i + 1] * 80, yvals[i + 1] * 800, 555/2),
#'                             radius = 10,
#'                             material = metal())
#' }
#' scene_segments = do.call(rbind,scene_list)
#' \donttest{
#' generate_cornell() %>% 
#'   add_object(scene_segments) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#'
#' #Draw the outline of a cube:
#' 
#' cube_outline = segment(start = c(100, 100, 100), end = c(100, 100, 455), radius = 10) %>%
#'   add_object(segment(start = c(100, 100, 100), end = c(100, 455, 100), radius = 10)) %>%
#'   add_object(segment(start = c(100, 100, 100), end = c(455, 100, 100), radius = 10)) %>%
#'   add_object(segment(start = c(100, 100, 455), end = c(100, 455, 455), radius = 10)) %>%
#'   add_object(segment(start = c(100, 100, 455), end = c(455, 100, 455), radius = 10)) %>%
#'   add_object(segment(start = c(100, 455, 455), end = c(100, 455, 100), radius = 10)) %>%
#'   add_object(segment(start = c(100, 455, 455), end = c(455, 455, 455), radius = 10)) %>%
#'   add_object(segment(start = c(455, 455, 100), end = c(455, 100, 100), radius = 10)) %>%
#'   add_object(segment(start = c(455, 455, 100), end = c(455, 455, 455), radius = 10)) %>%
#'   add_object(segment(start = c(455, 100, 100), end = c(455, 100, 455), radius = 10)) %>%
#'   add_object(segment(start = c(455, 100, 455), end = c(455, 455, 455), radius = 10)) %>%
#'   add_object(segment(start = c(100, 455, 100), end = c(455, 455, 100), radius = 10))
#' 
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cube_outline) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Shrink and rotate the cube
#' \donttest{
#' generate_cornell() %>%
#'   add_object(group_objects(cube_outline, pivot_point = c(555/2, 555/2, 555/2),
#'                            group_angle = c(45,45,45), group_scale = c(0.5,0.5,0.5))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
segment = function(start = c(0, -1, 0), end = c(0, 1, 0), radius = 1, 
                   phi_min = 0, phi_max = 360,
                   material = diffuse(), 
                   velocity = c(0, 0, 0), flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  assertthat::assert_that(phi_max > phi_min)
  x = (start[1] + end[1])/2
  y = (start[2] + end[2])/2
  z = (start[3] + end[3])/2
  order_rotation = c(3, 2, 1)
  phi =  atan2(end[1]-start[1], end[3]-start[3])/pi*180 + 90
  
  length_xy = sqrt((end[1]-start[1])^2 + (end[3]-start[3])^2)
  if(end[1] == start[1] && end[3] == start[3]) {
    theta = 0
  } else {
    theta = atan2(-length_xy, (end[2]-start[2]))/pi*180
  }
  fulllength = sqrt(sum((end-start)^2))
  angle = c(0, phi, theta)
  info = c(unlist(material$properties), fulllength, phi_min * pi / 180, phi_max * pi / 180)
  tibble::tibble(x = x, y = y, z = z, radius = radius, type = material$type, shape = "cylinder",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}

#' Ellipsoid Object
#'
#' Note: light importance sampling for this shape is currently approximated by a sphere. This will fail
#' for ellipsoids with large differences between axes.
#'
#' @param x Default `0`. x-coordinate of the center of the ellipsoid.
#' @param y Default `0`. y-coordinate of the center of the ellipsoid.
#' @param z Default `0`. z-coordinate of the center of the ellipsoid.
#' @param a Default `1`. Principal x-axis of the ellipsoid.
#' @param b Default `1`. Principal y-axis of the ellipsoid.
#' @param c Default `1`. Principal z-axis of the ellipsoid.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param velocity Default `c(0, 0, 0)`. Velocity of the segment.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly. Note: emissive objects may not currently function correctly when scaled.
#' 
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the ellipsoid in the scene.
#' @export
#' @examples
#' #Generate an ellipsoid in a Cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 50, c = 50)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Change the axes to make it taller rather than wide:
#' \donttest{
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 200, c = 100, material = metal())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Rotate it and make it dielectric:
#' \donttest{
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 200, c = 100, angle = c(0, 0, 45),
#'                        material = dielectric())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 400, parallel = TRUE, clamp_value = 5)
#' }
ellipsoid = function(x = 0, y = 0, z = 0, a = 1, b = 1, c = 1,
                  material = diffuse(), 
                  angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                  velocity = c(0, 0, 0), flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  radius = 1
  info = c(unlist(material$properties), a, b, c)
  tibble::tibble(x = x, y = y, z = z, radius = radius, type = material$type, shape = "ellipsoid",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA))
}
