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
  new_tibble_row(list(x = x, y = y, z = z, radius = radius, type = material$type, shape = "sphere",
                 properties = material$properties, velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image,  image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample, sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
#'                ambient_light = FALSE, samples = 500, parallel = TRUE, clamp_value = 5)
#' }
#' #Generate a gold cube in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 100, z = 555/2, 
#'                   xwidth = 200, ywidth = 200, zwidth = 200, angle = c(0, 30, 0),
#'                   material = metal(color = "gold", fuzz = 0.2))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 500, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a rotated dielectric box in the cornell box
#' \donttest{
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 200, z = 555/2, 
#'                   xwidth = 200, ywidth = 100, zwidth = 200, angle = c(30, 30, 30),
#'                   material = dielectric())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 500, parallel = TRUE, clamp_value = 5)
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
  new_tibble_row(list(x = x, y = y, z = z, radius = NA, type = material$type, shape = "box",
                 properties = list(boxinfo), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image,  image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation), 
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
  new_tibble_row(list(x = x, y = y, z = z, radius = NA, type = material$type, shape = "xy_rect",
                 properties = list(rectinfo), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
  new_tibble_row(list(x = x, y = y, z = z, radius = NA, type = material$type, shape = "yz_rect",
                 properties = list(rectinfo), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
  new_tibble_row(list(x = x, y = y, z = z, radius = NA, 
                 type = material$type, shape = "xz_rect",
                 properties = list(rectinfo), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
#' @param reversed Default `FALSE`. Similar to the `flipped` argument, but this reverses the handedness of the 
#' triangle so it will be oriented in the opposite direction.
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
                    flipped = FALSE, reversed = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  if(!reversed) {
    vertex_vec = c(v1, v2, v3)
    normal_vec = c(n1, n2, n3)
  } else {
    vertex_vec = c(v3, v2, v1)
    normal_vec = c(n3, n2, n1)
  }
  info = c(unlist(material$properties), vertex_vec, normal_vec)
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
  new_tibble_row(list(x = 0, y = 0, z = 0, radius = NA, 
                 type = material$type, shape = "triangle",
                 properties = list(info), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(colorvec), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
  new_tibble_row(list(x = x, y = y, z = z, radius = radius, type = material$type, shape = "disk",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
#' @param vertex_colors Default `FALSE`. Set to `TRUE` if the OBJ file has vertex colors to apply them
#' to the model.
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
#'   render_scene(parallel = TRUE, samples = 500, 
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
#'   render_scene(parallel = TRUE, samples = 500, ambient = TRUE, 
#'                backgroundhigh="blue", backgroundlow="red",
#'                aperture = 0.05, fov = 32, lookfrom = c(0, 2, 10),
#'                lookat = c(0,1,0)) 
#' }
obj_model = function(filename, x = 0, y = 0, z = 0, scale_obj = 1, 
                     texture = FALSE, vertex_colors = FALSE,
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
  if(vertex_colors) {
    shape = "objvertexcolor"
  }
  new_tibble_row(list(x = x, y = y, z = z, radius = NA, 
                 type = material$type, shape = shape,
                 properties = list(info), velocity = list(c(0, 0, 0)),
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = filename, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
  new_tibble_row(list(x = x, y = y, z = z, radius = radius, type = material$type, shape = "cylinder",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
#' @param direction Default `NA`. Alternative to `start` and `end`, specify the direction (via 
#' a length-3 vector) of the segment. Segment will be centered at `start`, and the length will be
#' determined by the magnitude of the direction vector.
#' @param from_center Default `TRUE`. If orientation specified via `direction`, setting this argument
#' to `FALSE` will make `start` specify the bottom of the segment, instead of the middle.
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
                   phi_min = 0, phi_max = 360, from_center = TRUE, direction = NA,
                   material = diffuse(), 
                   velocity = c(0, 0, 0), flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  assertthat::assert_that(phi_max > phi_min)
  if(all(!is.na(direction)) && length(direction) == 3) {
    if(from_center) {
      new_start = start - direction/2
      new_end = start + direction/2
    } else {
      new_start = start
      new_end = start + direction
    }
    start = new_start
    end = new_end
  } 
  x = (start[1] + end[1])/2
  y = (start[2] + end[2])/2
  z = (start[3] + end[3])/2
  order_rotation = c(3, 2, 1)
  phi =  atan2( as.numeric(end[1]-start[1]), as.numeric(end[3]-start[3]))/pi*180 + 90
  
  length_xy = sqrt((end[1]-start[1])^2 + (end[3]-start[3])^2)
  if(end[1] == start[1] && end[3] == start[3]) {
    theta = 0
  } else {
    theta = atan2(-length_xy, (end[2]-start[2]))/pi*180
  }
  fulllength = sqrt(sum((end-start)^2))
  angle = c(0, phi, theta)
  info = c(unlist(material$properties), fulllength, phi_min * pi / 180, phi_max * pi / 180)
  new_tibble_row(list(x = x, y = y, z = z, radius = radius, type = material$type, shape = "cylinder",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
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
#'                ambient_light = FALSE, samples = 500, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Change the axes to make it taller rather than wide:
#' \donttest{
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 200, c = 100, material = metal())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 500, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Rotate it and make it dielectric:
#' \donttest{
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 200, c = 100, angle = c(0, 0, 45),
#'                        material = dielectric())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 500, parallel = TRUE, clamp_value = 5)
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
  new_tibble_row(list(x = x, y = y, z = z, radius = radius, type = material$type, shape = "ellipsoid",
                 properties = list(info), velocity = list(velocity), 
                 checkercolor = material$checkercolor, 
                 gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                 world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                 gradient_type = material$gradient_type,
                 noise = material$noise, noisephase = material$noisephase, 
                 noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                 angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                 alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                 bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                 flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                 implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                 order_rotation = list(order_rotation),
                 pivot_point = list(NA), group_translate = list(NA),
                 group_angle = list(NA), group_order_rotation = list(NA),
                 tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                 material_id = NA))
}

#' Extruded Polygon Object
#'
#' @param polygon `sf` object, "SpatialPolygon" `sp` object,  or xy coordinates
#'   of polygon represented in a way that can be processed by `xy.coords()`.  If
#'   xy-coordinate based polygons are open, they will be closed by adding an
#'   edge from the last point to the first.
#' @param x Default `0`. x-coordinate to offset the extruded model.
#' @param y Default `0`. y-coordinate to offset the extruded model.
#' @param z Default `0`. z-coordinate to offset the extruded model.
#' @param plane Default `xz`. The plane the polygon is drawn in. All possibile orientations
#'  are `xz`, `zx`, `xy`, `yx`, `yz`, and `zy`.
#' @param top Default `1`. Extruded top distance. If this equals `bottom`, the polygon will not be
#' extruded and just the one side will be rendered.
#' @param bottom Default `0`. Extruded bottom distance. If this equals `top`, the polygon will not be
#' extruded and just the one side will be rendered.
#' @param holes Default `0`. If passing in a polygon directly, this specifies which index represents
#' the holes in the polygon. See the `earcut` function in the `decido` package for more information.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}. 
#' @param center Default `FALSE`. Whether to center the polygon at the origin.
#' @param flip_horizontal Default `FALSE`. Flip polygon horizontally in the plane defined by `plane`.
#' @param flip_vertical Default `FALSE`. Flip polygon vertically in the plane defined by `plane`.
#' @param data_column_top Default `NULL`. A string indicating the column in the `sf` object to use 
#' to specify the top of the extruded polygon.
#' @param data_column_bottom Default `NULL`. A string indicating the column in the `sf` object to use 
#' to specify the bottom of the extruded polygon.
#' @param scale_data Default `1`. If specifying `data_column_top` or `data_column_bottom`, how
#' much to scale that value when rendering.
#' @param material_id Default `NA`. A unique label/number to ensure the material is shared between
#' all triangles that make up the extruded polygon. Required if the material is `dielectric()`.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param pivot_point Default `c(0,0,0)`. Point at which to rotate the polygon around.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Multiple row tibble describing the extruded polygon in the scene.
#' @export
#'
#' @examples
#' #Manually create a polygon object, here a star:
#' 
#' angles = seq(0,360,by=36)
#' xx = rev(c(rep(c(1,0.5),5),1) * sinpi(angles/180))
#' yy = rev(c(rep(c(1,0.5),5),1) * cospi(angles/180))
#' star_polygon = data.frame(x=xx,y=yy)
#' 
#' \donttest{
#' generate_ground(depth=0,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(star_polygon,top=0.5,bottom=0,
#'                               material=diffuse(color="red",sigma=90))) %>%
#'   add_object(sphere(y=4,x=-3,z=-3,material=light(intensity=30))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,2,3),samples=400,lookat=c(0,0.5,0),fov=60)
#' }
#' 
#' #Now, let's add a hole to the center of the polygon. We'll make the polygon
#' #hollow by shrinking it, combining it with the normal size polygon,
#' #and specify with the `holes` argument that everything after `nrow(star_polygon)`
#' #in the following should be used to draw a hole:
#' 
#' hollow_star = rbind(star_polygon,0.8*star_polygon)
#' 
#' \donttest{
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon) + 1,
#'                               material=diffuse(color="red",sigma=90))) %>%
#'   add_object(sphere(y=4,x=-3,z=-3,material=light(intensity=30))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,2,4),samples=400,lookat=c(0,0,0),fov=30)
#' }
#' 
#' # Render one in the y-x plane as well by changing the `plane` argument,
#' # as well as offset it slightly.
#' \donttest{
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon),
#'                               material=diffuse(color="red",sigma=90))) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, y=1.2, z=-1.2, 
#'                               holes = nrow(star_polygon) + 1, plane = "yx", 
#'                               material=diffuse(color="green",sigma=90))) %>%
#'   add_object(sphere(y=4,x=-3,material=light(intensity=30))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,2,4),samples=400,lookat=c(0,0.9,0),fov=40)
#' }
#' 
#' # Now add the zy plane:
#' \donttest{
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon) + 1,
#'                               material=diffuse(color="red",sigma=90))) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, y=1.2, z=-1.2, 
#'                               holes = nrow(star_polygon) + 1, plane = "yx", 
#'                               material=diffuse(color="green",sigma=90))) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, y=1.2, x=1.2, 
#'                               holes = nrow(star_polygon) + 1, plane = "zy", 
#'                               material=diffuse(color="blue",sigma=90))) %>%
#'   add_object(sphere(y=4,x=-3,material=light(intensity=30))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(-4,2,4),samples=400,lookat=c(0,0.9,0),fov=40)
#' }
#' 
#' #We can also directly pass in sf polygons:
#' if("spData" %in% rownames(utils::installed.packages())) {
#'   us_states = spData::us_states
#'   texas = us_states[us_states$NAME == "Texas",]
#'   #Fix no sfc class in us_states geometry data
#'   class(texas$geometry) = c("list","sfc")
#' }
#' 
#' #This uses the raw coordinates, unless `center = TRUE`, which centers the bounding box
#' #of the polygon at the origin.
#' \donttest{
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(texas, center = TRUE,
#'                               material=diffuse(color="#ff2222",sigma=90))) %>%
#'   add_object(sphere(y=30,x=-30,radius=10,
#'                     material=light(color="lightblue",intensity=40))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,10,-10),samples=400,fov=60)
#' }
#' 
#' #Here we use the raw coordinates, but offset the polygon manually.
#' \donttest{
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(us_states, x=-96,z=-40, top=2,
#'                               material=diffuse(color="#ff2222",sigma=90))) %>%
#'   add_object(sphere(y=30,x=-100,radius=10,
#'                     material=light(color="lightblue",intensity=200))) %>%
#'   add_object(sphere(y=30,x=100,radius=10,
#'                     material=light(color="orange",intensity=200))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,120,-120),samples=400,fov=20)
#' }
#' 
#' #We can also set the map the height of each polygon to a column in the sf object,
#' #scaling it down by the maximum population state.
#' 
#' \donttest{
#' generate_ground(depth=0,
#'                 material = diffuse(color="grey50",checkercolor="grey20",sigma=90)) %>%
#'   add_object(extruded_polygon(us_states, x=-96,z=-45, data_column_top = "total_pop_15",
#'                               scale_data = 1/max(us_states$total_pop_15)*5,
#'                               material=diffuse(color="#ff2222",sigma=90))) %>%
#'   add_object(sphere(y=30,x=-100,z=60,radius=10,
#'                     material=light(color="lightblue",intensity=250))) %>%
#'   add_object(sphere(y=30,x=100,z=-60,radius=10,
#'                     material=light(color="orange",intensity=250))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(-60,50,-40),lookat=c(0,-5,0),samples=400,fov=30)
#' }
#' 
extruded_polygon = function(polygon = NULL, x = 0, y = 0, z = 0, plane = "xz",
                   top = 1, bottom = 0, holes = NULL, 
                   angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   pivot_point = c(0,0,0), material = diffuse(),
                   center = FALSE, flip_horizontal = FALSE, flip_vertical = FALSE,
                   data_column_top = NULL, data_column_bottom = NULL, scale_data = 1,
                   scale = c(1,1,1), material_id = NA) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  reversed = FALSE
  if(!is.character(plane) || length(plane) != 1L || is.na(plane)) {
    stop("plane must be scalar character and not NA")
  }
  if(tolower(plane) %in% c("xy","yx","xz","zx","zy","yz")) {
    if(plane == "xz") {
      planeval = 1
      reversed = TRUE
    } else if(plane == "zx") {
      planeval = 2
    } else if(plane == "xy") {
      planeval = 3
      reversed = TRUE
    } else if(plane == "yx") {
      planeval = 4
    } else if(plane == "zy") {
      planeval = 5
      reversed = TRUE
    } else if(plane == "yz") {
      planeval = 6
    }
  } else {
    stop("Plane ", plane, " not recognized.")
  }
  if(material$type == "dielectric" && is.na(material_id)) {
    stop("If an extruded polygon has a dielectric material, user must supply a unique material_id")
  }
  rot_coords = function(x1,x2,theta) {
    cos_theta = cospi(theta/180)
    sin_theta = sinpi(theta/180)
    return(c(cos_theta*x1 + sin_theta*x2,-sin_theta*x1 + cos_theta*x2))
  }
  permute_axes = function(x,plane_val) {
    if(plane_val == 1) {
      return(x)
    }
    if(plane_val == 2) {
      return(x[c(3,2,1)])
    }
    if(plane_val == 3) {
      return(x[c(3,1,2)])
    }
    if(plane_val == 4) {
      return(x[c(1,3,2)])
    }
    if(plane_val == 5) {
      return(x[c(2,3,1)])
    }
    if(plane_val == 6) {
      return(x[c(2,1,3)])
    }
  }
  proplen = length(material$properties[[1]])
  
  if(top == bottom) {
    extruded = FALSE
  } else {
    extruded = TRUE
  }
  x_off = x
  y_off = y
  z_off = z
  poly_list = list()
  vertex_list = list()
  height_list = list()
  bottom_list = list()
  holes_start_i_list = list()
  base_poly = FALSE
  counter = 1
  if(inherits(polygon,"sf")) {
    if(!"sf" %in% rownames(utils::installed.packages())) {
      stop("sf package required when handling sf objects")
    }
    poly_info = sf::st_drop_geometry(polygon)
    polygon = sf::as_Spatial(polygon)
    if(!is.null(data_column_top)) {
      if(data_column_top %in% colnames(poly_info)) {
        data_vals_top = poly_info[[data_column_top]] * scale_data
      } else {
        warning("Was not able to find data_column_top `", data_column_top, "` in sf object.")
        data_vals_top = rep(top, nrow(poly_info))
      }
    } else {
      data_vals_top = rep(top, nrow(poly_info))
    }
    if(!is.null(data_column_bottom)) {
      if(data_column_bottom %in% colnames(poly_info)) {
        data_vals_bottom = poly_info[[data_column_bottom]] * scale_data
      } else {
        warning("Was not able to find data_column_top `", data_column_bottom, "` in sf object.")
        data_vals_bottom = rep(bottom, nrow(poly_info))
      }
    } else {
      data_vals_bottom = rep(bottom, nrow(poly_info))
    }
  } else if(inherits(polygon,"SpatialPolygonsDataFrame") || inherits(polygon,"SpatialPolygons")) {
    data_vals_top = rep(top, nrow(polygon))
    data_vals_bottom = rep(bottom, nrow(polygon))
  } else {
    xylist = grDevices::xy.coords(polygon)
    data_vals_top = top[1]
    data_vals_bottom = bottom[1]
  }
  if(inherits(polygon,"SpatialPolygonsDataFrame") || inherits(polygon,"SpatialPolygons")) {
    if(!is.null(holes)) {
      warning("holes is not NULL, but is unused when input is sf or Spatial")
    }
    coord_data = raster::geom(polygon) # See ?raster::geom for col meaning
    coord_data = coord_data[order(coord_data[,1],coord_data[,2] ),] # just in case

    # Holes are only holes if every hole value of `part` column (col 2) is a
    # hole to match original code (seems like it should always be true?).

    opart = interaction(coord_data[, 1], coord_data[, 2], drop=TRUE)
    coord_data[, 4] = stats::ave(coord_data[, 4], opart, FUN=min)

    # Every `part` is considered a polygon, unless it is a hole, in which case
    # it is assumed to be a hole in the nearest preceding non-hole polygon.

    object = coord_data[, 1]
    object_new = coord_data[, 3]
    object_new[coord_data[, 4] == 1] = min(coord_data[,3]) - 1L
    object_new = cummax(object_new)
    coord_data[, 1] = object_new
    coord_data[, 5] = -coord_data[, 5]  # due to xz defalt plane need to negate x

    # give distinct hole ids instead of just 0,1 flag

    coord_data[, 4] = coord_data[, 4] * coord_data[, 3]

    # remap the new object id to old to get height data

    obj_new_len = rle(object_new)[['lengths']]
    obj_map_id = object[c(1L, cumsum(obj_new_len)[-length(obj_new_len)] + 1L)]

    height_list = as.list(data_vals_top[obj_map_id])    # list for compatibility
    bottom_list = as.list(data_vals_bottom[obj_map_id])
    poly_list = unname(split.data.frame(coord_data[, c(5,6,4,2)], coord_data[,1]))

    # assume that all holes are after all outer polygons

    holes_start_i_list = lapply(
      poly_list,
      function(x) {
        if(!any(x[,3] > 0)) {              # no hole
          0
        } else {
          which(
            c(TRUE, diff(x[,4]) > 0) &     # start position of each "part"
            x[,3] > 0                      # that correspond to holes
          )
    } } )
  } else {
    base_poly = TRUE
    xylist = grDevices::xy.coords(polygon)
    x = xylist$x
    y = xylist$y
    if(is.null(holes)) {
      holes = 0
    } else if (
      !is.numeric(holes) || anyNA(holes <- as.integer(holes))
    ) {
      stop("holes must be integer")
    } else if (
      any(holes < 0L) || any(holes) > length(x) ||
      (any(holes == 0L) && length(holes) != 1L)
    ) {
      stop(
        "holes must be zero, or contain indices to polygon vertices"
      )
    } else if (any(holes < 4)) {
      stop("holes cannot begin before vertex 4. Hole index here starts at: ", min(holes))
    }
    # label each vertex with a hole id

    if(isTRUE(holes == 0)) {
      xy_dat = data.frame(x, y, holes=0L)
    } else {
      xy_dat = data.frame(x, y, holes=cumsum(seq_along(x) %in% holes))
    }
    # close polygons if not closed, must do so for outer and each hole;
    # sf and Spatial polygons should be closed

    xy_dat_split = split(xy_dat, xy_dat[['holes']])
    close_poly = function(dat) {
      if(!all(dat[1L,] == dat[nrow(dat),])) {
        dat[c(seq_len(nrow(dat)), 1L),]
      } else dat
    }
    xy_dat_closed = lapply(xy_dat_split, close_poly)
    xy_dat_len = vapply(xy_dat_closed, nrow, 0)
    holes = if(length(xy_dat_closed) > 1) {
      cumsum(xy_dat_len[-length(xy_dat_closed)]) + 1L
    } else 0L
    xy_dat_fin = do.call(rbind, xy_dat_closed)
    rownames(xy_dat_fin) = NULL

    holes_start_i_list[[1]] = holes  # hole indices for decido::earcut
    poly_list[[1]] = as.matrix(xy_dat_fin)
    height_list[[1]] = data_vals_top
    bottom_list[[1]] = data_vals_bottom
  }
  # Processing common to base and SF: flip/earcut/center

  for(i in seq_along(poly_list)) {
    if(flip_horizontal) {
      poly_list[[i]][,1] = -poly_list[[i]][,1]
    }
    if(flip_vertical) {
      poly_list[[i]][,2] = -poly_list[[i]][,2]
    }
    vertex_list[[i]] = matrix(
      decido::earcut(poly_list[[i]][,1:2],holes = holes_start_i_list[[i]]),
      ncol=3, byrow=TRUE
    )
  }
  if(center) {
    all_vertices = do.call(rbind,poly_list)
    middle_x = (max(all_vertices[,1]) + min(all_vertices[,1]))/2
    middle_y = (max(all_vertices[,2]) + min(all_vertices[,2]))/2
    for(i in seq_along(poly_list)) {
      poly_list[[i]][,1] = poly_list[[i]][,1] - middle_x
      poly_list[[i]][,2] = poly_list[[i]][,2] - middle_y
    }
  }
  scenelist= list()
  counter = 1

  for(poly in 1:length(poly_list)) {
    x=poly_list[[poly]][,1]
    y=poly_list[[poly]][,2]
    vertices = vertex_list[[poly]]
    height_poly = height_list[[poly]]
    bottom_poly = bottom_list[[poly]]
    
    for(i in 1:nrow(vertices)) {
      scenelist[[counter]] = triangle(v1=scale*permute_axes(c(x[vertices[i,3]],bottom_poly,y[vertices[i,3]]),planeval),
                                      v2=scale*permute_axes(c(x[vertices[i,2]],bottom_poly,y[vertices[i,2]]),planeval),
                                      v3=scale*permute_axes(c(x[vertices[i,1]],bottom_poly,y[vertices[i,1]]),planeval),
                                      material = material, reversed = reversed)
      counter = counter + 1
      if(extruded) {
        scenelist[[counter]] = triangle(v1=scale*permute_axes(c(x[vertices[i,1]],height_poly,y[vertices[i,1]]),planeval),
                                        v2=scale*permute_axes(c(x[vertices[i,2]],height_poly,y[vertices[i,2]]),planeval),
                                        v3=scale*permute_axes(c(x[vertices[i,3]],height_poly,y[vertices[i,3]]),planeval),
                                        material = material, reversed = reversed)
        counter = counter + 1
      }
    }
    holes = poly_list[[poly]][,3]
    hole_n = table(holes)
    if(hole_n[1L] < 4L) {
      stop(
        "outer polygon for polygon ", poly, " only has ", hole_n[1L],
        " vertices including the closing vertex; needs at least 4."
      )
    } else if (any(hole_n[-1L] < 4L)) {
      stop(
        "hole polygon(s) ",
        paste0(names(hole_n[-1L])[hole_n[-1L] < 4L], collapse=", "),
        " for poly ", poly, " do not have at least 4 vertices ",
        "including the closing vertex."
      )
    }

    if(extruded) {
      for(polyv in split(seq_along(x), holes)) {
        # Find which direction polygon is wound by computing signed area,
        # assumes non-intersecting polygon (side is a closed polygon).  CW
        # outer polygon need to flip sides, as do CCW holes

        i = polyv[seq_len(length(polyv) - 1L)]
        ii = polyv[seq_len(length(polyv) - 1L) + 1L]   # i + 1
        area_s = sum(x[i] * y[ii] - x[ii] * y[i]) / 2

        ccw = area_s >= 0  # treat degenerates as counter-clockwise
        side_rev = (
          (holes[polyv[1L]] == 0L && !ccw) ||
          (holes[polyv[1L]] != 0L && ccw)
        )
        for(i in seq_len(length(polyv) - 1L)) {  # polygons are closed
          xi = x[polyv[i]]        # vertex i
          yi = y[polyv[i]]
          xii = x[polyv[i + 1L]]  # vertex i + 1
          yii = y[polyv[i + 1L]]

          scenelist[[counter]] = triangle(v1=scale*permute_axes(c(xi,height_poly,yi),planeval),
                                          v2=scale*permute_axes(c(xi,bottom_poly,yi),planeval),
                                          v3=scale*permute_axes(c(xii,bottom_poly,yii),planeval),
                                          material = material, reversed = xor(reversed, side_rev))
          counter = counter + 1
          scenelist[[counter]] = triangle(v1=scale*permute_axes(c(xi,height_poly,yi),planeval),
                                          v2=scale*permute_axes(c(xii,bottom_poly,yii),planeval),
                                          v3=scale*permute_axes(c(xii,height_poly,yii),planeval),
                                          material = material, reversed = xor(reversed, side_rev))
          counter = counter + 1
        }
      }
    }
  }
  if("dplyr" %in% rownames(utils::installed.packages())) {
    scenefull = dplyr::bind_rows(scenelist)
  } else {
    scenefull = do.call(rbind, scenelist)
  }
  if(any(angle != 0)) {
    if(any(pivot_point != 0)) {
      sceneprop = scenefull$properties
      add_at_indices = function(x, indices, off) {
        x[indices] = x[indices] + off
        x
      }
      sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 1,proplen + 4, proplen + 7), off = -pivot_point[1])
      sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 2,proplen + 5, proplen + 8), off = -pivot_point[2])
      sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 3,proplen + 6, proplen + 9), off = -pivot_point[3])
      scenefull$properties = sceneprop
    }
    rot_at_indices = function(x, indices, angle) {
      x[indices[1:2]] = rot_coords(x[indices[1]], x[indices[2]], angle)
      x[indices[3:4]] = rot_coords(x[indices[3]], x[indices[4]], angle)
      x[indices[5:6]] = rot_coords(x[indices[5]], x[indices[6]], angle)
      x
    }
    for(i in 1:3) {
      if(order_rotation[i] == 1) {
        sceneprop = scenefull$properties
        sceneprop = lapply(sceneprop, rot_at_indices, indices = proplen + c(2, 3, 5, 6, 8, 9), angle = angle[1])
        scenefull$properties = sceneprop
      }
      if(order_rotation[i] == 2) {
        sceneprop = scenefull$properties
        sceneprop = lapply(sceneprop, rot_at_indices, indices = proplen + c(1, 3, 4, 6, 7, 9), angle = angle[2])
        scenefull$properties = sceneprop
      }
      if(order_rotation[i] == 3) {
        sceneprop = scenefull$properties
        sceneprop = lapply(sceneprop, rot_at_indices, indices = proplen + c(1, 2, 4, 5, 7, 8), angle = angle[3])
        scenefull$properties = sceneprop
      }
    }
    if(any(pivot_point != 0)) {
      sceneprop = scenefull$properties
      add_at_indices = function(x, indices, off) {
        x[indices] = x[indices] + off
        x
      }
      sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 1,proplen + 4, proplen + 7), off = pivot_point[1])
      sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 2,proplen + 5, proplen + 8), off = pivot_point[2])
      sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 3,proplen + 6, proplen + 9), off = pivot_point[3])
      scenefull$properties = sceneprop
    }
  }
  if(any(x_off != 0 || y_off != 0 || z_off != 0)) {
    sceneprop = scenefull$properties
    add_at_indices = function(x, indices, off) {
      x[indices] = x[indices] + off
      x
    }
    sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 1,proplen + 4, proplen + 7), off = x_off)
    sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 2,proplen + 5, proplen + 8), off = y_off)
    sceneprop = lapply(sceneprop,add_at_indices,indices=c(proplen + 3,proplen + 6, proplen + 9), off = z_off)
    scenefull$properties = sceneprop
  }
  scenefull$material_id = rep(material_id, nrow(scenefull))
  tibble::validate_tibble(scenefull)
  return(scenefull)
}

#' Cone Object
#' 
#'
#' @param start Default `c(0, 0, 0)`. Base of the cone, specifying `x`, `y`, `z`.
#' @param end Default `c(0, 1, 0)`. Tip of the cone, specifying `x`, `y`, `z`.
#' @param radius Default `1`. Radius of the bottom of the cone.
#' @param direction Default `NA`. Alternative to `start` and `end`, specify the direction (via 
#' a length-3 vector) of the cone. Cone will be centered at `start`, and the length will be
#' determined by the magnitude of the direction vector.
#' @param from_center Default `TRUE`. If orientation specified via `direction`, setting this argument
#' to `FALSE` will make `start` specify the bottom of the cone, instead of the middle.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Rotation angle. Note: This will change the `start` and `end` coordinates.
#' @param velocity Default `c(0, 0, 0)`. Velocity of the segment.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly. Notes: this will change the stated start/end position of the cone. 
#' Emissive objects may not currently function correctly when scaled.
#' 
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the cone in the scene.
#' @export
#'
#' @examples
#' #Generate a cone in a studio, pointing upwards:
#' \donttest{
#' generate_studio() %>% 
#'  add_object(cone(start=c(0,-1,0), end=c(0,1,0), radius=1,material=diffuse(color="red"))) %>% 
#'  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'  render_scene(samples=400,clamp_value=10)
#'  
#'  #Change the radius, length, and direction
#' generate_studio() %>% 
#'  add_object(cone(start=c(0,0,0), end=c(0,-1,0), radius=0.5,material=diffuse(color="red"))) %>% 
#'  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'  render_scene(samples=400,clamp_value=10)
#'  
#' #Give custom start and end points (and customize the color/texture)
#' generate_studio() %>% 
#'  add_object(cone(start=c(-1,0.5,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="red"))) %>%
#'  add_object(cone(start=c(1,0.5,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="green"))) %>%
#'  add_object(cone(start=c(0,1,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="orange"))) %>% 
#'  add_object(cone(start=c(-1,-0.5,0), end=c(1,-0.5,0), radius=0.25,
#'    material = diffuse(color="red",gradient_color="green"))) %>% 
#'  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'  render_scene(samples=400,clamp_value=10)
#'  
#'  
#' #Specify cone via direction and location, instead of start and end positions
#' #Length is derived from the magnitude of the direction.
#' gold_mat = microfacet(roughness=0.1,eta=c(0.216,0.42833,1.3184), kappa=c(3.239,2.4599,1.8661))
#' generate_studio() %>% 
#'   add_object(cone(start = c(-1,0,0), direction = c(-0.5,0.5,0), material = gold_mat)) %>% 
#'   add_object(cone(start = c(1,0,0), direction = c(0.5,0.5,0), material = gold_mat)) %>% 
#'   add_object(cone(start = c(0,0,-1), direction = c(0,0.5,-0.5), material = gold_mat)) %>% 
#'   add_object(cone(start = c(0,0,1), direction = c(0,0.5,0.5), material = gold_mat)) %>% 
#'   add_object(sphere(y=5,material=light())) %>% 
#'   add_object(sphere(y=3,x=-3,z=-3,material=light(color="red"))) %>% 
#'   add_object(sphere(y=3,x=3,z=-3,material=light(color="green"))) %>% 
#'   render_scene(lookfrom=c(0,4,10), clamp_value=10, samples=400)
#'   
#'  #Render the position from the base, instead of the center of the cone:
#'  noise_mat = material = glossy(color="purple",noisecolor="blue", noise=5)
#'  generate_studio() %>% 
#'   add_object(cone(start = c(0,-1,0), from_center = FALSE, radius=1, direction = c(0,2,0), 
#'     material = noise_mat)) %>% 
#'   add_object(cone(start = c(-1.5,-1,0), from_center = FALSE, radius=0.5, direction = c(0,1,0), 
#'     material = noise_mat)) %>% 
#'   add_object(cone(start = c(1.5,-1,0), from_center = FALSE, radius=0.5, direction = c(0,1,0), 
#'     material = noise_mat)) %>% 
#'   add_object(cone(start = c(0,-1,1.5), from_center = FALSE, radius=0.5, direction = c(0,1,0), 
#'     material = noise_mat)) %>% 
#'   add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'   render_scene(lookfrom=c(0,4,10), clamp_value=10,fov=25, samples=400)
#'   
#' }
#' 
cone = function(start = c(0, 0, 0), end = c(0, 1, 0), radius = 0.5, 
                direction = NA, from_center = TRUE,
                material = diffuse(), angle = c(0,0,0),
                velocity = c(0, 0, 0), flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  if(!is.na(material$lightintensity) && material$implicit_sample) {
    message("Importance sampling not supported with cone lights: turning off")
    material$implicit_sample = FALSE
  }
  if(!is.na(material$alphaimage)) {
    message("Alpha transparency currently not supported with cones: turning off")
    material$alphaimage = NA
  }
  if(all(!is.na(direction)) && length(direction) == 3) {
    if(from_center) {
      new_start = start - direction/2
      new_end = start + direction/2
    } else {
      new_start = start
      new_end = start + direction
    }
    start = new_start
    end = new_end
  }
  x = start[1] 
  y = start[2] 
  z = start[3] 
  order_rotation = c(3, 2, 1)
  phi =  atan2( as.numeric(end[1]-start[1]), as.numeric(end[3]-start[3]))/pi*180 + 90
  
  length_xy = sqrt((end[1]-start[1])^2 + (end[3]-start[3])^2)
  if(end[1] == start[1] && end[3] == start[3]) {
    if((end[2]-start[2]) < 0) {
      theta = 180
    } else {
      theta = 0
    }
  } else {
    theta = atan2(-length_xy, (end[2]-start[2]))/pi*180
  }
  fulllength = sqrt(sum((end-start)^2))
  angle = c(0, phi, theta) + angle
  info = c(unlist(material$properties), fulllength)
  new_tibble_row(list(x = x, y = y, z = z, radius = radius, type = material$type, shape = "cone",
                      properties = list(info), velocity = list(velocity), 
                      checkercolor = material$checkercolor, 
                      gradient_color = material$gradient_color, gradient_transpose = material$gradient_transpose, 
                      world_gradient = material$world_gradient, gradient_point_info = material$gradient_point_info,
                      gradient_type = material$gradient_type,
                      noise = material$noise, noisephase = material$noisephase, 
                      noiseintensity = material$noiseintensity, noisecolor = material$noisecolor,
                      angle = list(angle), image = material$image, image_repeat = material$image_repeat,
                      alphaimage = list(material$alphaimage), bump_texture = list(material$bump_texture),
                      bump_intensity = material$bump_intensity, lightintensity = material$lightintensity,
                      flipped = flipped, fog = material$fog, fogdensity = material$fogdensity,
                      implicit_sample = material$implicit_sample,  sigma = material$sigma, glossyinfo = material$glossyinfo,
                      order_rotation = list(order_rotation),
                      pivot_point = list(NA), group_translate = list(NA),
                      group_angle = list(NA), group_order_rotation = list(NA),
                      tricolorinfo = list(NA), fileinfo = NA, scale_factor = list(scale), group_scale = list(NA),
                      material_id = NA))
}

#' Arrow Object
#' 
#' Composite object (cone + segment)
#'
#' @param start Default `c(0, 0, 0)`. Base of the arrow, specifying `x`, `y`, `z`.
#' @param end Default `c(0, 1, 0)`. Tip of the arrow, specifying `x`, `y`, `z`.
#' @param radius_top Default `0.5`. Radius of the top of the arrow.
#' @param radius_tail Default `0.2`.  Radius of the tail of the arrow.
#' @param tail_proportion Default `0.5`. Proportion of the arrow that is the tail.
#' @param direction Default `NA`. Alternative to `start` and `end`, specify the direction (via 
#' a length-3 vector) of the arrow. Arrow will be centered at `start`, and the length will be
#' determined by the magnitude of the direction vector.
#' @param from_center Default `TRUE`. If orientation specified via `direction`, setting this argument
#' to `FALSE` will make `start` specify the bottom of the cone, instead of the middle.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param velocity Default `c(0, 0, 0)`. Velocity of the segment.
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly. Notes: this will change the stated start/end position of the cone. 
#' Emissive objects may not currently function correctly when scaled.
#' 
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the cone in the scene.
#' @export
#'
#' @examples
#' #Draw a simple arrow from x = -1 to x = 1
#' \donttest{
#' generate_studio() %>% 
#'   add_object(arrow(start = c(-1,0,0), end = c(1,0,0), material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,material=light(intensity=20))) %>% 
#'   render_scene(clamp_value=10,  samples=400)
#' 
#' #Change the proportion of tail to top
#' generate_studio(depth=-2) %>% 
#'   add_object(arrow(start = c(-1,-1,0), end = c(1,-1,0), tail_proportion = 0.5,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(arrow(start = c(-1,0,0), end = c(1,0,0), tail_proportion = 0.75,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(arrow(start = c(-1,1,0), end = c(1,1,0), tail_proportion = 0.9,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,z=5,x=2,material=light(intensity=30))) %>% 
#'   render_scene(clamp_value=10, fov=25,  samples=400)
#'   
#' #Change the radius of the tail/top segments
#' generate_studio(depth=-1.5) %>% 
#'   add_object(arrow(start = c(-1,-1,0), end = c(1,-1,0), tail_proportion = 0.75,
#'                    radius_top = 0.1, radius_tail=0.03,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(arrow(start = c(-1,0,0), end = c(1,0,0), tail_proportion = 0.75,
#'                    radius_top = 0.2, radius_tail=0.1,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(arrow(start = c(-1,1,0), end = c(1,1,0), tail_proportion = 0.75,
#'                    radius_top = 0.3, radius_tail=0.2,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,z=5,x=2,material=light(intensity=30))) %>% 
#'   render_scene(clamp_value=10, samples=400)
#'   
#'   
#' #We can also specify arrows via a midpoint and direction:
#' generate_studio(depth=-1) %>% 
#'   add_object(arrow(start = c(-1,-0.5,0), direction = c(0,0,1),
#'                    material=glossy(color="green"))) %>% 
#'   add_object(arrow(start = c(1,-0.5,0), direction = c(0,0,-1),
#'                    material=glossy(color="red"))) %>% 
#'   add_object(arrow(start = c(0,-0.5,1), direction = c(1,0,0),
#'                    material=glossy(color="yellow"))) %>% 
#'   add_object(arrow(start = c(0,-0.5,-1), direction = c(-1,0,0),
#'                    material=glossy(color="purple"))) %>% 
#'   add_object(sphere(y=5,z=5,x=2,material=light(intensity=30))) %>% 
#'   render_scene(clamp_value=10, samples=400, 
#'                lookfrom=c(0,5,10), lookat=c(0,-0.5,0), fov=16)
#' 
#' #Plot a 3D vector field for a gravitational well:
#' 
#' r = 1.5
#' theta_vals = seq(0,2*pi,length.out = 16)[-16]
#' phi_vals = seq(0,pi,length.out = 16)[-16][-1]
#' arrow_list = list()
#' counter = 1
#' for(theta in theta_vals) {
#'   for(phi in phi_vals) {
#'     rval = c(r*sin(phi)*cos(theta),r*cos(phi),r*sin(phi)*sin(theta)) 
#'     arrow_list[[counter]] = arrow(rval, direction = -1/2*rval/sqrt(sum(rval*rval))^3,
#'                                   tail_proportion = 0.66, radius_top=0.03, radius_tail=0.01,
#'                                   material = diffuse(color="red"))
#'     counter = counter + 1
#'   }
#' }
#' vector_field = do.call(rbind,arrow_list)
#' sphere(material=diffuse(noise=1,color="blue",noisecolor="darkgreen")) %>% 
#'   add_object(vector_field) %>% 
#'   add_object(sphere(y=0,x=10,z=5,material=light(intensity=200))) %>% 
#'   render_scene(fov=20, ambient=TRUE, samples=400,
#'                backgroundlow="black",backgroundhigh="white")
#' }
arrow = function(start = c(0,0,0), end = c(0,1,0), 
                 radius_top = 0.2, radius_tail = 0.1, tail_proportion = 0.5,
                 direction = NA,  from_center = TRUE, 
                 material = diffuse(), velocity = c(0,0,0),
                 flipped = FALSE, scale = c(1,1,1)) {
  stopifnot(tail_proportion > 0 && tail_proportion < 1)
  head_proportion = 1 - tail_proportion
  from_direction = TRUE
  if(any(is.na(direction)) || length(direction) != 3) {
    direction = end - start
    from_direction = FALSE
  }
  if(from_center && from_direction) {
    start = start - direction/2
  }
  start_segment = start
  end_segment = start + direction * tail_proportion 
  start_tip = end_segment - direction/100
  end_tip = start + direction + direction/100
  segment(start = start_segment, 
          end = end_segment,
          radius = radius_tail, 
          material = material, 
          flipped = flipped, scale = scale, velocity = velocity) %>%
    add_object(cone(start = start_tip, 
                    end = end_tip,
                    radius = radius_top, 
                    material = material, 
                    flipped = flipped, scale = scale, velocity = velocity))
}