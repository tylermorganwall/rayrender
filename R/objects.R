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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(sphere(x = 555/2, y = 555/2, z = 555/2, radius = 100)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, clamp_value = 5)
#' }
#' 
#' #Generate a gold sphere in the cornell box
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(sphere(x = 555/2, y = 100, z = 555/2, radius = 100, 
#'                     material = microfacet(color = "gold"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, clamp_value = 5)
#' }
sphere = function(x = 0, y = 0, z = 0, radius = 1, material = diffuse(), 
                  angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                  flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  new_tibble_row(list(x = x, y = y, z = z, shape = "sphere",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(radius = radius),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1))
                 )
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 100, z = 555/2, 
#'                   xwidth = 200, ywidth = 200, zwidth = 200, angle = c(0, 30, 0))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' #Generate a gold cube in the cornell box
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 100, z = 555/2, 
#'                   xwidth = 200, ywidth = 200, zwidth = 200, angle = c(0, 30, 0),
#'                   material = metal(color = "gold", fuzz = 0.2))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a rotated dielectric box in the cornell box
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(cube(x = 555/2, y = 200, z = 555/2, 
#'                   xwidth = 200, ywidth = 100, zwidth = 200, angle = c(-30, 30, -30),
#'                   material = dielectric())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40,  
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5) 
#' }
cube = function(x = 0, y = 0, z = 0, width = 1, xwidth = 1, ywidth = 1, zwidth = 1, 
                material = diffuse(), angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  xwidth = ifelse(missing(xwidth), width, xwidth)
  ywidth = ifelse(missing(ywidth), width, ywidth)
  zwidth = ifelse(missing(zwidth), width, zwidth)
  boxinfo = c(xwidth, ywidth, zwidth)

  new_tibble_row(list(x = x, y = y, z = z, shape = "box",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(boxinfo = boxinfo),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(xy_rect(x = 555/2, y = 100, z = 555/2, xwidth = 200, ywidth = 200,
#'              material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800), lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a gold plane in the cornell box
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(xy_rect(x = 555/2, y = 100, z = 555/2, 
#'                      xwidth = 200, ywidth = 200, angle = c(0, 30, 0),
#'                      material = metal(color = "gold"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
xy_rect = function(x = 0, y = 0, z = 0, xwidth = 1, ywidth = 1,  
                   material = diffuse(), angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  rectinfo = c(xwidth, ywidth)

  new_tibble_row(list(x = x, y = y, z = z, shape = "xy_rect",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(rectinfo = rectinfo),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(yz_rect(x = 100, y = 100, z = 555/2, ywidth = 200, zwidth = 200,
#'                      material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' #Generate a gold plane in the cornell box
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(yz_rect(x = 100, y = 100, z = 555/2, 
#'                      ywidth = 200, zwidth = 200, angle = c(0, 30, 0),
#'                      material = metal(color = "gold"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
yz_rect = function(x = 0, y = 0, z = 0, ywidth = 1, zwidth = 1, material = diffuse(), 
                   angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  rectinfo = c(ywidth, zwidth)
  
  new_tibble_row(list(x = x, y = y, z = z, shape = "yz_rect",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(rectinfo = rectinfo),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(xz_rect(x = 555/2, y = 100, z = 555/2, xwidth = 200, zwidth = 200,
#'              material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Generate a gold plane in the cornell box
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(xz_rect(x = 555/2, y = 100, z = 555/2, 
#'              xwidth = 200, zwidth = 200, angle = c(0, 30, 0),
#'              material = metal(color = "gold"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
xz_rect = function(x = 0, xwidth = 1, z = 0, zwidth = 1, y = 0, material = diffuse(), 
                   angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  rectinfo = c(xwidth, zwidth)
  new_tibble_row(list(x = x, y = y, z = z, 
                 shape = "xz_rect",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(rectinfo = rectinfo),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(triangle(v1 = c(100, 100, 100), v2 = c(555/2, 455, 455), v3 = c(455, 100, 100),
#'                       material = diffuse(color = "purple"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' #Pass individual colors to each vertex: 
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(triangle(v1 = c(100, 100, 100), v2 = c(555/2, 455, 455), v3 = c(455, 100, 100),
#'                       color1 = "green", color2 = "yellow", color3 = "red")) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
triangle = function(v1 = c(1, 0, 0), v2 = c(0, 1, 0), v3 = c(-1, 0, 0), 
                    n1 = rep(NA, 3), n2 = rep(NA, 3), n3 = rep(NA, 3),
                    color1 = rep(NA, 3), color2 = rep(NA, 3), color3 = rep(NA, 3),
                    material = diffuse(), 
                    angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                    flipped = FALSE, reversed = FALSE, scale = c(1,1,1)) {
  if(!reversed) {
    vertex_vec = c(v1, v2, v3)
    normal_vec = c(n1, n2, n3)
  } else {
    vertex_vec = c(v3, v2, v1)
    normal_vec = c(n3, n2, n1)
  }
  vb = matrix(vertex_vec,nrow=3,ncol=3)
  it = matrix(c(1,2,3),nrow=3,ncol=1)
  if(all(!is.na(normal_vec))) {
    normals = matrix(normal_vec,nrow=3,ncol=3)
  } else {
    normals = matrix(0,nrow=3,ncol=0)
  }
  
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  vertex_colors = FALSE
  if(all(!is.na(color1))) {
    vertex_colors = TRUE
    color1 = convert_color(color1)
  }
  if(all(!is.na(color2))) {
    vertex_colors = TRUE
    color2 = convert_color(color2)
  }
  if(all(!is.na(color3))) {
    vertex_colors = TRUE
    color3 = convert_color(color3)
  }
  if(any(is.na(color1)) && any(!is.na(c(color2, color3)))) {
    color1 = unlist(material[[1]]$properties)[1:3]
  }
  if(any(is.na(color2)) && any(!is.na(c(color1, color3)))) {
    color2 = unlist(material[[1]]$properties)[1:3]
  }
  if(any(is.na(color3)) && any(!is.na(c(color1, color2)))) {
    color3 = unlist(material[[1]]$properties)[1:3]
  }
  color_matrix = matrix(c(color1, color2, color3),nrow=3,ncol=3,byrow=TRUE)
  mesh_obj = list(vb=vb,it=it,normals=normals)
  if(vertex_colors) {
    mesh_obj$material$color = color_matrix 
  } else {
    mesh_obj$material$color = matrix(nrow=0,ncol=0)
  }
  class(mesh_obj) = "mesh3d"
  mesh_obj$meshColor = ifelse(vertex_colors, "default", "vertex")
  mesh3d_model(mesh_obj,material = material, 
               angle = angle, order_rotation = order_rotation,
               scale = scale, 
               flipped = flipped)
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(disk(x = 555/2, y = 50, z = 555/2, radius = 150, 
#'                   material = diffuse(color = "orange"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' #Rotate the disk.
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(disk(x = 555/2, y = 555/2, z = 555/2, radius = 150, angle = c(-45, 0, 0), 
#'                   material = diffuse(color = "orange"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) , lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' #Pass a value for the inner radius.
#' if(run_documentation()) {
#' generate_cornell() %>% 
#'   add_object(disk(x = 555/2, y = 555/2, z = 555/2, 
#'                   radius = 150, inner_radius = 75, angle = c(-45, 0, 0), 
#'                   material = diffuse(color = "orange"))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
disk = function(x = 0, y = 0, z = 0, radius = 1, inner_radius = 0, material = diffuse(), 
                angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  info = c(inner_radius)

  new_tibble_row(list(x = x, y = y, z = z, shape = "disk",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(radius = radius,
                                                                     inner_radius = inner_radius),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' @param load_material Default `TRUE`. Whether to load the obj file material (MTL file). If material for faces
#' aren't specified, the default material will be used (specified by the user in `material`).
#' @param load_textures Default `TRUE`. If `load_material = TRUE`, whether to load textures in the MTL file (versus
#' just using the colors specified for each material).
#' @param load_normals Default `TRUE`. Whether to load the vertex normals if they exist in the OBJ file.
#' @param calculate_consistent_normals Default `TRUE`. Whether to calculate consistent vertex normals to prevent energy 
#' loss at edges.
#' @param vertex_colors Default `FALSE`. Set to `TRUE` if the OBJ file has vertex colors to apply them
#' to the model.
#' @param importance_sample_lights Default `TRUE`. Whether to importance sample lights specified in the OBJ material
#' (objects with a non-zero Ke MTL material).
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
#' if(run_documentation()) {
#' generate_ground(material = diffuse(checkercolor = "grey50")) %>%
#'   add_object(obj_model(y = -0.8, filename = r_obj(),
#'                        material = microfacet(color = "gold", roughness = 0.05))) %>%
#'   add_object(obj_model(x = 1.8, y = -0.8, filename = r_obj(), 
#'                        material = diffuse(color = "dodgerblue"))) %>%
#'   add_object(obj_model(x = -1.8, y = -0.8, filename = r_obj() , 
#'                        material = dielectric(attenuation = c(1,0.3,1)*2))) %>%
#'   add_object(sphere(z = 20, x = 20, y = 20, radius = 10,
#'                     material = light(intensity = 10))) %>%
#'   render_scene(parallel = TRUE, samples = 128, aperture = 0.05, 
#'                fov = 32, lookfrom = c(0, 2, 10))
#' 
#' }
#' 
#' #Use scale_obj to make objects bigger--this is more robust than the generic scale argument.
#' if(run_documentation()) {
#' generate_ground(material = diffuse(checkercolor = "grey50")) %>%
#'   add_object(obj_model(y = -0.8, filename = r_obj(), scale_obj = 2,
#'                        material = diffuse(noise = TRUE, noiseintensity = 10,noisephase=45))) %>%
#'   add_object(sphere(z = 20, x = 20, y = 20, radius = 10,
#'                     material = light(intensity = 10))) %>%
#'   render_scene(parallel = TRUE, samples = 128, ambient = TRUE, 
#'                backgroundhigh="blue", backgroundlow="red",
#'                aperture = 0.05, fov = 32, lookfrom = c(0, 2, 10),
#'                lookat = c(0,1,0)) 
#' }
obj_model = function(filename, x = 0, y = 0, z = 0, scale_obj = 1, 
                     load_material = TRUE, load_textures = TRUE, load_normals = TRUE,
                     vertex_colors = FALSE, calculate_consistent_normals = TRUE,
                     importance_sample_lights = TRUE,
                     material = diffuse(), 
                     angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                     flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  if(!load_material) {
    load_textures = FALSE
  }
  base_dir = function(x) {
    dirname_processed = dirname(x)
    if(dirname_processed == ".") {
      return("")
    } else {
      return(dirname_processed)
    }
  }
  new_tibble_row(list(x = x, y = y, z = z, 
                 shape = "obj",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(scale_obj = scale_obj,
                                                                     load_textures = load_textures,
                                                                     load_material = load_material,
                                                                     vertex_colors = vertex_colors,
                                                                     importance_sample_lights = importance_sample_lights,
                                                                     load_normals = load_normals,
                                                                     calculate_consistent_normals = calculate_consistent_normals,
                                                                     basename = base_dir(filename)),
                                             tricolorinfo = list(NA), 
                                             fileinfo = filename,
                                             material_id = NA_integer_, 
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' @param capped Default `TRUE`. Whether to add caps to the segment. Turned off when using the `light()` material.
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
#'                       length = 300, radius = 100, material = metal())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' #Rotate the cylinder
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(cylinder(x = 555/2, y = 250, z = 555/2, 
#'                       length = 300, radius = 100, angle = c(0, 0, 45),
#'                       material = diffuse())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' # Only render a subtended arc of the cylinder, flipping the normals.
#' if(run_documentation()) {
#' generate_cornell(lightintensity=3) %>%
#'   add_object(cylinder(x = 555/2, y = 250, z = 555/2, capped = FALSE,
#'                       length = 300, radius = 100, angle = c(45, 0, 0), phi_min = 0, phi_max = 180,
#'                       material = diffuse(), flipped = TRUE)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
cylinder = function(x = 0, y = 0, z = 0, radius = 1, length = 1, 
                    phi_min = 0, phi_max = 360, material = diffuse(), 
                    angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                    flipped = FALSE, scale = c(1,1,1), capped = TRUE) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  stopifnot(phi_max > phi_min)
  new_tibble_row(list(x = x, y = y, z = z,shape = "cylinder",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(radius = radius,
                                                                     length = length,
                                                                     phi_min = phi_min * pi / 180,
                                                                     phi_max = phi_max * pi / 180,
                                                                     has_cap = capped),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly. Notes: this will change the stated start/end position of the segment. 
#' Emissive objects may not currently function correctly when scaled.
#' @param capped Default `TRUE`. Whether to add caps to the segment. Turned off when using the `light()` material.
#' 
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the segment in the scene.
#' @export
#'
#' @examples
#' #Generate a segment in the cornell box. 
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(segment(start = c(100, 100, 100), end = c(455, 455, 455), radius = 50)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
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
#' if(run_documentation()) {
#' generate_cornell() %>% 
#'   add_object(scene_segments) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(cube_outline) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Shrink and rotate the cube
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(group_objects(cube_outline, pivot_point = c(555/2, 555/2, 555/2),
#'                            angle = c(45,45,45), scale = c(0.5,0.5,0.5))) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
segment = function(start = c(0, -1, 0), end = c(0, 1, 0), radius = 0.1, 
                   phi_min = 0, phi_max = 360, from_center = TRUE, direction = NA,
                   material = diffuse(), capped = TRUE, 
                   flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  stopifnot(phi_max > phi_min)
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
  
  # cap_int = ifelse(capped, 1, 0)
  length_xy = sqrt((end[1]-start[1])^2 + (end[3]-start[3])^2)
  if(end[1] == start[1] && end[3] == start[3]) {
    theta = 0
  } else {
    theta = atan2(length_xy, (end[2]-start[2]))/pi*180
  }
  fulllength = sqrt(sum((end-start)^2))
  angle = c(0, phi, theta)
  new_tibble_row(list(x = x, y = y, z = z, shape = "cylinder",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(radius = radius,
                                                                     length = fulllength,
                                                                     phi_min = phi_min * pi / 180,
                                                                     phi_max = phi_max * pi / 180,
                                                                     has_cap = capped),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
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
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 50, c = 50)) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Change the axes to make it taller rather than wide:
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 200, c = 100, material = metal())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
#' 
#' #Rotate it and make it dielectric:
#' if(run_documentation()) {
#' generate_cornell() %>%
#'   add_object(ellipsoid(x = 555/2, y = 555/2, z = 555/2, 
#'                        a = 100, b = 200, c = 100, angle = c(0, 0, 45),
#'                        material = dielectric())) %>%
#'   render_scene(lookfrom = c(278, 278, -800) ,lookat = c(278, 278, 0), fov = 40, 
#'                ambient_light = FALSE, samples = 128, parallel = TRUE, clamp_value = 5)
#' }
ellipsoid = function(x = 0, y = 0, z = 0, a = 1, b = 1, c = 1,
                  material = diffuse(), 
                  angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                  flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  new_tibble_row(list(x = x, y = y, z = z, shape = "ellipsoid",
                 material = material,
                 shape_info = ray_shape_info(shape_properties = list(a = a, 
                                                                     b = b,
                                                                     c = c),
                                             tricolorinfo = list(NA), 
                                             fileinfo = NA,
                                             material_id = NA_integer_,  
                                             csg_object = list(NA), 
                                             mesh_info = list(NA),
                                             flipped = flipped),
                 transforms = ray_transform(angle = list(angle),
                                            order_rotation = list(order_rotation),
                                            scale = list(scale),
                                            group_transform = list(matrix(NA_real_))),
                 animation_info = ray_animated_transform(
                   start_transform_animation = list(matrix(NA_real_)), 
                   end_transform_animation = list(matrix(NA_real_)),
                   start_time = 0, end_time = 1)
                 ))
}

#' Extruded Polygon Object
#'
#' @param polygon `sf` object, "SpatialPolygon" `sp` object,  or xy coordinates
#'   of polygon represented in a way that can be processed by `xy.coords()`.  If
#'   xy-coordinate based polygons are open, they will be closed by adding an
#'   edge from the last point to the first. If the `sf` object contains MULTIPOLYGONZ data, it will
#'   flattened.
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
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
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
#' if(run_documentation()) {
#' angles = seq(0,360,by=36)
#' xx = rev(c(rep(c(1,0.5),5),1) * sinpi(angles/180))
#' yy = rev(c(rep(c(1,0.5),5),1) * cospi(angles/180))
#' star_polygon = data.frame(x=xx,y=yy)
#' }
#' 
#' if(run_documentation()) {
#' generate_ground(depth=0,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(star_polygon,top=0.5,bottom=0,
#'                               material=diffuse(color="red",sigma=90))) %>%
#'   add_object(sphere(y=4,x=-3,z=-3,material=light(intensity=30))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,2,3),samples=128,lookat=c(0,0.5,0),fov=60)
#' }
#' 
#' #Now, let's add a hole to the center of the polygon. We'll make the polygon
#' #hollow by shrinking it, combining it with the normal size polygon,
#' #and specify with the `holes` argument that everything after `nrow(star_polygon)`
#' #in the following should be used to draw a hole:
#' 
#' if(run_documentation()) {
#' hollow_star = rbind(star_polygon,0.8*star_polygon)
#' }
#' 
#' if(run_documentation()) {
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon) + 1,
#'                               material=diffuse(color="red",sigma=90))) %>%
#'   add_object(sphere(y=4,x=-3,z=-3,material=light(intensity=30))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,2,4),samples=128,lookat=c(0,0,0),fov=30)
#' }
#' 
#' # Render one in the y-x plane as well by changing the `plane` argument,
#' # as well as offset it slightly.
#' if(run_documentation()) {
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, holes = nrow(star_polygon),
#'                               material=diffuse(color="red",sigma=90))) %>%
#'   add_object(extruded_polygon(hollow_star,top=0.25,bottom=0, y=1.2, z=-1.2, 
#'                               holes = nrow(star_polygon) + 1, plane = "yx", 
#'                               material=diffuse(color="green",sigma=90))) %>%
#'   add_object(sphere(y=4,x=-3,material=light(intensity=30))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,2,4),samples=128,lookat=c(0,0.9,0),fov=40)
#' }
#' 
#' # Now add the zy plane:
#' if(run_documentation()) {
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
#'   render_scene(parallel=TRUE,lookfrom = c(-4,2,4),samples=128,lookat=c(0,0.9,0),fov=40)
#' }
#' 
#' #We can also directly pass in sf polygons:
#' if(run_documentation()) {
#' if(length(find.package("spData",quiet=TRUE)) > 0) {
#'   us_states = spData::us_states
#'   texas = us_states[us_states$NAME == "Texas",]
#'   #Fix no sfc class in us_states geometry data
#'   class(texas$geometry) = c("list","sfc")
#' }
#' }
#' 
#' #This uses the raw coordinates, unless `center = TRUE`, which centers the bounding box
#' #of the polygon at the origin.
#' if(run_documentation()) {
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(texas, center = TRUE,
#'                               material=diffuse(color="#ff2222",sigma=90))) %>%
#'   add_object(sphere(y=30,x=-30,radius=10,
#'                     material=light(color="lightblue",intensity=40))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,10,-10),samples=128,fov=60)
#' }
#' 
#' #Here we use the raw coordinates, but offset the polygon manually.
#' if(run_documentation()) {
#' generate_ground(depth=-0.01,
#'                 material = diffuse(color="grey50",checkercolor="grey20")) %>%
#'   add_object(extruded_polygon(us_states, x=-96,z=-40, top=2,
#'                               material=diffuse(color="#ff2222",sigma=90))) %>%
#'   add_object(sphere(y=30,x=-100,radius=10,
#'                     material=light(color="lightblue",intensity=200))) %>%
#'   add_object(sphere(y=30,x=100,radius=10,
#'                     material=light(color="orange",intensity=200))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(0,120,-120),samples=128,fov=20)
#' }
#' 
#' #We can also set the map the height of each polygon to a column in the sf object,
#' #scaling it down by the maximum population state.
#' 
#' if(run_documentation()) {
#' generate_ground(depth=0,
#'                 material = diffuse(color="grey50",checkercolor="grey20",sigma=90)) %>%
#'   add_object(extruded_polygon(us_states, x=-96,z=-45, data_column_top = "total_pop_15",
#'                               scale_data = 1/max(us_states$total_pop_15)*5,
#'                               material=diffuse(color="#ff2222",sigma=90))) %>%
#'   add_object(sphere(y=30,x=-100,z=60,radius=10,
#'                     material=light(color="lightblue",intensity=250))) %>%
#'   add_object(sphere(y=30,x=100,z=-60,radius=10,
#'                     material=light(color="orange",intensity=250))) %>%
#'   render_scene(parallel=TRUE,lookfrom = c(-60,50,-40),lookat=c(0,-5,0),samples=128,fov=30)
#' }
#' 
extruded_polygon = function(polygon = NULL, x = 0, y = 0, z = 0, plane = "xz",
                   top = 1, bottom = 0, holes = NULL, 
                   angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                   material = diffuse(),
                   center = FALSE, flip_horizontal = FALSE, flip_vertical = FALSE,
                   data_column_top = NULL, data_column_bottom = NULL, scale_data = 1,
                   scale = c(1,1,1)) {
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
    if(length(find.package("sf",quiet=TRUE)) == 0) {
      stop("sf package required when handling sf objects")
    }
    poly_info = sf::st_drop_geometry(polygon)
    #Remove z dimension from multipolygon z geometry
    if(ncol(as.matrix(sf::st_geometry(polygon)[[1]])) == 3) {
      polygon = sf::as_Spatial(sf::st_zm(polygon))
    } else {
      polygon = sf::as_Spatial(polygon)
    }
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
    if(is.null(holes) || holes == 0) {
      holes = 0
    } else if (!is.numeric(holes) || anyNA(as.integer(holes))) {
      stop("holes must be integer")
    } else if (any(holes < 0L) || any(holes) > length(x) ||
               (any(holes == 0L) && length(holes) != 1L)) {
      stop("holes must be zero, or contain indices to polygon vertices")
    } else if (any(holes < 4)) {
      stop("holes cannot begin before vertex 4. Hole index here starts at: ", min(holes))
    }
    holes = as.integer(holes)
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
  counter = 1
  
  vert_list = list()
  idx_list = list()
  idx_counter = rev(c(1,2,3))
  
  for(poly in 1:length(poly_list)) {
    x=poly_list[[poly]][,1]
    y=poly_list[[poly]][,2]
    vertices = vertex_list[[poly]]
    height_poly = height_list[[poly]]
    bottom_poly = bottom_list[[poly]]
    
    
    for(i in 1:nrow(vertices)) {
      vert_list[[counter]] = matrix(c(scale*permute_axes(c(x[vertices[i,3]],bottom_poly,y[vertices[i,3]]),planeval),
                                      scale*permute_axes(c(x[vertices[i,2]],bottom_poly,y[vertices[i,2]]),planeval),
                                      scale*permute_axes(c(x[vertices[i,1]],bottom_poly,y[vertices[i,1]]),planeval)),
                                    ncol=3,nrow=3,byrow=TRUE)
      idx_list[[counter]] = matrix(idx_counter,nrow=1,ncol=3)
      idx_counter = idx_counter + 3
      counter = counter + 1
      if(extruded) {
        vert_list[[counter]] = matrix(c(scale*permute_axes(c(x[vertices[i,1]],height_poly,y[vertices[i,1]]),planeval),
                                        scale*permute_axes(c(x[vertices[i,2]],height_poly,y[vertices[i,2]]),planeval),
                                        scale*permute_axes(c(x[vertices[i,3]],height_poly,y[vertices[i,3]]),planeval)),
                                      ncol=3,nrow=3,byrow=TRUE)
        idx_list[[counter]] = matrix(idx_counter,nrow=1,ncol=3)
        idx_counter = idx_counter + 3
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
          vert_list[[counter]] = matrix(c(scale*permute_axes(c(xi,height_poly ,yi),planeval),
                                          scale*permute_axes(c(xi,bottom_poly ,yi),planeval),
                                          scale*permute_axes(c(xii,bottom_poly,yii),planeval)),
                                        ncol=3,nrow=3,byrow=TRUE)
          if(!xor(reversed, side_rev)) {
            idx_list[[counter]] = matrix(rev(idx_counter),nrow=1,ncol=3)
          } else {
            idx_list[[counter]] = matrix(idx_counter,nrow=1,ncol=3)
          }
          idx_counter = idx_counter + 3
          counter = counter + 1
          vert_list[[counter]] = matrix(c(scale*permute_axes(c(xi,height_poly ,yi),planeval),
                                          scale*permute_axes(c(xii,bottom_poly ,yii),planeval),
                                          scale*permute_axes(c(xii,height_poly,yii),planeval)),
                                        ncol=3,nrow=3,byrow=TRUE)
          if(!xor(reversed, side_rev)) {
            idx_list[[counter]] = matrix(rev(idx_counter),nrow=1,ncol=3)
          } else {
            idx_list[[counter]] = matrix(idx_counter,nrow=1,ncol=3)
          }
          idx_counter = idx_counter + 3
          counter = counter + 1
        }
      }
    }
  }
  
  v_mat = do.call(rbind,vert_list)
  i_mat = do.call(rbind,idx_list)
  mesh = list()
  mesh$vb = t(v_mat)
  mesh$it = t(i_mat)
  class(mesh) = "mesh3d"
  return(mesh3d_model(mesh, x=x_off,y=y_off,z=z_off,
                      angle = angle, order_rotation = order_rotation,
                      override_material=TRUE, material=material))
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
#' if(run_documentation()) {
#' generate_studio() %>% 
#'  add_object(cone(start=c(0,-1,0), end=c(0,1,0), radius=1,material=diffuse(color="red"))) %>% 
#'  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'  render_scene(samples=128,clamp_value=10)
#' }
#' if(run_documentation()) {
#'  #Change the radius, length, and direction
#' generate_studio() %>% 
#'  add_object(cone(start=c(0,0,0), end=c(0,-1,0), radius=0.5,material=diffuse(color="red"))) %>% 
#'  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'  render_scene(samples=128,clamp_value=10)
#' }
#' if(run_documentation()) {
#' #Give custom start and end points (and customize the color/texture)
#' generate_studio() %>% 
#'  add_object(cone(start=c(-1,0.5,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="red"))) %>%
#'  add_object(cone(start=c(1,0.5,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="green"))) %>%
#'  add_object(cone(start=c(0,1,-1), end=c(0,0,0), radius=0.5,material=diffuse(color="orange"))) %>% 
#'  add_object(cone(start=c(-1,-0.5,0), end=c(1,-0.5,0), radius=0.25,
#'    material = diffuse(color="red",gradient_color="green"))) %>% 
#'  add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'  render_scene(samples=128,clamp_value=10)
#' }
#' if(run_documentation()) {  
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
#'   render_scene(lookfrom=c(0,4,10), clamp_value=10, samples=128)
#' }
#' if(run_documentation()) {
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
#'   render_scene(lookfrom=c(0,4,10), clamp_value=10,fov=25, samples=128)
#' }
cone = function(start = c(0, 0, 0), end = c(0, 1, 0), radius = 0.5, 
                direction = NA, from_center = TRUE,
                material = diffuse(), angle = c(0,0,0),
                flipped = FALSE, scale = c(1,1,1)) {
  return(raymesh_model(rayvertex::cone_mesh(start = start,
                                            end = end, radius = radius, direction = direction,
                                            from_center = from_center),
                       material = material, angle = angle, flipped = flipped, scale = scale))
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
#' if(run_documentation()) {
#' generate_studio() %>% 
#'   add_object(arrow(start = c(-1,0,0), end = c(1,0,0), material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,material=light(intensity=20))) %>% 
#'   render_scene(clamp_value=10,  samples=128)
#' }
#' if(run_documentation()) {
#' #Change the proportion of tail to top
#' generate_studio(depth=-2) %>% 
#'   add_object(arrow(start = c(-1,-1,0), end = c(1,-1,0), tail_proportion = 0.5,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(arrow(start = c(-1,0,0), end = c(1,0,0), tail_proportion = 0.75,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(arrow(start = c(-1,1,0), end = c(1,1,0), tail_proportion = 0.9,
#'                    material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,z=5,x=2,material=light(intensity=30))) %>% 
#'   render_scene(clamp_value=10, fov=25,  samples=128)
#' }
#' if(run_documentation()) {
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
#'   render_scene(clamp_value=10, samples=128)
#' }
#' if(run_documentation()) {
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
#'   render_scene(clamp_value=10, samples=128, 
#'                lookfrom=c(0,5,10), lookat=c(0,-0.5,0), fov=16)
#' }
#' if(run_documentation()) {
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
#'   render_scene(fov=20, ambient=TRUE, samples=128,
#'                backgroundlow="black",backgroundhigh="white") 
#' }
arrow = function(start = c(0,0,0), end = c(0,1,0), 
                 radius_top = 0.2, radius_tail = 0.1, tail_proportion = 0.5,
                 direction = NA,  from_center = TRUE, 
                 material = diffuse(), 
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
          flipped = flipped, scale = scale) %>%
    add_object(cone(start = start_tip, 
                    end = end_tip,
                    radius = radius_top, 
                    material = material, 
                    flipped = flipped, scale = scale))
}

#' Bezier Curve Object
#' 
#' Bezier curve, defined by 4 control points.
#'
#' @param p1 Default `c(0,0,0)`. First control point. Can also be a list of 4 length-3 numeric vectors 
#' or 4x3 matrix/data.frame specifying the x/y/z control points.
#' @param p2 Default `c(-1,0.33,0)`. Second control point.
#' @param p3 Default `c(1,0.66,0)`. Third control point.
#' @param p4 Default `c(0,1,0)`. Fourth control point.
#' @param x Default `0`. x-coordinate offset for the curve.
#' @param y Default `0`. y-coordinate offset for the curve.
#' @param z Default `0`. z-coordinate offset for the curve.
#' @param width Default `0.1`. Curve width.
#' @param width_end Default `NA`. Width at end of path. Same as `width`, unless specified.
#' @param u_min Default `0`. Minimum parametric coordinate for the curve.
#' @param u_max Default `1`. Maximum parametric coordinate for the curve.
#' @param type Default `cylinder`. Other options are `flat` and `ribbon`.
#' @param normal Default `c(0,0,-1)`. Orientation surface normal for the start of ribbon curves.
#' @param normal_end Default `NA`. Orientation surface normal for the start of ribbon curves. If not
#' specified, same as `normal`.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
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
#' #Generate the default curve:
#' if(run_documentation()) {
#' generate_studio(depth=-0.2) %>%
#'   add_object(bezier_curve(material=diffuse(color="red"))) %>%
#'   add_object(sphere(y=3,z=5,x=2,radius=0.3,
#'                     material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
#'   render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=13,
#'                samples=128)
#' }
#' 
#' if(run_documentation()) {
#' #Change the control points to change the direction of the curve. Here, we place spheres
#' #at the control point locations.
#' generate_studio(depth=-0.2) %>%
#'   add_object(bezier_curve(material=diffuse(color="red"))) %>%
#'   add_object(sphere(radius=0.075,material=glossy(color="green"))) %>% 
#'   add_object(sphere(radius=0.075,x=-1,y=0.33,material=glossy(color="green"))) %>% 
#'   add_object(sphere(radius=0.075,x=1,y=0.66,material=glossy(color="green"))) %>% 
#'   add_object(sphere(radius=0.075,y=1,material=glossy(color="green"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,radius=0.3,
#'                     material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
#'   render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=15,
#'                samples=128)
#' }
#' if(run_documentation()) {
#' #We can make the curve flat (always facing the camera) by setting the type to `flat`
#' generate_studio(depth=-0.2) %>%
#'   add_object(bezier_curve(type="flat", material=glossy(color="red"))) %>%
#'   add_object(sphere(y=3,z=5,x=2,radius=0.3,
#'                     material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
#'   render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=13,
#'                samples=128)
#' }
#' if(run_documentation()) {
#' #We can also plot a ribbon, which is further specified by a start and end orientation with
#' #two surface normals.
#' generate_studio(depth=-0.2) %>%
#'   add_object(bezier_curve(type="ribbon", width=0.2,
#'                    p1 = c(0,0,0), p2 = c(0,0.33,0), p3 = c(0,0.66,0), p4 = c(0.3,1,0),
#'                    normal_end = c(0,0,1),
#'                    material=glossy(color="red"))) %>%
#'   add_object(sphere(y=3,z=5,x=2,radius=0.3,
#'                     material=light(intensity=200, spotlight_focus = c(0,0.5,0)))) %>%
#'   render_scene(clamp_value = 10, lookat = c(0,0.5,0), fov=13,
#'                samples=128)
#' }
#' if(run_documentation()) {
#' #Create a single curve and copy and rotate it around the y-axis to create a wavy fountain effect:
#' scene_curves = list()
#' for(i in 1:90) {
#'   scene_curves[[i]] = bezier_curve(p1 = c(0,0,0),p2 = c(0,5-sinpi(i*16/180),2),
#'                             p3 = c(0,5-0.5 * sinpi(i*16/180),4),p4 = c(0,0,6),
#'                             angle=c(0,i*4,0), type="cylinder",
#'                             width = 0.1, width_end =0.1,material=glossy(color="red"))
#' }
#' all_curves = do.call(rbind, scene_curves)
#' generate_ground(depth=0,material=diffuse(checkercolor="grey20")) %>%
#'   add_object(all_curves) %>%
#'   add_object(sphere(y=7,z=0,x=0,material=light(intensity=100))) %>% 
#'   render_scene(lookfrom = c(12,20,50),samples=100,
#'                lookat=c(0,1,0), fov=15, clamp_value = 10)
#' 
#' }
bezier_curve = function(p1 = c(0,0,0), p2 = c(-1,0.33,0), p3 = c(1,0.66,0), p4=c(0,1,0), 
                        x=0, y=0, z=0,
                        width = 0.1, width_end = NA, u_min = 0, u_max = 1, type = "cylinder",
                        normal = c(0,0,-1), normal_end = NA, 
                        material = diffuse(), angle = c(0, 0, 0),
                        order_rotation = c(1, 2, 3), 
                        flipped = FALSE, scale = c(1,1,1)) {
  if(inherits(p1,"list")) {
    stopifnot(length(p1) == 4 && all(lapply(p1, (function(x) length(x) == 3))))
    p1 = do.call(rbind, p1)
  }
  if(inherits(p1,"matrix") || inherits(p1,"data.frame")) {
    if(dim(p1)[1] == 4 && dim(p1)[2] == 3) {
      p2 = p1[2,]
      p3 = p1[3,]
      p4 = p1[4,]
      p1 = p1[1,]
    }
  }
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  if(is.na(width_end)) {
    width_end = width
  }
  stopifnot(u_min < u_max)
  stopifnot(length(width) == 1 && is.numeric(width))
  stopifnot(length(width_end) == 1 && is.numeric(width_end))
  
  if(all(!is.na(normal)) && all(is.na(normal_end))) {
    normal_end = normal
  }
  if(all(!is.na(normal)) || all(!is.na(normal_end))) {
    stopifnot(length(normal) == 3 && is.numeric(normal))
    stopifnot(length(normal_end) == 3 && is.numeric(normal_end))
  }
  if(material[[1]]$type == "hair") {
    type = "flat"
  }
  if(all(is.na(normal)) || type == "cylinder" || type == "flat") {
    normal = c(0,0,0)
    normal_end = c(0,0,0)
  }
  curvetype = unlist(lapply(tolower(type),switch,
                           "flat" = 1,"cylinder" = 2, "ribbon" = 3))
  new_tibble_row(list(x = x, y = y, z = z, shape = "curve",
                      material = material,
                      shape_info = ray_shape_info(shape_properties = list(p1 = p1,
                                                                          p2 = p2,
                                                                          p3 = p3,
                                                                          p4 = p4,
                                                                          width = width,
                                                                          width_end = width_end,
                                                                          u_min = u_min,
                                                                          u_max = u_max,
                                                                          curvetype = curvetype,
                                                                          normal = normal,
                                                                          normal_end = normal_end),
                                                  tricolorinfo = list(NA), 
                                                  fileinfo = NA,
                                                  material_id = NA_integer_,  
                                                  csg_object = list(NA), 
                                                  mesh_info = list(NA),
                                                  flipped = flipped),
                      transforms = ray_transform(angle = list(angle),
                                                 order_rotation = list(order_rotation),
                                                 scale = list(scale),
                                                 group_transform = list(matrix(NA_real_))),
                      animation_info = ray_animated_transform(
                        start_transform_animation = list(matrix(NA_real_)), 
                        end_transform_animation = list(matrix(NA_real_)),
                        start_time = 0, end_time = 1)
                      ))
}

#' Path Object
#' 
#' Either a closed or open path made up of bezier curves that go through the specified points 
#' (with continuous first and second derivatives), or straight line segments.
#'
#' @param points Either a list of length-3 numeric vectors or 3-column matrix/data.frame specifying
#' the x/y/z points that the path should go through.
#' @param x Default `0`. x-coordinate offset for the path.
#' @param y Default `0`. y-coordinate offset for the path.
#' @param z Default `0`. z-coordinate offset for the path.
#' @param closed Default `FALSE`. If `TRUE`, the path will be closed by smoothly connecting the first
#' and last points.
#' @param closed_smooth Default `TRUE`. If `closed = TRUE`, this will ensure C2 (second derivative) 
#' continuity between the ends. If `closed = FALSE`, the curve will only have C1 (first derivative)
#' continuity between the ends.
#' @param straight Default `FALSE`. If `TRUE`, straight lines will be used to connect the points instead
#' of bezier curves. 
#' @param precomputed_control_points Default `FALSE`. If `TRUE`, `points` argument will expect
#' a list of control points calculated with the internal rayrender function `rayrender:::calculate_control_points()`.
#' @param width Default `0.1`. Curve width.
#' @param width_end Default `NA`. Width at end of path. Same as `width`, unless specified.
#' @param u_min Default `0`. Minimum parametric coordinate for the path.
#' @param u_max Default `1`. Maximum parametric coordinate for the path.
#' @param type Default `cylinder`. Other options are `flat` and `ribbon`.
#' @param normal Default `c(0,0,-1)`. Orientation surface normal for the start of ribbon curves.
#' @param normal_end Default `NA`. Orientation surface normal for the start of ribbon curves. If not
#' specified, same as `normal`.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
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
#' if(run_documentation()) {
#' #Generate a wavy line, showing the line goes through the specified points:
#' wave = list(c(-2,1,0),c(-1,-1,0),c(0,1,0),c(1,-1,0),c(2,1,0))
#' point_mat = glossy(color="green")
#' generate_studio(depth=-1.5) %>% 
#'   add_object(path(points = wave,material=glossy(color="red"))) %>% 
#'   add_object(sphere(x=-2,y=1,radius=0.1,material=point_mat)) %>% 
#'   add_object(sphere(x=-1,y=-1,radius=0.1,material=point_mat)) %>% 
#'   add_object(sphere(x=0,y=1,radius=0.1,material=point_mat)) %>% 
#'   add_object(sphere(x=1,y=-1,radius=0.1,material=point_mat)) %>% 
#'   add_object(sphere(x=2,y=1,radius=0.1,material=point_mat)) %>% 
#'   add_object(sphere(z=5,x=5,y=5,radius=2,material=light(intensity=15))) %>% 
#'   render_scene(samples=128, clamp_value=10,fov=30)
#' }
#' if(run_documentation()) {
#' #Here we use straight lines by setting `straight = TRUE`:
#' generate_studio(depth=-1.5) %>% 
#'   add_object(path(points = wave,straight = TRUE, material=glossy(color="red"))) %>% 
#'   add_object(sphere(z=5,x=5,y=5,radius=2,material=light(intensity=15))) %>% 
#'   render_scene(samples=128, clamp_value=10,fov=30)
#' }
#' if(run_documentation()) {
#' #We can also pass a matrix of values, specifying the x/y/z coordinates. Here,
#' #we'll create a random curve:
#' set.seed(21)
#' random_mat = matrix(runif(3*9)*2-1, ncol=3)
#' generate_studio(depth=-1.5) %>% 
#'   add_object(path(points=random_mat, material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,radius=1,material=light(intensity=30))) %>% 
#'   render_scene(samples=128, clamp_value=10)
#' }
#' if(run_documentation()) {
#' #We can ensure the curve is closed by setting `closed = TRUE`
#' generate_studio(depth=-1.5) %>% 
#'   add_object(path(points=random_mat, closed = TRUE, material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,radius=1,material=light(intensity=30))) %>% 
#'   render_scene(samples=128, clamp_value=10)
#' }
#' if(run_documentation()) {
#' #Finally, let's render a pretzel to show how you can render just a subset of the curve:
#' pretzel = list(c(-0.8,-0.5,0.1),c(0,-0.2,-0.1),c(0,0.3,0.1),c(-0.5,0.5,0.1), c(-0.6,-0.5,-0.1),
#'                c(0,-0.8,-0.1),
#'                c(0.6,-0.5,-0.1),c(0.5,0.5,-0.1), c(0,0.3,-0.1),c(-0,-0.2,0.1), c(0.8,-0.5,0.1))
#'                
#' #Render the full pretzel:
#' generate_studio(depth = -1.1) %>% 
#'   add_object(path(pretzel, width=0.17,  material = glossy(color="#db5b00"))) %>% 
#'   add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
#'   render_scene(samples=128, clamp_value=10)
#' }
#' if(run_documentation()) {
#' #Here, we'll render only the first third of the pretzel by setting `u_max = 0.33`
#' generate_studio(depth = -1.1) %>% 
#'   add_object(path(pretzel, width=0.17, u_max=0.33, material = glossy(color="#db5b00"))) %>% 
#'   add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
#'   render_scene(samples=128, clamp_value=10)
#' }
#' if(run_documentation()) {
#' #Here's the last third, by setting `u_min = 0.66`
#' generate_studio(depth = -1.1) %>% 
#'   add_object(path(pretzel, width=0.17, u_min=0.66, material = glossy(color="#db5b00"))) %>% 
#'   add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
#'   render_scene(samples=128, clamp_value=10)
#' }
#' if(run_documentation()) {
#' #Here's the full pretzel, decomposed into thirds using the u_min and u_max coordinates
#' generate_studio(depth = -1.1) %>% 
#'   add_object(path(pretzel, width=0.17, u_max=0.33, x = -0.8, y =0.6,
#'                   material = glossy(color="#db5b00"))) %>% 
#'   add_object(path(pretzel, width=0.17, u_min=0.66, x = 0.8, y =0.6,
#'                   material = glossy(color="#db5b00"))) %>% 
#'   add_object(path(pretzel, width=0.17, u_min=0.33, u_max=0.66, x=0,
#'                   material = glossy(color="#db5b00"))) %>% 
#'   add_object(sphere(y=5,x=2,z=4,material=light(intensity=20,spotlight_focus = c(0,0,0)))) %>% 
#'   render_scene(samples=128, clamp_value=10, lookfrom=c(0,3,10))
#' }
path = function(points,
                x=0,y=0,z=0, closed = FALSE, closed_smooth = TRUE, 
                straight = FALSE, precomputed_control_points = FALSE,
                width = 0.1, width_end = NA, u_min = 0, u_max = 1, type = "cylinder",
                normal = c(0,0,-1), normal_end = NA, 
                material = diffuse(), angle = c(0, 0, 0),
                order_rotation = c(1, 2, 3),
                flipped = FALSE, scale = c(1,1,1)) {
  if(is.na(width_end)) {
    width_end = width
  }
  if(u_min == u_max) {
    return()
  }
  stopifnot(u_min <= u_max)
  stopifnot(length(width) == 1 && is.numeric(width))
  stopifnot(length(width_end) == 1 && is.numeric(width_end))
  
  if(all(!is.na(normal)) && all(is.na(normal_end))) {
    normal_end = normal
  }
  if(all(!is.na(normal)) || all(!is.na(normal_end))) {
    stopifnot(length(normal) == 3 && is.numeric(normal))
    stopifnot(length(normal_end) == 3 && is.numeric(normal_end))
  }
  if(all(is.na(normal)) || type == "cylinder" || type == "flat") {
    normal = c(0,0,0)
    normal_end = c(0,0,0)
  }
  if(inherits(points,"numeric")) {
    stop("Input must either be list, matrix, or data.frame, not numeric.")
  }
  if(inherits(points,"list") && !precomputed_control_points) {
    if(any(unlist(lapply(points,(function(x) length(x) != 3))))) {
      stop("If `points` is a list, each entry must be a length-3 vector")
    }
    points = do.call(rbind,points)
  }
  if(inherits(points,"data.frame")) {
    points = as.matrix(points)
  }
  if(is.array(points)) {
    if(nrow(points) == 1) {
      stop("Only one point passed, no path specified.")
    }
    if(nrow(points) == 2 && closed) {
      closed=FALSE
    }
  }
  if(!precomputed_control_points) {
    if(inherits(points,"matrix")) {
      if(ncol(points) == 3) {
        if(closed && closed_smooth) {
          points = rbind(points[c(nrow(points)-2,nrow(points)-1,nrow(points)),], points, points[1:3,])
        }
        if(!straight) {
          full_control_points = calculate_control_points(points)
        } else {
          full_control_points = calculate_control_points_straight(points)
        }
        if(closed && closed_smooth) {
          full_control_points[[length(full_control_points)]] = NULL
          full_control_points[[length(full_control_points)]] = NULL
          full_control_points[[1]] = NULL
          full_control_points[[1]] = NULL
          full_control_points[[1]] = NULL
        }
      } else {
        stop("If points a matrix or data.frame, must have 3 columns")
      }
    } else {
      stop("points not of supported type (function expects matrix/data.frame/list, got ", class(points),")")
    }
  } else {
    full_control_points = points
  }
  if(closed && !closed_smooth) {
    first_point = full_control_points[[1]]
    last_point = full_control_points[[length(full_control_points)]]
    full_control_points[[length(full_control_points) + 1]] = last_point
    full_control_points[[length(full_control_points)]][4,] = first_point[1,]
    full_control_points[[length(full_control_points)]][3,] = 2*first_point[1,] - first_point[2,]
    full_control_points[[length(full_control_points)]][2,] = 2*last_point[4,]  - last_point[3,]
    full_control_points[[length(full_control_points)]][1,] = last_point[4,]
  }
  u_min_segment = length(full_control_points) * u_min 
  u_max_segment = length(full_control_points) * u_max 
  
  seg_begin = floor(u_min_segment) + 1
  seg_end = floor(u_max_segment) + 1
  curve_list = list()
  for(i in seq_len(length(full_control_points))) {
    u_min_temp = 0
    u_max_temp = 1
    if(u_min_segment <= i && u_min_segment > i-1) {
      u_min_temp = u_min_segment - (i - 1)
    }
    if(u_max_segment >= i-1 && u_max_segment < i) {
      u_max_temp = u_max_segment - (i - 1)
    }
    if(u_min_temp == u_max_temp) {
      next
    }
    if(round(u_max_temp,8) == 0) {
      break
    }
    temp_width_start = width + (width_end - width) * (u_min_temp + i - 1) / length(full_control_points)
    temp_width_end = width + (width_end - width) * (u_max_temp + i - 1) / length(full_control_points)
    temp_normal_start = normal + (normal_end - normal) * (u_min_temp + i - 1) / length(full_control_points)
    temp_normal_end = normal + (normal_end - normal) * (u_max_temp + i - 1) / length(full_control_points)

    if(seg_begin <= i  && seg_end >= i) {
      curve_list[[i]] = bezier_curve(x=x,y=y,z=z,
                                     p1=full_control_points[[i]], type = type,
                                     u_min = u_min_temp,
                                     u_max = u_max_temp,
                                     width = temp_width_start, width_end = temp_width_end,
                                     normal = temp_normal_start, normal_end = temp_normal_end,
                                     material = material, angle = angle,
                                     order_rotation = order_rotation, 
                                     flipped = flipped, scale = scale)
    }
  }
  if(length(find.package("dplyr",quiet=TRUE)) > 0) {
    return(dplyr::bind_rows(curve_list))
  } else {
    return(do.call(rbind, curve_list))
  }
}

#' Text Object
#'
#' @param label Text string.
#' @param x Default `0`. x-coordinate of the center of the label.
#' @param y Default `0`. y-coordinate of the center of the label.
#' @param z Default `0`. z-coordinate of the center of the label.
#' @param text_height Default `1`. Height of the text.
#' @param orientation Default `xy`. Orientation of the plane. Other options are `yz` and `xz`.
#' @param material Default  \code{\link{diffuse}}. The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the text in the scene.
#' @export
#'
#' @examples
#' #Generate a label in the cornell box.
#' if(run_documentation()) {
#' generate_cornell() %>% 
#'   add_object(text3d(label="Cornell Box", x=555/2,y=555/2,z=555/2,text_height=60,
#'                     material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
#'   render_scene(samples=128, clamp_value=10)
#' }
#' if(run_documentation()) {
#' #Change the orientation
#' generate_cornell() %>% 
#'   add_object(text3d(label="YZ Plane", x=550,y=555/2,z=555/2,text_height=100,
#'                     orientation = "yz",
#'                     material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
#'  add_object(text3d(label="XY Plane", z=550,y=555/2,x=555/2,text_height=100,
#'                     orientation = "xy",
#'                     material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
#'  add_object(text3d(label="XZ Plane", z=555/2,y=5,x=555/2,text_height=100,
#'                     orientation = "xz",
#'                     material=diffuse(color="grey10"))) %>% 
#'   render_scene(samples=128, clamp_value=10)
#' }
#' if(run_documentation()) {
#' #Add an label in front of a sphere
#' generate_cornell() %>% 
#'   add_object(text3d(label="Cornell Box", x=555/2,y=555/2,z=555/2,text_height=60,
#'                     material=diffuse(color="grey10"), angle=c(0,180,0))) %>% 
#'   add_object(text3d(label="Sphere", x=555/2,y=100,z=100,text_height=30,
#'                     material=diffuse(color="white"), angle=c(0,180,0))) %>% 
#'   add_object(sphere(y=100,radius=100,z=555/2,x=555/2,
#'                     material=glossy(color="purple"))) %>% 
#'   add_object(sphere(y=555,radius=100,z=-1000,x=555/2,
#'                     material=light(intensity=100,
#'                                    spotlight_focus=c(555/2,100,100)))) %>%                   
#'   render_scene(samples=128, clamp_value=10)
#' }
#'   
#' if(run_documentation()) {
#' #A room full of bees
#' bee_list = list()
#' for(i in 1:100) {
#' bee_list[[i]] = text3d("B", x=20+runif(1)*525, y=20+runif(1)*525, z=20+runif(1)*525, 
#'                        text_height = 50, angle=c(0,180,0))
#' }
#' bees = do.call(rbind,bee_list)
#' generate_cornell() %>% 
#'   add_object(bees) %>%                   
#'   render_scene(samples=128, clamp_value=10)
#' }
text3d = function(label, x = 0, y = 0, z = 0, text_height = 1, orientation = "xy",
                  material = diffuse(), 
                  angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                  flipped = FALSE, scale = c(1,1,1)) {
  labelfile = tempfile(fileext = ".png")
  rayimage::add_title(matrix(0,ncol = nchar(label)*60, nrow=60*1.2), 
                      title_size  = 60,
                      title_offset = c(0,0),title_text = label, title_color = "white",
                      title_position = "center", filename = labelfile)
  material[[1]]$alphaimage = list(labelfile)
  if(orientation == "xy" || orientation == "yx") {
    rayrender::xy_rect(x=x,y=y,z=z, angle = angle,
                       xwidth = nchar(label)*text_height, ywidth = text_height,
                       material = material)
  } else if (orientation == "yz" || orientation == "zy") {
    rayrender::yz_rect(x=x,y=y,z=z, angle = angle,
                       zwidth = nchar(label)*text_height, ywidth = text_height,
                       material = material)
  } else if (orientation == "xz" || orientation == "zx") {
    rayrender::xz_rect(x=x,y=y,z=z, angle = angle,
                       xwidth = nchar(label)*text_height, zwidth = text_height,
                       material = material)
  } else {
    stop("Orientation ", orientation, " not recognized")
  }
}

#' `ply` File Object
#' 
#' Load an PLY file via a filepath. 
#' Note: light importance sampling currently not supported for this shape.
#'
#' @param filename Filename and path to the `ply` file. Can also be a `txt` file, if it's in the correct `ply` internally.
#' @param x Default `0`. x-coordinate to offset the model.
#' @param y Default `0`. y-coordinate to offset the model.
#' @param z Default `0`. z-coordinate to offset the model.
#' @param scale_ply Default `1`. Amount to scale the model. Use this to scale the object up or down on all axes, as it is
#' more robust to numerical precision errors than the generic scale option.
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
#' #See the documentation for `obj_model()`--no example PLY models are included with this package,
#' #but the process of loading a model is the same (without support for vertex colors).
ply_model = function(filename, x = 0, y = 0, z = 0, scale_ply = 1, 
                     material = diffuse(), 
                     angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                     flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  tempcon = file(filename, open="rt")
  on.exit(close(tempcon))
  is_ply = scan(tempcon,what=character(),n=1, quiet=TRUE) == "ply"
  if(!is_ply) {
    stop(filename, " does not appear to be PLY file.")
  }
  tokenval = scan(tempcon,what=character(),n=1, quiet=TRUE)
  while(tokenval != "end_header") {
    if(tokenval == "vertex") {
      if(scan(tempcon,what=integer(),n=1, quiet=TRUE) == 0) {
        warning(filename, " contains no vertices, skipping.")
        return()
      }
    }
    if(tokenval == "face") {
      if(scan(tempcon,what=integer(),n=1, quiet=TRUE) == 0) {
        warning(filename, " contains no polygon faces, skipping.")
        return()
      }
    }
    tokenval = scan(tempcon,what=character(),n=1, quiet=TRUE)
  }
  base_dir = function(x) {
    dirname_processed = dirname(x)
    if(dirname_processed == ".") {
      return("")
    } else {
      return(dirname_processed)
    }
  }
  new_tibble_row(list(x = x, y = y, z = z, 
                      shape = "ply",
                      material = material,
                      shape_info = ray_shape_info(shape_properties = list(scale_ply = scale_ply,
                                                                          basename = base_dir(filename)),
                                                  tricolorinfo = list(NA), 
                                                  fileinfo = filename,
                                                  material_id = NA_integer_,  
                                                  csg_object = list(NA), 
                                                  mesh_info = list(NA),
                                                  flipped = flipped),
                      transforms = ray_transform(angle = list(angle),
                                                 order_rotation = list(order_rotation),
                                                 scale = list(scale),
                                                 group_transform = list(matrix(NA_real_))),
                      animation_info = ray_animated_transform(
                        start_transform_animation = list(matrix(NA_real_)), 
                        end_transform_animation = list(matrix(NA_real_)),
                        start_time = 0, end_time = 1)
                      ))
}

#' `mesh3d` model
#' 
#' Load an `mesh3d` (or `shapelist3d`) object, as specified in the `rgl` package. 
#'
#' @param mesh A `mesh3d` or `shapelist3d` object. Pulls the vertex, index, texture coordinates, 
#' normals, and material information. If the material references an image texture, the 
#' `mesh$material$texture` argument should be set to the image filename. The `mesh3d` format
#' only supports one image texture per mesh. All quads will be triangulated.
#' @param x Default `0`. x-coordinate to offset the model.
#' @param y Default `0`. y-coordinate to offset the model.
#' @param z Default `0`. z-coordinate to offset the model.
#' @param swap_yz Default `FALSE`. Swap the Y and Z coordinates.
#' @param reverse Default `FALSE`. Reverse the orientation of the indices, flipping their normals.
#' @param scale_mesh Default `1`. Amount to scale the size of the mesh in all directions.
#' @param verbose Default `FALSE`. If `TRUE`, prints information about the mesh to the console.
#' @param override_material Default `FALSE`. If `TRUE`, overrides the material specified in the 
#' `mesh3d` object with the one specified in `material`.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}. 
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Single row of a tibble describing the mesh3d model in the scene.
#' @export
#'
#' @examples
#' #Load a mesh3d object (from the Rvcg) and render it:
#' if(length(find.package("Rcvg", quiet=TRUE)) > 0) {
#'   library(Rvcg)
#'   data(humface)
#'   
#'   generate_studio() %>% 
#'     add_object(mesh3d_model(humface,y=-0.3,x=0,z=0,
#'                           material=glossy(color="dodgerblue4"), scale_mesh = 1/70)) %>%
#'     add_object(sphere(y=5,x=5,z=5,material=light(intensity=50))) %>% 
#'     render_scene(samples=128,width=800,height=800,
#'                  lookat = c(0,0.5,1), aperture=0.0)
#' }
mesh3d_model = function(mesh, x = 0, y = 0, z = 0, swap_yz = FALSE, reverse = FALSE,
                        scale_mesh = 1, verbose = FALSE, 
                        override_material = FALSE, material = diffuse(), 
                        angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                        flipped = FALSE, scale = c(1,1,1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  if(inherits(mesh,"shapelist3d")) {
    shapes = list()
    for(shape in seq_len(length(mesh))) {
      shapes[[shape]] = mesh3d_model(mesh[[shape]], x=x,y=y,z=z,
                                     swap_yz=swap_yz,reverse=reverse,scale_mesh = scale_mesh,
                                     override_material=override_material,material=material,
                                     angle=angle,order_rotation=order_rotation,flipped=flipped,
                                     scale=scale,verbose=verbose)
    }
    return(do.call(rbind,shapes))
  }
  if(!inherits(mesh,"mesh3d")) {
    stop("mesh must be of class 'mesh3d': actual class is ", class(mesh))
  }
  vertices = t(mesh$vb)
  if(ncol(vertices) == 4) {
    vertices = vertices[,1:3]
  }
  if(swap_yz) {
    vertices = vertices[,c(1,3,2)]
  }
  ## there might be triangles, quads, or both
  indices = NULL
  if (!is.null(mesh$it)) {
    indices = t(mesh$it)-1
  }
  if(!is.null(mesh$ib)) {
    quads = mesh$ib
    tri_ind = t(matrix(rbind(quads[c(1L, 2L, 4L),], 
                        quads[c(2L, 3L, 4L),]), 3L)) - 1
    indices = rbind(indices,tri_ind)
  }
  normals = mesh$normals
  if(is.null(normals)) {
    normals = matrix(nrow=0,ncol=3)
  } else {
    normals = t(normals)
  }
  texcoords = mesh$texcoords
  if(is.null(texcoords)) {
    texcoords = matrix(nrow=0,ncol=2)
  } else {
    texcoords = t(texcoords)
  }
  texture = mesh$material$texture
  if(!is.null(texture)) {
    texture = path.expand(texture)
  } else {
    texture = ""
  }
  bump_texture = mesh$material$bump_texture
  bump_intensity = 1
  if(!is.null(bump_texture)) {
    bump_texture = path.expand(bump_texture)
    if(!is.null(mesh$material$bump_intensity)) {
      bump_intensity = mesh$material$bump_intensity
    } else {
      bump_intensity = 1
    }
  } else {
    bump_texture = ""
  }
  face_color_vals = mesh$material$color
  if(!is.null(face_color_vals)) {
    if(length(face_color_vals) == 1 && texture == "") {
      face_color_vals = rep(face_color_vals, nrow(indices))
      mesh$meshColor = "faces"
    }
    if(!is.matrix(face_color_vals)) {
      color_vals = matrix(convert_color(face_color_vals), ncol=3, byrow=TRUE)
    } else {
      color_vals = face_color_vals
    }
  } else {
    color_vals = matrix(nrow=0,ncol=0)
  }
  if(!is.null(mesh$meshColor)) {
    color_type = switch(mesh$meshColor,"vertices" = 1, "faces" = 2, 3)
    if(color_type == 1) {
      if(is.null(texcoords) && nrow(color_vals) == nrow(vertices)) {
        color_type = 4
      } else {
        if(texture == "" && bump_texture == "") {
          warning("material set as vertex color but no texture or bump map passed--ignoring mesh3d material.")
          color_type = 3
        }
      }
    }
  } else {
    color_type = 3
  }
  if(override_material) {
    color_type = 1
  }
  if(reverse) {
    indices = indices[,c(3,2,1)]
  }
  material_type = switch(material[[1]]$type,
                         "diffuse" = 1,"metal" = 2,"dielectric" = 3, 
                         "oren-nayar" = 4, "light" = 5, "mf" = 6, 
                         "glossy" = 7, "spotlight" = 8, "hair" = 9, "mf-t" = 10)

  mesh_info = list(vertices=vertices,indices=indices,
                   normals=normals,texcoords=texcoords,
                   texture=texture,bump_texture=bump_texture,
                   bump_intensity=bump_intensity,
                   color_vals=color_vals,
                   color_type=color_type,scale_mesh=scale_mesh,
                   material_type = material_type)
  if(verbose) {
    bbox = apply(vertices,2,range)
    message(sprintf("mesh3d Bounding Box: %0.1f-%0.1f x %0.1f-%0.1f x %0.1f-%0.1f", 
                    bbox[1,1],bbox[2,1],bbox[1,2],bbox[2,2],bbox[1,3],bbox[2,3]))
  }
  new_tibble_row(list(x = x, y = y, z = z, 
                      shape = "mesh3d",
                      material = material,
                      shape_info = ray_shape_info(shape_properties = list(NA),
                                                  tricolorinfo = list(NA), 
                                                  fileinfo = NA,
                                                  material_id = NA_integer_,  
                                                  csg_object = list(NA), 
                                                  mesh_info = list(mesh_info),
                                                  flipped = flipped),
                      transforms = ray_transform(angle = list(angle),
                                                 order_rotation = list(order_rotation),
                                                 scale = list(scale),
                                                 group_transform = list(matrix(NA_real_))),
                      animation_info = ray_animated_transform(
                        start_transform_animation = list(matrix(NA_real_)), 
                        end_transform_animation = list(matrix(NA_real_)),
                        start_time = 0, end_time = 1)
                      ))
}

#' Extruded Path Object
#' 
#' Note: Bump mapping with non-diffuse materials does not work correctly, and smoothed normals will be flat when
#' using a bump map.
#'
#' @param points Either a list of length-3 numeric vectors or 3-column matrix/data.frame specifying
#' the x/y/z points that the path should go through.
#' @param x Default `0`. x-coordinate offset for the path.
#' @param y Default `0`. y-coordinate offset for the path.
#' @param z Default `0`. z-coordinate offset for the path.
#' @param polygon Defaults to a circle. A polygon with no holes, specified by a data.frame() parsable by `xy.coords()`. Vertices
#' are taken as sequential rows. If the polygon isn't closed (the last vertex equal to the first), it will be closed automatically.
#' @param polygon_end Defaults to `polygon`. If specified, the number of vertices should equal the to the number of vertices 
#' of the polygon set in `polygon`. Vertices are taken as sequential rows. If the polygon isn't closed (the last vertex equal to the first), it will be closed automatically.
#' @param breaks Defaults to `20` times the number of control points in the bezier curve.
#' @param closed Default `FALSE`. If `TRUE`, the path will be closed by smoothly connecting the first
#' and last points, also ensuring the final polygon is aligned to the first.
#' @param closed_smooth Default `TRUE`. If `closed = TRUE`, this will ensure C2 (second derivative) 
#' continuity between the ends. If `closed = FALSE`, the curve will only have C1 (first derivative)
#' continuity between the ends.
#' @param polygon_add_points Default `0`. Positive integer specifying the number of points to fill in between polygon
#' vertices. Higher numbers can give smoother results (especially when combined with `smooth_normals = TRUE`.
#' @param twists Default `0`. Number of twists in the polygon from one end to another.
#' @param straight Default `FALSE`. If `TRUE`, straight lines will be used to connect the points instead
#' of bezier curves. 
#' @param precomputed_control_points Default `FALSE`. If `TRUE`, `points` argument will expect
#' a list of control points calculated with the internal rayrender function `rayrender:::calculate_control_points()`.
#' @param width Default `0.1`. Curve width. If `width_ease == "spline"`, `width` is specified in a format that can be read by 
#' `xy.coords()` (with `y` as the width), and the `x` coordinate is between `0` and `1`, this can also specify the exact 
#' positions along the curve for the corresponding width values. If a numeric vector, specifies the different values of the width evenly along the curve.
#' If not a single value, `width_end` will be ignored.
#' @param width_end Default `NA`. Width at end of path. Same as `width`, unless specified. Ignored if multiple width values 
#' specified in `width`.
#' @param width_ease Default `spline`. Ease function between width values. Other options: `linear`, `quad`, `cubic`, `exp`.
#' @param u_min Default `0`. Minimum parametric coordinate for the path. If `closed = TRUE`, values greater than one will refer to the beginning 
#' of the loop (but the path will be generated as two objects).
#' @param u_max Default `1`. Maximum parametric coordinate for the path. If `closed = TRUE`, values greater than one will refer to the beginning 
#' of the loop (but the path will be generated as two objects).
#' @param texture_repeats Default `1`. Number of times to repeat the texture along the length of the path.
#' @param smooth_normals Default `FALSE`. Whether to smooth the normals of the polygon to remove sharp angles.
#' @param linear_step Default `FALSE`. Whether the polygon intervals should be set at linear intervals,
#' rather than intervals based on the underlying bezier curve parameterization.
#' @param end_caps Default `c(TRUE, TRUE)`. Specifies whether to add an end cap to the beginning and end of a path.
#' @param material Default  \code{\link{diffuse}}. The material, called from one of the material 
#' functions.
#' @param material_caps Defaults to the same material set in `material`.
#' Note: emissive objects may not currently function correctly when scaled.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' @importFrom  grDevices col2rgb
#'
#' @return Single row of a tibble describing the cube in the scene.
#' @export
#'
#' @examples
#' if(run_documentation()) {
#' #Specify the points for the path to travel though and the ground material
#' points = list(c(0,0,1),c(-0.5,0,-1),c(0,1,-1),c(1,0.5,0),c(0.6,0.3,1))
#' ground_mat = material=diffuse(color="grey50",
#'                               checkercolor = "grey20",checkerperiod = 1.5)
#' }
#' if(run_documentation()) {
#' #Default path shape is a circle
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, width=0.25, 
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0.5),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Change the width evenly along the tube
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, width=0.25, 
#'                            width_end = 0.5,
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0.5),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Change the width along the full length of the tube
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, 
#'                            width=0.25*sinpi(0:72*20/180),
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0.5),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Specify the exact parametric x positions for the width values:
#' custom_width = data.frame(x=c(0,0.2,0.5,0.8,1), y=c(0.25,0.5,0,0.5,0.25))
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, 
#'                            width=custom_width,
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0.5),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Generate a star polygon
#' angles = seq(360,0,length.out=21)
#' xx = c(rep(c(1,0.75,0.5,0.75),5),1) * sinpi(angles/180)/4
#' yy = c(rep(c(1,0.75,0.5,0.75),5),1) * cospi(angles/180)/4
#' star_polygon = data.frame(x=xx,y=yy)
#' 
#' #Extrude a path using a star polygon
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, width=0.5, 
#'                            polygon = star_polygon,
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,1),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Specify a circle polygon
#' angles = seq(360,0,length.out=21)
#' xx = sinpi(angles/180)/4
#' yy = cospi(angles/180)/4
#' circ_polygon = data.frame(x=xx,y=yy)
#' 
#' #Transform from the circle polygon to the star polygon and change the end cap material
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, width=0.5, 
#'                            polygon=circ_polygon, polygon_end = star_polygon,
#'                            material_cap  = diffuse(color="white"),
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0.5),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Add three and a half twists along the path, and make sure the breaks are evenly spaced
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, width=0.5, twists = 3.5,
#'                            polygon=star_polygon, linear_step = TRUE, breaks=360,
#'                            material_cap  = diffuse(color="white"),
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Smooth the normals for a less sharp appearance:
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, width=0.5, twists = 3.5,
#'                            polygon=star_polygon, 
#'                            linear_step = TRUE, breaks=360,
#'                            smooth_normals = TRUE,
#'                            material_cap  = diffuse(color="white"),
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Only generate part of the curve, specified by the u_min and u_max arguments
#' generate_studio(depth=-0.4,material=ground_mat) %>%
#'   add_object(extruded_path(points = points, width=0.5, twists = 3.5,
#'                            u_min = 0.2, u_max = 0.8,
#'                            polygon=star_polygon, linear_step = TRUE, breaks=360,
#'                            material_cap  = diffuse(color="white"),
#'                            material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=3,z=5,x=2,material=light(intensity=15))) %>% 
#'   render_scene(lookat=c(0.3,0.5,0),fov=12, width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
#' if(run_documentation()) {
#' #Render a Mobius strip with 1.5 turns 
#' points = list(c(0,0,0),c(0.5,0.5,0),c(0,1,0),c(-0.5,0.5,0))
#' square_polygon = matrix(c(-1, -0.1, 0,
#'                            1, -0.1, 0,
#'                            1,  0.1, 0,
#'                           -1,  0.1, 0)/10, ncol=3,byrow = T)
#'
#' generate_studio(depth=-0.2,
#'                material=diffuse(color = "dodgerblue4", checkercolor = "#002a61",
#'                                 checkerperiod = 1)) %>%
#'  add_object(extruded_path(points = points,  polygon=square_polygon, closed = TRUE,
#'                           linear_step = TRUE, twists = 1.5, breaks = 720, 
#'                           material = diffuse(noisecolor = "black", noise = 10, 
#'                                              noiseintensity = 10))) %>%
#'  add_object(sphere(y=20,x=0,z=21,material=light(intensity = 1000))) %>% 
#'  render_scene(lookat=c(0,0.5,0), fov=10, samples=128, sample_method = "sobol_blue",
#'               width = 800, height=800)
#' }
#' if(run_documentation()) {
#' #Create a green glass tube with the dielectric priority interface
#' #and fill it with a purple neon tube light
#' generate_ground(depth=-0.4,material=diffuse(color="grey50",
#'                                             checkercolor = "grey20",checkerperiod = 1.5)) %>%
#'   add_object(extruded_path(points = points, width=0.7, linear_step = TRUE, 
#'                            polygon = star_polygon, twists = 2, closed = TRUE,
#'                            polygon_end = star_polygon, breaks=500,
#'                            material=dielectric(priority = 1, refraction = 1.2, 
#'                                                attenuation=c(1,0.3,1),
#'                                                attenuation_intensity=20))) %>% 
#'   add_object(extruded_path(points = points, width=0.4, linear_step = TRUE,
#'                            polygon = star_polygon,twists = 2, closed = TRUE,
#'                            polygon_end = star_polygon, breaks=500,
#'                            material=dielectric(priority = 0,refraction = 1))) %>%  
#'   add_object(extruded_path(points = points, width=0.05, closed = TRUE,
#'                            material=light(color="purple", intensity = 5,
#'                                           importance_sample = FALSE))) %>%
#'   add_object(sphere(y=10,z=-5,x=0,radius=5,material=light(color = "white",intensity = 5))) %>%
#'   render_scene(lookat=c(0,0.5,1),fov=10, 
#'                width=800,height=800, clamp_value = 10,
#'                aperture=0.025, samples=128, sample_method="sobol_blue")
#' }
extruded_path = function(points, x = 0, y = 0, z = 0, 
                         polygon = NA, polygon_end = NA, breaks=NA,
                         closed = FALSE, closed_smooth = TRUE, 
                         polygon_add_points = 0,
                         twists = 0, 
                         texture_repeats = 1, 
                         straight = FALSE, precomputed_control_points = FALSE, 
                         width = 1, width_end = NA, width_ease = "spline",
                         smooth_normals = FALSE,
                         u_min = 0, u_max = 1, linear_step = FALSE, end_caps = c(TRUE, TRUE),
                         material = diffuse(), material_caps = NA, angle = c(0, 0, 0),
                         order_rotation = c(1, 2, 3),
                         flipped = FALSE, scale = c(1,1,1)) {
  if(closed) {
    if((u_max - u_min) >= 1) {
      u_min = 0
      u_max = 1
      end_caps = c(FALSE, FALSE)
    } else {
      while(u_min >= 1) {
        u_min = u_min - 1
      }
      while(u_max > 1) {
        u_max = u_max - 1
      }
      if(u_max == 0 && u_min > 0) {
        u_max = 1
      }
    }
    if(u_max < u_min) {
      combined_return = add_object(
        extruded_path(points=points,x=x,y=y,z=z,polygon=polygon,polygon_end=polygon_end,
                                                 breaks=breaks,closed = closed, closed_smooth = closed_smooth,
                                                 polygon_add_points = polygon_add_points,
                                                 twists = twists,
                                                 texture_repeats = texture_repeats,
                                                 straight = straight, precomputed_control_points = precomputed_control_points,
                                                 width = width, width_end = width_end, width_ease = width_ease,
                                                 smooth_normals = smooth_normals,
                                                 u_min = 0, u_max = u_max, linear_step = linear_step, end_caps = c(FALSE, TRUE),
                                                 material = material, material_caps = material_caps, angle = angle,
                                                 order_rotation = order_rotation,
                                                 flipped = flipped, scale = scale),
                                   extruded_path(points=points,x=x,y=y,z=z,polygon=polygon,polygon_end=polygon_end,
                                                 breaks=breaks,closed = closed, closed_smooth = closed_smooth,
                                                 polygon_add_points = polygon_add_points,
                                                 twists = twists,
                                                 texture_repeats = texture_repeats,
                                                 straight = straight, precomputed_control_points = precomputed_control_points,
                                                 width = width, width_end = width_end, width_ease = width_ease,
                                                 smooth_normals = smooth_normals,
                                                 u_min = u_min, u_max = 1, linear_step = linear_step, end_caps = c(TRUE, FALSE),
                                                 material = material, material_caps = material_caps, angle = angle,
                                                 order_rotation = order_rotation,
                                                 flipped = flipped, scale = scale)
                                   )
      material_id_old = get("max_material_id", envir = ray_environment)
      material_id_new = material_id_old + 1L
      assign("max_material_id", material_id_new, envir = ray_environment)
      combined_return$shape_info[[1]]$material_id = material_id_new
      combined_return$shape_info[[2]]$material_id = material_id_new
      
      return(combined_return)
    }
  }
  if(is.null(dim(polygon))) {
    angles = seq(0,360,length.out=31)
    xx = 1 / 2 * sinpi(angles/180)
    yy = 1 / 2 * cospi(angles/180)
    polygon = as.matrix(data.frame(x=xx,y=yy,z=0))
  } else {
    polygon = grDevices::xy.coords(polygon)
    polygon = data.frame(x=polygon$x,y=polygon$y,z=0)
  }
  same_polygon = TRUE
  if(is.null(dim(polygon_end))) {
    polygon_end = polygon
  } else {
    if(nrow(polygon_end) != nrow(polygon)) {
      stop("`polygon` and `polygon_end` must have same number of vertices")
    }
    same_polygon = FALSE
    polygon_end = grDevices::xy.coords(polygon_end)
    polygon_end = data.frame(x=polygon_end$x,y=polygon_end$y,z=0)
  }
  if(ncol(polygon) == 2) {
    polygon = cbind(polygon,rep(0,nrow(polygon)))
  }
  if(ncol(polygon_end) == 2) {
    polygon_end = cbind(polygon_end,rep(0,nrow(polygon_end)))
  }
  #Process polygon orientation
  i = seq_len(nrow(polygon) - 1)
  reverse_poly = sum(polygon[i,1] * polygon[i+1,2] - 
                     polygon[i+1,1] * polygon[i,2]) > 0
  if(reverse_poly) {
    polygon = polygon[nrow(polygon):1,]
  }
  i = seq_len(nrow(polygon_end) - 1)
  reverse_poly = sum(polygon_end[i,1] * polygon_end[i+1,2] - 
                     polygon_end[i+1,1] * polygon_end[i,2]) > 0
  if(reverse_poly) {
    polygon_end = polygon_end[nrow(polygon_end):1,]
  }
  if(any(polygon[1,] != polygon[nrow(polygon),])) {
    polygon = rbind(polygon,polygon[1,])
  }
  if(any(polygon_end[1,] != polygon_end[nrow(polygon_end),])) {
    polygon_end = rbind(polygon_end,polygon_end[1,])
  }
  if(polygon_add_points > 0) {
    polygon = add_points_polygon(polygon,polygon_add_points)
    polygon_end = add_points_polygon(polygon_end,polygon_add_points)
  }
  normal_poly = list()
  normal_poly_end = list()
  if(smooth_normals) {
    revisit = rep(FALSE, nrow(polygon))
    for(i in seq_len(nrow(polygon))) {
      if(i == 1 || i == nrow(polygon)) {
        vert1 = nrow(polygon)-1
        vert2 = 1
        vert3 = 2
      } else if (i < nrow(polygon)) {
        vert1 = i-1
        vert2 = i
        vert3 = i+1
      }
      edge1 = polygon[vert1,] - polygon[vert2,]
      edge2 = polygon[vert3,] - polygon[vert2,]
      edge1 = edge1 / sqrt(sum(edge1^2))
      edge2 = edge2 / sqrt(sum(edge2^2))
      
      if(abs(abs(sum(edge1*edge2))-1) < 1e-8) {
        temp_norm = NA
        revisit[i] = TRUE
      } else {
        temp_norm = as.numeric(edge1 + edge2)
      }
      if(!revisit[i] && cross_prod(edge1,temp_norm)[3] > 0) {
        temp_norm = -temp_norm
      } 
      normal_poly[[i]] = matrix(temp_norm, nrow=1,ncol=3)
      if(!revisit[i]) {
        normal_poly[[i]] = normal_poly[[i]] / sqrt(sum(normal_poly[[i]]^2))
        if(any(is.nan(normal_poly[[i]]))) {
          stop("NaN coordinates found in smoothed normals--is polygon correct?")
        }
      }
    }
    if(any(revisit)) {
      for(i in seq_len(nrow(polygon))) {
        if(revisit[i]) {
          idx1 = i
          idx2 = i
          n = 0L

          inds = c(i)
          while(revisit[idx1]) {
            idx1 = idx1 - 1L
            if(idx1 == 0L) {
              idx1 = nrow(polygon)
            }
            inds = c(idx1,inds)
            n = n + 1L
          }
          while(revisit[idx2]) {
            idx2 = idx2 + 1L
            if(idx2 == nrow(polygon)+1L) {
              idx2 = 1L
            }
            inds = c(inds,idx2)
            n = n + 1L
          }
          norm_prev = normal_poly[[idx1]]
          norm_next = normal_poly[[idx2]]
          tween_norms = slerp(norm_prev, norm_next, n = n+1)
          for(idx in seq_len(length(inds)-1)[-1]) {
            normal_poly[[inds[idx] ]] = tween_norms[[idx]]
            revisit[inds[idx]] = FALSE
          }
        }
      }
    }
    revisit = rep(FALSE, nrow(polygon_end))
    for(i in seq_len(nrow(polygon_end))) {
      if(i == 1 || i == nrow(polygon_end)) {
        vert1 = nrow(polygon_end)-1
        vert2 = 1
        vert3 = 2
      } else if (i < nrow(polygon_end)) {
        vert1 = i-1
        vert2 = i
        vert3 = i+1
      }
      edge1 = polygon_end[vert1,] - polygon_end[vert2,]
      edge2 = polygon_end[vert3,] - polygon_end[vert2,]
      edge1 = edge1 / sqrt(sum(edge1^2))
      edge2 = edge2 / sqrt(sum(edge2^2))
      
      if(abs(abs(sum(edge1*edge2))-1) < 1e-8) {
        temp_norm = NA
        revisit[i] = TRUE
      } else {
        temp_norm = as.numeric(edge1 + edge2)
      }
      if(!revisit[i] && cross_prod(edge1,temp_norm)[3] > 0) {
        temp_norm = -temp_norm
      } 
      normal_poly_end[[i]] = matrix(temp_norm, nrow=1,ncol=3)
      if(!revisit[i]) {
        normal_poly_end[[i]] = normal_poly_end[[i]] / sqrt(sum(normal_poly_end[[i]]^2))
        if(any(is.nan(normal_poly_end[[i]]))) {
          stop("NaN coordinates found in smoothed normals--is polygon_end correct?")
        }
      }
    }
    if(any(revisit)) {
      for(i in seq_len(nrow(polygon_end))) {
        if(revisit[i]) {
          idx1 = i
          idx2 = i
          n = 0L
          
          inds = c(i)
          while(revisit[idx1]) {
            idx1 = idx1 - 1L
            if(idx1 == 0L) {
              idx1 = nrow(polygon_end)
            }
            inds = c(idx1,inds)
            n = n + 1L
          }
          while(revisit[idx2]) {
            idx2 = idx2 + 1L
            if(idx2 == nrow(polygon_end)+1L) {
              idx2 = 1L
            }
            inds = c(inds,idx2)
            n = n + 1L
          }
          norm_prev = normal_poly_end[[idx1]]
          norm_next = normal_poly_end[[idx2]]
          tween_norms = slerp(norm_prev, norm_next, n = n+1)
          for(idx in seq_len(length(inds)-1)[-1]) {
            normal_poly_end[[inds[idx] ]] = tween_norms[[idx]]
            revisit[inds[idx]] = FALSE
          }
        }
      }
    }
    normal_polys = do.call(rbind,normal_poly)
    normal_polys_end = do.call(rbind,normal_poly_end)
  }
  if(u_min == u_max) {
    return()
  }
  
  if(inherits(points,"numeric")) {
    stop("Input must either be list, matrix, or data.frame, not numeric.")
  }
  if(inherits(points,"list") && !precomputed_control_points) {
    if(any(unlist(lapply(points,(function(x) length(x) != 3))))) {
      stop("If `points` is a list, each entry must be a length-3 vector")
    }
    points = do.call(rbind,points)
  }
  if(inherits(points,"data.frame")) {
    points = as.matrix(points)
  }
  if(is.array(points)) {
    if(nrow(points) == 1) {
      stop("Only one point passed, no path specified.")
    }
    if(nrow(points) == 2 && closed) {
      closed=FALSE
    }
  }
  if(!precomputed_control_points) {
    if(inherits(points,"matrix")) {
      if(ncol(points) == 3) {
        if(closed && closed_smooth) {
          points = rbind(points[c(nrow(points)-2,nrow(points)-1,nrow(points)),], points, points[1:3,])
        }
        if(!straight) {
          full_control_points = calculate_control_points(points)
        } else {
          full_control_points = calculate_control_points_straight(points)
        }
        if(closed && closed_smooth) {
          full_control_points[[length(full_control_points)]] = NULL
          full_control_points[[length(full_control_points)]] = NULL
          full_control_points[[1]] = NULL
          full_control_points[[1]] = NULL
          full_control_points[[1]] = NULL
        }
      } else {
        stop("If points a matrix or data.frame, must have 3 columns")
      }
    } else {
      stop("points not of supported type (function expects matrix/data.frame/list, got ", class(points),")")
    }
  } else {
    full_control_points = points
  }
  if(closed && !closed_smooth) {
    first_point = full_control_points[[1]]
    last_point = full_control_points[[length(full_control_points)]]
    full_control_points[[length(full_control_points) + 1]] = last_point
    full_control_points[[length(full_control_points)]][4,] = first_point[1,]
    full_control_points[[length(full_control_points)]][3,] = 2*first_point[1,] - first_point[2,]
    full_control_points[[length(full_control_points)]][2,] = 2*last_point[4,]  - last_point[3,]
    full_control_points[[length(full_control_points)]][1,] = last_point[4,]
  }
  if(is.na(breaks)) {
    breaks = length(full_control_points) * 20
  }
  u_min_segment = floor(breaks * u_min)
  u_max_segment = floor(breaks * u_max)
  
  seg_begin = u_min_segment + 1
  seg_end = u_max_segment + 1
  if(seg_end > breaks-1) {
    seg_end = breaks-1
  }

  if(is.na(width_end) && is.numeric(width)) {
    width_end = width[1]
  }
  if(length(width) == 1 && is.numeric(width)) {
    width =  c(width,width_end)
  }
  if(is.null(dim(width))) {
    width = data.frame(x=seq(0,1,length.out=length(width)),y=width)
  }
  width = grDevices::xy.coords(width)

  if(linear_step) {
    dist_df = calculate_distance_along_bezier_curve(full_control_points,20)
    t_vals = stats::predict(stats::smooth.spline(dist_df$total_dist,dist_df$t),
                            seq(0,max(dist_df$total_dist),length.out=breaks))$y * length(full_control_points)
    #Numerical precision fix
    t_vals[length(t_vals)] = length(full_control_points)
  } else {
    t_vals = seq(0, length(full_control_points), length.out=breaks)
    #Numerical precision fix
    t_vals[length(t_vals)] = length(full_control_points)
  }
  if(width_ease != "spline") {
    width_vals = tween(width$y, n = breaks,  ease = width_ease)
  } else {
    width_vals = stats::spline(width, n = breaks)$y
  }
  width_vals = abs(width_vals)
  morph_vals = seq(0, 1, length.out=breaks)
  t_init = t_vals[seg_begin]
  if(t_init < 0 || closed) {
    t_init = 0
  }
  initial_deriv = eval_bezier_deriv(full_control_points[[1]],t_init)
  initial_2nd_deriv = eval_bezier_2nd_deriv(full_control_points[[1]],t_init)
  if(all(abs(initial_2nd_deriv) < 1e-6)) {
    if(any(cross_prod(c(1,0,0),initial_deriv) != 0)) {
      initial_2nd_deriv = cross_prod(c(1,0,0),initial_deriv)
    } else {
      initial_2nd_deriv = cross_prod(c(0,1,0),initial_deriv)
    }
  } 
  t_vec = initial_deriv/sqrt(sum(initial_deriv*initial_deriv))
  s_vec = cross_prod(initial_deriv,initial_2nd_deriv)
  s_vec = s_vec/sqrt(sum(s_vec*s_vec))
  r_vec = cross_prod(s_vec,t_vec)
  t_vec0 = t_vec
  s_vec0 = s_vec
  r_vec0 = r_vec
  
  vertices = list()
  texcoords = list()
  normals = list()
  poly_tex = seq(0,1,length.out=nrow(polygon))
  counter = 1
  if(closed && !closed_smooth) {
    twist_amount = calculate_final_twist(full_control_points,
                                         breaks, t_vals,
                                         t_vec, s_vec, r_vec)
    end_angle_r = twists*2*pi + twist_amount[1]
    # end_angle_s = twist_amount[2]
    # end_angle_t = twist_amount[3]
    end_angle_s = 0
    end_angle_t = 0
  } else {
    end_angle_r = twists*2*pi
    end_angle_s = 0
    end_angle_t = 0
  }
  
  for(i in seq(1,seg_end,by=1)) {
    t_val0 = t_vals[i]
    if(t_val0 < 0) {
      t_val0 = 0
    }
    t_val1 = t_vals[i+1]
    if(t_val1 < 0) {
      t_val1 = 0
    }
    width_temp = width_vals[i]
    
    temp_poly = morph_vals[i] * polygon_end + (1-morph_vals[i]) * polygon
    temp_angle_r = morph_vals[i] * end_angle_r
    temp_angle_s = morph_vals[i] * end_angle_s
    temp_angle_t = morph_vals[i] * end_angle_t
    
    ca = cos(temp_angle_t)
    sa = sin(temp_angle_t)
    cb = cos(temp_angle_s)
    sb = sin(temp_angle_s)
    cc = cos(temp_angle_r)
    sc = sin(temp_angle_r)
    
    twist_mat = matrix(c(cb*cc, sa*sb*cc-ca*sc, ca*sb*cc + sa*sc,
                         cb*sc, sa*sb*sc+ca*cc, ca*sb*sc - sa*cc,
                         -sb,              sa*cb,            ca*cb), nrow=3,ncol=3,byrow=TRUE)
    
    rot_mat = matrix(c(s_vec,r_vec,t_vec),3,3)
    t_temp0 = t_val0-floor(t_val0)
    if(i != breaks-1) {
      t_temp1 = t_val1-floor(t_val1)
    } else {
      t_temp1 = 1
    }

    i0 = floor(t_val0) + 1
    if(i != breaks-1) {
      i1 = floor(t_val1) + 1
    } else {
      i1 = max(c(1,floor(t_val1+1e-8)))
    }


    cp0 = full_control_points[[i0]]
    if(i1 <= length(full_control_points)) {
      cp1 = full_control_points[[i1]]
    } else {
      cp1 = cp0
    }
    
    x0 = eval_bezier(cp0,t_temp0)
    x1 = eval_bezier(cp1,t_temp1)
    
    if(i < length(width_vals)) {
      width_norm = -(width_vals[i+1]-width_vals[i])/sqrt(sum((x1-x0)^2))
    } else if (i == length(width_vals)) {
      width_norm = -(width_vals[i]-width_vals[i-1])/sqrt(sum((x1-x0)^2))
    }
    if(i >= seg_begin) {
      vertices[[counter]] = matrix(x0,ncol=3,nrow=nrow(polygon), byrow=TRUE) + 
        t((rot_mat %*% twist_mat %*% t(temp_poly*width_temp)))
      texcoords[[counter]] = matrix(c(poly_tex,rep(morph_vals[i] * texture_repeats,nrow(polygon))), 
                                    ncol=2,nrow=nrow(polygon))
      if(smooth_normals) {
        temp_norm = morph_vals[i] * normal_polys_end + (1-morph_vals[i]) * normal_polys
        temp_norm[,3] = width_norm
        norm_transform = t(solve(rot_mat %*% twist_mat))
        normals[[counter]] = t((norm_transform %*% t(temp_norm)))
      }
      
      counter = counter + 1
    }
    #Evaluate next set of vectors
    v1 = x1-x0
    c1 = sum(v1*v1)
    rl = r_vec - (2 / c1) * sum(v1*r_vec) * v1
    tl = t_vec - (2 / c1) * sum(v1*t_vec) * v1
    
    next_deriv = eval_bezier_deriv(cp1,t_temp1)
    t_vec_prev = next_deriv/sqrt(sum(next_deriv*next_deriv))

    v2 = t_vec_prev - tl
    c2 = sum(v2*v2)
    if(c2 != 0) {
      t_vec = t_vec_prev
      r_vec = rl - (2 / c2) * sum(v2*rl) * v2
      s_vec = cross_prod(t_vec,r_vec)
    }
  }
  rot_mat = matrix(c(s_vec,r_vec,t_vec),3,3)
  
  temp_angle_r = morph_vals[seg_end] * end_angle_r
  temp_angle_s = morph_vals[seg_end] * end_angle_s
  temp_angle_t = morph_vals[seg_end] * end_angle_t
  ca = cos(temp_angle_t)
  sa = sin(temp_angle_t)
  cb = cos(temp_angle_s)
  sb = sin(temp_angle_s)
  cc = cos(temp_angle_r)
  sc = sin(temp_angle_r)
  
  twist_mat = matrix(c(cb*cc, sa*sb*cc-ca*sc, ca*sb*cc + sa*sc,
                       cb*sc, sa*sb*sc+ca*cc, ca*sb*sc - sa*cc,
                       -sb,            sa*cb,            ca*cb), nrow=3,ncol=3,byrow=TRUE)
  vertices[[counter]] = matrix(x1,ncol=3,nrow=nrow(polygon), byrow=TRUE) +
    t((rot_mat %*% twist_mat %*% t(polygon_end*width_vals[seg_end+1])))
  texcoords[[counter]] = matrix(c(poly_tex,rep(1 * texture_repeats,nrow(polygon))),
                                ncol=2,nrow=nrow(polygon))
  if(smooth_normals) {
    norm_transform = t(solve(rot_mat %*% twist_mat))
    normals[[counter]] = t((norm_transform %*% t(normal_polys_end)))
  }
  mesh = list()
  mesh_caps = list()
  end_caps[1] = end_caps[1] & width_vals[1] > 0
  end_caps[2] = end_caps[2] & width_vals[length(width_vals)] > 0
  
  vb = do.call(rbind,vertices)
  tex = do.call(rbind,texcoords)
  if(smooth_normals) {
    normal_mat = do.call(rbind,normals)
  }
  
  faces = (nrow(polygon)-1)*2*(length(vertices)-1)
  if(faces == 0) {
    return()
  }
  band_faces  = (nrow(polygon)-1)*2
  it = matrix(0,nrow=3,ncol=faces)
  polyadd = c(0,0,nrow(polygon),nrow(polygon),0,nrow(polygon))
  single_structure = c(1,2,1,1,2,2)
  single_band = matrix(0,nrow=6,ncol=band_faces/2)
  for(i in seq_len((nrow(polygon)-1))) {
    single_band[,i] = polyadd + single_structure + (i-1)
  }
  single_band = matrix(as.vector(single_band),nrow=3)
  for(i in seq_len(length(vertices)-1)) {
    it[,(1+(i-1)*band_faces):(i*band_faces)] = single_band + (i-1) * nrow(polygon)
  }
  cap_it_start = matrix(decido::earcut(polygon[,1:2]),nrow=3)
  end_poly_triangulated = decido::earcut(polygon_end[,1:2])
  if(all(end_caps)) {
    cap_it_end_full = matrix(
      rev(end_poly_triangulated) + nrow(polygon)*(length(vertices)-1),
      nrow=3)
    cap_it_end = matrix(rev(end_poly_triangulated) + max(cap_it_start)-min(cap_it_start)+1,
                        nrow=3)
    vb_ends = rbind(vb[seq(1,max(cap_it_start)),], vb[seq(min(cap_it_end_full), max(cap_it_end_full)),]) 
    it_ends = cbind(cap_it_start, cap_it_end)
  } else if (end_caps[1]) {
    vb_ends = vb[seq(1,max(cap_it_start)),]
    it_ends = cap_it_start
  } else if (end_caps[2]) {
    cap_it_end_full = matrix(
      rev(end_poly_triangulated) + nrow(polygon)*(length(vertices)-1),
      nrow=3)
    cap_it_end = matrix(rev(end_poly_triangulated),
                        nrow=3)
    if(min(cap_it_end) == 2) {
      cap_it_end = cap_it_end - 1
    }
    vb_ends = vb[seq(min(cap_it_end_full), max(cap_it_end_full)),]
    it_ends = cap_it_end
  }
  if(any(end_caps)) {
    mesh_caps$vb = t(cbind(vb_ends,rep(1,nrow(vb_ends))))          
    mesh_caps$it = it_ends
  }
  
  mesh$vb = t(cbind(vb,rep(1,nrow(vb))))
  if(any(is.nan(mesh$vb))) {
    stop("NaN coordinates in mesh generated.")
  }
  mesh$it = it
  if(smooth_normals) {
    mesh$normals = t(normal_mat)
  }
  mesh$texcoords = t(tex)
  if(!is.na(material[[1]]$image[[1]])) {
    mesh$material$texture = material[[1]]$image[[1]]
    mesh$meshColor = "vertices"
  }
  if(!is.na(material[[1]]$bump_texture[[1]])) {
    mesh$material$bump_texture = material[[1]]$bump_texture[[1]]
    mesh$material$bump_intensity =  material[[1]]$bump_intensity
    mesh$meshColor = "vertices"
  } else {
    mesh$material$bump_texture = ""
    mesh$material$bump_intensity =  1
  }
  class(mesh) = "mesh3d"
  class(mesh_caps) = "mesh3d"
  same_material = FALSE
  if(is.null(dim(material_caps))) {
    material_caps = material
    same_material = TRUE
  }
  material_id_old = get("max_material_id", envir = ray_environment)
  material_id_new = material_id_old + 1L
  assign("max_material_id", material_id_new, envir = ray_environment)
  if(any(end_caps)) {
    return_scene = add_object(mesh3d_model(mesh,x=x,y=y,z=z, override_material = TRUE,
                                angle=angle, order_rotation = order_rotation, flipped=flipped,
                                scale=scale, material = material),
                              mesh3d_model(mesh_caps,x=x,y=y,z=z, material = material_caps,
                                           override_material = TRUE,
                                           angle=angle, order_rotation = order_rotation, flipped=flipped,
                                           scale=scale))
    if(same_material) {
      for(i in seq_len(nrow(return_scene))) {
        return_scene$shape_info[[i]]$material_id = material_id_new
      }
    }
  } else {
    return_scene = mesh3d_model(mesh,x=x,y=y,z=z, override_material = TRUE,
                                           angle=angle, order_rotation = order_rotation, flipped=flipped,
                                           scale=scale, material = material)
    if(same_material) {
      for(i in seq_len(nrow(return_scene))) {
        return_scene$shape_info[[i]]$material_id = material_id_new
      }
    }
  }
  return(return_scene)
}

#' `raymesh` model
#' 
#' Load an `raymesh` object, as specified in the `rayvertex` package. 
#'
#' @param mesh A `raymesh` object. Pulls the vertex, index, texture coordinates, 
#' normals, and material information. 
#' @param x Default `0`. x-coordinate to offset the model.
#' @param y Default `0`. y-coordinate to offset the model.
#' @param z Default `0`. z-coordinate to offset the model.
#' @param verbose Default `FALSE`. If `TRUE`, prints information about the mesh to the console.
#' @param override_material Default `TRUE`. If `TRUE`, overrides the material specified in the 
#' `raymesh` object with the one specified in `material`.
#' @param flip_transmittance Default `TRUE`. Flips `(1-t)` the transmittance values to match the way the colors
#' would be interpreted in a rasterizer (where it specifies the transmitted color). Turn off to specify
#' the attenuation values directly.
#' @param calculate_consistent_normals Default `TRUE`. Whether to calculate consistent vertex normals to prevent energy 
#' loss at edges.
#' @param importance_sample_lights Default `TRUE`. Whether to importance sample lights specified in the OBJ material
#' (objects with a non-zero Ke MTL material).
#' @param validate_mesh Default `TRUE`. Validates the `raymesh` object using `rayvertex::validate_mesh()` 
#' before parsing to ensure correct parsing. Set to `FALSE` to speed up scene construction if `raymesh_model()` 
#' is taking a long time (Note: this does not affect rendering time).
#' @param material Default  \code{\link{diffuse}}, but ignored unless `override_material = TRUE`. The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}. 
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param flipped Default `FALSE`. Whether to flip the normals.
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#' Note: emissive objects may not currently function correctly when scaled.
#' 
#' @return Single row of a tibble describing the raymesh model in the scene.
#' @export
#'
#' @examples
#' #Render a simple raymesh object
#' library(rayvertex)
#' if(run_documentation()) {
#' raymesh_model(sphere_mesh(position = c(-1, 0, 0),
#'               material = material_list(transmittance = "red"))) %>%
#'   add_object(generate_ground(material = diffuse(checkercolor="grey20"))) %>%
#'   render_scene(fov = 30, samples=128, sample_method="sobol_blue")
#' }
#' 
#' # We create a complex rayvertex mesh, using the `rayvertex::add_shape` function which
#' # creates a new `raymesh` object out of individual `raymesh` objects
#' rm_scene = sphere_mesh(position = c(-1, 0, 0),
#'             material = material_list(transmittance = "red")) %>% 
#'     add_shape(sphere_mesh(position = c(1, 0, 0),
#'             material = material_list(transmittance = "green", ior = 1.5)))
#'
#' # Pass the single raymesh object to `raymesh_model()`
#' # `raymesh_model()`
#' if(run_documentation()) {
#' raymesh_model(rm_scene) %>%
#'   add_object(generate_ground(material = diffuse(checkercolor="grey20"))) %>%
#'   render_scene(fov = 30, samples=128, sample_method="sobol_blue")
#' }
#' 
#' # Set `flip_transmittance = FALSE` argument to specify attenuation coefficients directly
#' # (as specified in the `dielectric()` material). We change the material's numerical attenuation
#' # constants using `rayvertex::change_material`
#' rm_scene_new= change_material(rm_scene, transmittance = c(1,2,0.3), id = 1) %>% 
#'   change_material(transmittance = c(3,1,2), id = 2)
#' if(run_documentation()) {
#' raymesh_model(rm_scene_new, flip_transmittance = FALSE) %>%
#'   add_object(generate_ground(material = diffuse(checkercolor="grey20"))) %>%
#'   render_scene(fov = 30, samples=128, sample_method="sobol_blue")
#' }
#'
#' # Override the material specified in the `raymesh` object and render the scene
#' if(run_documentation()) {
#' raymesh_model(rm_scene,
#'               material = dielectric(attenuation = "dodgerblue2", attenuation_intensity = 4), 
#'   override_material = TRUE) %>%
#'   add_object(generate_ground(material = diffuse(checkercolor="grey20"))) %>%
#'   render_scene(fov = 30, samples=128, sample_method="sobol_blue")
#' }
#'
#' # Adjusting the scale, position, and rotation parameters of the `raymesh` model
#' if(run_documentation()) {
#' raymesh_model(rm_scene,
#'               x = 0, y = 0.5, z = -1, angle = c(0, 0, 20)) %>%
#'   add_object(generate_ground(material = diffuse(checkercolor="grey20"))) %>%
#'   render_scene(fov = 30,lookat=c(0,0.5,0), samples=128, sample_method="sobol_blue")
#' }
raymesh_model = function(mesh, x = 0, y = 0, z = 0, 
                         flip_transmittance = TRUE, verbose = FALSE, 
                         importance_sample_lights = FALSE,
                         calculate_consistent_normals = TRUE,
                         override_material = TRUE, material = diffuse(), 
                         angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
                         flipped = FALSE, scale = c(1,1,1), validate_mesh = TRUE) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
  if(!inherits(mesh,"ray_mesh")) {
    stop("mesh must be of class 'ray_mesh': actual class is ", class(mesh))
  }
  if(validate_mesh) {
    raymesh = rayvertex::validate_mesh(mesh)
  }
  new_tibble_row(list(x = x, y = y, z = z, 
                      shape = "raymesh", 
                      material = material,
                      shape_info = ray_shape_info(shape_properties = list(importance_sample_lights = importance_sample_lights,
                                                                          calculate_consistent_normals = calculate_consistent_normals,
                                                                          override_material = override_material,
                                                                          flip_transmittance = flip_transmittance),
                                                  tricolorinfo = list(NA), 
                                                  fileinfo = NA,
                                                  material_id = NA_integer_,  
                                                  csg_object = list(NA), 
                                                  mesh_info = list(raymesh),
                                                  flipped = flipped),
                      transforms = ray_transform(angle = list(angle),
                                                 order_rotation = list(order_rotation),
                                                 scale = list(scale),
                                                 group_transform = list(matrix(NA_real_))),
                      animation_info = ray_animated_transform(
                        start_transform_animation = list(matrix(NA_real_)), 
                        end_transform_animation = list(matrix(NA_real_)),
                        start_time = 0, end_time = 1)
                      ))
}
