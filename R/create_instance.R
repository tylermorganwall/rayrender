#' Create Instances of an Object
#' 
#' This creates multiple instances of the `ray_scene` passed, each with it's own transformation applied (measured
#' from the origin of the ray_scene). This means the scene only uses the memory of the object once and each 
#' copy only requires a 4x4 matrix in memory.
#'
#' @param ray_scene A `ray_scene` object to be copied at the specified transformed coordinates. 
#' @param x Default `0`. A vector of x-coordinates to offset the instances. Note that this can also be a 3 column matrix or
#' `data.frame()` parsable by `xyz.coords()`: if so, the other axes will be ignored.
#' @param y Default `0`. A vector of y-coordinates to offset the instances.
#' @param z Default `0`. A vector of z-coordinates to offset the instances.
#' @param angle_x Default `0`. A vector of angles around the x axis to rotate the instances.
#' @param angle_y Default `0`. A vector of angles around the y axis to rotate the instances.
#' @param angle_z Default `0`. A vector of angles around the z axis to rotate the instances.
#' @param scale_x Default `0`. A vector of values around the scale the instances on the x-axis.
#' @param scale_y Default `0`. A vector of values around the scale the instances on the y-axis.
#' @param scale_z Default `0`. A vector of values around the scale the instances on the z-axis.
#' @param material Default  \code{\link{diffuse}}.The material, called from one of the material 
#' functions \code{\link{diffuse}}, \code{\link{metal}}, or \code{\link{dielectric}}. 
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z" axes. 
#' @param flipped Default `FALSE`. Whether to flip the normals.
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
create_instances = function(ray_scene, 
                            x = 0, y = 0, z = 0, 
                            angle_x = 0, angle_y = 0, angle_z = 0,
                            scale_x = 1, scale_y = 1, scale_z = 1,
                            material = diffuse(), 
                            order_rotation = c(1, 2, 3), 
                            flipped = FALSE) {
  ray_scene_processed = process_scene(ray_scene, process_material_ids = FALSE)

  # Use data.frame for recycling
  x_is_df = inherits(x, "data.frame") || inherits(x, "matrix")
  
  if(x_is_df && ncol(x) == 3) {
    xyz = xyz.coords(x=x)
  } else {
    xyz = data.frame(x=x,y=y,z=z)
  }
  x = xyz$x
  y = xyz$y
  z = xyz$z
  angle_is_df = inherits(angle_x, "data.frame") || inherits(angle_x, "matrix")
  if(angle_is_df && ncol(angle_x) == 3) {
    angles = xyz.coords(angle_x)
    stopifnot(length(angles$x) == length(x))
  } else {
    angles = data.frame(x_val=x,x=angle_x,y=angle_y,z=angle_z)
  }
  angle_x = angles$x
  angle_y = angles$y
  angle_z = angles$z
  
  scale_is_df = inherits(scale_x, "data.frame") || inherits(scale_x, "matrix")
  if(scale_is_df && ncol(scale_x) == 3) {
    scales = xyz.coords(scale_x)
    stopifnot(length(scales$x) == length(x))
  } else {
    scales = data.frame(x_val=x,x=scale_x,y=scale_y,z=scale_z)
  }
  scale_x = scales$x
  scale_y = scales$y
  scale_z = scales$z
  # Add detecting importance sampling and error
  # 
  # 
  new_tibble_row(list(x = 0, y = 0, z = 0, 
                      shape = "instance",
                      material = material,
                      shape_info = ray_shape_info(shape_properties = list(original_scene = list(ray_scene_processed$scene), 
                                                                          x_values = x,
                                                                          y_values = y,
                                                                          z_values = z,
                                                                          angle_x = angle_x,
                                                                          angle_y = angle_y,
                                                                          angle_z = angle_z,
                                                                          scale_x = scale_x,
                                                                          scale_y = scale_y,
                                                                          scale_z = scale_z,
                                                                          any_light = ray_scene_processed$any_light),
                                                  tricolorinfo = list(NA), 
                                                  fileinfo = NA,
                                                  material_id = NA_integer_,  
                                                  csg_object = list(NA), 
                                                  mesh_info = list(NA),
                                                  flipped = flipped),
                      transforms = ray_transform(angle = list(c(0,0,0)),
                                                 order_rotation = list(order_rotation),
                                                 scale = list(c(1,1,1)),
                                                 group_transform = list(matrix(NA_real_))),
                      animation_info = ray_animated_transform(
                        start_transform_animation = list(matrix(NA_real_)), 
                        end_transform_animation = list(matrix(NA_real_)),
                        start_time = 0, end_time = 1)
  ))
}
