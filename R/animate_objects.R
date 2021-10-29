#' Animate Objects
#'
#' Animate objects together. 
#'
#' @param scene Tibble of pre-existing object locations and properties to group together.
#' @param start_time Default `0`.
#' @param end_time Default `1`.
#' @param start_pivot_point Default `c(0,0,0)`. The point about which to pivot, scale, and move the group.
#' @param start_position Default `c(0,0,0)`. Vector indicating where to offset the group.
#' @param start_angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, 
#' applied in the order specified in `order_rotation`.
#' @param start_order_rotation Default `c(1,2,3)`. The order to apply the rotations, 
#' referring to "x", "y", and "z".
#' @param start_scale Default `c(1,1,1)`. Scaling factor for x, y, and z directions for all objects in group.
#' @param start_axis_rotation Default `NA`. Provide an axis of rotation and a single angle (via `angle`) of rotation
#' @param end_pivot_point Default `c(0,0,0)`. The point about which to pivot, scale, and move the group.
#' @param end_position Default `c(0,0,0)`. Vector indicating where to offset the group.
#' @param end_angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, 
#' applied in the order specified in `order_rotation`.
#' @param end_order_rotation Default `c(1,2,3)`. The order to apply the rotations, 
#' referring to "x", "y", and "z".
#' @param end_scale Default `c(1,1,1)`. Scaling factor for x, y, and z directions for all objects in group.
#' @param end_axis_rotation Default `NA`. Provide an axis of rotation and a single angle (via `angle`) of rotation
#' around that axis.
#'
#' @return Tibble of grouped object locations and properties.
#' @export
#' @examples
#' #Generate the ground and add some objects
#' 
animate_objects = function(scene, start_time = 0, end_time = 1,
                           start_pivot_point = c(0,0,0), start_position = c(0,0,0),
                           start_angle = c(0,0,0), start_order_rotation = c(1,2,3),
                           start_scale = c(1,1,1), start_axis_rotation = NA,
                           end_pivot_point = c(0,0,0), end_position = c(0,0,0),
                           end_angle = c(0,0,0), end_order_rotation = c(1,2,3),
                           end_scale = c(1,1,1), end_axis_rotation = NA) {
  stopifnot(length(start_pivot_point) == 3)
  stopifnot(length(start_position) == 3)
  stopifnot(length(start_order_rotation) == 3)
  
  PivotTranslateStart = generate_translation_matrix(-start_pivot_point)
  Scale = diag(4) * c(start_scale,1)
  PivotTranslateEnd = generate_translation_matrix(start_pivot_point)
  Translation = generate_translation_matrix(start_position)
  if(any(is.na(start_axis_rotation))) {
    stopifnot(length(start_angle) == 3)
    Rotation = generate_rotation_matrix(start_angle,start_order_rotation)
  } else {
    stopifnot(length(start_angle) == 1)
    stopifnot(length(start_axis_rotation) == 3)
    stopifnot(sum(start_axis_rotation*start_axis_rotation) > 0)
    Rotation = RotateAxis(start_angle,start_axis_rotation)
  }
  for(i in seq_len(nrow(scene))) {
    scene$start_transform_animation[i] = list(Translation %*% 
                                      PivotTranslateEnd %*% 
                                      Rotation %*% 
                                      Scale %*% 
                                      PivotTranslateStart)
    scene$start_time[i] = start_time
  }
  stopifnot(length(end_pivot_point) == 3)
  stopifnot(length(end_position) == 3)
  stopifnot(length(end_order_rotation) == 3)
  
  PivotTranslateStart = generate_translation_matrix(-end_pivot_point)
  Scale = diag(4) * c(end_scale,1)
  PivotTranslateEnd = generate_translation_matrix(end_pivot_point)
  Translation = generate_translation_matrix(end_position)
  if(any(is.na(end_axis_rotation))) {
    stopifnot(length(end_angle) == 3)
    Rotation = generate_rotation_matrix(end_angle,end_order_rotation)
  } else {
    stopifnot(length(end_angle) == 1)
    stopifnot(length(end_axis_rotation) == 3)
    stopifnot(sum(end_axis_rotation*end_axis_rotation) > 0)
    Rotation = RotateAxis(end_angle,end_axis_rotation)
  }
  for(i in seq_len(nrow(scene))) {
    scene$end_transform_animation[i] = list(Translation %*% 
                                                PivotTranslateEnd %*% 
                                                Rotation %*% 
                                                Scale %*% 
                                                PivotTranslateStart)
    scene$end_time[i] = end_time
  }
  return(scene)
}
