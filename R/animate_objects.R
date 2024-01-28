#' Animate Objects
#'
#' This function animates an object between two states. This animates objects separately from the transformations set
#' in `group_objects()` and in the object transformations themselves. This creates motion blur, controlled by the shutter
#' open/close options in `render_scene()`.
#' 
#'
#' @param scene Tibble of pre-existing object locations.
#' @param start_time Default `0`. Start time of movement.
#' @param end_time Default `1`. End time of movement.
#' @param start_pivot_point Default `c(0,0,0)`. The point about which to pivot, scale, and move the objects.
#' @param start_position Default `c(0,0,0)`. Vector indicating where to offset the objects.
#' @param start_angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, 
#' applied in the order specified in `order_rotation`.
#' @param start_order_rotation Default `c(1,2,3)`. The order to apply the rotations, 
#' referring to "x", "y", and "z".
#' @param start_scale Default `c(1,1,1)`. Scaling factor for x, y, and z directions for all objects.
#' @param start_axis_rotation Default `NA`. Provide an axis of rotation and a single angle (via `angle`) of rotation
#' @param end_pivot_point Default `c(0,0,0)`. The point about which to pivot, scale, and move the group.
#' @param end_position Default `c(0,0,0)`. Vector indicating where to offset the objects.
#' @param end_angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, 
#' applied in the order specified in `order_rotation`.
#' @param end_order_rotation Default `c(1,2,3)`. The order to apply the rotations, 
#' referring to "x", "y", and "z".
#' @param end_scale Default `c(1,1,1)`. Scaling factor for x, y, and z directions for all objects.
#' @param end_axis_rotation Default `NA`. Provide an axis of rotation and a single angle (via `angle`) of rotation
#' around that axis.
#'
#' @return Tibble of animated object.
#' @export
#' @examples
#' #Render a pig
#' if(run_documentation()) {
#' generate_studio() %>% 
#'   add_object(pig(y=-1.2,scale=0.5,angle=c(0,-70,0)))%>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=2,material=light())) %>% 
#'   render_scene(samples=128,sample_method = "sobol_blue")
#' }
#' if(run_documentation()) {
#' #Render a moving pig
#' generate_studio() %>% 
#'   add_object(
#'     animate_objects(
#'       pig(y=-1.2,scale=0.5,angle=c(0,-70,0)),
#'       start_position = c(-0.1,0,0), end_position = c(0.1,0.2,0))
#'   ) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=2,material=light())) %>% 
#'   render_scene(samples=128,sample_method = "sobol_blue",clamp_value = 10)
#' }
#' if(run_documentation()) {
#' 
#' #Render a shrinking pig
#' generate_studio() %>% 
#'   add_object(
#'     animate_objects(
#'       pig(y=-1.2,scale=0.5,angle=c(0,-70,0)),
#'       start_scale = c(1,1,1), end_scale = c(0.5,0.5,0.5))
#'   ) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=2,material=light())) %>% 
#'   render_scene(samples=128,sample_method = "sobol_blue",clamp_value = 10)
#' }
#' if(run_documentation()) {
#' #Render a spinning pig
#' generate_studio() %>% 
#'   add_object(
#'     animate_objects(
#'       pig(y=-1.2,scale=0.5,angle=c(0,-70,0)),
#'       start_angle = c(0,-30,0), end_angle = c(0,30,0))
#'   ) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=2,material=light())) %>% 
#'   render_scene(samples=128,sample_method = "sobol_blue",clamp_value = 10)
#' }
#' if(run_documentation()) {
#' 
#' #Shorten the open shutter time frame
#' generate_studio() %>% 
#'   add_object(
#'     animate_objects(
#'       pig(y=-1.2,scale=0.5,angle=c(0,-70,0)),
#'       start_angle = c(0,-30,0), end_angle = c(0,30,0))
#'   ) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=2,material=light())) %>% 
#'   render_scene(samples=128,sample_method = "sobol_blue",clamp_value = 10, 
#'                shutteropen=0.4, shutterclose = 0.6)
#' }
#' if(run_documentation()) {
#' #Change the time frame when the shutter is open
#' generate_studio() %>% 
#'   add_object(
#'     animate_objects(
#'       pig(y=-1.2,scale=0.5,angle=c(0,-70,0)),
#'       start_angle = c(0,-30,0), end_angle = c(0,30,0))
#'   ) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=2,material=light())) %>% 
#'   render_scene(samples=128,sample_method = "sobol_blue",clamp_value = 10, 
#'                shutteropen=0, shutterclose = 0.1)
#' }
#' if(run_documentation()) {    
#' #Shorten the time span in which the movement occurs (which, in effect, 
#' #increases the speed of the transition).
#' generate_studio() %>% 
#'   add_object(
#'     animate_objects(start_time = 0, end_time=0.1,
#'       pig(y=-1.2,scale=0.5,angle=c(0,-70,0)),
#'       start_angle = c(0,-30,0), end_angle = c(0,30,0))
#'   ) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=2,material=light())) %>% 
#'   render_scene(samples=128,sample_method = "sobol_blue",clamp_value = 10, 
#'                shutteropen=0, shutterclose = 0.1)
#' }
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
    scene$animation_info[[i]]$start_transform_animation[[1]] = 
      (Translation %*% 
       PivotTranslateEnd %*% 
       Rotation %*% 
       Scale %*% 
       PivotTranslateStart)
    scene$animation_info[[i]]$start_time = start_time
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
    scene$animation_info[[i]]$end_transform_animation[[1]] = 
      (Translation %*% 
       PivotTranslateEnd %*% 
       Rotation %*% 
       Scale %*% 
       PivotTranslateStart)
    scene$animation_info[[i]]$end_time = end_time
  }
  return(scene)
}
