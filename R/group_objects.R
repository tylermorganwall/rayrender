#' Group Object
#'
#' @param scene Tibble of pre-existing object locations and properties to group together.
#' @param pivot_point Defaults to the mean location of all the objects. 
#' The point about which to pivot and move the group.
#' @param group_angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, 
#' applied in the order specified in `order_rotation`.
#' @param group_order_rotation Default `c(1,2,3)`. The order to apply the rotations, 
#' referring to "x", "y", and "z".
#'
#' @return Tibble of grouped object locations and properties.
#' @export
#'
#' @examples
#' #Generate the ground and add some objects
group_objects = function(scene, pivot_point=c(0,0,0), group_angle = c(0,0,0), group_order_rotation = c(1,2,3)) {
  if(missing(pivot_point)) {
    pivot_point = c(mean(scene$x), mean(scene$y), mean(scene$z))
  }
  assertthat::assert_that(length(pivot_point) == 3)
  scene$pivot_point = list(pivot_point)
  scene$group_angle = list(group_angle)
  scene$group_order_rotation = list(group_order_rotation)
  return(scene)
}