#' Group Objects
#'
#' Group and transform objects together. 
#'
#' @param scene Tibble of pre-existing object locations and properties to group together.
#' @param pivot_point Default `c(0,0,0)`. The point about which to pivot, scale, and move the group.
#' @param translate Default `c(0,0,0)`. Vector indicating where to offset the group.
#' @param angle Default `c(0,0,0)`. Angle of rotation around the x, y, and z axes, 
#' applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1,2,3)`. The order to apply the rotations, 
#' referring to "x", "y", and "z".
#' @param scale Default `c(1,1,1)`. Scaling factor for x, y, and z directions for all objects in group.
#' @param axis_rotation Default `NA`. Provide an axis of rotation and a single angle (via `angle`) of rotation
#' around that axis.
#'
#' @return Tibble of grouped object locations and properties.
#' @export
#' @examples
#' #Generate the ground and add some objects
#' scene = generate_cornell() %>%
#'         add_object(cube(x=555/2,y=555/8,z=555/2,width=555/4)) %>%
#'         add_object(cube(x=555/2,y=555/4+555/16,z=555/2,width=555/8))
#' \donttest{
#' render_scene(scene,lookfrom=c(278,278,-800),lookat = c(278,278,0), aperture=0,
#'              samples=500, fov=50, parallel=TRUE, clamp_value=5)
#' }
#' 
#' #Group the entire room and rotate around its center, but keep the cubes in the same place.
#' scene2 = group_objects(generate_cornell(), 
#'                        pivot_point=c(555/2,555/2,555/2),
#'                        angle=c(0,30,0)) %>%
#'          add_object(cube(x=555/2,y=555/8,z=555/2,width=555/4)) %>%
#'         add_object(cube(x=555/2,y=555/4+555/16,z=555/2,width=555/8))
#'                        
#' \donttest{
#' render_scene(scene2,lookfrom=c(278,278,-800),lookat = c(278,278,0), aperture=0,
#'              samples=500, fov=50, parallel=TRUE, clamp_value=5)
#' }.
#' 
#' #Now group the cubes instead of the Cornell box, and rotate/translate them together
#' twocubes = cube(x=555/2,y=555/8,z=555/2,width=555/4) %>%
#'            add_object(cube(x=555/2, y=555/4 + 555/16, z=555/2, width=555/8))
#' scene3 = generate_cornell() %>%
#'          add_object(group_objects(twocubes, translate = c(0,50,0),angle = c(0,45,0), 
#'          pivot_point = c(555/2,0,555/2)))
#' \donttest{
#' render_scene(scene3,lookfrom=c(278,278,-800),lookat = c(278,278,0), aperture=0,
#'              samples=500, fov=50, parallel=TRUE, clamp_value=5)
#' }
#' 
#' #Flatten and stretch the cubes together on two axes
#' scene4 = generate_cornell() %>%
#'          add_object(group_objects(twocubes, translate = c(0,-40,0), 
#'                                   angle = c(0,45,0), scale = c(2,0.5,1), 
#'                                   pivot_point = c(555/2,0,555/2)))
#' \donttest{
#' render_scene(scene4,lookfrom=c(278,278,-800),lookat = c(278,278,0), aperture=0,
#'              samples=500, fov=50, parallel=TRUE, clamp_value=5)
#' }
#' 
#' #Add another layer of grouping, including the Cornell box
#' \donttest{
#' scene4 %>% 
#'   group_objects(pivot_point = c(555/2,555/2,555/2),scale=c(1.5,0.5,0.3), angle=c(-20,0,20)) %>% 
#'   render_scene(lookfrom=c(278,278,-800),lookat = c(278,278,0), aperture=0,
#'              samples=509, fov=50, parallel=TRUE, clamp_value=5)
#' }
group_objects = function(scene, pivot_point=c(0,0,0), translate = c(0,0,0),
                         angle = c(0,0,0), order_rotation = c(1,2,3),
                         scale = c(1,1,1), axis_rotation = NA) {
  if(missing(pivot_point)) {
    pivot_point = c(0,0,0)
  }
  stopifnot(length(pivot_point) == 3)
  stopifnot(length(translate) == 3)
  stopifnot(length(order_rotation) == 3)
  
  PivotTranslateStart = generate_translation_matrix(-pivot_point)
  Scale = diag(4) * c(scale,1)
  PivotTranslateEnd = generate_translation_matrix(pivot_point)
  Translation = generate_translation_matrix(translate)
  if(any(is.na(axis_rotation))) {
    stopifnot(length(angle) == 3)
    Rotation = generate_rotation_matrix(angle,order_rotation)
  } else {
    stopifnot(length(angle) == 1)
    stopifnot(length(axis_rotation) == 3)
    stopifnot(sum(axis_rotation*axis_rotation) > 0)
    Rotation = RotateAxis(angle,axis_rotation)
  }
  for(i in seq_len(nrow(scene))) {
    if(is.na(scene$group_transform[i])) {
      scene$group_transform[i] = list(Translation %*% 
                                      PivotTranslateEnd %*% 
                                      Rotation %*% 
                                      Scale %*% 
                                      PivotTranslateStart)
    } else {
      scene$group_transform[i] = list(Translation %*% 
                                      PivotTranslateEnd %*% 
                                      Rotation %*% 
                                      Scale %*% 
                                      PivotTranslateStart %*% scene$group_transform[[i]])
    }
  }
  return(scene)
}
