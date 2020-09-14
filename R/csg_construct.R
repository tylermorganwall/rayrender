#' Constructive Solid Geometry Object
#' 
#' This object takes an object constructed using the `csg_*` functions. The object is drawn using
#' ray marching/sphere tracing.
#' 
#' Note: For dielectric objects, any other objects not included in the CSG object and
#' nested inside will be ignored.
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
csg_object = function(object, x = 0, y = 0, z = 0, radius = 1, material = diffuse(), 
                      angle = c(0, 0, 0), order_rotation = c(1, 2, 3), velocity = c(0, 0, 0), 
                      flipped = FALSE, scale = c(1,1,1)) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  new_tibble_row(list(x = x, y = y, z = z, radius = radius, type = material$type, shape = "csg_object",
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
                      material_id = NA, csg_object = list(object)))
}

#' CSG Sphere
#'
#' @param x Default `0`. x-coordinate of the center of the sphere.
#' @param y Default `0`. y-coordinate of the center of the sphere.
#' @param z Default `0`. z-coordinate of the center of the sphere.
#' @param radius Default `1`. Radius of the sphere.
#' @return List describing the sphere in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
#' generate_cornell() %>% 
#'   add_object(csg_object(csg_sphere(x=555/2,y=555/2,z=555/2,radius=150),
#'                         material=glossy(color="purple"))) %>% 
#'   render_scene()
csg_sphere = function(x=0,y=0,z=0, radius=1) {
  csg_obj = list(csg_type = 2, x=x,y=y,z=z,radius=radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Plane
#'
#' @param x Default `0`. An x-coordinate on the plane.
#' @param y Default `0`. A y-coordinate on the plane.
#' @param z Default `0`. A z-coordinate on the plane.
#' @param normal Default `c(0,1,0)`. Surface normal of the plane.
#' @param width_x Default `1`.
#' @param width_z Default `1`.
#' @return List describing the plane in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_plane = function(x=0,y=0,z=0, normal=c(0,1,0)) {
  csg_obj = list(csg_type = 3,x=x,y=y,z=z,normal=normal)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Box
#'
#' @param x Default `0`. An x-coordinate on the box.
#' @param y Default `0`. A y-coordinate on the box.
#' @param z Default `0`. A z-coordinate on the box
#' @param normal Default `c(1,1,1)`. Length-3 vector describing the x/y/z widths of the box
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_box = function(x=0,y=0,z=0, width=c(1,1,1), corner_radius = 0) {
  if(corner_radius == 0) {
    csg_obj = list(csg_type = 4,x=x,y=y,z=z,width=width/2)
    class(csg_obj) = c("list", "ray_csg")
    csg_obj
  } else {
    csg_obj = list(csg_type = 5,x=x,y=y,z=z,width=width/2, radius=corner_radius)
    class(csg_obj) = c("list", "ray_csg")
    csg_obj
  }
}

#' CSG Torus
#'
#' @param x Default `0`. x-coordinate on the torus.
#' @param y Default `0`. y-coordinate on the torus.
#' @param z Default `0`. z-coordinate on the torus.
#' @param radius Default `1`. Torus radius.
#' @param cross_radius Default `0.5`. Cross section radius of the torus.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_torus = function(x=0,y=0,z=0, radius=1, minor_radius=0.5) {
  csg_obj = list(csg_type = 7,x=x,y=y,z=z,ring_radius=radius, cross_radius=minor_radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Capsule
#'
#' @param x Default `0`. x-coordinate on the torus.
#' @param y Default `0`. y-coordinate on the torus.
#' @param z Default `0`. z-coordinate on the torus.
#' @param radius Default `1`. Torus radius.
#' @param cross_radius Default `0.5`. Cross section radius of the torus.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_capsule = function(start = c(0,0,0), end = c(0,1,0), radius=1) {
  csg_obj = list(csg_type = 8,start=start, end=end,radius=radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Cylinder
#'
#' @param x Default `0`. x-coordinate on the torus.
#' @param y Default `0`. y-coordinate on the torus.
#' @param z Default `0`. z-coordinate on the torus.
#' @param radius Default `1`. Torus radius.
#' @param cross_radius Default `0.5`. Cross section radius of the torus.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_cylinder = function(start = c(0,0,0), end = c(0,1,0), radius=1, corner_radius=0) {
  csg_obj = list(csg_type = 9,start=start, end=end,radius=radius,corner_radius=corner_radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Ellipsoid
#'
#' @param x Default `0`. x-coordinate on the torus.
#' @param y Default `0`. y-coordinate on the torus.
#' @param z Default `0`. z-coordinate on the torus.
#' @param axes Default `c(0.5,1,0.5)`. Ellipsoid principle axes.
#' @param cross_radius Default `0.5`. Cross section radius of the torus.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_ellipsoid = function(x=0,y=0,z=0,axes=c(0.5,1,0.5)) {
  csg_obj = list(csg_type = 10,x=x, y=y,z=z,end=end,axes=axes)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Rounded Cone
#'
#' @param x Default `0`. x-coordinate on the torus.
#' @param y Default `0`. y-coordinate on the torus.
#' @param z Default `0`. z-coordinate on the torus.
#' @param axes Default `c(0.5,1,0.5)`. Ellipsoid principle axes.
#' @param cross_radius Default `0.5`. Cross section radius of the torus.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_rounded_cone = function(start = c(0,0,0), end = c(0,1,0), radius=0.5, upper_radius = 0.2) {
  csg_obj = list(csg_type = 11, start = start, end=end,
                 radius=radius,upper_radius=upper_radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Cone
#'
#' @param x Default `0`. x-coordinate on the torus.
#' @param y Default `0`. y-coordinate on the torus.
#' @param z Default `0`. z-coordinate on the torus.
#' @param axes Default `c(0.5,1,0.5)`. Ellipsoid principle axes.
#' @param cross_radius Default `0.5`. Cross section radius of the torus.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_cone = function(start = c(0,0,0), end = c(0,1,0), radius=0.5, upper_radius = 0.2) {
  csg_obj = list(csg_type = 12, start = start, end=end,
                 radius=radius,upper_radius=upper_radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Pyramid
#'
#' @param x Default `0`. x-coordinate on the torus.
#' @param y Default `0`. y-coordinate on the torus.
#' @param z Default `0`. z-coordinate on the torus.
#' @param height Default `1`. Ellipsoid principle axes.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_pyramid = function(x=0,y=0,z=0,height=1,base=1) {
  csg_obj = list(csg_type = 13, x=x,y=y,z=z,
                 height=height,base=base)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Triangle
#'
#' @param v1 Default `c(0,1,0)`. First vertex.
#' @param v2 Default `c(1,0,0)`. Second vertex.
#' @param v3 Default `c(-1,0,0)`. Third vertex.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_triangle = function(object, v1=c(0,1,0),v2=c(1,0,0),v3=c(-1,0,0)) {
  csg_obj = list(csg_type = 14, v1=v1,v2=v2,v3=v3)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Elongate
#'
#' @param object CSG object.
#' @param x Default `0`. x-elongation.
#' @param y Default `0`. y-elongation.
#' @param z Default `0`. z-elongation.
#' @param robust Default `TRUE`. `FALSE` switches to a faster (but less robust in 2D) method.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_elongate = function(object, x=0,y=0,z=0, robust = TRUE) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  if(!robust) {
    csg_obj = list(csg_type = 15, object = object, x=x,y=y,z=z)
  } else {
    csg_obj = list(csg_type = 16, object = object, x=x,y=y,z=z)
  }
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Round
#'
#' @param object CSG object.
#' @param radius Default `0.1`. Rounding distance.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_round = function(object, radius=0.1) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  csg_obj = list(csg_type = 17, object = object, radius=radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Onion
#'
#' @param object CSG object.
#' @param thickness Default `0.1`. Onioning distance.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_onion = function(object, thickness=0.1) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  csg_obj = list(csg_type = 18, object = object, thickness=thickness)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Scale
#'
#' @param object CSG object.
#' @param scale Default `1`.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_scale = function(object, scale=1) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  if(scale == 1) {
    return(object)
  }
  csg_obj = list(csg_type = 19, object = object, scale=scale)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Rotate
#'
#' @param object CSG object.
#' @param up Default `c(0,1,0)`. New up vector.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_rotate = function(object, up = c(0,1,0)) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  csg_obj = list(csg_type = 20, object = object, angles=up)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Translate
#'
#' @param object CSG object.
#' @param x Default `0`. x translation.
#' @param y Default `0`. y translation.
#' @param z Default `0`. z translation.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_translate = function(object, x=0,y=0,z=0) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  csg_obj = list(csg_type = 21, object = object, translate = c(x,y,z))
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Construction
#'
#' @param object1 First CSG object
#' @param object2 Second CSG object
#' @param operation Default `union`. Can be `union`, `subtract`, `intersection`, `blend`, `subtractblend`, or `mix`.
#' @param blend Default `0.5`. Blending radius.
#'
#' @return List describing the combined csg object in the scene.
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_combine = function(object1, object2, operation = "union", blend = 0.5) {
  if(!inherits(object1,"ray_csg")) {
    stop("`object1` must be constructed with rayrender csg_* functions")
  }
  if(!inherits(object2,"ray_csg")) {
    stop("`object2` must be constructed with rayrender csg_* functions")
  }
  combine_type = unlist(lapply(tolower(operation),switch,
                          "union" = 1,"subtract" = 2,"intersection" = 3, 
                          "blend" = 4, "mix" = 5,"subtractblend" = 6, 1))
  csg_obj = list(csg_type = 1, object1 = object1, object2 = object2,
                 operation = combine_type, blend = blend)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Group
#'
#' @param object_list List of objects created with the csg_* functions. This will make all further operations
#' be applied to this object as a group.
#'
#' @return List describing the group in the scene. 
#' @export
#'
#' @examples
#' #Generate a sphere in the cornell box.
csg_group = function(object_list) {
  if(inherits(object_list, "ray_csg")) {
    return(object_list)
  }
  if(!inherits(object_list, "list")) {
    stop("`object_list` must be list() of CSG objects")
  }
  if(length(object_list) == 1) {
    return(object_list[[1]])
  }
  all_ray = all(unlist(lapply(object_list, inherits, "ray_csg")))
  if(!all_ray) {
    stop("All objects in `object_list` must come from rayrender CSG functions.")
  }
  csg_obj = list(csg_type = 6, shapes = object_list)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}