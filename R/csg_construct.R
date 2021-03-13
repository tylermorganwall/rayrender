#' Constructive Solid Geometry Object
#' 
#' This object takes an object constructed using the `csg_*` functions. The object is drawn using
#' ray marching/sphere tracing.
#' 
#' Note: For dielectric objects, any other objects not included in the CSG object and
#' nested inside will be ignored.
#'
#' @param object Object created with CSG interface.
#' @param x Default `0`. x-offset of the center of the object.
#' @param y Default `0`. y-offset of the center of the object.
#' @param z Default `0`. z-offset of the center of the object.
#' @param radius Default `1`. Radius of the sphere..
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
#' \donttest{
#' #We will combine these three objects:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>%
#'   add_object(csg_object(csg_box(), material=glossy(color="red"))) %>% 
#'   add_object(csg_object(csg_sphere(radius=0.707), material=glossy(color="green"))) %>% 
#'   add_object(csg_object(csg_group(list(csg_cylinder(start=c(-1,0,0), end=c(1,0,0), radius=0.4),
#'                    csg_cylinder(start=c(0,-1,0), end=c(0,1,0), radius=0.4),
#'                    csg_cylinder(start=c(0,0,-1), end=c(0,0,1), radius=0.4))),
#'                    material=glossy(color="blue"))) %>% 
#'   add_object(sphere(y=5,x=3,radius=1,material=light(intensity=30))) %>%
#'   render_scene(clamp_value=10, fov=15,lookfrom=c(5,5,10))
#' 
#' #Standard CSG sphere + box - crossed cylinder combination:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>%
#'   add_object(csg_object(csg_combine(
#'     csg_combine(
#'       csg_box(),
#'       csg_sphere(radius=0.707),
#'       operation="intersection"),
#'     csg_group(list(csg_cylinder(start=c(-1,0,0), end=c(1,0,0), radius=0.4),
#'                    csg_cylinder(start=c(0,-1,0), end=c(0,1,0), radius=0.4),
#'                    csg_cylinder(start=c(0,0,-1), end=c(0,0,1), radius=0.4))),
#'     operation="subtract"),
#'     material=glossy(color="red"))) %>%
#'   add_object(sphere(y=5,x=3,radius=1,material=light(intensity=30))) %>%
#'   render_scene(clamp_value=10, fov=10,lookfrom=c(5,5,10))
#'   
#' #Blend them all instead:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>%
#'   add_object(csg_object(csg_combine(
#'     csg_combine(
#'       csg_box(),
#'       csg_sphere(radius=0.707),
#'       operation="blend"),
#'     csg_group(list(csg_cylinder(start=c(-1,0,0), end=c(1,0,0), radius=0.4),
#'                    csg_cylinder(start=c(0,-1,0), end=c(0,1,0), radius=0.4),
#'                    csg_cylinder(start=c(0,0,-1), end=c(0,0,1), radius=0.4))),
#'     operation="blend"),
#'     material=glossy(color="purple"))) %>%
#'   add_object(sphere(y=5,x=3,radius=1,material=light(intensity=30))) %>%
#'   render_scene(clamp_value=10, fov=15,lookfrom=c(5,5,10))
#' }
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
                      material_id = NA, csg_object = list(object), mesh_info = list(NA)))
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
#' \donttest{
#' #Generate a simple sphere:
#' generate_ground() %>% 
#'   add_object(csg_object(csg_sphere(),
#'                         material=glossy(color="purple"))) %>% 
#'   render_scene(clamp_value=10)
#' 
#' #Generate a bigger sphere in the cornell box.
#' generate_cornell() %>% 
#'   add_object(csg_object(csg_sphere(x=555/2,y=555/2,z=555/2,radius=100),
#'                         material=glossy(checkercolor="purple", checkerperiod=100))) %>% 
#'   render_scene(clamp_value=10)
#'   
#' #Combine two spheres of different sizes
#' generate_cornell() %>% 
#'   add_object(csg_object(
#'     csg_combine(
#'       csg_sphere(x=555/2,y=555/2-50,z=555/2,radius=100),
#'       csg_sphere(x=555/2,y=555/2+50,z=555/2,radius=80)),
#'     material=glossy(color="purple"))) %>% 
#'   render_scene(clamp_value=10)
#'
#'#Subtract two spheres to create an indented region
#' generate_cornell() %>% 
#'   add_object(csg_object(
#'     csg_combine(
#'       csg_sphere(x=555/2,y=555/2-50,z=555/2,radius=100),
#'       csg_sphere(x=555/2+30,y=555/2+20,z=555/2-90,radius=40),
#'       operation="subtract"),
#'     material=glossy(color="grey20"))) %>% 
#'   render_scene(clamp_value=10)
#'   
#'#Use csg_combine(operation="blend") to melt the two together
#' generate_cornell() %>% 
#'   add_object(csg_object(
#'     csg_combine(
#'       csg_sphere(x=555/2,y=555/2-50,z=555/2,radius=100),
#'       csg_sphere(x=555/2,y=555/2+50,z=555/2,radius=80),
#'       operation="blend", radius=20),
#'     material=glossy(color="purple"))) %>% 
#'   render_scene(clamp_value=10)
#' }
csg_sphere = function(x=0,y=0,z=0, radius=1) {
  csg_obj = list(csg_type = 2, x=x,y=y,z=z,radius=radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Plane
#' 
#' Note: This shape isn't closed, so there may be odd lighting issues if it's oriented the wrong
#' way.
#'
#' @param x Default `0`. An x-coordinate on the plane.
#' @param y Default `0`. A y-coordinate on the plane.
#' @param z Default `0`. A z-coordinate on the plane.
#' @param normal Default `c(0,1,0)`. Surface normal of the plane.
#' @param width_x Default `10`.
#' @param width_z Default `10`.
#' @return List describing the plane in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a plane
#' csg_object(csg_plane(width_x=4, width_z=4), material=diffuse(checkercolor="purple")) %>% 
#'   add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'   render_scene(clamp_value=10)
#'  
#' #Combine the plane with a sphere
#' csg_object(csg_combine(
#'     csg_sphere(radius=0.5),
#'     csg_plane(width_x=4, width_z=4,y=-0.5), 
#'     operation="blend"),material=diffuse(checkercolor="purple")) %>% 
#'   add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'   render_scene(clamp_value=10)
#'   
#' #Re-orient the plane using the normal and 
#' csg_object(csg_combine(
#'     csg_sphere(radius=0.5),
#'     csg_plane(normal = c(1,1,0),width_x=4, width_z=4,y=-0.5), 
#'     operation="blend"),material=diffuse(checkercolor="purple")) %>% 
#'   add_object(sphere(y=5,x=5,material=light(intensity=40))) %>% 
#'   render_scene(clamp_value=10)
#' }
csg_plane = function(x=0,y=0,z=0, normal=c(0,1,0),width_x=4, width_z=4) {
  csg_obj = list(csg_type = 3,x=x,y=y,z=z,normal=normal,width_x=width_x, width_z=width_z)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Box
#'
#' @param x Default `0`. An x-coordinate on the box.
#' @param y Default `0`. A y-coordinate on the box.
#' @param z Default `0`. A z-coordinate on the box
#' @param width Default `c(1,1,1)`. Length-3 vector describing the x/y/z widths of the box
#' @param corner_radius Default `0`. Radius if rounded box.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a box
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_box(), material=glossy(color="#FF69B4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=5))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(7,3,7))
#'   
#' #Change the width
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_box(width = c(2,1,0.5)), material=glossy(color="#FF69B4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=5))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(7,3,7))
#'
#' #Subtract two boxes to make stairs
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'     csg_box(),
#'     csg_box(x=0.5,y=0.5,width=c(1,1,1.1)),operation="subtract"),
#'    material=glossy(color="#FF69B4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=5))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(7,3,7),fov=13)
#'   }
csg_box = function(x=0,y=0,z=0, width=c(1,1,1), corner_radius = 0) {
  if(corner_radius == 0) {
    csg_obj = list(csg_type = 4,x=x,y=y,z=z,width=width)
    class(csg_obj) = c("list", "ray_csg")
    csg_obj
  } else {
    csg_obj = list(csg_type = 5,x=x,y=y,z=z,width=width, radius=corner_radius)
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
#' @param minor_radius Default `0.5`. Cross section radius of the torus.
#' @return List describing the torus in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a torus:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_torus(), material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(0,5,10),fov=30)
#'   
#' #Change the radius of the torus:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_torus(radius=2), material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(0,5,10),fov=30)
#'
#'#Change the minor radius of the torus:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_torus(radius=2, minor_radius=0.25), 
#'                         material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(0,5,10),fov=30)
#'   
#'#Generate a rotated torus in the Cornell Box
#' generate_cornell() %>% 
#'   add_object(csg_object(csg_rotate(
#'     csg_torus(x=555/2,y=555/2,z=555/2,radius=100, minor_radius=50), 
#'     pivot_point = c(555/2,555/2,555/2), up =c(0,1,-1)), 
#'                         material=glossy(color="dodgerblue4"))) %>%
#'   render_scene(clamp_value=10)
#' }
csg_torus = function(x=0,y=0,z=0, radius=1, minor_radius=0.5) {
  csg_obj = list(csg_type = 7,x=x,y=y,z=z,ring_radius=radius, cross_radius=minor_radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Capsule
#'
#' @param start Default `c(0, 0, 0)`. Start point of the capsule, specifying `x`, `y`, `z`.
#' @param end Default `c(0, 1, 0)`. End point of the capsule, specifying `x`, `y`, `z`.
#' @param radius Default `1`. Capsule radius.
#' @return List describing the capsule in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a basic capsule:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_capsule(radius=0.5),material=glossy(color="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Change the orientation by specifying a start and end
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(csg_capsule(start = c(-1,0.5,-2), end = c(1,0.5,-2),
#'   radius=0.5),material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(3,3,10))
#'  
#' #Show the effect of changing the radius
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(
#'     csg_combine(
#'     csg_capsule(start = c(-1,0.5,-2), end = c(1,0.5,-2), radius=0.5),
#'     csg_capsule(start = c(-0.5,1.5,-2), end = c(0.5,1.5,-2), radius=0.25)),
#'     material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(-3,3,10))
#'                
#' #Render a capsule in a Cornell box
#' generate_cornell() %>% 
#'   add_object(csg_object(
#'     csg_capsule(start = c(555/2-100,555/2,555/2), end = c(555/2+100,555/2,555/2), radius=100),
#'     material=glossy(color="dodgerblue4"))) %>% 
#'   render_scene(clamp_value=10)
#'}
csg_capsule = function(start = c(0,0,0), end = c(0,1,0), radius=1) {
  csg_obj = list(csg_type = 8,start=start, end=end,radius=radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Cylinder
#'
#' @param start Default `c(0, 0, 0)`. Start point of the cylinder, specifing `x`, `y`, `z`.
#' @param end Default `c(0, 1, 0)`. End point of the cylinder, specifing `x`, `y`, `z`.
#' @param radius Default `1`. Cylinder radius.
#' @param corner_radius Default `0`. Radius if rounded cylinder.
#' @return List describing the cylinder in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a basic cylinder:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_cylinder(radius=0.25),material=glossy(color="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Change the orientation by specifying a start and end
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(csg_cylinder(start = c(-1,0.5,-2), end = c(1,0.5,-2),
#'     radius=0.5),material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(3,3,10))
#'  
#' #Show the effect of changing the radius
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(
#'     csg_combine(
#'     csg_cylinder(start = c(-1,0.5,-2), end = c(1,0.5,-2), radius=0.5),
#'     csg_cylinder(start = c(-0.5,1.5,-2), end = c(0.5,1.5,-2), radius=0.25)),
#'     material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(-3,3,10))
#'                
#' #Render a red marble cylinder in a Cornell box
#' generate_cornell(light=FALSE) %>% 
#'   add_object(csg_object(
#'     csg_cylinder(start = c(555/2,0,555/2), end = c(555/2,350,555/2), radius=100),
#'     material=glossy(color="darkred",noisecolor="white",noise=0.03))) %>% 
#'     add_object(sphere(y=555,x=5,z=5, radius=5,
#'                material=light(intensity=10000,
#'                               spotlight_focus = c(555/2,555/2,555/2),spotlight_width = 45))) %>% 
#'   render_scene(clamp_value=4)
#'}
csg_cylinder = function(start = c(0,0,0), end = c(0,1,0), radius=1, corner_radius=0) {
  csg_obj = list(csg_type = 9,start=start, end=end,radius=radius,corner_radius=corner_radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Ellipsoid
#'
#' @param x Default `0`. x-coordinate on the ellipsoid.
#' @param y Default `0`. y-coordinate on the ellipsoid.
#' @param z Default `0`. z-coordinate on the ellipsoid.
#' @param axes Default `c(0.5,1,0.5)`. Ellipsoid principle axes.
#' @return List describing the ellipsoid in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a basic ellipsoid:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_ellipsoid(),material=glossy(color="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Three different ellipsoids:
#'generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'  add_object(csg_object(csg_group(list(
#'    csg_ellipsoid(x=-1.2, axes = c(0.2,0.5,0.5)),
#'    csg_ellipsoid(x=0, axes = c(0.5,0.2,0.5)),
#'    csg_ellipsoid(x=1.2, axes = c(0.5,0.5,0.2)))),
#'    material=glossy(color="red"))) %>% 
#'  render_scene(clamp_value=10,fov=20,lookfrom=c(0,5,10))
#'  
#' #Generate a glass ellipsoid:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_ellipsoid(),material=dielectric(attenuation = c(1,1,0.3)))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Generate a glass ellipsoid in a Cornell box:
#' generate_cornell() %>% 
#'   add_object(csg_object(csg_ellipsoid(x=555/2,y=555/2,z=555/2,axes=c(100,150,200)),
#'     material=dielectric(attenuation = c(1,0.3,1)/200))) %>% 
#'   render_scene(clamp_value=10)
#'}
csg_ellipsoid = function(x=0,y=0,z=0,axes=c(0.5,1,0.5)) {
  csg_obj = list(csg_type = 10,x=x, y=y,z=z,axes=axes)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Rounded Cone
#'
#' @param start Default `c(0, 0, 0)`. Start point of the cone, specifing `x`, `y`, `z`.
#' @param end Default `c(0, 1, 0)`. End point of the cone, specifing `x`, `y`, `z`.
#' @param radius Default `0.5`. Radius of the bottom of the cone.
#' @param upper_radius Default `0.2`. Radius from the top of the cone.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a basic rounded cone:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_rounded_cone(),material=glossy(color="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Change the orientation by specifying a start and end
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(csg_rounded_cone(start = c(-1,0.5,-2), end = c(1,0.5,-2),
#'   radius=0.5),material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(3,3,10))
#'  
#' #Show the effect of changing the radius
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(
#'     csg_combine(
#'     csg_rounded_cone(start = c(-1,0.5,-2), end = c(1,0.5,-2), radius=0.5),
#'     csg_rounded_cone(start = c(-0.5,1.5,-2), end = c(0.5,1.5,-2), radius=0.2,upper_radius = 0.5)),
#'     material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(-3,3,10))
#'                
#' #Render a glass rounded cone in a Cornell box
#' generate_cornell() %>% 
#'   add_object(csg_object(
#'     csg_rounded_cone(start = c(555/2,555/2-100,555/2), end = c(555/2,555/2+100,555/2), radius=100),
#'     material=dielectric(attenuation=c(1,1,0.3)/100))) %>% 
#'   render_scene(clamp_value=10)
#'}
csg_rounded_cone = function(start = c(0,0,0), end = c(0,1,0), radius=0.5, upper_radius = 0.2) {
  csg_obj = list(csg_type = 11, start = start, end=end,
                 radius=radius,upper_radius=upper_radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Cone
#'
#' @param start Default `c(0, 0, 0)`. Start point of the cone, specifing `x`, `y`, `z`.
#' @param end Default `c(0, 1, 0)`. End point of the cone, specifing `x`, `y`, `z`.
#' @param radius Default `1`. Radius of the bottom of the cone.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a basic cone:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_cone(),material=glossy(color="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Change the orientation by specifying a start and end
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(csg_cone(start = c(-1,0.5,-2), end = c(1,0.5,-2),
#'   radius=0.5),material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(3,3,10))
#'  
#' #Show the effect of changing the radius
#' generate_ground(material=diffuse(color="dodgerblue4",checkercolor="grey10")) %>% 
#'   add_object(csg_object(
#'     csg_combine(
#'     csg_cone(start = c(-1,0.5,-2), end = c(1,0.5,-2), radius=0.5),
#'     csg_cone(start = c(-0.5,1.5,-2), end = c(0.5,1.5,-2), radius=0.2)),
#'     material=glossy(checkercolor="red"))) %>% 
#'   render_scene(clamp_value=10,fov=20,
#'                lookat=c(0,0.5,-2),lookfrom=c(-3,3,10))
#'                
#' #Render a glass cone in a Cornell box
#' generate_cornell() %>% 
#'   add_object(csg_object(
#'     csg_cone(start = c(555/2,0,555/2), end = c(555/2,555/2+100,555/2), radius=100),
#'     material=dielectric(attenuation=c(1,1,0.3)/100))) %>% 
#'   render_scene(clamp_value=10)
#'}
csg_cone = function(start = c(0,0,0), end = c(0,1,0), radius=0.5) {
  csg_obj = list(csg_type = 12, start = start, end=end,
                 radius=radius)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Pyramid
#' 
#' Note: This primitive slows down immensely for large values of base and height. Try using csg_scale()
#' with this object for large pyramids instead.
#'
#' @param x Default `0`. x-coordinate on the pyramid.
#' @param y Default `0`. y-coordinate on the pyramid.
#' @param z Default `0`. z-coordinate on the pyramid.
#' @param height Default `1`. Pyramid height.
#' @param base Default `1`. Pyramid base width.
#' @return List describing the box in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Generate a simple pyramid:
#' generate_ground() %>% 
#'   add_object(csg_object(csg_pyramid(y=-0.99),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,material=light(intensity=20))) %>% 
#'   render_scene(clamp_value=10,lookfrom=c(-3,1,10), 
#'                fov=15, lookat=c(0,-0.5,0))
#' 
#' #Make a taller pyramid
#' generate_ground() %>% 
#'   add_object(csg_object(csg_pyramid(y=-0.95, height=1.5),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,material=light(intensity=20))) %>% 
#'   render_scene(clamp_value=10,lookfrom=c(-3,1,10), 
#'                fov=15, lookat=c(0,-0.5,0))
#'   
#' #Make a wider pyramid
#' generate_ground() %>% 
#'   add_object(csg_object(csg_pyramid(y=-0.95, base=1.5),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,material=light(intensity=20))) %>% 
#'   render_scene(clamp_value=10,lookfrom=c(-3,1,10), 
#'                fov=15, lookat=c(0,-0.5,0))
#'}
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
#' \donttest{
#' #Generate a basic triangle:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_triangle(),material=diffuse(color="red"))) %>% 
#'   add_object(sphere(y=5,z=3,material=light(intensity=30))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Change a vertex:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_triangle(v1 = c(1,1,0)),material=diffuse(color="green"))) %>% 
#'   add_object(sphere(y=5,z=3,material=light(intensity=30))) %>% 
#'   render_scene(clamp_value=10,fov=20)
#'   
#' #Change all three vertices:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_triangle(v1 = c(0.5,1,0), v2 = c(1,-0.5,0), v3 = c(-1,0.5,0)),
#'                         material=diffuse(color="blue"))) %>% 
#'   add_object(sphere(y=5,z=3,material=light(intensity=30))) %>% 
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(0,5,10))
#'}
csg_triangle = function(v1=c(0,1,0),v2=c(1,0,0),v3=c(-1,0,0)) {
  csg_obj = list(csg_type = 14, v1=v1,v2=v2,v3=v3)
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Elongate
#' 
#' This operation elongates an existing CSG object in a direction.
#'
#' @param object CSG object.
#' @param x Default `0`. Center of x-elongation.
#' @param y Default `0`. Center of y-elongation.
#' @param z Default `0`. Center of z-elongation.
#' @param elongate Default `c(0,0,0)` (no elongation). Elongation amount. 
#' @param robust Default `TRUE`. `FALSE` switches to a faster (but less robust in 2D) method.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Elongate a sphere to create a capsule in 1D or a rounded rectangle in 2D:
#' generate_ground(material=diffuse(checkercolor="grey20",color="dodgerblue4")) %>% 
#'  add_object(csg_object(csg_sphere(z=-3,x=-3),
#'                         material=glossy(color="purple"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_sphere(z=-3,x=3),x=3,z=-3, elongate = c(0.8,0,0)),
#'                         material=glossy(color="red"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_sphere(z=2),z=2, elongate = c(0.8,0,0.8)),
#'                         material=glossy(color="white"))) %>% 
#'  add_object(sphere(y=10,radius=3,material=light(intensity=8))) %>% 
#'  render_scene(clamp_value=10,fov=40,lookfrom=c(0,10,10))
#'   
#' #Elongate a torus:
#' generate_ground(material=diffuse(checkercolor="grey20",color="dodgerblue4")) %>% 
#'  add_object(csg_object(csg_torus(z=-3,x=-3),
#'                         material=glossy(color="purple"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_torus(z=-3,x=3),x=3,z=-3, elongate = c(0.8,0,0)),
#'                         material=glossy(color="red"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_torus(z=2),z=2, elongate = c(0.8,0,0.8)),
#'                         material=glossy(color="white"))) %>% 
#'  add_object(sphere(y=10,radius=3,material=light(intensity=8))) %>% 
#'  render_scene(clamp_value=10,fov=40,lookfrom=c(0,10,10))
#'  
#' #Elongate a cylinder:
#' generate_ground(material=diffuse(checkercolor="grey20",color="dodgerblue4")) %>% 
#'  add_object(csg_object(csg_cylinder(start=c(-3,0,-3), end = c(-3,1,-3)),
#'                         material=glossy(color="purple"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_cylinder(start=c(3,0,-3), end = c(3,1,-3)), x=3, z=-3, 
#'                        elongate = c(0.8,0,0)),
#'                        material=glossy(color="red"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_cylinder(start=c(0,0,3), end = c(0,1,3)), z=3, 
#'                        elongate = c(0.8,0,0.8)),
#'                        material=glossy(color="white"))) %>% 
#'  add_object(sphere(y=10,radius=3,material=light(intensity=8))) %>% 
#'  render_scene(clamp_value=10,fov=40,lookfrom=c(0,10,10))
#'  
#' #Elongate a pyramid:
#' generate_ground(material=diffuse(checkercolor="grey20",color="dodgerblue4")) %>% 
#'  add_object(csg_object(csg_pyramid(z=-3,x=-3),
#'                         material=glossy(color="purple"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_pyramid(z=-3,x=3),x=3,z=-3, elongate = c(0.8,0,0)),
#'                         material=glossy(color="red"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_pyramid(z=2),z=2, elongate = c(0.8,0,0.8)),
#'                         material=glossy(color="white"))) %>% 
#'  add_object(sphere(y=10,radius=3,material=light(intensity=8))) %>% 
#'  render_scene(clamp_value=10,fov=40,lookfrom=c(0,10,10))
#'
#' #Change the elongation point to start the elongation on the side of the pyramid:
#' generate_ground(material=diffuse(checkercolor="grey20",color="dodgerblue4")) %>% 
#'  add_object(csg_object(csg_pyramid(z=-3,x=-3),
#'                         material=glossy(color="purple"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_pyramid(z=-3,x=3),x=2.75,z=-2.75, elongate = c(0.8,0,0)),
#'                         material=glossy(color="red"))) %>% 
#'  add_object(csg_object(csg_elongate(csg_pyramid(z=2),z=2.25, elongate = c(0.8,0,0.8)),
#'                         material=glossy(color="white"))) %>% 
#'  add_object(sphere(y=10,radius=3,material=light(intensity=8))) %>% 
#'  render_scene(clamp_value=10,fov=40,
#'               lookfrom=c(5,5,10),lookat=c(0,0,-1.5))
#' }
csg_elongate = function(object, x=0,y=0,z=0, elongate = c(0,0,0), robust = TRUE) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  if(!robust) {
    csg_obj = list(csg_type = 15, object = object, x=x,y=y,z=z, elongate=elongate)
  } else {
    csg_obj = list(csg_type = 16, object = object, x=x,y=y,z=z, elongate=elongate)
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
#' \donttest{
#' #Generate a rounded pyramid:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_pyramid(x=-1,y=-0.99,z=1),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(csg_object(csg_round(csg_pyramid(x=1,y=-0.89)),
#'                         material=glossy(color="blue"))) %>% 
#'   add_object(csg_object(csg_round(csg_pyramid(x=0,z=-2,y=-0.5), radius=0.5),
#'                         material=glossy(color="green"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=1,material=light(intensity=50))) %>% 
#'   render_scene(lookfrom=c(-3,4,10), fov=22, 
#'                lookat=c(0,-0.5,0),clamp_value=10)
#'
#' #Round a blend of two objects
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_round(csg_combine(
#'     csg_pyramid(x=-0.5,y=-0.99,z=1.5),
#'     csg_pyramid(x=0.5,y=-0.99,z=2), operation="blend"), radius=0),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(csg_object(csg_round(csg_combine(
#'     csg_pyramid(x=-0.5,y=-0.79,z=-1.5),
#'     csg_pyramid(x=0.5,y=-0.79,z=-1), operation="blend"), radius=0.2),
#'                         material=glossy(color="green"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,radius=1,material=light(intensity=50))) %>% 
#'   render_scene(lookfrom=c(-3,5,10), fov=22, 
#'                lookat=c(0,-0.5,0),clamp_value=10)
#'}
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
#' Note: This operation has no overt effect on the external appearance of an object--it carves
#' regions on the interior. Thus, you will only see an effect with a transparent material or when
#' you carve into the object.
#'
#' @param object CSG object.
#' @param thickness Default `0.1`. Onioning distance.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Cut and onion a sphere:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>%
#'   add_object(csg_object(csg_combine(
#'     csg_onion(csg_sphere(z=2,x=2,radius=1), thickness = 0.2),
#'     csg_box(y=1,width=c(10,2,10)), operation = "subtract"),
#'     material=glossy(color="red"))) %>%
#'     add_object(csg_object(csg_combine(
#'       csg_onion(csg_sphere(radius=1), thickness = 0.4),
#'       csg_box(y=1,width=c(10,2,10)), operation = "subtract"),
#'       material=glossy(color="purple"))) %>%
#'     add_object(csg_object(csg_combine(
#'       csg_onion(csg_sphere(z=-2.5,x=-2.5,radius=1), thickness = 0.6),
#'       csg_box(y=1,width=c(10,2,10)), operation = "subtract"),
#'       material=glossy(color="green"))) %>%
#'  add_object(sphere(y=5,x=5,radius=2,material=light())) %>% 
#'  render_scene(clamp_value=10,lookat=c(0,-0.5,0),
#'               lookfrom=c(3,5,10),fov=35)
#'
#'#Multiple onion layers:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>%
#'   add_object(csg_object(csg_combine(
#'     csg_onion(csg_onion(csg_onion(csg_sphere(radius=1), 0.4), 0.2),0.1),
#'     csg_box(y=1,width=c(10,2,10)), operation = "subtract"),
#'     material=glossy(color="purple"))) %>%
#'   add_object(sphere(y=5,x=5,radius=2,material=light())) %>% 
#'   render_scene(clamp_value=10,lookat=c(0,-0.5,0),
#'                lookfrom=c(3,5,10),fov=20)
#'   
#'#Onion with dielectric sphere to make a bubble:
#' generate_cornell() %>%
#'   add_object(csg_object(
#'     csg_onion(csg_sphere(x=555/2,y=555/2,z=555/2, radius=150), 5),
#'     material=dielectric(attenuation=c(1,1,0.3)/100))) %>%
#'   render_scene(clamp_value=10)
#'   
#'#Multiple onion operations to make a bubble within a bubble:
#' generate_cornell() %>%
#'   add_object(csg_object(
#'     csg_onion(csg_onion(csg_sphere(x=555/2,y=555/2,z=555/2, radius=150), 10),5),
#'     material=dielectric(attenuation=c(1,1,0.3)/100))) %>%
#'   render_scene(clamp_value=10)
#'}
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
#' \donttest{
#' #Scale a pyramid (translating it upwards because the object is scaled from the center):
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_pyramid(z=1,y=-0.99),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(csg_object(csg_scale(csg_pyramid(z=-1,y=-0.5),2),
#'                         material=glossy(color="green"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,material=light(intensity=40))) %>% 
#'   render_scene(lookfrom=c(-3,4,10), fov=20, 
#'                lookat=c(0,-0.5,-0.5),clamp_value=10)
#'}
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
#' @param pivot_point Default `c(0,0,0)`. Pivot point for the rotation.
#' @param angles Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param up Default `c(0,1,0). Alternative method for specifying rotation--change the new "up" vector.
#' @param axis_x Default `NULL`, computed automatically if not passed. Given the `up` vector as the y-axis, this is the x vector.
#' @param axis_z Default `NULL`, computed automatically if not passed. Given the `up` vector as the y-axis, this is the z vector.
#' @return List describing the triangle in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Rotate a pyramid (translating it upwards because the object is scaled from the center):
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_pyramid(z=1,y=-0.99),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(csg_object(csg_rotate(csg_pyramid(z=-1.5,y=-0.99),
#'                         pivot_point = c(0,-0.99,-1.5),angle=c(0,45,0)),
#'                         material=glossy(color="green"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,material=light(intensity=40))) %>% 
#'   render_scene(lookfrom=c(-3,4,10), fov=15, 
#'                lookat=c(0,-0.5,0),clamp_value=10)
#'   
#' #Rotate by specifying a new up vector:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_pyramid(z=1,y=-0.99),
#'                         material=glossy(color="red"))) %>% 
#'   add_object(csg_object(csg_rotate(csg_pyramid(z=-1.5,y=-0.49),
#'                         pivot_point = c(0,-0.49,-1.5), up =c(1,1,0)),
#'                         material=glossy(color="green"))) %>% 
#'   add_object(sphere(y=5,x=5,z=5,material=light(intensity=40))) %>% 
#'   render_scene(lookfrom=c(-3,4,10), fov=15, 
#'                lookat=c(0,-0.5,0),clamp_value=10)
#'}
csg_rotate = function(object, pivot_point = c(0,0,0), 
                      angles = c(0,0,0), order_rotation = c(1,2,3),
                      up = c(0,1,0), axis_x = NULL, axis_z = NULL) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  if(any(angles != 0)) {
    an = angles/180*pi
    or = order_rotation
    mat1 = matrix(c(1,0,0,0,cos(an[or[1]]), sin(an[or[1]]), 0, -sin(an[or[1]]),  cos(an[or[1]])), nrow=3,ncol=3)
    mat2 = matrix(c(cos(an[or[2]]),0,-sin(an[or[2]]),0,1, 0, sin(an[or[2]]), 0,  cos(an[or[2]])), nrow=3,ncol=3)
    mat3 = matrix(c(cos(an[or[3]]), sin(an[or[3]]), 0, -sin(an[or[3]]),  cos(an[or[3]]),0,0,0,1), nrow=3,ncol=3)
    matfull = mat3 %*% mat2 %*% mat1
    axis_x = matfull[,1]
    up = matfull[,2]
    axis_z = matfull[,3]
  } else {
    if(is.null(axis_x) ) {
      axis_x = c(0,0,0)
    }
    if(is.null(axis_z)) {
      axis_z = c(0,0,0)
    }
  }
  csg_obj = list(csg_type = 20, object = object, pivot_point=pivot_point,up=up,
                 axis_x = axis_x, axis_z = axis_z)
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
#' \donttest{
#' #Translate a simple object:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_torus(), material=glossy(color="dodgerblue4"))) %>%
#'   add_object(csg_object(csg_translate(csg_torus(),x=-2,y=1,z=-2), 
#'                         material=glossy(color="red"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(0,5,10),fov=30,
#'                lookat=c(-1,0.5,-1))
#' 
#' #Translate a blended object:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'     csg_torus(),
#'     csg_torus(y=1, radius=0.8), operation="blend"), material=glossy(color="dodgerblue4"))) %>%
#'   add_object(csg_object(csg_translate(
#'     csg_combine(
#'       csg_torus(),
#'       csg_torus(y=1, radius=0.8), operation="blend"),
#'     x=-3,y=1,z=-3),
#'     material=glossy(color="red"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,lookfrom=c(0,5,10),fov=30,
#'                lookat=c(-1.5,0.5,-1.5))
#'}
csg_translate = function(object, x=0,y=0,z=0) {
  if(!inherits(object,"ray_csg")) {
    stop("`object` must be constructed with rayrender csg_* functions")
  }
  csg_obj = list(csg_type = 21, object = object, translate = c(x,y,z))
  class(csg_obj) = c("list", "ray_csg")
  csg_obj
}

#' CSG Combine
#' 
#' Note: Subtract operations aren't commutative: the second object is subtracted from the first. 
#'
#' @param object1 First CSG object
#' @param object2 Second CSG object
#' @param operation Default `union`. Can be `union`, `subtract`, `intersection`, `blend`, `subtractblend`, or `mix`.
#' @param radius Default `0.5`. Blending radius.
#'
#' @return List describing the combined csg object in the scene.
#' @export
#'
#' @examples
#' \donttest{
#' #Combine two spheres:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(x=-0.4,z=-0.4),
#'      csg_sphere(x=0.4,z=0.4), operation="union"),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-3,5,10))
#'   
#' #Subtract one sphere from another:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(x=-0.4,z=-0.4),
#'      csg_sphere(x=0.4,z=0.4), operation="subtract"),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-3,5,10))
#'   
#' #Get the intersection of two spheres:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(x=-0.4,z=-0.4),
#'      csg_sphere(x=0.4,z=0.4), operation="intersection"),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-3,5,10))
#'   
#' #Get the blended union of two spheres:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(x=-0.4,z=-0.4),
#'      csg_sphere(x=0.4,z=0.4), operation="blend"),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-3,5,10))
#'   
#' #Get the blended subtraction of two spheres:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(x=-0.4,z=-0.4),
#'      csg_sphere(x=0.4,z=0.4), operation="subtractblend"),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-3,5,10))
#'   
#' #Change the blending radius:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(x=-0.4,z=-0.4),
#'      csg_sphere(x=0.4,z=0.4), operation="blend", radius=0.2),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-3,5,10))
#'   
#' #Change the subtract blending radius:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(x=-0.4,z=-0.4),
#'      csg_sphere(x=0.4,z=0.4), operation="subtractblend", radius=0.2),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(sphere(y=5,x=5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-3,5,10))
#'   
#' #Get the mixture of various objects:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'      csg_sphere(),
#'      csg_box(), operation="mix"),
#'   material=glossy(color="dodgerblue4"))) %>%
#'   add_object(csg_object(csg_translate(csg_combine(
#'      csg_box(),
#'      csg_torus(), operation="mix"),z=-2.5),
#'   material=glossy(color="red"))) %>%
#'   add_object(csg_object(csg_translate(csg_combine(
#'      csg_pyramid(),
#'      csg_box(), operation="mix"),z=2.5),
#'   material=glossy(color="green"))) %>%
#'   add_object(sphere(y=10,x=-5,radius=3,material=light(intensity=10))) %>%  
#'   render_scene(clamp_value=10,fov=20,lookfrom=c(-15,10,10))
#'}
csg_combine = function(object1, object2, operation = "union", radius = 0.5) {
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
                 operation = combine_type, blend = radius)
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
#' \donttest{
#' #Group four spheres together and merge them with a box:
#' generate_ground(material=diffuse(checkercolor="grey20")) %>% 
#'   add_object(csg_object(csg_combine(
#'   csg_group(list(csg_sphere(x=1,z=1, radius=0.5),csg_sphere(x=-1,z=1, radius=0.5),
#'                  csg_sphere(x=1,z=-1, radius=0.5),csg_sphere(x=-1,z=-1, radius=0.5))),
#'   csg_box(y=0.5, width=c(2,0.2,2)), operation="blend"), material=glossy(color="red"))) %>%
#'   add_object(sphere(y=10,x=-5,radius=3,material=light(intensity=10))) %>% 
#'   render_scene(clamp_value=10,lookfrom=c(5,5,10)) 
#'}
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
