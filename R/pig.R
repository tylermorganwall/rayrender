#' Pig Object
#'
#' @param x Default `0`. x-coordinate of the center of the pig.
#' @param y Default `0`. y-coordinate of the center of the pig.
#' @param z Default `0`. z-coordinate of the center of the pig.
#' @param material Default  \code{\link{lambertian}}. The material, called from one of the material 
#' functions \code{\link{lambertian}}, \code{\link{metal}}, or \code{\link{dielectric}}.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".s
#'
#' @return Single row of a tibble describing the sphere in the scene.
#' @export
#'
#' @examples
#' #Generate a pig in the cornell box.
pig = function(x = 0, y = 0, z = 0, angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
               scale = c(1, 1, 1), emotion = "neutral") {
  spiral = list()
  tail_angles = seq(0,500,length.out = 33)
  for(i in 1:32) {
    spiral[[i]] = segment(start = c(-1.9 - 0.02 * i, 
                                    1 - 0.2 * cospi(tail_angles[i]/180- pi/2), 
                                    -0.2 * sinpi(tail_angles[i]/180- pi/2)),
                          end = c(-1.9 - 0.02 * (i+1), 
                                  1 - 0.2 * cospi(tail_angles[i+1]/180 - pi/2), 
                                  -0.2 * sinpi(tail_angles[i+1]/180- pi/2)),
                          radius=0.05,
                          material = lambertian(color="#f09089"))
  }
  spiralscene = do.call(rbind,spiral)
  if(emotion == "skeptical") {
    eyebrow_offset_left = c(0.1, -0.1)
    eyebrow_offset_right = c(0, 0)
  } else if (emotion == "worried") {
    eyebrow_offset_left = c(0.1, -0.1)
    eyebrow_offset_right = c(-0.1, 0.1)
  } else if (emotion == "angry") {
    eyebrow_offset_left = c(-0.1, 0.1)
    eyebrow_offset_right = c(0.1, -0.1)
  } else {
    eyebrow_offset_left = c(0, 0)
    eyebrow_offset_right = c(0, 0)
  }
  pig = group_objects(
    ellipsoid(y = 1, a = 2, b = 1, c = 1, material = lambertian(color="#f09089")) %>%
    add_object(sphere(x = 1.5, y = 2, z = 0, radius=0.8,material = lambertian(color="#f09089"))) %>%
    add_object(cylinder(x = 1, y = 0.4, z = 0.5, length=2,radius=0.25,material = lambertian(color="#f09089"))) %>%
    add_object(cylinder(x = 1, y = 0.4, z = -0.5, length=2,radius=0.25,material = lambertian(color="#f09089"))) %>%
    add_object(cylinder(x = -1, y = 0.4, z = 0.5, length=2,radius=0.25,material = lambertian(color="#f09089"))) %>%
    add_object(cylinder(x = -1, y = 0.4, z = -0.5, length=2,radius=0.25,material = lambertian(color="#f09089"))) %>%
    add_object(segment(start = c(1.5,2,0), end = c(2.5,2,0),radius=0.3,material = lambertian(color="#f09089"))) %>%
    add_object(sphere(x = 2, y = 2.5, z = 0.3, radius=0.25,material = lambertian(color="white"))) %>%
    add_object(sphere(x = 2, y = 2.5, z = -0.3, radius=0.25,material = lambertian(color="white"))) %>%
    add_object(sphere(x = 2.2, y = 2.5, z = 0.3, radius=0.1,material = lambertian(color="black"))) %>%
    add_object(sphere(x = 2.2, y = 2.5, z = -0.3, radius=0.1,material = lambertian(color="black"))) %>%
    add_object(spiralscene) %>%
    add_object(ellipsoid(x = 1.98+0.02, y = 1.6, z = 0, a= 0.2, b= 0.2, c= 0.3,material = lambertian(color="black"))) %>%
    add_object(ellipsoid(x = 2.5, y = 2, z = 0.1, a= 0.05, b= 0.1, c= 0.05,material = lambertian(color="black"))) %>%
    add_object(ellipsoid(x = 2.5, y = 2, z = -0.1, a= 0.05, b= 0.1, c= 0.05,material = lambertian(color="black"))) %>%
    add_object(disk(x = 2.5, y = 2, z = 0, radius=0.3, angle = c(0,0,90),material = lambertian(color="#f09089"))) %>%
    add_object(segment(start = c(2,2.8 + eyebrow_offset_left[2],-0.5), end = c(2,2.8+ eyebrow_offset_left[1],-0.1),
                       radius=0.05,material = lambertian(color="black"))) %>%
    add_object(segment(start = c(2,2.8 + eyebrow_offset_right[1],0.5), end = c(2,2.8+eyebrow_offset_right[2],0.1),
                       radius=0.05,material = lambertian(color="black"))), group_scale = scale,
    group_translate = c(x,y,z), group_angle = angle, group_order_rotation = order_rotation, pivot_point = c(0,1,0))
  return(pig)
}