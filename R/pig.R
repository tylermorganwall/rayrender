#' Pig Object
#'
#' @param x Default `0`. x-coordinate of the center of the pig.
#' @param y Default `0`. y-coordinate of the center of the pig.
#' @param z Default `0`. z-coordinate of the center of the pig.
#' @param emotion Default `neutral`. Other options include `skeptical`, `worried`, and `angry`.
#' @param angle Default `c(0, 0, 0)`. Angle of rotation around the x, y, and z axes, applied in the order specified in `order_rotation`.
#' @param order_rotation Default `c(1, 2, 3)`. The order to apply the rotations, referring to "x", "y", and "z".
#' @param scale Default `c(1, 1, 1)`. Scale transformation in the x, y, and z directions. If this is a single value,
#' number, the object will be scaled uniformly.
#'
#' @return Single row of a tibble describing the pig in the scene.
#' @export
#'
#' @examples
#' #Generate a pig in the cornell box.
#' 
#' \donttest{
#' generate_cornell() %>%
#'   add_object(pig(x=555/2,z=555/2,y=120,scale=c(80,80,80), angle = c(0,135,0))) %>%
#'   render_scene(parallel=TRUE, samples=400,clamp_value=10)
#' }
#' 
#' # Show the pig staring into a mirror, worried 
#' \donttest{
#' generate_cornell() %>%
#'   add_object(pig(x=555/2-70,z=555/2+50,y=120,scale=c(80,80,80), 
#'                  angle = c(0,-40,0), emotion = "worried")) %>%
#'   add_object(xy_rect(x=450,z=450,y=250, ywidth=500, xwidth=200,  
#'                  angle = c(0,45,0), material = metal())) %>%
#'   render_scene(parallel=TRUE, samples=500,clamp_value=10)
#' }
#' 
#' # Render many small pigs facing random directions, with an evil pig overlord
#' set.seed(1)
#' lots_of_pigs = list() 
#' 
#' for(i in 1:10) {
#'   lots_of_pigs[[i]] = pig(x=50 + 450 * runif(1), z = 50 + 450 * runif(1), y=50, 
#'                              scale = c(30,30,30), angle = c(0,360*runif(1),0), emotion = "worried")
#' }
#' 
#' many_pigs_scene = do.call(rbind, lots_of_pigs) %>%
#'  add_object(generate_cornell(lightintensity=20)) %>%
#'  add_object(pig(z=500,x=555/2,y=400, emotion = "angry",
#'             scale=c(100,100,100),angle=c(30,90,0), order_rotation=c(2,1,3)))
#' \donttest{
#' render_scene(many_pigs_scene,parallel=TRUE,clamp_value=10, samples=500, tonemap = "reinhold")
#' }
pig = function(x = 0, y = 0, z = 0, emotion = "neutral",
               angle = c(0, 0, 0), order_rotation = c(1, 2, 3), 
               scale = c(1, 1, 1)) {
  if(length(scale) == 1) {
    scale = c(scale, scale, scale)
  }
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
