#' Process Points to Control Poinst
#'
#' @param points
#' @param closed 
#' @param straight
#' @keywords internal
process_point_series = function(points, closed=FALSE, straight=FALSE) {
  if(inherits(points,"numeric")) {
    stop("Input must either be list, matrix, or data.frame, not numeric.")
  }
  if(inherits(points,"list")) {
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
    if(closed && all(points[1,] != points[nrow(points),])) {
      points = rbind(points,points[1,])
    }
  }
  if(inherits(points,"matrix")) {
    if(ncol(points) == 3) {
      if(!straight) {
        full_control_points = rayrender:::calculate_control_points(points)
      } else {
        full_control_points = rayrender:::calculate_control_points_straight(points)
      }
    } else {
      stop("If points a matrix or data.frame, must have 3 columns")
    }
  } else {
    stop("points not of supported type (function expects matrix/data.frame/list, got ", class(points),")")
  }
  if(closed) {
    first_point = full_control_points[[1]]
    last_point = full_control_points[[length(full_control_points)-1]]
    full_control_points[[length(full_control_points)]][4,] = first_point[1,]
    full_control_points[[length(full_control_points)]][3,] = 2*first_point[1,] - first_point[2,]
    full_control_points[[length(full_control_points)]][2,] = 2*last_point[4,]  - last_point[3,]
  }
  return(full_control_points)
}

#' Process Points to Control Points
#'
#' @param points
#' @param closed 
#' @param straight
#' @keywords internal
process_point_series_1d = function(values, closed=FALSE, straight=FALSE) {
  mat_values = matrix(0,ncol=3,nrow=length(values))
  mat_values[,1] = values
  if(!straight) {
    full_control_points = rayrender:::calculate_control_points(mat_values)
  } else {
    full_control_points = rayrender:::calculate_control_points_straight(mat_values)
  }
  return(full_control_points)
}

#' Process Points to Control Points
#'
#' @param points
#' @param closed 
#' @param straight
#' @keywords internal
process_point_series_2d = function(values, closed=FALSE, straight=FALSE) {
  mat_values = matrix(0,ncol=3,nrow=length(values))
  mat_values[,1] = values[1]
  mat_values[,2] = values[2]
  if(!straight) {
    full_control_points = rayrender:::calculate_control_points(mat_values)
  } else {
    full_control_points = rayrender:::calculate_control_points_straight(mat_values)
  }
  return(full_control_points)
}

#' Lerp
#'
#' @param t
#' @param v1 
#' @param v2
#' @keywords internal
lerp = function(t, v1, v2) {
  return((1-t) * v1 + t * v2)
}

#' Evaluate Bezier
#'
#' @param cp Matrix (4x3)
#' @param t
#' @keywords internal
eval_bezier = function(cp,t) {
  return(cp[1,] * (1-t)^3 + cp[2,] * 3*t*(1-t)^2 + cp[3,] * 3*t^2*(1-t) + cp[4,] * t^3)
}

#' Get Distance Along Bezier Curve
#'
#' @param cps Control points
#' @param breaks
#' @keywords internal
calculate_distance_along_bezier_curve = function(cps,breaks=20) {
  distance_list = list()
  prev_pos = eval_bezier(cps[[1]], 0)
  distance_list[[1]] = data.frame(dist=0,segment=1,t=0, x=prev_pos[1],y=prev_pos[2],z=prev_pos[3])
  counter = 2
  for(seg in seq_len(length(cps))) {
    for(t in seq(0,1,length.out = breaks+1)[-1]) {
      temp_pos = eval_bezier(cps[[seg]], t)
      distance_list[[counter]] = data.frame(dist=sqrt(sum((temp_pos-prev_pos)^2)),
                                            segment=seg,t=(seg+t-1)/length(cps), x=temp_pos[1],y=temp_pos[2],z=temp_pos[3])
      prev_pos = temp_pos
      counter = counter + 1
    }
  }
  final_values = do.call(rbind,distance_list)
  final_values$total_dist = cumsum(final_values$dist)
  return(final_values)
}

#' Linearize and Calculate Final Points (with constant stepsize)
#'
#' @param linearized_cp Matrix (4x3)
#' @param steps
#' @keywords internal
calculate_final_path = function(linearized_cp, steps, constant_step = TRUE, offset = 0) {
  if(constant_step) {
    stepsize = max(linearized_cp$total_dist)/steps
    final_points = list()
    current_dist = offset
    for(i in 1:steps) {
      row = which.min(abs(floor(linearized_cp$total_dist - current_dist)))
      if(row+1 > nrow(linearized_cp)) {
        row = nrow(linearized_cp) - 1
      }
      tval = (current_dist - linearized_cp$total_dist[row])/(linearized_cp$total_dist[row+1] - linearized_cp$total_dist[row])
      
      final_points[[i]] = lerp(tval,
                               linearized_cp[row,c("x","y","z")],
                               linearized_cp[row+1,c("x","y","z")])
      current_dist = current_dist + stepsize
    }
  } else {
    rowvals = seq(1,nrow(linearized_cp), length.out = steps)
    final_points = list()
    for(i in 1:(length(rowvals)-1)) {
      final_points[[i]] = lerp(rowvals[i]-floor(rowvals[i]),
                               linearized_cp[floor(rowvals[i]),c("x","y","z")],
                               linearized_cp[max(c(floor(rowvals[i+1]),rowvals[i]+1)),c("x","y","z")])
    }
    final_points[[length(rowvals)]] = linearized_cp[nrow(linearized_cp),c("x","y","z")]
  }
  return(do.call(rbind,final_points))
}


#' Generate Camera Movement
#' 
#' Takes the scene description and renders an image, either to the device or to a filename. 
#'
#' @param positions
#' @param lookats
#' @param apertures
#' @param frames Default `30`. Number of frames between each key frame.
#' @param closed Default `FALSE`. Whether to close the camera curve.
#' @param straight Default `FALSE`. Whether the camera movement should follow straight lines.
#' @param constant_step Default `FALSE`. If `TRUE`, 
#' @param smooth Default `TRUE`.
#' @export
#' @return Nothing
#'
#' @examples
#' #Generate a camera moving through space
#' \donttest{
#' set.seed(3)
#' elliplist = list()
#' ellip_colors = rainbow(8)
#' for(i in 1:1200) {
#'   elliplist[[i]] = ellipsoid(x=10*runif(1)-5,y=10*runif(1)-5,z=10*runif(1)-5,
#'                              angle = 360*runif(3), a=0.1,b=0.05,c=0.1,
#'                              material=glossy(color=sample(ellip_colors,1)))
#' }
#' ellip_scene = do.call(rbind, elliplist)
#' 
#' camera_pos = list(c(0,1,15),c(5,-5,5),c(-5,5,-5),c(0,1,-15))
#' 
#' #Plot the camera path and render from above using the path object:
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(0,20,0),camera_up=c(0,0,1), width=800,height=800,samples=4,
#'                  fov=80)
#'             
#' #Side view     
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(20,0,0),width=800,height=800,samples=4,
#'                  fov=80)
#'  
#' #View from the start        
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(0,1.5,16),width=800,height=800,samples=4,
#'                  fov=80)
#'                  
#' #Generate Camera movement, setting the lookat position to be same as camera position, but offset
#' #slightly in front. We'll render 12 frames.
#' 
#' camera_motion =  generate_camera_motion(positions = camera_pos, lookats = camera_pos, 
#'                                         offset_lookat = 0.1, fovs=80, frames=12) 
#'                                         
#' #This returns a data frame of individual camera positions, interpolated by cubic bezier curves.
#' camera_motion
#' 
#' #Pass NA filename to plot to the device. We'll keep the path and offset it slightly to see
#' #where we're going. This results in a "roller coaster" effect.
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(obj_model(r_obj(),x=10,y=-10,scale_obj=3, angle=c(0,-45,0),
#'                        material=dielectric(attenuation=c(1,1,0.3)))) %>% 
#'   add_object(pig(x=-7,y=10,z=-5,scale=1,angle=c(0,-45,80),emotion="angry")) %>% 
#'   add_object(pig(x=0,y=-0.25,z=-15,scale=1,angle=c(0,225,-20), order_rotation=c(3,2,1),
#'                  emotion="angry", spider=TRUE)) %>% 
#'   add_object(path(camera_pos, y=-0.2,material=diffuse(color="red"))) %>% 
#'   render_animation(filename = NA, camera_motion = camera_motion, samples=100,
#'                    sample_method="sobol_blue",  
#'                    clamp_value=10, width=800, height=800)
#' }
generate_camera_motion = function(positions, lookats = NULL, apertures = 0, fovs = 40,
                                  focal_distances = NULL, ortho_dims = NULL,
                                  frames = 30, closed=FALSE, straight=FALSE, constant_step = TRUE,
                                  aperture_linear = TRUE,  fovs_linear = TRUE, focal_linear = TRUE,
                                  offset_lookat = 0) {
  position_control_points = process_point_series(positions,closed=closed,straight=straight)
  points_dist = calculate_distance_along_bezier_curve(position_control_points, breaks = 50)
  points_final = calculate_final_path(points_dist, steps = frames, constant_step = constant_step)
  if(!is.null(lookats)) {
    position_control_lookats = process_point_series(lookats,closed=closed,straight=straight)
    lookats_dist = calculate_distance_along_bezier_curve(position_control_lookats, breaks = 50)
    lookats_final = calculate_final_path(lookats_dist, steps = frames, constant_step = constant_step, offset = offset_lookat)
  } else {
    lookats_final = matrix(c(0,0,0), nrow = nrow(points_final), ncol=3, byrow=TRUE)
  }
  
  if(length(apertures) != 1) {
    apertures_control_lookats = process_point_series_1d(apertures,closed=closed,straight=aperture_linear)
    apertures_dist = calculate_distance_along_bezier_curve(apertures_control_lookats, breaks = 50)
    apertures_final = calculate_final_path(apertures_dist, steps = frames, constant_step = constant_step)
    apertures_final = apertures_final[,1]
  } else {
    apertures_final = matrix(apertures,  nrow = nrow(points_final), ncol=1)
  }
  if(length(fovs) != 1) {
    fovs_control_lookats = process_point_series_1d(fovs,closed=closed,straight=fovs_linear)
    fovs_dist = calculate_distance_along_bezier_curve(fovs_control_lookats, breaks = 50)
    fovs_final = calculate_final_path(fovs_dist, steps = frames, constant_step = constant_step)
    fovs_final = fovs_final[,1]
  } else {
    fovs_final = matrix(fovs,  nrow = nrow(points_final), ncol=1)
  }
  if(!is.null(focal_distances)) {
    focal_control_lookats = process_point_series_1d(focal_distances,closed=closed,straight=focal_linear)
    focal_dist = calculate_distance_along_bezier_curve(focal_control_lookats, breaks = 50)
    focal_final = calculate_final_path(focal_dist, steps = frames, constant_step = constant_step)
    focal_final = focal_final[,1]
  } else {
    focal_final = sqrt(apply((points_final-lookats_final)^2,1,sum))
  }
  if(!is.null(ortho_dims)) {
    ortho_control_lookats = process_point_series_2d(ortho_dims,closed=closed,straight=focal_linear)
    ortho_dist = calculate_distance_along_bezier_curve(ortho_control_lookats, breaks = 50)
    ortho_final = calculate_final_path(ortho_dist, steps = frames, constant_step = constant_step)
    ortho_final = ortho_final[,1:2]
  } else {
    ortho_final = matrix(c(1,1),  nrow = nrow(points_final), ncol=2)
  }
  return_values = as.data.frame(cbind(points_final,lookats_final,apertures_final,fovs_final,focal_final,ortho_final))
  colnames(return_values) = c("x","y","z","dx","dy","dz","aperture","fov","focal","orthox","orthoy")
  return(return_values)
}
