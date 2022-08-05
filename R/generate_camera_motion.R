#' Generate Camera Movement
#' 
#' Takes a series of key frame camera positions and smoothly interpolates between them. Generates
#' a data.frame that can be passed to `render_animation()`.
#' 
#' @param positions A list or 3-column XYZ matrix of camera positions. 
#' These will serve as key frames for the camera position. Alternatively, this can also be the a 
#' dataframe of the keyframe output from an interactive rayrender session (`ray_keyframes`).
#' @param lookats Default `NULL`, which sets the camera lookat to the origin `c(0,0,0)` 
#' for the animation. A list or 3-column XYZ matrix
#' of `lookat` points. Must be the same number of points as `positions`.
#' @param apertures Default `0`. A numeric vector of aperture values.
#' @param fovs Default `40`. A numeric vector of field of view values.
#' @param focal_distances Default `NULL`, automatically the distance between positions and lookats. 
#' Numeric vector of focal distances.
#' @param ortho_dims Default `NULL`, which results in `c(1,1)` orthographic dimensions.  A list or 2-column matrix
#' of orthographic dimensions.
#' @param camera_ups Default `NULL`, which gives at up vector of `c(0,1,0)`. Camera up orientation.
#' @param type Default `cubic`. Type of transition between keyframes. 
#' Other options are `linear`, `quad`, `bezier`, `exp`, and `manual`. `manual` just returns the values 
#' passed in, properly formatted to be passed to `render_animation()`.
#' @param frames Default `30`. Total number of frames.
#' @param closed Default `FALSE`. Whether to close the camera curve so the first position matches the last. Set this to `TRUE` for perfect loops.
#' @param constant_step Default `TRUE`. This will make the camera travel at a constant speed. 
#' @param aperture_linear Default `TRUE`. This linearly interpolates focal distances, rather than using a smooth Bezier curve  or easing function.
#' @param fov_linear Default `TRUE`. This linearly interpolates focal distances, rather than using a smooth Bezier curve  or easing function.
#' @param focal_linear Default `TRUE`. This linearly interpolates focal distances, rather than using a smooth Bezier curve or easing function. 
#' @param ortho_linear Default `TRUE`. This linearly interpolates orthographic dimensions, rather than using a smooth Bezier curve or easing function. 
#' @param curvature_adjust Default `none`. Other options are `position`, `lookat`, and `both`. Whether to slow down the camera at areas of high curvature
#' to prevent fast swings. Only used for curve `type = bezier`. This does not preserve key frame positions.
#' Note: This feature will likely result in the `lookat` and `position` diverging if they do not 
#' have similar curvatures at each point. This feature is best used when passing the same set of points to `positions` and `lookats` 
#' and providing an `offset_lookat` value, which ensures the curvature will be the same.
#' @param curvature_scale Default `30`. Constant dividing factor for curvature. Higher values will subdivide the
#' path more, potentially finding a smoother path, but increasing the calculation time. Only used for curve `type = bezier`.
#' Increasing this value after a certain point will not increase the quality of the path, but it is scene-dependent.
#' @param offset_lookat Default `0`. Amount to offset the lookat position, either along the path (if `constant_step = TRUE`)
#' or towards the derivative of the Bezier curve. 
#' @param progress Default `TRUE`. Whether to display a progress bar.
#' 
#' @export
#' @return Data frame of camera positions, orientations, apertures, focal distances, and field of views
#'
#' @examples
#' #Create and animate flying through a scene on a simulated roller coaster
#' if(rayrender:::run_documentation()) {
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
#'   render_scene(lookfrom=c(0,20,0), width=800,height=800,samples=32,
#'                camera_up = c(0,0,1),
#'                fov=80)
#' }
#' if(rayrender:::run_documentation()) { 
#' #Side view     
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(20,0,0),width=800,height=800,samples=32,
#'                  fov=80)
#'  }
#' if(rayrender:::run_documentation()) {
#' #View from the start        
#' generate_ground(material=diffuse(checkercolor="grey20"),depth=-10) %>% 
#'   add_object(ellip_scene) %>% 
#'   add_object(sphere(y=50,radius=10,material=light(intensity=30))) %>% 
#'   add_object(path(camera_pos, material=diffuse(color="red"))) %>% 
#'   render_scene(lookfrom=c(0,1.5,16),width=800,height=800,samples=32,
#'                  fov=80)
#'  }
#' if(rayrender:::run_documentation()) {  
#' #Generate Camera movement, setting the lookat position to be same as camera position, but offset
#' #slightly in front. We'll render 12 frames, but you'd likely want more in a real animation.
#' 
#' camera_motion =  generate_camera_motion(positions = camera_pos, lookats = camera_pos, 
#'                                         offset_lookat = 1, fovs=80, frames=12,
#'                                         type="bezier") 
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
#'   add_object(pig(x=0,y=-0.25,z=-15,scale=1,angle=c(0,225,-20),
#'                  emotion="angry", spider=TRUE)) %>% 
#'   add_object(path(camera_pos, y=-0.2,material=diffuse(color="red"))) %>% 
#'   render_animation(filename = NA, camera_motion = camera_motion, samples=100,
#'                    sample_method="sobol_blue",  
#'                    clamp_value=10, width=400, height=400)
#' 
#' }
generate_camera_motion = function(positions, 
                                  lookats = NULL, 
                                  apertures = 0, fovs = 40,
                                  focal_distances = NULL, 
                                  ortho_dims = NULL, 
                                  camera_ups = NULL,
                                  type = "cubic",
                                  frames = 30, 
                                  closed = FALSE, 
                                  aperture_linear = TRUE,  
                                  fov_linear = TRUE, 
                                  focal_linear = TRUE,
                                  ortho_linear = TRUE,
                                  constant_step = TRUE, 
                                  curvature_adjust = "none", 
                                  curvature_scale = 30,
                                  offset_lookat = 0, 
                                  progress = TRUE) {
  curve_adj_type = unlist(lapply(tolower(curvature_adjust),switch,
                                 "none" = 0,"position" = 1,"lookats" = 2, "both" = 3, 0))
  curvature_adjust_pos = curve_adj_type %in% c(1,3)
  curvature_adjust_look = curve_adj_type %in% c(2,3)
  if(!is.null(dim(positions)) && ncol(positions) == 14) {
    temp = positions
    positions = temp[,1:3]
    lookats = temp[,4:6]
    apertures = temp[,7]
    fovs = temp[,8]
    focal_distances = temp[,9]
    ortho_dims = temp[,10:11]
    camera_ups = temp[,12:14]
  }
  if(type=="bezier") {
    position_control_points = process_point_series(positions,closed=closed,straight=FALSE)
    points_dist = calculate_distance_along_bezier_curve(position_control_points, breaks = 50)
    points_final = calculate_final_path(points_dist, steps = frames, constant_step = constant_step,
                                        curvature_adjust = curvature_adjust_pos, 
                                        curvature_scale=curvature_scale,
                                        progress = progress, string = "Position   ")
    if(!is.null(lookats)) {
      position_control_lookats = process_point_series(lookats,closed=closed,straight=FALSE)
      lookats_dist = calculate_distance_along_bezier_curve(position_control_lookats, breaks = 50)
      lookats_final = calculate_final_path(lookats_dist, steps = frames, 
                                           constant_step = constant_step, offset = offset_lookat,
                                           curvature_adjust = curvature_adjust_look, 
                                           curvature_scale=curvature_scale,
                                           progress = progress, string = "Orientation")
    } else {
      lookats_final = matrix(c(0,0,0), nrow = nrow(points_final), ncol=3, byrow=TRUE)
    }
    if(curvature_adjust_pos || curvature_adjust_look) {
      temp_points  = points_final[1:frames,]
      temp_lookats = lookats_final[1:frames,]
      
      rowvals = seq(1,nrow(points_final),length.out = frames)
      for(i in 1:frames) {
        rv = rowvals[i]
        frv = floor(rv)
        if(frv == rv) {
          temp_points[i,] = points_final[rv,]
          temp_lookats[i,] = lookats_final[rv,]
          next
        }
        temp_points[i,] = lerp(rv-frv,
                               points_final[rv,],
                               points_final[frv+1,])
        temp_lookats[i,] = lerp(rv-frv,
                                lookats_final[rv,],
                                lookats_final[frv+1,])
      }
      if(curvature_adjust_pos) {
        points_final = temp_points
      }
      if(curvature_adjust_look) {
        lookats_final = temp_lookats
      }
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
      fovs_control_lookats = process_point_series_1d(fovs,closed=closed,straight=fov_linear)
      fovs_dist = calculate_distance_along_bezier_curve(fovs_control_lookats, breaks = 50)
      fovs_final = calculate_final_path(fovs_dist, steps = frames, constant_step = constant_step)
      fovs_final = fovs_final[,1]
    } else {
      fovs_final = matrix(fovs,  nrow = nrow(points_final), ncol=1)
    }
    if(!is.null(focal_distances)) {
      if(is.numeric(focal_distances)) {
        if(length(focal_distances) == 1) {
          focal_final = rep(focal_distances,nrow(points_final))
        } else {
          focal_final = tween(focal_distances,n=frames,ease="linear")
        }
      } else {
        focal_control_lookats = process_point_series_1d(focal_distances,closed=closed,straight=focal_linear)
        focal_dist = calculate_distance_along_bezier_curve(focal_control_lookats, breaks = 50)
        focal_final = calculate_final_path(focal_dist, steps = frames, constant_step = constant_step)
        focal_final = focal_final[,1]
      }
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
    if(!is.null(camera_ups)) {
      ups_control = process_point_series_2d(camera_ups,closed=closed,straight=ortho_linear)
      up_dist = calculate_distance_along_bezier_curve(ups_control, breaks = 50)
      up_final = calculate_final_path(up_dist, steps = frames, constant_step = constant_step)
    } else {
      up_final = matrix(c(0,1,0), nrow = nrow(points_final), ncol=3, byrow=TRUE)
    }
    return_values = as.data.frame(cbind(points_final,lookats_final,apertures_final,fovs_final,focal_final,ortho_final,
                                        up_final))
    rownames(return_values) = NULL
    colnames(return_values) = c("x","y","z",
                                "dx","dy","dz",
                                "aperture","fov","focal","orthox","orthoy",
                                "upx","upy","upz")
    return(return_values)
  } else if (type %in% c("exp", "quad",  "cubic", "linear")) {
    if(inherits(positions,"list")) {
      positions = do.call(rbind,positions)
    }
    if(inherits(lookats,"list")) {
      lookats = do.call(rbind,lookats)
    }
    if(inherits(camera_ups,"list")) {
      camera_ups = do.call(rbind,camera_ups)
    }
    if(inherits(ortho_dims,"list")) {
      ortho_dims = do.call(rbind,ortho_dims)
    }
    positions = grDevices::xyz.coords(positions)
    
    if(is.null(lookats)) {
      lookats = matrix(0,ncol=3,nrow=length(positions$x))
    }
    lookats = grDevices::xyz.coords(lookats)
    if(is.null(camera_ups)) {
      camera_ups = matrix(c(0,1,0),ncol=3,nrow=length(positions$x),byrow=TRUE)
    }
    camera_ups = grDevices::xyz.coords(camera_ups)
    
    if(is.null(focal_distances)) {
      focal_distances = sqrt((positions$x - lookats$x)^2 + 
                             (positions$y - lookats$y)^2 + 
                             (positions$z - lookats$z)^2)
    } else if (length(focal_distances) == 1) {
      focal_distances = rep(focal_distances,length(positions$x))
    }
    if(is.null(ortho_dims)) {
      ortho_dims = matrix(c(1,1),ncol=2,nrow=length(positions$x),byrow=TRUE)
    }
    ortho = grDevices::xy.coords(ortho_dims)
    if(length(apertures) == 1) {
      apertures = rep(apertures,length(positions$x))
    }
    if(length(apertures) == 1) {
      fovs = rep(fovs,length(positions$x))
    }
    tween_df = data.frame(x=positions$x, y=positions$y, z=positions$z,
                          dx=lookats$x, dy=lookats$y, dz=lookats$z,
                          aperture = apertures,
                          fov = fovs,
                          focal = focal_distances,
                          orthox = ortho$x, orthoy = ortho$y,
                          upx =camera_ups$x, upy = camera_ups$y ,upz = camera_ups$z)
    if(closed) {
      tween_df = rbind(tween_df,tween_df[1,])
      apertures = c(apertures,apertures[1])
      fovs = c(fovs,fovs[1])
      focal_distances = c(focal_distances,focal_distances[1])
    }
    return_values = as.data.frame(apply(tween_df,2,tween, n = frames,ease=type))
    rownames(return_values) = NULL
    if(aperture_linear) {
      return_values$aperture =tween(apertures,n=frames,ease="linear")
    }
    if(fov_linear) {
      return_values$fov = tween(fovs,n=frames,ease="linear")
    }
    if(focal_linear) {
      return_values$focal = tween(focal_distances,n=frames,ease="linear")
    }
    if(ortho_linear) {
      return_values$orthox = tween(ortho$x,n=frames,ease="linear")
      return_values$orthoy = tween(ortho$y,n=frames,ease="linear")
    }
    return(return_values)
  } else if (type == "manual") {
    if(inherits(positions,"list")) {
      positions = do.call(rbind,positions)
    }
    if(inherits(lookats,"list")) {
      lookats = do.call(rbind,lookats)
    }
    if(inherits(camera_ups,"list")) {
      camera_ups = do.call(rbind,camera_ups)
    }
    if(inherits(ortho_dims,"list")) {
      ortho_dims = do.call(rbind,ortho_dims)
    }
    if (is.numeric(ortho_dims) && length(ortho_dims) == 2) {
      ortho_dims = matrix(ortho_dims,ncol=2,nrow=length(positions$x),byrow=TRUE)
    }
    
    positions = grDevices::xyz.coords(positions)
    
    if(is.null(lookats)) {
      lookats = matrix(0,ncol=3,nrow=length(positions$x))
    } else if (is.numeric(lookats) && length(lookats) == 3) {
      lookats = matrix(lookats,ncol=3,nrow=length(positions$x),byrow=TRUE)
    }
    lookats = grDevices::xyz.coords(lookats)
    if(is.null(camera_ups)) {
      camera_ups = matrix(c(0,1,0),ncol=3,nrow=length(positions$x),byrow=TRUE)
    }
    camera_ups = grDevices::xyz.coords(camera_ups)
    
    if(is.null(focal_distances)) {
      focal_distances = sqrt((positions$x - lookats$x)^2 + 
                               (positions$y - lookats$y)^2 + 
                               (positions$z - lookats$z)^2)
    } else if (length(focal_distances) == 1) {
      focal_distances = rep(focal_distances,length(positions$x))
    }
    if(is.null(ortho_dims)) {
      ortho_dims = matrix(c(1,1),ncol=2,nrow=length(positions$x),byrow=TRUE)
    }
    ortho = grDevices::xy.coords(ortho_dims)
    return(data.frame(x=positions$x, y=positions$y, z=positions$z,
                      dx=lookats$x, dy=lookats$y, dz=lookats$z,
                      aperture = rep(apertures,length(positions$x)),
                      fov = rep(fovs,length(positions$x)),
                      focal = focal_distances,
                      orthox = ortho$x, orthoy = ortho$y,
                      upx =camera_ups$x, upy = camera_ups$y ,upz = camera_ups$z))
  } else {
    stop("type '", type, "' not recognized")
  }
}

#' Process Points to Control Points
#'
#' @param points Points
#' @param closed Whether to be closed
#' @param straight Whether to be straight
#' @return Matrix of control points
#' 
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
  }
  if(inherits(points,"matrix")) {
    if(ncol(points) == 3) {
      if(!straight) {
        full_control_points = calculate_control_points(points)
      } else {
        full_control_points = calculate_control_points_straight(points)
      }
    } else {
      stop("If points a matrix or data.frame, must have 3 columns")
    }
  } else {
    stop("points not of supported type (function expects matrix/data.frame/list, got ", class(points),")")
  }
  if(closed) {
    first_point = full_control_points[[1]]
    last_point = full_control_points[[length(full_control_points)]]
    full_control_points[[length(full_control_points) + 1]] = last_point
    full_control_points[[length(full_control_points)]][4,] = first_point[1,]
    full_control_points[[length(full_control_points)]][3,] = 2*first_point[1,] - first_point[2,]
    full_control_points[[length(full_control_points)]][2,] = 2*last_point[4,]  - last_point[3,]
    full_control_points[[length(full_control_points)]][1,] = last_point[4,]
  }
  return(full_control_points)
}

#' Process Points to Control Points
#'
#' @param points Points
#' @param closed Whether to be closed
#' @param straight Whether to be straight
#' @return Matrix of control points
#' 
#' @keywords internal
process_point_series_1d = function(values, closed=FALSE, straight=FALSE) {
  mat_values = matrix(0,ncol=3,nrow=length(values))
  mat_values[,1] = values
  if(!straight) {
    full_control_points = calculate_control_points(mat_values)
  } else {
    full_control_points = calculate_control_points_straight(mat_values)
  }
  return(full_control_points)
}

#' Process Points to Control Points
#'
#' @param points Points
#' @param closed Whether to be closed
#' @param straight Whether to be straight
#' @return Matrix of control points
#' 
#' @keywords internal
process_point_series_2d = function(values, closed=FALSE, straight=FALSE) {
  mat_values = matrix(0,ncol=3,nrow=nrow(values))
  mat_values[,1] = values[,1]
  mat_values[,2] = values[,2]
  if(!straight) {
    full_control_points = calculate_control_points(mat_values)
  } else {
    full_control_points = calculate_control_points_straight(mat_values)
  }
  return(full_control_points)
}

#' Evaluate Bezier
#'
#' @param cp Control point matrix (4x3)
#' @param t Interpolation distance
#' @return 3D Numeric value for point in space
#' 
#' @keywords internal
eval_bezier = function(cp,t) {
  return(cp[1,] * (1-t)^3 + cp[2,] * 3*t*(1-t)^2 + cp[3,] * 3*t^2*(1-t) + cp[4,] * t^3)
}

#' Evaluate Deriv Bezier
#'
#' @param cp Control point matrix (4x3)
#' @param t Interpolation distance
#' @return 3D Numeric value for point in space
#' 
#' @keywords internal
eval_bezier_deriv = function(cp,t) {
  return(-3*cp[1,] * (1-t)^2 + 
           cp[2,] * 3*(1-t)^2 - cp[2,] * 6*t*(1-t) +
           cp[3,] * 6*t*(1-t) - cp[3,] * 3*t^2 + 
           3*cp[4,] * t^2)
}

#' Evaluate Deriv Bezier
#'
#' @param cp Control point matrix (4x3)
#' @param t Interpolation distance
#' @return 3D Numeric value for point in space
#' 
#' @keywords internal
eval_bezier_2nd_deriv = function(cp,t) {
  return(6*cp[1,] * (1-t) - 
           cp[2,] * 6 * (1 - t) - cp[2,] * 6 *(1-t) + cp[2,] * 6*t + 
           cp[3,] * 6 * (1 - t) - cp[3,] * 6 * t    - cp[3,] * 6*t + 
         6*cp[4,] * t)
}

#' Cross product (vec)
#'
#' @param x vec1
#' @param y vec2
#' @return 3D Numeric value for vector in space
#' 
#' @keywords internal
cross_prod = function(x,y) {
  return(c(x[2]*y[3]-x[3]*y[2],
           -(x[1]*y[3]-x[3]*y[1]),
           x[1]*y[2]-x[2]*y[1]))
}

#' Get Distance Along Bezier Curve
#'
#' @param cps Control points
#' @param breaks Number of interpolation breaks
#' @return Data frame of points along curve, along with distances
#' 
#' @keywords internal
calculate_distance_along_bezier_curve = function(cps,breaks=20) {
  distance_list = list()
  prev_pos = eval_bezier(cps[[1]], 0)
  temp_deriv = eval_bezier_deriv(cps[[1]], 0)
  temp_second_deriv = eval_bezier_2nd_deriv(cps[[1]], 0)
  if(any(temp_deriv != 0)) {
    curve = sqrt(sum(cross_prod(temp_deriv,temp_second_deriv)^2))/sum(temp_deriv^2)^(3/2)
  } else {
    curve = Inf
  }
  
  distance_list[[1]] = data.frame(dist=0,segment=1,t=0, x=prev_pos[1],y=prev_pos[2],z=prev_pos[3],
                                  curvature = curve,
                                  dx=temp_deriv[1],dy=temp_deriv[2],dz=temp_deriv[3])
  counter = 2
  for(seg in seq_len(length(cps))) {
    for(t in seq(0,1,length.out = breaks+1)[-1]) {
      temp_pos = eval_bezier(cps[[seg]], t)
      temp_deriv = eval_bezier_deriv(cps[[seg]], t)
      temp_second_deriv = eval_bezier_2nd_deriv(cps[[seg]], t)
      if(any(temp_deriv != 0)) {
        curve = sqrt(sum(cross_prod(temp_deriv,temp_second_deriv)^2))/sum(temp_deriv^2)^(3/2)
        distance_list[[counter]] = data.frame(dist=sqrt(sum((temp_pos-prev_pos)^2)),
                                              segment=seg,t=(seg+t-1)/length(cps), x=temp_pos[1],y=temp_pos[2],z=temp_pos[3],
                                              curvature = curve,
                                              dx=temp_deriv[1],dy=temp_deriv[2],dz=temp_deriv[3])
      } else {
        curve = Inf
        distance_list[[counter]] = data.frame(dist=sqrt(sum((temp_pos-prev_pos)^2)),
                                              segment=seg,t=(seg+t-1)/length(cps), x=temp_pos[1],y=temp_pos[2],z=temp_pos[3],
                                              curvature = curve,
                                              dx=temp_deriv[1],dy=temp_deriv[2],dz=temp_deriv[3])
      }
      
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
#' @param steps Number of steps
#' @param constant_step Whether to step at constant steps or not
#' @param offset Offset along curve
#' @return Matrix of points along curve
#' 
#' @keywords internal
calculate_final_path = function(linearized_cp, steps, constant_step = TRUE, 
                                curvature_adjust = FALSE, curvature_scale = 1, offset = 0,
                                progress = FALSE, string = "") {
  if(max(linearized_cp$total_dist) < 1e-7) {
    single_row = linearized_cp[1,c("x","y","z")]
    single_df = as.data.frame(matrix(as.numeric(single_row),nrow=steps,ncol=3, byrow=T))
    colnames(single_df) = c("x","y","z")
    return(single_df)
  }
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
    if(!curvature_adjust) {
      rowvals = seq(1,nrow(linearized_cp), length.out = steps)
      final_points = list()
      for(i in 1:(length(rowvals)-1)) {
        final_points[[i]] = lerp(rowvals[i]-floor(rowvals[i]),
                                 linearized_cp[floor(rowvals[i]),c("x","y","z")],
                                 linearized_cp[floor(rowvals[i])+1,c("x","y","z")])
        if(offset != 0) {
          direction = lerp(rowvals[i]-floor(rowvals[i]),
                           linearized_cp[floor(rowvals[i]),c("dx","dy","dz")],
                           linearized_cp[floor(rowvals[i])+1,c("dx","dy","dz")]) 
          if(all(direction != 0)) {
            direction = direction / sqrt(sum(direction^2))
          } else {
            direction = c(0,0,0)
          }
          final_points[[i]] = final_points[[i]] + direction * offset
        }
      }
      if(offset != 0) {
        direction = linearized_cp[nrow(linearized_cp),c("dx","dy","dz")]
        direction = direction / sqrt(sum(direction^2))
        final_points[[length(rowvals)]] = linearized_cp[nrow(linearized_cp),c("x","y","z")] + direction * offset
      } else {
        final_points[[length(rowvals)]] = linearized_cp[nrow(linearized_cp),c("x","y","z")]
      }
    } else {
      maxdist = max(linearized_cp$total_dist)
      default_stepsize = maxdist/steps
      final_points = list()
      current_dist = 0
      i = 1
      if(progress) {
        pb = progress::progress_bar$new(format = sprintf("Generating Camera %s [:bar] :percent",string),
                                        total = maxdist, width=70)
      }
      while(current_dist < maxdist) {
        if(progress) {
          pb$update(current_dist/maxdist)
        }
        row = which.min(abs(linearized_cp$total_dist - current_dist))
        if(linearized_cp$total_dist[row] - current_dist > 0) {
          row = row - 1
        }
        if(row+1 > nrow(linearized_cp)) {
          row = nrow(linearized_cp) - 1
        }
        
        tval = (current_dist - linearized_cp$total_dist[row])/(linearized_cp$total_dist[row+1] - linearized_cp$total_dist[row])
        
        final_points[[i]] = lerp(tval,
                                 linearized_cp[row,c("x","y","z")],
                                 linearized_cp[row+1,c("x","y","z")])
        direction = lerp(tval,
                         linearized_cp[row,c("dx","dy","dz")],
                         linearized_cp[row+1,c("dx","dy","dz")])

        direction = direction/sqrt(sum(direction^2))
        final_points[[i]] = final_points[[i]] + direction * offset
        temp_curve = lerp(tval,
                          linearized_cp[row,c("curvature")],
                          linearized_cp[row+1,c("curvature")])
        step = min(c(1/temp_curve/curvature_scale, default_stepsize))
        current_dist = current_dist + step
        i = i + 1
      }
    }
  }
  return(do.call(rbind,final_points))
}

#' Get Saved Keyframes
#'
#' @description Get a dataframe of the saved keyframes (using the interactive renderer) to pass to `generate_camera_motion()`
#' @return Data frame of keyframes
#' 
#' @export
#' @examples 
#' #This will return an empty data frame if no keyframes have been set.
#' get_saved_keyframes()
get_saved_keyframes = function() {
  keyframes = get("keyframes",  envir = ray_environment)
  if(nrow(keyframes) == 0) {
    message("No keyframes saved: press K when using interactive preview to save keyframe")
    return(keyframes)
  }
  return(keyframes)
}
