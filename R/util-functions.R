#' Calculate Control Points
#'
#' @keywords internal
calculate_control_points = function(s_mat) {
  nr = nrow(s_mat)
  nr2 = nr-2
  if(nr == 1) {
    stop("Only one point passed, unable to draw curve.")
  }
  if(nr == 2) {
    vec = s_mat[2,] - s_mat[1,]
    return(list(matrix(c(s_mat[1,],
                         s_mat[1,] + 1/3*vec,
                         s_mat[1,] + 2/3*vec,
                         s_mat[2,]),ncol=3,byrow=TRUE)))
  }
  spline_matrix = diag(nr2) * 4 + diag(nr2+1)[-1,-(nr2+1)] + t(diag(nr2+1)[-1,-(nr2+1)])
  inv_spline_matrix = solve(spline_matrix)
  if(nr == 3) {
    new_b = 1/4 * (6 * s_mat[2,] - s_mat[1,] - s_mat[3,])
    b_vec = rbind(s_mat[1,],new_b,s_mat[3,])
    return_points = list()
    for(i in seq_len(nrow(s_mat)-1)) {
      vec = b_vec[i+1,] - b_vec[i,]
      return_points[[i]] = matrix(c(s_mat[i,],
                                    b_vec[i,] + 1/3*vec,
                                    b_vec[i,] + 2/3*vec,
                                    s_mat[i+1,]),ncol=3,byrow=TRUE)
    }
    return(return_points)
  }
  s_vec = matrix(0,nrow=nr-2,ncol=3)
  s_vec[1,] = 6 * s_mat[2,] - s_mat[1,]
  for(i in seq_len(nr-4)) {
    s_vec[i+1,] = 6 * s_mat[i+2,]
  }
  s_vec[nr-2,] = 6 * s_mat[nr-1,] - s_mat[nr,]
  b_vec = inv_spline_matrix %*% s_vec
  b_vec = rbind(s_mat[1,],b_vec,s_mat[nr,])
  return_points = list()
  for(i in seq_len(nrow(s_mat)-1)) {
    vec = b_vec[i+1,] - b_vec[i,]
    return_points[[i]] = matrix(c(s_mat[i,],
                                  b_vec[i,] + 1/3*vec,
                                  b_vec[i,] + 2/3*vec,
                                  s_mat[i+1,]),ncol=3,byrow=TRUE)
  }
  return(return_points)
}

#' Calculate Control Points (straight)
#'
#' @keywords internal
calculate_control_points_straight = function(s_mat) {
  nr = nrow(s_mat)
  if(nr == 1) {
    stop("Only one point passed, unable to draw curve.")
  }
  return_points = list()
  for(i in seq_len(nrow(s_mat)-1)) {
    vec = s_mat[i+1,] - s_mat[i,]
    return_points[[i]] = matrix(c(s_mat[i,],
                                  s_mat[i,] + 1/3*vec,
                                  s_mat[i,] + 2/3*vec,
                                  s_mat[i+1,]),ncol=3,byrow=TRUE)
  }
  return(return_points)
}

#' Clamp Values
#'
#' @keywords internal
clamp = function(v, min=0, max=Inf) {
  v[v < min] = min
  v[v > max] = max
  v
}

#' Lerp
#'
#' @param t Interpolation distance
#' @param v1 Value 1
#' @param v2 Value 2
#' @return Linearly interpolated value
#' 
#' @keywords internal
lerp = function(t, v1, v2) {
  return((1-t) * v1 + t * v2)
}

#' Quad-in-out
#' 
#' @param t Value
#' @return number
#'
#' @keywords internal
quadInOut = function(t) {
  ifelse(t * 2 <= 1, (2*t) ^ 2/2, (2-(2*t-2)^2)/2)
}

#' Cubic-in-out
#' 
#' @param t Value
#' @return number
#'
#' @keywords internal
cubicInOut = function(t) {
  ifelse(t * 2 <= 1, (2*t) ^ 3/2, ((2*t-2)^3+2)/2)
}

#' Cubic-in-out
#' 
#' @param t Value
#' @return number
#'
#' @keywords internal
expInOut = function(t) {
  ifelse(t * 2 <= 1, 2^(-10 * (1-2*t))/2, (2 - 2^(-10 * (2*t-1)))/2)
}

#' Cubic-in-out
#' 
#' @param vals Numeric values
#' @param n Frames
#' @param type typenae
#' @return number
#'
#' @keywords internal
tween = function(vals, n, ease = "cubic") {
  len_vals = rep(0, length(vals)-1)
  free_vals = n - length(vals)
  counter = 1
  for(i in seq_len(free_vals)) {
    len_vals[counter] = len_vals[counter] + 1
    counter = counter + 1
    if(counter > (length(vals)-1)) {
      counter = 1
    }
  }
  tlist = list()
  for(i in seq_len((length(vals)-1))) {
    if(ease == "cubic") {
      tlist[[i]] = cubicInOut(seq(0,1,length.out = len_vals[i]+2)[-(len_vals[i]+2)])
    } else if(ease == "exp") {
      tlist[[i]] = expInOut(seq(0,1,length.out = len_vals[i]+2)[-(len_vals[i]+2)])
    } else if(ease == "quad") {
      tlist[[i]] = quadInOut(seq(0,1,length.out = len_vals[i]+2)[-(len_vals[i]+2)])
    } else {
      tlist[[i]] = seq(0,1,length.out = len_vals[i]+2)[-(len_vals[i]+2)]
    }
  }
  final_vals = list()
  counter = 1
  for(i in seq_len(length(tlist))) {
    seg = tlist[[i]]
    for(j in seq_len(length(seg))) {
      final_vals[[counter]] = lerp(seg[j], vals[i], vals[i+1])
      counter = counter + 1
    }
  }
  final_vals[[counter]] = vals[length(vals)]
  if(length(final_vals) != n) {
    stop("Length of interpolated sequence (", length(final_vals) ,") not equal to frames (", n, ")")
  }
  return(unlist(final_vals))
}

#' Generate Translation Matrix
#' 
#' @param delta Distance 
#' @return number
#'
#' @keywords internal
#' 
generate_translation_matrix = function(delta) {
  m = matrix(c(1, 0, 0, delta[1], 0, 1, 0, delta[2], 0, 0, 1, delta[3], 0, 0, 0,
              1),4,4,byrow=T)
  return(m)
}

#' Generate Rotation Matrix (order)
#' 
#' @param angles Angles
#' @param order_rotation Order of Rotation
#' @return number
#'
#' @keywords internal
#' 
generate_rotation_matrix = function(angles, order_rotation) {
  M = diag(4)
  for(i in 1:3) {
    if(order_rotation[i] == 1) {
      if(angles[1] != 0) {
        M = RotateX(angles[1]) %*% M
      }
    }
    if(order_rotation[i] == 2) {
      if(angles[2] != 0) {
        M = RotateY(angles[2])  %*% M
      }
    }
    if(order_rotation[i] == 3) {
      if(angles[3] != 0) {
        M = RotateZ(angles[3])  %*% M
      }
    }
  }
  return(M)
}

#' Generate Rotation Matrix X
#' 
#' @param theta Angle
#' @return number
#'
#' @keywords internal
#' 
RotateX = function(theta) {
  sinTheta = sinpi(theta/180)
  cosTheta = cospi(theta/180)
  M = matrix(c(1, 0, 0, 0, 0, cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0,
              0, 0, 0, 1),4,4,byrow=T)
  return(M)
}

#' Generate Rotation Matrix Y
#' 
#' @param theta Angle
#' @return number
#'
#' @keywords internal
#' 
RotateY = function(theta) {
  sinTheta = sinpi(theta/180)
  cosTheta = cospi(theta/180)
  M = matrix(c(cosTheta, 0, sinTheta, 0, 0, 1, 0, 0, -sinTheta, 0, cosTheta, 0,
              0, 0, 0, 1),4,4,byrow=T)
  return(M)
}

#' Generate Rotation Matrix Z
#' 
#' @param theta Angle
#' @return number
#'
#' @keywords internal
#' 
RotateZ = function(theta) {
  sinTheta = sinpi(theta/180)
  cosTheta = cospi(theta/180)
  M = matrix(c(cosTheta, -sinTheta, 0, 0, sinTheta, cosTheta, 0, 0, 0, 0, 1, 0,
              0, 0, 0, 1),4,4)
  return(M)
}

#' Generate Rotation Matrix Axis
#' 
#' @param theta Angle
#' @param axis The rotation axis
#' @return matrix
#'
#' @keywords internal
#' 
RotateAxis = function(theta, axis) {
  a = axis/sqrt(sum(axis*axis))
  sinTheta = sinpi(theta/180)
  cosTheta = cospi(theta/180)
  m = diag(4)
  # Compute rotation of first basis vector
  m[1,1] = a[1] * a[1] + (1 - a[1] * a[1]) * cosTheta
  m[1,2] = a[1] * a[2] * (1 - cosTheta) - a[3] * sinTheta
  m[1,3] = a[1] * a[3] * (1 - cosTheta) + a[2] * sinTheta
  m[1,4] = 0
  
  # Compute rotations of second and third basis vectors
  m[2,1] = a[1] * a[2] * (1 - cosTheta) + a[3] * sinTheta
  m[2,2] = a[2] * a[2] + (1 - a[2] * a[2]) * cosTheta
  m[2,3] = a[2] * a[3] * (1 - cosTheta) - a[1] * sinTheta
  m[2,4] = 0
  
  m[3,1] = a[1] * a[3] * (1 - cosTheta) - a[2] * sinTheta
  m[3,2] = a[2] * a[3] * (1 - cosTheta) + a[1] * sinTheta
  m[3,3] = a[3] * a[3] + (1 - a[3] * a[3]) * cosTheta
  m[3,4] = 0
  return(m)
}
