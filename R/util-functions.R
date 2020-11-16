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