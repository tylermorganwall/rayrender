#' Title
#'
#' @return DUM
#' @export
#' @import raster
#'
#' @examples
#' #fake example
runcode = function(width = 200, height = 100, aliasing = 100) {
  rgb_mat = generate_initial(nx = width, ny = height, ns = aliasing)
  full_array = array(0,c(ncol(rgb_mat$r),nrow(rgb_mat$r),3))
  full_array[,,1] = t(rgb_mat$r)
  full_array[,,2] = t(rgb_mat$g)
  full_array[,,3] = t(rgb_mat$b)
  array_from_mat = array(full_array,dim=c(nrow(full_array),ncol(full_array),3))
  suppressWarnings(raster::plotRGB(raster::brick(flipud(array_from_mat), xmn = 0.5, xmx = dim(array_from_mat)[2] + 0.5,ymn =  0.5, ymx = dim(array_from_mat)[1] +  0.5),scale=1, maxpixels=nrow(rgb_mat$r)*ncol(rgb_mat$r)))
  
}