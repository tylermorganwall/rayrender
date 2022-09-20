#' Post-process Frame
#'
#' @return Nothing
#'
#' @keywords internal
post_process_frame = function(rgb_mat, debug_channel, filename, toneval, bloom = TRUE) {
  full_array = array(0,c(ncol(rgb_mat$r),nrow(rgb_mat$r),3))
  full_array[,,1] = flipud(t(rgb_mat$r))
  full_array[,,2] = flipud(t(rgb_mat$g))
  full_array[,,3] = flipud(t(rgb_mat$b))
  if(debug_channel == 1) {
    returnmat = full_array[,,1]
    returnmat[is.infinite(returnmat)] = NA
    save_png((full_array-min(full_array,na.rm=TRUE))/(max(full_array,na.rm=TRUE) - min(full_array,na.rm=TRUE)),
             filename)
  } else if (debug_channel %in% c(2,3,4,5,6,7,8,9)) {
    save_png(full_array,filename)
  } 
  if(bloom) {
    kernel = rayimage::generate_2d_exponential(0.1,11,3)
    full_array = rayimage::render_convolution(image = full_array, kernel = kernel, min_value = 1, preview=FALSE)
  }
  tonemapped_channels = tonemap_image(full_array[,,1],full_array[,,2],full_array[,,3],toneval)
  full_array = array(0,c(nrow(tonemapped_channels$r),ncol(tonemapped_channels$r),3))
  full_array[,,1] = tonemapped_channels$r
  full_array[,,2] = tonemapped_channels$g
  full_array[,,3] = tonemapped_channels$b
  
  array_from_mat = array(full_array,dim=c(nrow(full_array),ncol(full_array),3))
  if(any(is.na(array_from_mat ))) {
    array_from_mat[is.na(array_from_mat)] = 0
  }
  if(any(array_from_mat > 1 | array_from_mat < 0,na.rm = TRUE)) {
    array_from_mat[array_from_mat > 1] = 1
    array_from_mat[array_from_mat < 0] = 0
  }
  if(filename != "") {
    save_png(array_from_mat,filename)
  } else {
    rayimage::plot_image(array_from_mat, new_page = TRUE)
  }
}