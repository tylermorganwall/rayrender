#'@title Post-process the scene
#'
#'@keywords internal
#'@examples
#'#internal
post_process_scene = function(rgb_mat, iso, tonemap, debug_channel, filename, return_raw_array, bloom) {
  toneval = switch(tonemap, "gamma" = 1,"reinhold" = 2,"uncharted" = 3,"hbd" = 4, "raw" = 5)
  
  if(!is.numeric(debug_channel)) {
    debug_channel = unlist(lapply(tolower(debug_channel),switch,
                                  "none" = 0,"depth" = 1,"normals" = 2, "uv" = 3, "bvh" = 4,
                                  "variance" = 5, "normal" = 2, "dpdu" = 6, "dpdv" = 7, "color" = 8, 
                                  "position" = 10, "direction" = 11, "time" = 12, "shape" = 13,
                                  "pdf" = 14, "error" = 15, "bounces" = 16, "camera" = 17,
                                  0))
    light_direction = c(0,1,0)
  } else {
    light_direction = debug_channel
    debug_channel = 9
  }
  
  full_array = array(0,c(ncol(rgb_mat$r),nrow(rgb_mat$r),3))
  full_array[,,1] = flipud(t(rgb_mat$r))
  full_array[,,2] = flipud(t(rgb_mat$g))
  full_array[,,3] = flipud(t(rgb_mat$b))
  if(debug_channel == 1) {
    returnmat = fliplr(t(full_array[,,1]))
    returnmat[is.infinite(returnmat)] = NA
    if(is.null(filename)) {
      plot_map((returnmat-min(returnmat,na.rm=TRUE))/(max(returnmat,na.rm=TRUE) - min(returnmat,na.rm=TRUE)))
      return(invisible(returnmat))
    } else {
      save_png((returnmat-min(returnmat,na.rm=TRUE))/(max(returnmat,na.rm=TRUE) - min(returnmat,na.rm=TRUE)),
               filename)
      return(invisible(returnmat))
    }
  } else if (debug_channel %in% c(2,3,4,5,17)) {
    if(is.null(filename)) {
      if(!return_raw_array) {
        if(debug_channel == 4) {
          plot_map(full_array/(max(full_array,na.rm=TRUE)))
        } else {
          plot_map(full_array)
        }
      }
      return(invisible(full_array))
    } else {
      save_png(full_array,filename)
      return(invisible(full_array))
    }
  } else if (debug_channel %in% c(10,13)) {
    full_array_ret = full_array
    full_array[,,1][is.infinite(full_array[,,1])] = max(full_array[,,1][!is.infinite(full_array[,,1])])
    full_array[,,2][is.infinite(full_array[,,2])] = max(full_array[,,2][!is.infinite(full_array[,,2])])
    full_array[,,3][is.infinite(full_array[,,3])] = max(full_array[,,3][!is.infinite(full_array[,,3])])
    
    full_array[,,1] = (full_array[,,1] - min(full_array[,,1]))/(max(full_array[,,1]) - min(full_array[,,1]))
    full_array[,,2] = (full_array[,,2] - min(full_array[,,2]))/(max(full_array[,,2]) - min(full_array[,,2]))
    full_array[,,3] = (full_array[,,3] - min(full_array[,,3]))/(max(full_array[,,3]) - min(full_array[,,3]))
    if(is.null(filename)) {
      plot_map(full_array)
    } else {
      save_png(full_array,filename)
    }
    return(invisible(full_array_ret))
  } else if (debug_channel == 11) {
    full_array_ret = full_array
    
    full_array[,,1] = (full_array[,,1]+1)/2
    full_array[,,2] = (full_array[,,2]+1)/2
    full_array[,,3] = (full_array[,,3]+1)/2
    if(is.null(filename)) {
      plot_map(full_array)
    } else {
      save_png(full_array,filename)
    }
    return(invisible(full_array_ret))
  } else if (debug_channel %in% c(12,14,15,16)) {
    full_array_ret = full_array
    full_array[is.infinite(full_array)] = max(full_array[!is.infinite(full_array)])
    
    full_array = (full_array - min(full_array))/(max(full_array) - min(full_array))
    if(is.null(filename)) {
      plot_map(full_array)
    } else {
      save_png(full_array,filename)
    }
    return(invisible(full_array_ret))
  }
  if(!is.matrix(bloom)) {
    if(is.numeric(bloom) && length(bloom) == 1) {
      kernel = rayimage::generate_2d_exponential(0.1,11,3*1/bloom)
      full_array = rayimage::render_convolution(image = full_array, kernel = kernel, min_value = 1, preview=FALSE)
    } else {
      if(bloom) {
        kernel = rayimage::generate_2d_exponential(0.1,11,3)
        full_array = rayimage::render_convolution(image = full_array, kernel = kernel, min_value = 1, preview=FALSE)
      }
    }
  } else {
    kernel = bloom
    if(ncol(kernel) %% 2 == 0) {
      newkernel = matrix(0, ncol = ncol(kernel) + 1, nrow = nrow(kernel))
      newkernel[,1:ncol(kernel)] = kernel
      kernel = newkernel
    }
    if(nrow(kernel) %% 2 == 0) {
      newkernel = matrix(0, ncol = ncol(kernel), nrow = nrow(kernel) + 1)
      newkernel[1:nrow(kernel),] = kernel
      kernel = newkernel
    }
    full_array = rayimage::render_convolution(image = full_array, kernel = kernel,  min_value = 1, preview=FALSE)
  }
  tonemapped_channels = tonemap_image(full_array[,,1],full_array[,,2],full_array[,,3],toneval)
  full_array = array(0,c(nrow(tonemapped_channels$r),ncol(tonemapped_channels$r),3))
  full_array[,,1] = tonemapped_channels$r
  full_array[,,2] = tonemapped_channels$g
  full_array[,,3] = tonemapped_channels$b
  if(toneval == 5) {
    return(full_array)
  }
  
  array_from_mat = array(full_array,dim=c(nrow(full_array),ncol(full_array),3))
  if(any(is.na(array_from_mat ))) {
    array_from_mat[is.na(array_from_mat)] = 0
  }
  if(any(array_from_mat > 1 | array_from_mat < 0,na.rm = TRUE)) {
    array_from_mat[array_from_mat > 1] = 1
    array_from_mat[array_from_mat < 0] = 0
  }
  if(is.null(filename)) {
    if(!return_raw_array) {
      plot_map(array_from_mat)
    }
  } else {
    save_png(array_from_mat,filename)
  }
  return(invisible(array_from_mat))
}