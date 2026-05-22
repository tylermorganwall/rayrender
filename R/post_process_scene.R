#'@title Post-process the scene
#'
#'@keywords internal
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#'#internal
post_process_scene = function(
  rgb_mat,
  iso,
  use_iso,
  tonemap,
  debug_channel,
  filename,
  plot_scene,
  bloom,
  new_page = TRUE,
  transparent_background = FALSE,
  auto_exposure = FALSE,
  verbose = FALSE,
  screen_text = NULL,
  screen_line = NULL,
  camera_info = NULL,
  screen_text_visible = NULL,
  screen_text_overlay = NULL,
  screen_line_visible = NULL,
  screen_line_overlay = NULL,
  exposure_adjustment = 1
) {
  if (!is.numeric(debug_channel)) {
    debug_channel = unlist(lapply(
      tolower(debug_channel),
      switch,
      "none" = 0,
      "depth" = 1,
      "normals" = 2,
      "uv" = 3,
      "bvh" = 4,
      "variance" = 5,
      "normal" = 2,
      "dpdu" = 6,
      "dpdv" = 7,
      "color" = 8,
      "position" = 10,
      "direction" = 11,
      "time" = 12,
      "shape" = 13,
      "pdf" = 14,
      "error" = 15,
      "bounces" = 16,
      "camera" = 17,
      0
    ))
    light_direction = c(0, 1, 0)
  } else {
    light_direction = debug_channel
    debug_channel = 9
  }
  if (!transparent_background) {
    full_array = array(0, c(ncol(rgb_mat$r), nrow(rgb_mat$r), 3))
  } else {
    full_array = array(0, c(ncol(rgb_mat$r), nrow(rgb_mat$r), 4))
  }
  if (!use_iso || debug_channel != 0) {
    iso = 1
  }
  full_array[,, 1] = fliplr(flipud(t(rgb_mat$r))) * iso
  full_array[,, 2] = fliplr(flipud(t(rgb_mat$g))) * iso
  full_array[,, 3] = fliplr(flipud(t(rgb_mat$b))) * iso
  if (transparent_background) {
    full_array[,, 4] = fliplr(flipud(t(rgb_mat$a)))
  }
  if (debug_channel == 1) {
    # returnmat = fliplr(t(full_array[,,1]))
    returnmat = full_array[,, 1]
    returnmat[is.infinite(returnmat)] = NA
    if (is.na(filename)) {
      if (plot_scene) {
        rayimage::plot_image(
          (returnmat - min(returnmat, na.rm = TRUE)) /
            (max(returnmat, na.rm = TRUE) - min(returnmat, na.rm = TRUE)),
          new_page = new_page
        )
      }
      return(invisible(returnmat))
    } else {
      rayimage::ray_write_image(
        (returnmat - min(returnmat, na.rm = TRUE)) /
          (max(returnmat, na.rm = TRUE) - min(returnmat, na.rm = TRUE)),
        filename
      )
      return(invisible(returnmat))
    }
  } else if (debug_channel %in% c(2, 3, 4, 5, 17)) {
    if (is.na(filename)) {
      if (plot_scene) {
        if (debug_channel == 4) {
          rayimage::plot_image(
            full_array / (max(full_array, na.rm = TRUE)),
            new_page = new_page
          )
        } else {
          rayimage::plot_image(
            full_array,
            new_page = new_page
          )
        }
      }
      return(invisible(full_array))
    } else {
      rayimage::ray_write_image(full_array, filename)
      return(invisible(full_array))
    }
  } else if (debug_channel %in% c(10, 13)) {
    full_array_ret = full_array
    full_array[,, 1][is.infinite(full_array[,, 1])] = max(full_array[,, 1][
      !is.infinite(full_array[,, 1])
    ])
    full_array[,, 2][is.infinite(full_array[,, 2])] = max(full_array[,, 2][
      !is.infinite(full_array[,, 2])
    ])
    full_array[,, 3][is.infinite(full_array[,, 3])] = max(full_array[,, 3][
      !is.infinite(full_array[,, 3])
    ])

    full_array[,, 1] = (full_array[,, 1] - min(full_array[,, 1])) /
      (max(full_array[,, 1]) - min(full_array[,, 1]))
    full_array[,, 2] = (full_array[,, 2] - min(full_array[,, 2])) /
      (max(full_array[,, 2]) - min(full_array[,, 2]))
    full_array[,, 3] = (full_array[,, 3] - min(full_array[,, 3])) /
      (max(full_array[,, 3]) - min(full_array[,, 3]))
    if (is.na(filename)) {
      if (plot_scene) {
        rayimage::plot_image(full_array, new_page = new_page)
      }
    } else {
      rayimage::ray_write_image(full_array, filename)
    }
    return(invisible(full_array_ret))
  } else if (debug_channel == 11) {
    full_array_ret = full_array

    full_array[,, 1] = (full_array[,, 1] + 1) / 2
    full_array[,, 2] = (full_array[,, 2] + 1) / 2
    full_array[,, 3] = (full_array[,, 3] + 1) / 2
    if (is.na(filename)) {
      if (plot_scene) {
        rayimage::plot_image(full_array, new_page = new_page)
      }
    } else {
      rayimage::ray_write_image(full_array, filename)
    }
    return(invisible(full_array_ret))
  } else if (debug_channel %in% c(12, 14, 15, 16)) {
    full_array_ret = full_array
    full_array[is.infinite(full_array)] = max(full_array[
      !is.infinite(full_array)
    ])

    full_array = (full_array - min(full_array)) /
      (max(full_array) - min(full_array))
    if (is.na(filename)) {
      if (plot_scene) {
        rayimage::plot_image(full_array, new_page = new_page)
      }
    } else {
      rayimage::ray_write_image(full_array, filename)
    }
    return(invisible(full_array_ret))
  }
  if (debug_channel == 0) {
    if (!is.matrix(bloom)) {
      if (is.numeric(bloom) && length(bloom) == 1) {
        kernel = rayimage::generate_2d_exponential(0.1, 11, 3 * 1 / bloom)
        full_array = rayimage::render_convolution(
          image = full_array,
          kernel = kernel,
          min_value = 1,
          preview = FALSE
        )
      } else {
        if (bloom) {
          kernel = rayimage::generate_2d_exponential(0.1, 11, 3)
          full_array = rayimage::render_convolution(
            image = full_array,
            kernel = kernel,
            min_value = 1,
            preview = FALSE
          )
        }
      }
    } else {
      kernel = bloom
      if (ncol(kernel) %% 2 == 0) {
        newkernel = matrix(0, ncol = ncol(kernel) + 1, nrow = nrow(kernel))
        newkernel[, 1:ncol(kernel)] = kernel
        kernel = newkernel
      }
      if (nrow(kernel) %% 2 == 0) {
        newkernel = matrix(0, ncol = ncol(kernel), nrow = nrow(kernel) + 1)
        newkernel[1:nrow(kernel), ] = kernel
        kernel = newkernel
      }
      full_array = rayimage::render_convolution(
        image = full_array,
        kernel = kernel,
        min_value = 1,
        preview = FALSE
      )
    }
  }
  if (isTRUE(auto_exposure)) {
    full_array = rayimage::render_exposure(
      full_array,
      auto = TRUE,
      preview = FALSE,
      verbose = verbose
    )
  }
  if (
    is.numeric(exposure_adjustment) &&
      length(exposure_adjustment) == 1 &&
      is.finite(exposure_adjustment)
  ) {
    full_array[,, 1:3] = full_array[,, 1:3] * exposure_adjustment
  }
  full_array = full_array |>
    rayimage::ray_read_image(
      convert_to_array = TRUE,
      assume_colorspace = rayimage::CS_SRGB,
      assume_white = "D65"
    ) |>
    rayimage::render_tonemap(method = tonemap)

  if (any(is.na(full_array))) {
    full_array[is.na(full_array)] = 0
  }
  if (!is.null(screen_text_overlay)) {
    full_array = composite_screen_text_image(
      full_array,
      screen_text_overlay,
      x = 1,
      y = 1,
      hjust = 0,
      vjust = 0
    )
  }
  if (!is.null(screen_line_overlay)) {
    full_array = composite_screen_text_image(
      full_array,
      screen_line_overlay,
      x = 1,
      y = 1,
      hjust = 0,
      vjust = 0
    )
  }
  if (!is.null(screen_line)) {
    screen_line_spec = normalize_screen_line(screen_line)
    if (!is.null(screen_line_overlay)) {
      draw_screen_line = !(screen_line_spec$occlusion &
        screen_line_spec$occlusion_mode == "line")
      if (!is.null(screen_line_visible)) {
        screen_line_visible = screen_line_visible[draw_screen_line]
      }
      screen_line_spec = screen_line_spec[
        draw_screen_line,
        ,
        drop = FALSE
      ]
    }
    full_array = add_screen_line(
      full_array,
      screen_line_spec,
      camera_info,
      screen_line_visible
    )
  }
  if (!is.null(screen_text)) {
    screen_text_spec = normalize_screen_text(screen_text)
    if (!is.null(screen_text_overlay)) {
      draw_screen_text = !(screen_text_spec$occlusion &
        screen_text_spec$occlusion_mode == "label")
      if (!is.null(screen_text_visible)) {
        screen_text_visible = screen_text_visible[draw_screen_text]
      }
      screen_text_spec = screen_text_spec[
        draw_screen_text,
        ,
        drop = FALSE
      ]
    }
    full_array = add_screen_text(
      full_array,
      screen_text_spec,
      camera_info,
      screen_text_visible
    )
  }
  if (is.na(filename)) {
    if (plot_scene) {
      rayimage::plot_image(
        full_array,
        new_page = new_page
      )
    } else {
      return(full_array)
    }
  } else {
    rayimage::ray_write_image(full_array, filename)
    if (plot_scene) {
      rayimage::plot_image(
        full_array,
        new_page = new_page
      )
    }
  }
  return(invisible(full_array))
}
