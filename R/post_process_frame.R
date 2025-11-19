#' Post-process Frame
#'
#' @return Nothing
#'
#' @keywords internal
post_process_frame = function(
	rgb_mat,
	debug_channel,
	filename,
	tonemap,
	bloom = TRUE,
	transparent_background = FALSE,
	write_file = TRUE
) {
	if (!transparent_background) {
		full_array = array(0, c(ncol(rgb_mat$r), nrow(rgb_mat$r), 3))
	} else {
		full_array = array(0, c(ncol(rgb_mat$r), nrow(rgb_mat$r), 4))
	}
	full_array[,, 1] = flipud(t(rgb_mat$r))
	full_array[,, 2] = flipud(t(rgb_mat$g))
	full_array[,, 3] = flipud(t(rgb_mat$b))
	if (transparent_background) {
		full_array[,, 4] = flipud(t(rgb_mat$a))
	}
	if (debug_channel == 1) {
		returnmat = full_array[,, 1]
		returnmat[is.infinite(returnmat)] = NA
		rayimage::ray_write_image(
			(full_array - min(full_array, na.rm = TRUE)) /
				(max(full_array, na.rm = TRUE) - min(full_array, na.rm = TRUE)),
			filename
		)
	} else if (debug_channel %in% c(2, 3, 4, 5, 6, 7, 8, 9)) {
		rayimage::ray_write_image(full_array, filename)
	}
	if (bloom) {
		kernel = rayimage::generate_2d_exponential(0.1, 11, 3)
		full_array = rayimage::render_convolution(
			image = full_array,
			kernel = kernel,
			min_value = 1,
			preview = FALSE
		)
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
	if (write_file) {
		rayimage::ray_write_image(full_array, filename)
	} else {
		rayimage::plot_image(full_array, new_page = TRUE)
	}
}
