#' Screen-space Lines
#'
#' Creates line annotations for `render_scene()` that are anchored to 3D
#' world-space endpoints but drawn in 2D screen space after rendering.
#'
#' @param x,y,z Default `0`. World-space coordinates of the line start point.
#' @param xend,yend,zend Default `0`. World-space coordinates of the line end
#' point.
#' @param start,end Default `NULL`. Optional 3-column matrix/data frame of
#' world-space start and end points. If supplied, overrides `x`, `y`, `z` and
#' `xend`, `yend`, `zend`.
#' @param offset Default `c(0,0)`. Pixel offset for the start point, as
#' `c(x,y)`, with positive y moving down the image.
#' @param end_offset Default `offset`. Pixel offset for the end point.
#' @param width Default `2`. Line width in pixels.
#' @param color Default `"black"`. Line color.
#' @param alpha Default `1`. Line alpha.
#' @param lineend Default `"round"`. Line end style. Options are `"round"`,
#' `"butt"`, and `"square"`.
#' @param clip Default `TRUE`. If `TRUE`, lines whose endpoints are both behind
#' the camera are skipped.
#' @param occlusion Default `FALSE`. If `TRUE`, use scene geometry to hide the
#' line.
#' @param occlusion_mode Default `"anchor"`. If `"anchor"`, occlusion skips the
#' entire line when its midpoint is blocked. If `"line"` or `"partial"`, each
#' line pixel is hidden only when the scene depth at that pixel is closer than
#' the interpolated screen-space line depth.
#' @param occlusion_tolerance Default `0.001`. Endpoint tolerance for occlusion.
#'
#' @return A data frame describing screen-space line annotations.
#' @export
screen_line = function(
  x = 0,
  y = 0,
  z = 0,
  xend = 0,
  yend = 0,
  zend = 0,
  start = NULL,
  end = NULL,
  offset = c(0, 0),
  end_offset = offset,
  width = 2,
  color = "black",
  alpha = 1,
  lineend = "round",
  clip = TRUE,
  occlusion = FALSE,
  occlusion_mode = "anchor",
  occlusion_tolerance = 0.001
) {
  if (!is.null(start)) {
    start = as.matrix(start)
    if (ncol(start) != 3) {
      stop("`start` must have three columns")
    }
    x = start[, 1]
    y = start[, 2]
    z = start[, 3]
  }
  if (!is.null(end)) {
    end = as.matrix(end)
    if (ncol(end) != 3) {
      stop("`end` must have three columns")
    }
    xend = end[, 1]
    yend = end[, 2]
    zend = end[, 3]
  }
  n = max(
    length(x),
    length(y),
    length(z),
    length(xend),
    length(yend),
    length(zend),
    length(width),
    length(color),
    length(alpha),
    length(lineend),
    length(clip),
    length(occlusion),
    length(occlusion_mode),
    length(occlusion_tolerance)
  )
  recycle_value = function(value, name) {
    if (length(value) == n) {
      return(value)
    }
    if (length(value) == 1) {
      return(rep(value, n))
    }
    stop("`", name, "` must have length 1 or ", n)
  }
  normalize_offset = function(value, name) {
    if (is.null(dim(value)) && length(value) == 2) {
      value = matrix(rep(as.numeric(value), n), ncol = 2, byrow = TRUE)
    } else {
      value = as.matrix(value)
    }
    if (ncol(value) != 2 || !(nrow(value) %in% c(1, n))) {
      stop("`", name, "` must be a length-2 vector or an n x 2 matrix")
    }
    if (nrow(value) == 1) {
      value = matrix(rep(as.numeric(value), n), ncol = 2, byrow = TRUE)
    }
    value
  }
  offset = normalize_offset(offset, "offset")
  end_offset = normalize_offset(end_offset, "end_offset")

  output = data.frame(
    x = as.numeric(recycle_value(x, "x")),
    y = as.numeric(recycle_value(y, "y")),
    z = as.numeric(recycle_value(z, "z")),
    xend = as.numeric(recycle_value(xend, "xend")),
    yend = as.numeric(recycle_value(yend, "yend")),
    zend = as.numeric(recycle_value(zend, "zend")),
    x_offset = as.numeric(offset[, 1]),
    y_offset = as.numeric(offset[, 2]),
    xend_offset = as.numeric(end_offset[, 1]),
    yend_offset = as.numeric(end_offset[, 2]),
    width = as.numeric(recycle_value(width, "width")),
    color = as.character(recycle_value(color, "color")),
    alpha = as.numeric(recycle_value(alpha, "alpha")),
    lineend = tolower(as.character(recycle_value(lineend, "lineend"))),
    clip = as.logical(recycle_value(clip, "clip")),
    occlusion = as.logical(recycle_value(occlusion, "occlusion")),
    occlusion_mode = tolower(as.character(recycle_value(
      occlusion_mode,
      "occlusion_mode"
    ))),
    occlusion_tolerance = as.numeric(recycle_value(
      occlusion_tolerance,
      "occlusion_tolerance"
    )),
    stringsAsFactors = FALSE
  )
  output$occlusion_mode[output$occlusion_mode == "partial"] = "line"
  if (!all(output$occlusion_mode %in% c("anchor", "line"))) {
    stop('`occlusion_mode` must be "anchor", "line", or "partial"')
  }
  if (!all(output$lineend %in% c("round", "butt", "square"))) {
    stop('`lineend` must be "round", "butt", or "square"')
  }
  output$width = pmax(output$width, .Machine$double.eps)
  output$alpha = pmin(pmax(output$alpha, 0), 1)
  class(output) = c("ray_screen_line", class(output))
  output
}

normalize_screen_line = function(screen_line_spec) {
  if (is.null(screen_line_spec)) {
    return(NULL)
  }
  if (is.data.frame(screen_line_spec)) {
    output = screen_line_spec
  } else if (
    is.list(screen_line_spec) &&
      length(screen_line_spec) > 0 &&
      all(vapply(screen_line_spec, is.data.frame, logical(1)))
  ) {
    output = do.call(rbind, screen_line_spec)
  } else if (is.list(screen_line_spec)) {
    output = do.call(screen_line, screen_line_spec)
  } else {
    stop("`screen_line` must be created by `screen_line()` or be a data frame")
  }
  required = c("x", "y", "z", "xend", "yend", "zend")
  missing_required = setdiff(required, names(output))
  if (length(missing_required) > 0) {
    stop(
      "`screen_line` is missing required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  }
  defaults = list(
    x_offset = 0,
    y_offset = 0,
    xend_offset = 0,
    yend_offset = 0,
    width = 2,
    color = "black",
    alpha = 1,
    lineend = "round",
    clip = TRUE,
    occlusion = FALSE,
    occlusion_mode = "anchor",
    occlusion_tolerance = 0.001
  )
  for (name in names(defaults)) {
    if (is.null(output[[name]])) {
      output[[name]] = defaults[[name]]
    }
  }
  output$occlusion_mode = tolower(output$occlusion_mode)
  output$occlusion_mode[output$occlusion_mode == "partial"] = "line"
  output
}

add_screen_line = function(
  image_array,
  screen_line_spec,
  camera_info,
  screen_line_visible = NULL
) {
  screen_line_spec = normalize_screen_line(screen_line_spec)
  if (is.null(screen_line_spec) || nrow(screen_line_spec) == 0) {
    return(image_array)
  }
  if (!is.null(screen_line_visible)) {
    if (length(screen_line_visible) != nrow(screen_line_spec)) {
      stop("`screen_line_visible` must match the number of screen lines")
    }
    screen_line_spec = screen_line_spec[screen_line_visible, , drop = FALSE]
    if (nrow(screen_line_spec) == 0) {
      return(image_array)
    }
  }
  start_projection = project_points_to_screen(
    screen_line_spec[, c("x", "y", "z"), drop = FALSE],
    camera_info
  )
  end_projection = project_points_to_screen(
    screen_line_spec[, c("xend", "yend", "zend"), drop = FALSE],
    camera_info
  )
  for (i in seq_len(nrow(screen_line_spec))) {
    if (!start_projection$in_front[i] || !end_projection$in_front[i]) {
      next
    }
    line_image = rasterize_screen_line(
      x0 = start_projection$x[i] + screen_line_spec$x_offset[i],
      y0 = start_projection$y[i] + screen_line_spec$y_offset[i],
      x1 = end_projection$x[i] + screen_line_spec$xend_offset[i],
      y1 = end_projection$y[i] + screen_line_spec$yend_offset[i],
      width = screen_line_spec$width[i],
      color = screen_line_spec$color[i],
      alpha = screen_line_spec$alpha[i],
      lineend = screen_line_spec$lineend[i]
    )
    image_array = composite_screen_text_image(
      image_array,
      line_image,
      x = attr(line_image, "left"),
      y = attr(line_image, "top"),
      hjust = 0,
      vjust = 0
    )
  }
  image_array
}

rasterize_screen_line = function(x0, y0, x1, y1, width, color, alpha, lineend) {
  radius = width / 2
  pad = ceiling(radius + 1)
  left = floor(min(x0, x1) - pad)
  right = ceiling(max(x0, x1) + pad)
  top = floor(min(y0, y1) - pad)
  bottom = ceiling(max(y0, y1) + pad)
  xs = seq(left, right)
  ys = seq(top, bottom)
  grid = expand.grid(x = xs, y = ys)
  dx = x1 - x0
  dy = y1 - y0
  len2 = dx * dx + dy * dy
  if (len2 <= .Machine$double.eps) {
    dist = sqrt((grid$x - x0)^2 + (grid$y - y0)^2)
  } else {
    if (lineend == "square") {
      len = sqrt(len2)
      ux = dx / len
      uy = dy / len
      x0 = x0 - ux * radius
      y0 = y0 - uy * radius
      x1 = x1 + ux * radius
      y1 = y1 + uy * radius
      dx = x1 - x0
      dy = y1 - y0
      len2 = dx * dx + dy * dy
    }
    t = ((grid$x - x0) * dx + (grid$y - y0) * dy) / len2
    if (lineend == "butt") {
      inside = t >= 0 & t <= 1
    } else {
      inside = rep(TRUE, length(t))
    }
    t = pmin(pmax(t, 0), 1)
    closest_x = x0 + t * dx
    closest_y = y0 + t * dy
    dist = sqrt((grid$x - closest_x)^2 + (grid$y - closest_y)^2)
    dist[!inside] = Inf
  }
  coverage = pmin(pmax(radius + 0.5 - dist, 0), 1)
  line_color = as.vector(grDevices::col2rgb(color)) / 255
  line_image = array(0, c(length(ys), length(xs), 4))
  alpha_channel = matrix(coverage * alpha, nrow = length(xs), byrow = FALSE)
  alpha_channel = t(alpha_channel)
  for (channel in 1:3) {
    line_image[,, channel] = line_color[channel]
  }
  line_image[,, 4] = alpha_channel
  attr(line_image, "left") = left
  attr(line_image, "top") = top
  line_image
}

prepare_screen_line_preview = function(screen_line_spec) {
  screen_line_spec = normalize_screen_line(screen_line_spec)
  if (is.null(screen_line_spec) || nrow(screen_line_spec) == 0) {
    return(list(active = FALSE))
  }
  lines = lapply(seq_len(nrow(screen_line_spec)), function(i) {
    row = screen_line_spec[i, , drop = FALSE]
    anchor_occlusion = isTRUE(row$occlusion) && row$occlusion_mode == "anchor"
    partial_occlusion = isTRUE(row$occlusion) && row$occlusion_mode == "line"
    rgb = as.vector(grDevices::col2rgb(row$color)) / 255
    list(
      x = row$x,
      y = row$y,
      z = row$z,
      xend = row$xend,
      yend = row$yend,
      zend = row$zend,
      x_offset = row$x_offset,
      y_offset = row$y_offset,
      xend_offset = row$xend_offset,
      yend_offset = row$yend_offset,
      width = row$width,
      red = rgb[1],
      green = rgb[2],
      blue = rgb[3],
      alpha = row$alpha,
      lineend = row$lineend,
      clip = row$clip,
      occlusion = anchor_occlusion,
      partial_occlusion = partial_occlusion,
      occlusion_tolerance = row$occlusion_tolerance
    )
  })
  list(active = length(lines) > 0, lines = lines)
}

screen_line_needs_native_overlay = function(screen_line_spec) {
  screen_line_spec = normalize_screen_line(screen_line_spec)
  if (is.null(screen_line_spec) || nrow(screen_line_spec) == 0) {
    return(FALSE)
  }
  any(screen_line_spec$occlusion & screen_line_spec$occlusion_mode == "line")
}

prepare_screen_line_occlusion = function(screen_line_spec) {
  screen_line_spec = normalize_screen_line(screen_line_spec)
  if (is.null(screen_line_spec) || nrow(screen_line_spec) == 0) {
    return(list(active = FALSE))
  }
  anchor_occlusion = screen_line_spec$occlusion &
    screen_line_spec$occlusion_mode == "anchor"
  list(
    active = any(anchor_occlusion),
    x = (screen_line_spec$x + screen_line_spec$xend) / 2,
    y = (screen_line_spec$y + screen_line_spec$yend) / 2,
    z = (screen_line_spec$z + screen_line_spec$zend) / 2,
    occlusion = anchor_occlusion,
    occlusion_tolerance = screen_line_spec$occlusion_tolerance
  )
}
