#' Screen-space Text
#'
#' Creates text annotations for `render_scene()` that are anchored to 3D
#' world-space points but drawn in 2D screen space after rendering.
#'
#' @param label Character vector of labels.
#' @param x Default `0`. World-space x-coordinate of the anchor point.
#' @param y Default `0`. World-space y-coordinate of the anchor point.
#' @param z Default `0`. World-space z-coordinate of the anchor point.
#' @param point Default `NULL`. Optional 3-column matrix/data frame of
#' world-space anchor points. If supplied, overrides `x`, `y`, and `z`.
#' @param offset Default `c(0,0)`. Pixel offset from the projected anchor point,
#' as `c(x,y)`, with positive y moving down the image.
#' @param hjust Default `0`. Horizontal text adjustment, where `0` places the
#' left edge at the anchor, `0.5` centers the label, and `1` right-aligns it.
#' @param vjust Default `0.5`. Vertical text adjustment, where `0` places the
#' top edge at the anchor, `0.5` centers the label, and `1` bottom-aligns it.
#' @param size Default `16`. Text size in pixels.
#' @param color Default `"black"`. Text color.
#' @param font Default `"sans"`. Font family.
#' @param lineheight Default `1`. Line height passed to
#' `rayimage::render_text_image()`.
#' @param background_color Default `"white"`. Background color passed to
#' `rayimage::render_text_image()`.
#' @param background_alpha Default `0`. Background alpha passed to
#' `rayimage::render_text_image()`.
#' @param just Default `"left"`. Text justification passed to
#' `rayimage::render_text_image()`.
#' @param clip Default `TRUE`. If `TRUE`, labels whose anchor point projects
#' outside the image are skipped.
#' @param halo_color Default `NA`, no halo. If a color is specified, the text
#' label will be surrounded by a halo of this color.
#' @param halo_expand Default `0`. Number of pixels to expand the halo.
#' @param halo_alpha Default `1`. Transparency of the halo.
#' @param halo_offset Default `c(0,0)`. Pixel offset to apply to the halo,
#' as `c(x,y)`, with positive y moving down the image.
#' @param halo_blur Default `0`. Amount of blur to apply to the halo.
#' @param halo_edge_softness Default `0.1`. Width of the softened halo edge
#' transition, in pixels.
#' @param halo_gap_fill Default `2`. Maximum alpha gap width, in pixels, to
#' bridge in the halo outline.
#' @param halo_gap_fill_alpha_threshold Default `0.25`. Alpha threshold used to
#' protect enclosed interior halo gaps from `halo_gap_fill`.
#' @param occlusion Default `FALSE`. If `TRUE`, trace a ray from the camera to
#' the anchor point and skip the label when another scene object blocks it.
#' @param occlusion_mode Default `"anchor"`. If `"anchor"`, occlusion skips the
#' entire label when the anchor point is blocked. If `"label"` or `"partial"`,
#' the label and halo are treated as a screen-space banner at the anchor depth
#' and each text pixel is hidden only when the scene depth at that pixel is
#' closer.
#' @param occlusion_tolerance Default `0.001`. Endpoint tolerance for the
#' occlusion ray. Values below `1` are treated as a fraction of the
#' camera-to-anchor distance, which helps avoid self-occlusion when the anchor
#' lies on a surface. Values of `1` or greater are treated as scene-unit
#' distances.
#'
#' @return A data frame describing screen-space text annotations.
#' @export
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' # Start with a Cornell box that has one bright object and one dark object.
#' scene = generate_cornell(lightwidth = 250, lightdepth = 250) |>
#'   add_object(sphere(
#'     x = 180, y = 90, z = 260, radius = 90,
#'     material = diffuse(color = "#f2c14e")
#'   )) |>
#'   add_object(cube(
#'     x = 385, y = 90, z = 330, xwidth = 120, ywidth = 180, zwidth = 120,
#'     angle = c(0, 25, 0), material = diffuse(color = "#243846")
#'   ))
#'
#' # Labels can be built in separate calls and passed as a list to render_scene().
#' sphere_label = screen_text(
#'   label = "matte sphere\nanchor test",
#'   x = 180, y = 190, z = 260,
#'   offset = c(18, -28),
#'   hjust = 0.5,
#'   vjust = 1,
#'   size = 18,
#'   color = "black",
#'   lineheight = 0.95,
#'   background_alpha = 0,
#'   just = "left",
#'   halo_color = "#f2c14e",
#'   halo_expand = 4,
#'   halo_alpha = 0.8,
#'   halo_offset = c(1, 1),
#'   halo_blur = 0.5,
#'   halo_edge_softness = 0.5,
#'   halo_gap_fill = 1,
#'   halo_gap_fill_alpha_threshold = 0.2,
#'   occlusion = TRUE,
#'   occlusion_mode = "anchor",
#'   occlusion_tolerance = 0.002
#' )
#'
#' # Dark text over a dark object remains legible with a subtle white halo.
#' box_label = screen_text(
#'   label = "dark text\nwhite halo",
#'   x = 385, y = 195, z = 20,
#'   offset = c(-16, -34),
#'   hjust = 0.5,
#'   vjust = 1,
#'   size = 16,
#'   color = "#111111",
#'   lineheight = 1.05,
#'   background_alpha = 0,
#'   just = "right",
#'   halo_color = "white",
#'   halo_expand = 6,
#'   halo_alpha = 0.72,
#'   halo_offset = c(-1, 2),
#'   halo_blur = 1.5,
#'   halo_edge_softness = 10,
#'   halo_gap_fill = 3,
#'   halo_gap_fill_alpha_threshold = 0.35,
#'   occlusion = TRUE,
#'   occlusion_mode = "label",
#'   occlusion_tolerance = 0.01
#' )
#'
#' # A floor label shows multi-line text, custom justification, and partial occlusion.
#' floor_label = screen_text(
#'   label = "floor label\n\nocclusion",
#'   point = matrix(c(280, 0, 50), ncol = 3),
#'   offset = c(0, 18),
#'   hjust = 0.5,
#'   vjust = 0.7,
#'   size = 15,
#'   color = "#222222",
#'   lineheight = 1,
#'   background_color = "#fff7cc",
#'   background_alpha = 0.65,
#'   just = "center",
#'   clip = FALSE,
#'   halo_color = "white",
#'   halo_expand = 3,
#'   halo_alpha = 0.7,
#'   halo_blur = 0,
#'   halo_edge_softness = 0.25,
#'   halo_gap_fill = 2,
#'   halo_gap_fill_alpha_threshold = 0.25,
#'   occlusion = TRUE,
#'   occlusion_mode = "partial",
#'   occlusion_tolerance = 0.005
#' )
#'
#' # Later list entries draw after earlier entries.
#' wall_label = screen_text(
#'   label = "back wall\nalways in front",
#'   x = 278, y = 330, z = 545,
#'   offset = c(0, -20),
#'   hjust = 0.5,
#'   vjust = 1,
#'   size = 12,
#'   color = "white",
#'   background_color = "black",
#'   background_alpha = 0.35,
#'   halo_color = "black",
#'   halo_expand = 2,
#'   clip = TRUE
#' )
#'
#' text_layers = list(sphere_label, box_label, floor_label, wall_label)
#'
#' render_scene(
#'   scene,
#'   samples = 32,
#'   clamp_value = 5,
#'   aperture = 0,
#'   fov = 50,
#'   ambient_light = FALSE,
#'   screen_text = text_layers
#' )
screen_text = function(
  label,
  x = 0,
  y = 0,
  z = 0,
  point = NULL,
  offset = c(0, 0),
  hjust = 0,
  vjust = 0.5,
  size = 16,
  color = "black",
  font = "sans",
  lineheight = 1,
  background_color = "white",
  background_alpha = 0,
  just = "left",
  clip = TRUE,
  halo_color = NA,
  halo_expand = 0,
  halo_alpha = 1,
  halo_offset = c(0, 0),
  halo_blur = 0,
  halo_edge_softness = 0.1,
  halo_gap_fill = 2,
  halo_gap_fill_alpha_threshold = 0.25,
  occlusion = FALSE,
  occlusion_mode = "anchor",
  occlusion_tolerance = 0.001
) {
  if (!is.null(point)) {
    point = as.matrix(point)
    if (ncol(point) != 3) {
      stop("`point` must have three columns")
    }
    x = point[, 1]
    y = point[, 2]
    z = point[, 3]
  }
  n = max(
    length(label),
    length(x),
    length(y),
    length(z),
    length(hjust),
    length(vjust),
    length(size),
    length(color),
    length(font),
    length(lineheight),
    length(background_color),
    length(background_alpha),
    length(just),
    length(clip),
    length(halo_color),
    length(halo_expand),
    length(halo_alpha),
    length(halo_blur),
    length(halo_edge_softness),
    length(halo_gap_fill),
    length(halo_gap_fill_alpha_threshold),
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
  if (is.null(dim(offset)) && length(offset) == 2) {
    offset = matrix(rep(as.numeric(offset), n), ncol = 2, byrow = TRUE)
  } else {
    offset = as.matrix(offset)
  }
  if (ncol(offset) != 2 || !(nrow(offset) %in% c(1, n))) {
    stop("`offset` must be a length-2 vector or an n x 2 matrix")
  }
  if (nrow(offset) == 1) {
    offset = matrix(rep(as.numeric(offset), n), ncol = 2, byrow = TRUE)
  }
  if (is.null(dim(halo_offset)) && length(halo_offset) == 2) {
    halo_offset = matrix(
      rep(as.numeric(halo_offset), n),
      ncol = 2,
      byrow = TRUE
    )
  } else {
    halo_offset = as.matrix(halo_offset)
  }
  if (ncol(halo_offset) != 2 || !(nrow(halo_offset) %in% c(1, n))) {
    stop("`halo_offset` must be a length-2 vector or an n x 2 matrix")
  }
  if (nrow(halo_offset) == 1) {
    halo_offset = matrix(
      rep(as.numeric(halo_offset), n),
      ncol = 2,
      byrow = TRUE
    )
  }
  output = data.frame(
    label = as.character(recycle_value(label, "label")),
    x = as.numeric(recycle_value(x, "x")),
    y = as.numeric(recycle_value(y, "y")),
    z = as.numeric(recycle_value(z, "z")),
    x_offset = as.numeric(offset[, 1]),
    y_offset = as.numeric(offset[, 2]),
    hjust = as.numeric(recycle_value(hjust, "hjust")),
    vjust = as.numeric(recycle_value(vjust, "vjust")),
    size = as.numeric(recycle_value(size, "size")),
    color = as.character(recycle_value(color, "color")),
    font = as.character(recycle_value(font, "font")),
    lineheight = as.numeric(recycle_value(lineheight, "lineheight")),
    background_color = as.character(recycle_value(
      background_color,
      "background_color"
    )),
    background_alpha = as.numeric(recycle_value(
      background_alpha,
      "background_alpha"
    )),
    just = as.character(recycle_value(just, "just")),
    clip = as.logical(recycle_value(clip, "clip")),
    halo_color = as.character(recycle_value(halo_color, "halo_color")),
    halo_expand = as.numeric(recycle_value(halo_expand, "halo_expand")),
    halo_alpha = as.numeric(recycle_value(halo_alpha, "halo_alpha")),
    halo_x_offset = as.numeric(halo_offset[, 1]),
    halo_y_offset = as.numeric(halo_offset[, 2]),
    halo_blur = as.numeric(recycle_value(halo_blur, "halo_blur")),
    halo_edge_softness = as.numeric(recycle_value(
      halo_edge_softness,
      "halo_edge_softness"
    )),
    halo_gap_fill = as.numeric(recycle_value(
      halo_gap_fill,
      "halo_gap_fill"
    )),
    halo_gap_fill_alpha_threshold = as.numeric(recycle_value(
      halo_gap_fill_alpha_threshold,
      "halo_gap_fill_alpha_threshold"
    )),
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
  output$occlusion_mode[output$occlusion_mode == "partial"] = "label"
  if (!all(output$occlusion_mode %in% c("anchor", "label"))) {
    stop('`occlusion_mode` must be "anchor", "label", or "partial"')
  }
  class(output) = c("ray_screen_text", class(output))
  output
}
#' @keywords internal
normalize_screen_text = function(screen_text_spec) {
  if (is.null(screen_text_spec)) {
    return(NULL)
  }
  if (is.data.frame(screen_text_spec)) {
    output = screen_text_spec
  } else if (
    is.list(screen_text_spec) &&
      length(screen_text_spec) > 0 &&
      all(vapply(screen_text_spec, is.data.frame, logical(1)))
  ) {
    output = do.call(rbind, screen_text_spec)
  } else if (is.list(screen_text_spec)) {
    output = do.call(screen_text, screen_text_spec)
  } else {
    stop("`screen_text` must be created by `screen_text()` or be a data frame")
  }
  required = c("label", "x", "y", "z")
  missing_required = setdiff(required, names(output))
  if (length(missing_required) > 0) {
    stop(
      "`screen_text` is missing required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  }
  defaults = list(
    x_offset = 0,
    y_offset = 0,
    hjust = 0,
    vjust = 0.5,
    size = 16,
    color = "black",
    font = "sans",
    lineheight = 1,
    background_color = "white",
    background_alpha = 0,
    just = "left",
    clip = TRUE,
    halo_color = NA,
    halo_expand = 0,
    halo_alpha = 1,
    halo_x_offset = 0,
    halo_y_offset = 0,
    halo_blur = 0,
    halo_edge_softness = 0.1,
    halo_gap_fill = 2,
    halo_gap_fill_alpha_threshold = 0.25,
    occlusion = FALSE,
    occlusion_mode = "anchor",
    occlusion_tolerance = 0.001
  )
  for (name in names(defaults)) {
    if (is.null(output[[name]])) {
      output[[name]] = defaults[[name]]
    }
  }
  output
}

#' @keywords internal
project_points_to_screen = function(points, camera_info) {
  points = as.matrix(points)
  if (ncol(points) != 3) {
    stop("`points` must have three columns")
  }
  if (!is.null(camera_info$screen_camera_origin)) {
    lookfrom = camera_info$screen_camera_origin
    right = camera_info$screen_camera_u
    up = camera_info$screen_camera_v
    forward = camera_info$screen_camera_w
    fov = camera_info$screen_camera_fov
    ortho_dimensions = camera_info$screen_camera_ortho_dimensions
  } else {
    lookfrom = camera_info$lookfrom
    lookat = camera_info$lookat
    camera_up = camera_info$camera_up
    forward = lookat - lookfrom
    forward = forward / sqrt(sum(forward^2))
    up = camera_up / sqrt(sum(camera_up^2))
    right = cross_prod(up, forward)
    right = right / sqrt(sum(right^2))
    up = cross_prod(forward, right)
    up = up / sqrt(sum(up^2))
    fov = camera_info$fov
    ortho_dimensions = camera_info$ortho_dimensions
  }

  relative_points = sweep(points, 2, lookfrom)
  x_camera = as.vector(relative_points %*% right)
  y_camera = as.vector(relative_points %*% up)
  z_camera = as.vector(relative_points %*% forward)

  if (fov == 0) {
    s = 0.5 + x_camera / ortho_dimensions[1]
    t = 0.5 + y_camera / ortho_dimensions[2]
    in_front = z_camera >= 0
  } else if (fov == 360) {
    direction = relative_points / sqrt(rowSums(relative_points^2))
    local_x = as.vector(direction %*% right)
    local_y = as.vector(direction %*% up)
    local_z = as.vector(direction %*% forward)
    theta = acos(pmin(pmax(local_y, -1), 1))
    phi = atan2(local_x, local_z)
    s = ((phi - pi) / (2 * pi)) %% 1
    t = 1 - theta / pi
    in_front = rowSums(relative_points^2) > 0
  } else if (fov > 0) {
    aspect = camera_info$nx / camera_info$ny
    half_height = tan(fov * pi / 360)
    half_width = aspect * half_height
    s = 0.5 + x_camera / (2 * z_camera * half_width)
    t = 0.5 + y_camera / (2 * z_camera * half_height)
    in_front = z_camera > 0
  } else {
    warning(
      "`screen_text` projection is not supported for realistic camera ",
      "description files; skipping labels."
    )
    return(data.frame(
      x = rep(NA_real_, nrow(points)),
      y = rep(NA_real_, nrow(points)),
      in_front = FALSE,
      in_frame = FALSE
    ))
  }
  data.frame(
    x = 1 + (1 - s) * (camera_info$nx - 1),
    y = 1 + (1 - t) * (camera_info$ny - 1),
    in_front = in_front,
    in_frame = in_front & s >= 0 & s <= 1 & t >= 0 & t <= 1
  )
}

#' @keywords internal
add_screen_text = function(
  image_array,
  screen_text_spec,
  camera_info,
  screen_text_visible = NULL
) {
  screen_text_spec = normalize_screen_text(screen_text_spec)
  if (is.null(screen_text_spec) || nrow(screen_text_spec) == 0) {
    return(image_array)
  }
  if (!is.null(screen_text_visible)) {
    if (length(screen_text_visible) != nrow(screen_text_spec)) {
      stop("`screen_text_visible` must match the number of screen text labels")
    }
    screen_text_spec = screen_text_spec[screen_text_visible, , drop = FALSE]
    if (nrow(screen_text_spec) == 0) {
      return(image_array)
    }
  }
  projection = project_points_to_screen(
    screen_text_spec[, c("x", "y", "z"), drop = FALSE],
    camera_info
  )
  for (i in seq_len(nrow(screen_text_spec))) {
    if (!projection$in_front[i]) {
      next
    }
    if (isTRUE(screen_text_spec$clip[i]) && !projection$in_frame[i]) {
      next
    }
    text_image = rayimage::render_text_image(
      screen_text_spec$label[i],
      lineheight = screen_text_spec$lineheight[i],
      color = screen_text_spec$color[i],
      size = screen_text_spec$size[i],
      font = screen_text_spec$font[i],
      just = screen_text_spec$just[i],
      background_color = screen_text_spec$background_color[i],
      background_alpha = screen_text_spec$background_alpha[i]
    )
    if (!is.na(screen_text_spec$halo_color[i])) {
      halo_image = generate_screen_text_halo(
        text_image = text_image,
        halo_color = screen_text_spec$halo_color[i],
        halo_expand = screen_text_spec$halo_expand[i],
        halo_alpha = screen_text_spec$halo_alpha[i],
        halo_blur = screen_text_spec$halo_blur[i],
        halo_edge_softness = screen_text_spec$halo_edge_softness[i],
        halo_gap_fill = screen_text_spec$halo_gap_fill[i],
        halo_gap_fill_alpha_threshold = screen_text_spec$halo_gap_fill_alpha_threshold[
          i
        ]
      )
      halo_offset = screen_text_padding_offset(
        attr(halo_image, "padding"),
        screen_text_spec$hjust[i],
        screen_text_spec$vjust[i]
      )
      image_array = composite_screen_text_image(
        image_array = image_array,
        text_image = halo_image,
        x = projection$x[i] +
          screen_text_spec$x_offset[i] +
          screen_text_spec$halo_x_offset[i] +
          halo_offset[["x"]],
        y = projection$y[i] +
          screen_text_spec$y_offset[i] +
          screen_text_spec$halo_y_offset[i] +
          halo_offset[["y"]],
        hjust = screen_text_spec$hjust[i],
        vjust = screen_text_spec$vjust[i]
      )
    }
    image_array = composite_screen_text_image(
      image_array = image_array,
      text_image = text_image,
      x = projection$x[i] + screen_text_spec$x_offset[i],
      y = projection$y[i] + screen_text_spec$y_offset[i],
      hjust = screen_text_spec$hjust[i],
      vjust = screen_text_spec$vjust[i]
    )
  }
  image_array
}

#' @keywords internal
prepare_screen_text_preview = function(screen_text_spec) {
  screen_text_spec = normalize_screen_text(screen_text_spec)
  if (is.null(screen_text_spec) || nrow(screen_text_spec) == 0) {
    return(list(active = FALSE))
  }
  overlays = list()
  for (i in seq_len(nrow(screen_text_spec))) {
    text_image = rayimage::render_text_image(
      screen_text_spec$label[i],
      lineheight = screen_text_spec$lineheight[i],
      color = screen_text_spec$color[i],
      size = screen_text_spec$size[i],
      font = screen_text_spec$font[i],
      just = screen_text_spec$just[i],
      background_color = screen_text_spec$background_color[i],
      background_alpha = screen_text_spec$background_alpha[i]
    )
    if (!is.na(screen_text_spec$halo_color[i])) {
      halo_image = generate_screen_text_halo(
        text_image = text_image,
        halo_color = screen_text_spec$halo_color[i],
        halo_expand = screen_text_spec$halo_expand[i],
        halo_alpha = screen_text_spec$halo_alpha[i],
        halo_blur = screen_text_spec$halo_blur[i],
        halo_edge_softness = screen_text_spec$halo_edge_softness[i],
        halo_gap_fill = screen_text_spec$halo_gap_fill[i],
        halo_gap_fill_alpha_threshold = screen_text_spec$halo_gap_fill_alpha_threshold[
          i
        ]
      )
      halo_offset = screen_text_padding_offset(
        attr(halo_image, "padding"),
        screen_text_spec$hjust[i],
        screen_text_spec$vjust[i]
      )
      overlays[[length(overlays) + 1]] = make_screen_text_preview_overlay(
        image = halo_image,
        row = screen_text_spec[i, , drop = FALSE],
        x_offset = screen_text_spec$x_offset[i] +
          screen_text_spec$halo_x_offset[i] +
          halo_offset[["x"]],
        y_offset = screen_text_spec$y_offset[i] +
          screen_text_spec$halo_y_offset[i] +
          halo_offset[["y"]]
      )
    }
    overlays[[length(overlays) + 1]] = make_screen_text_preview_overlay(
      image = text_image,
      row = screen_text_spec[i, , drop = FALSE],
      x_offset = screen_text_spec$x_offset[i],
      y_offset = screen_text_spec$y_offset[i]
    )
  }
  list(active = length(overlays) > 0, overlays = overlays)
}

#' @keywords internal
make_screen_text_preview_overlay = function(image, row, x_offset, y_offset) {
  image = rayimage::ray_read_image(image, convert_to_array = TRUE)
  image = unclass(image)
  image = array(image, dim = dim(image))
  anchor_occlusion = isTRUE(row$occlusion) && row$occlusion_mode == "anchor"
  partial_occlusion = isTRUE(row$occlusion) && row$occlusion_mode == "label"
  list(
    image = image,
    x = row$x,
    y = row$y,
    z = row$z,
    x_offset = x_offset,
    y_offset = y_offset,
    hjust = row$hjust,
    vjust = row$vjust,
    clip = row$clip,
    occlusion = anchor_occlusion,
    partial_occlusion = partial_occlusion,
    occlusion_tolerance = row$occlusion_tolerance
  )
}

#' @keywords internal
screen_text_needs_native_overlay = function(screen_text_spec) {
  screen_text_spec = normalize_screen_text(screen_text_spec)
  if (is.null(screen_text_spec) || nrow(screen_text_spec) == 0) {
    return(FALSE)
  }
  any(screen_text_spec$occlusion & screen_text_spec$occlusion_mode == "label")
}

prepare_screen_text_occlusion = function(screen_text_spec) {
  screen_text_spec = normalize_screen_text(screen_text_spec)
  if (is.null(screen_text_spec) || nrow(screen_text_spec) == 0) {
    return(list(active = FALSE))
  }
  anchor_occlusion = screen_text_spec$occlusion &
    screen_text_spec$occlusion_mode == "anchor"
  list(
    active = any(anchor_occlusion),
    x = screen_text_spec$x,
    y = screen_text_spec$y,
    z = screen_text_spec$z,
    occlusion = anchor_occlusion,
    occlusion_tolerance = screen_text_spec$occlusion_tolerance
  )
}

#' @keywords internal
generate_screen_text_halo = function(
  text_image,
  halo_color,
  halo_expand,
  halo_alpha,
  halo_blur,
  halo_edge_softness,
  halo_gap_fill,
  halo_gap_fill_alpha_threshold
) {
  halo_padding = ceiling(
    max(0, halo_expand) +
      max(.Machine$double.eps, halo_edge_softness) +
      3 * max(0, halo_blur)
  )
  rayimage::render_alpha_outline(
    image = text_image,
    expand = halo_expand,
    edge_softness = halo_edge_softness,
    blur = halo_blur,
    color = halo_color,
    alpha = halo_alpha,
    pad = halo_padding,
    gap_fill = halo_gap_fill,
    gap_fill_alpha_threshold = halo_gap_fill_alpha_threshold,
    composite = FALSE,
    preview = FALSE
  )
}

#' @keywords internal
screen_text_padding_offset = function(padding, hjust, vjust) {
  if (is.null(padding)) {
    padding = 0
  }
  if (length(padding) == 1) {
    padding = rep(padding, 4)
  }
  padding = as.numeric(padding[1:4])
  c(
    x = hjust * (padding[4] + padding[2]) - padding[4],
    y = vjust * (padding[1] + padding[3]) - padding[1]
  )
}

#' @keywords internal
composite_screen_text_image = function(
  image_array,
  text_image,
  x,
  y,
  hjust,
  vjust
) {
  output = rayimage::render_sprite_overlay(
    image = image_array,
    image_overlay = text_image,
    convert_overlay_colorspace = FALSE,
    overlay_coords = c(x, y),
    hjust = hjust,
    vjust = vjust,
    preserve_channels = TRUE,
    preview = FALSE
  )
  if (!inherits(image_array, "rayimg")) {
    output = unclass(output)
    output = array(output, dim = dim(output))
  }
  output
}
