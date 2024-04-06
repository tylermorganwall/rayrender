#'@title Determines if rendering in knitr
#'
#'@return boolean
#'@keywords internal
is_rendering_in_knitr = function() {
  isTRUE(getOption('knitr.in.progress'))
}

#'@title Constructor for ray_material
#'
#'@return ray
#'@keywords internal
ray_scene = function(...) {
  structure(list(...), class = c("ray_scene", "tbl_df", "tbl", "data.frame"))
}

#' @keywords internal
style_subtle = function(x) {
  sprintf("\033[38;5;246m%s\033[39m",x)
}

#'@export
print.ray_scene = function(x, ...) {
  if(getOption("ray.ANSI", TRUE) && !is_rendering_in_knitr()) {
    boldstart = "\u001b[1m"
    formatend = "\u001b[0m"
    color = "\u001b[0;31m"
    blue_color = "\u001b[34m"
    colorend = "\u001b[0m"
    bullet = "\u2022"
  } else {
    boldstart = ""
    formatend = ""
    color = ""
    blue_color = ""
    colorend = ""
    bullet = "*"
  }
  generate_text = function(label, output) {
    sprintf("%s%s %s - %s%s \n",boldstart, bullet, label, formatend, output)
  }
  red = function(text) {
    sprintf("%s%s%s%s",color, boldstart, text,colorend)
  }
  blue = function(text) {
    sprintf("%s%s%s%s",boldstart, blue_color,text,colorend)
  }
  # Count total objects and lights
  total_objects <- nrow(x)
  total_lights = sum(unlist(lapply(x$material, \(x) x$type == "light")))
  
  # Count each type of object
  shape_counts <- table(x$shape)
  shape_summary <- sprintf("Objects - %s", 
                           paste(cli::col_blue(names(shape_counts)), 
                                 cli::col_red(shape_counts), sep = ": ", collapse = " | "))
  
  # Calculate bounding box
  bbxmin <- c(min(x$x), min(x$y), min(x$z))
  bbxmax <- c(max(x$x), max(x$y), max(x$z))
  
  # Construct the print output
  line1 = sprintf("Summary - %s: %s | %s: %s", 
                  cli::col_blue("Objects"), 
                  cli::col_red(as.character(total_objects)),
                  cli::col_blue("Lights"), 
                  cli::col_red(as.character(total_lights)))
  min_bbox = sprintf("c(%0.2f, %0.2f, %0.2f)",bbxmin[1],bbxmin[2],bbxmin[3])
  max_bbox = sprintf("c(%0.2f, %0.2f, %0.2f)",bbxmax[1],bbxmax[2],bbxmax[3])
  bbox_text = sprintf("XYZ Bounds - %s: %s | %s: %s", cli::col_blue("Min"), cli::col_red(min_bbox), cli::col_blue("Max"), cli::col_red(max_bbox))
  
  cli_output = function() {
    cli::cli_rule(left = "Scene Description")
    cli::cli_bullets(c(
      "*" = line1,
      ">" = shape_summary,
      "i" = bbox_text
    ))
    # cli::cli_li(c("{.emph Shapes}" = shape_summary))
    # cli_li("{.strong Strong} importance")
    # cli_li("A piece of code: {.code sum(a) / length(a)}")
    # cli_li("A package name: {.pkg cli}")
    # cli_li("A function name: {.fn cli_text}")
    # cli_li("A keyboard key: press {.kbd ENTER}")
    # cli_li("A file name: {.file /usr/bin/env}")
    # cli_li("An email address: {.email bugs.bunny@acme.com}")
    # cli_li("A URL: {.url https://acme.com}")
    # cli_li("An environment variable: {.envvar R_LIBS}")
    # cli_li("Some {.field field}")
  }
  
  # bbox_text = sprintf("%s: %s | %s: %s", cli::col_blue("Min"), cli::col_red(min_bbox), cli::col_blue("Max"), cli::col_red(max_bbox))
  # cat(generate_text("Rayscene Info", line1))
  cli_output()
  # cat(generate_text("Scene Details", shape_summary))
  # cat(generate_text("Center Bounds", bbox_text ))
  if(length(find.package("tibble",quiet=TRUE)) > 0) {
    withr::local_options(list(pillar.advice = FALSE, pillar.print_max = 3, pillar.width = 50),
      code = {
        print(tibble::as_tibble(x), ...)
      }
    )
  }
}
