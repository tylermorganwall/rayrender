#' Internal vctrs methods
#'
#' @import pillar vctrs
#' @keywords internal
#' @name ray_animated_transform
NULL

#'@title Constructor for ray_transform
#'
#'@description list(start_transform_animation = list(matrix(4x4) OR ..NA_real..),
#'     end_transform_animation = list(matrix(4x4 OR ..NA_real..)),
#'     start_time  = numeric(1),
#'     end_time = numeric(1))
#'
#'@return ray_animated_transform
#'@keywords internal
ray_animated_transform = function(...) {
  vctrs::new_vctr(.data = list(list(...)), class = c("ray_animated_transform"))
}

#' @export
vec_ptype_abbr.ray_animated_transform = function(x, ...) {
  "ray_ani_tf"
}

#' @export
print.ray_animated_transform = function(x, ...) {
  num = length(x)
  cat(cli::col_grey(sprintf("ray_animated_transform list <%i>\n",num)))
  fmt = function(x) {
    cli::ansi_columns(
      text = paste0(paste(cli::col_red(names(x)), cli::col_grey(x), sep=": "),
                    collapse=" "),
      width = 100,
      fill = "rows",
      max_cols=6,
      align = "left",
      sep = "  "
    )
  }
  nprint = ifelse(num > 3, 3, num)
  for(i in seq_len(nprint)) {
    cat(sprintf("[[%i]] %s\n",i,fmt(x[[i]])))
  }
  if(num > 3) {
    plrl = ifelse(num > 4, "s", "")
    cat(cli::col_grey(sprintf("# %i more animated transform%s\n",num-3,plrl)))
  }
} 



#' @keywords internal
format_pillar_animated_transform = function(x) {
  format_material <- function(x) {
    output_desc = pillar::style_subtle(c("","|","", "|"))
    animated = FALSE
    if(nrow(x$start_transform_animation[[1]]) == 4) {
      animated = TRUE
      output_desc[1] = cli::style_bold(cli::col_green("Moving"))
    } else {
      output_desc = pillar::style_subtle("Not Moving")
    }
    if(animated) {
      output_desc[3] = cli::col_yellow(sprintf("%0.1f", x$start_time))
      output_desc[5] = pillar::style_subtle(cli::col_red(sprintf("%0.1f", x$end_time)))
    }
    sprintf("%s%s%s",
            pillar::style_subtle("<"), 
            paste0(output_desc, collapse=""),
            pillar::style_subtle(">"))
  }
  vapply(x, format_material, character(1))
}

#' @export
pillar_shaft.ray_animated_transform = function(x, ...) {
  pillar::new_pillar_shaft_simple(format_pillar_animated_transform(x),
                                  width = 10)
}
