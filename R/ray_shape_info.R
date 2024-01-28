#' Internal vctrs methods
#'
#' @import pillar vctrs
#' @keywords internal
#' @name ray_shape_info
NULL

#'@title Constructor for ray_transform
#'
#'@description list(angle = list(numeric(3)),
#'     order_rotation = list(numeric(3)),
#'     scale_factor  = list(numeric(3)),
#'     group_transform = list(matrix(4x4) ..OR.. NA_real))
#'
#'@return ray_transform
#'@keywords internal
ray_shape_info = function(...) {
  vctrs::new_vctr(.data = list(list(...)), class = c("ray_shape_info"))
}

#' @export
vec_ptype_abbr.ray_shape_info = function(x, ...) {
  "ray_shp"
}

#' @export
print.ray_shape_info = function(x, ...) {
  num = length(x)
  cat(cli::col_grey(sprintf("ray_shape_info list <%i>\n",num)))
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
    cat(cli::col_grey(sprintf("# %i more shapes%s\n",num-3,plrl)))
  }
} 



#' @keywords internal
format_pillar_shp = function(x) {
  format_material <- function(x) {
    sprintf("%s%s%s",
            pillar::style_subtle("<"), 
            pillar::style_subtle("shp_info"),
            pillar::style_subtle(">"))
  }
  vapply(x, format_material, character(1))
}

#' @export
pillar_shaft.ray_shape_info = function(x, ...) {
  pillar::new_pillar_shaft_simple(format_pillar_shp(x),
                                  width = 10)
}
