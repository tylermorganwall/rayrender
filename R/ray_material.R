#' Internal vctrs methods
#'
#' @import vctrs pillar
#' @keywords internal
#' @name ray_material
NULL

#'@title Constructor for ray_material
#'
#'@return ray_material
#'@keywords internal
ray_material = function(...) {
  vctrs::new_vctr(.data = list(...), class = c("ray_material"))
}

#' @export
vec_ptype_abbr.ray_material <- function(x, ...) {
  "ray_mat"
}

#' @export
print.ray_material = function(x, ...) {
  format_material = function(x) {
    mat_col = switch(x,
                     "diffuse" = cli::col_none,
                     "metal" = cli::col_grey,
                     "dielectric" = cli::col_blue, 
                     "oren-nayar" = cli::col_cyan, 
                     "light" = cli::col_yellow, 
                     "mf" = cli::col_green, 
                     "glossy" = cli::col_magenta, 
                     "spotlight" = cli::col_yellow, 
                     "hair" = cli::col_black, 
                     "mf-t" = cli::col_br_blue)
    mat_col(sprintf("%s",x))
  }
  num = length(x)
  cat(cli::col_grey(sprintf("ray_material list <%i>\n",num)))
  fmt = function(x) {
    cli::ansi_columns(
      text = paste0(c(paste(c(cli::col_red(names(x[1])), format_material(x[[1]])), collapse=": "),
                      paste(cli::col_red(names(x[-1])), cli::col_grey(x[-1]), sep=": ")),
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
    cat(cli::col_grey(sprintf("# %i more material%s\n",num-3,plrl)))
  }
} 

#' @keywords internal
format_pillar = function(x) {
  format_material <- function(x) {
    x_char = get_material_name(x$type)
    mat_col = switch(x$type,
                     "diffuse" = cli::col_none,
                     "metal" = cli::col_grey,
                     "dielectric" = cli::col_blue, 
                     "oren-nayar" = cli::col_cyan, 
                     "light" = cli::col_yellow, 
                     "mf" = cli::col_green, 
                     "glossy" = cli::col_magenta, 
                     "spotlight" = cli::col_yellow, 
                     "hair" = cli::col_black, 
                     "mf-t" = cli::col_br_blue)
    mat_col(sprintf("<%s>",x_char))
  }
  vapply(x, format_material, character(1))
}

#' @export
pillar_shaft.ray_material <- function(x, ...) {
  pillar::new_pillar_shaft_simple(format_pillar(x),
                                  width = 11)
}
