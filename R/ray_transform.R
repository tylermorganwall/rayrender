#' Internal vctrs methods
#'
#' @import pillar vctrs
#' @keywords internal
#' @name ray_transform
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
ray_transform = function(...) {
  vctrs::new_vctr(.data = list(list(...)), class = c("ray_transform"))
}

#' @export
vec_ptype_abbr.ray_transform = function(x, ...) {
  "ray_tf"
}

#' @export
print.ray_transform = function(x, ...) {
  num = length(x)
  cat(cli::col_grey(sprintf("ray_transform list <%i>\n",num)))
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
    cat(cli::col_grey(sprintf("# %i more transform%s\n",num-3,plrl)))
  }
} 



#' @keywords internal
format_pillar_transform = function(x) {
  format_material <- function(x) {
    ord = paste0(c("x","y","z")[x$order_rotation[[1]]], collapse="")
    if(ord != "xyz") {
      ord = cli::style_underline(ord)
    }
    output_desc = rep(pillar::style_subtle("|"),5)
    if(any(x$angle[[1]] != 0)) {
      output_desc[1] = cli::style_bold(cli::col_cyan(sprintf("R%s",ord)))
    } else {
      output_desc[1] = pillar::style_subtle(sprintf("R%s",ord))
    }
    if(any(x$scale[[1]] != 1)) {
      output_desc[3] = cli::style_bold(cli::col_blue("S"))
    } else {
      output_desc[3] = pillar::style_subtle("S")
    }
    if(nrow(x$group_transform[[1]]) == 4) {
      output_desc[5] = cli::style_bold(cli::col_br_magenta("Gr"))
    } else {
      output_desc[5] = pillar::style_subtle("Gr")
    }
    sprintf("%s%s%s",
            pillar::style_subtle("<"), 
            paste0(output_desc, collapse=""),
            pillar::style_subtle(">"))
  }
  vapply(x, format_material, character(1))
}

#' @export
pillar_shaft.ray_transform = function(x, ...) {
  pillar::new_pillar_shaft_simple(format_pillar_transform(x),
                                  width = 10)
}
