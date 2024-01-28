#' New Tibble Row
#' 
#' Creates a row of a tibble, without the parsing and checks in tibble::new_tibble(). Internal use only.
#'
#' @param x Named list.
#'
#' @return Tibble row.
#' @keywords internal
#'
#' @examples
#' #none
new_tibble_row = function(x) {
  x = unclass(x)
  attr(x, "row.names") = 1L
  attr(x, "names") = names(x)
  class(x) = c("ray_scene", "tbl_df", "tbl", "data.frame")
  return(x)
}