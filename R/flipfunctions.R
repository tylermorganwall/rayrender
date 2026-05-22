#' Flip Left-Right
#'
#' @param x Matrix
#'
#' @return Flipped matrix
#' @keywords internal
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' #Fake example
fliplr = function(x) {
  if (length(dim(x)) == 2) {
    x[, ncol(x):1]
  } else {
    x[, ncol(x):1, ]
  }
}


#' Flip Up-Down
#'
#' @param x Matrix
#'
#' @return Flipped matrix
#' @keywords internal
#'
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' #Fake example
flipud = function(x) {
  if (length(dim(x)) == 2) {
    x[nrow(x):1, ]
  } else {
    x[nrow(x):1, , ]
  }
}
