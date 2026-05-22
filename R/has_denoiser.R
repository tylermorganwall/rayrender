#' Check for Denoiser Support
#'
#' @description Returns TRUE if rayrender was compiled with Open Image Denoise (OIDN) support.
#'
#' @return Logical value.
#' @export
#'@examplesIf interactive() || identical(Sys.getenv("IN_PKGDOWN"), "true")
#' has_denoiser()
has_denoiser = function() {
  cppdef_HAS_OIDN()
}
