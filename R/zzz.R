ray_environment = new.env(parent = emptyenv())

.onLoad = function(libname, pkgname) {
  assign("keyframes", data.frame(), envir = ray_environment)
}