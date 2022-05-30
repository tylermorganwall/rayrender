ray_environment = new.env(parent = emptyenv())

.onLoad = function(libname, pkgname) {
  assign("keyframes", data.frame(), envir = ray_environment)
  assign("max_material_id", 1, envir = ray_environment)
}