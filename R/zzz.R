ray_environment = new.env(parent = emptyenv())

.onLoad = function(libname, pkgname) {
  assign("keyframes", data.frame(), envir = ray_environment)
  assign("max_material_id", 1, envir = ray_environment)
  assign("init_time", 0, envir = ray_environment)
  assign("prev_time", 0, envir = ray_environment)
}