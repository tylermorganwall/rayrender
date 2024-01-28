register_s3_method <- function(pkg, generic, class, fun = NULL) {
  stopifnot(is.character(pkg), length(pkg) == 1)
  stopifnot(is.character(generic), length(generic) == 1)
  stopifnot(is.character(class), length(class) == 1)
  
  if (is.null(fun)) {
    fun <- get(paste0(generic, ".", class), envir = parent.frame())
  } else {
    stopifnot(is.function(fun))
  }
  
  if (pkg %in% loadedNamespaces()) {
    registerS3method(generic, class, fun, envir = asNamespace(pkg))
  }
  
  # Always register hook in case package is later unloaded & reloaded
  setHook(
    packageEvent(pkg, "onLoad"),
    function(...) {
      registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
  )
}

.onLoad <- function(...) {
  register_s3_method("rayrender", "print", "ray_scene")
  register_s3_method("rayrender", "print", "ray_material")
  register_s3_method("rayrender", "format", "ray_material")
  register_s3_method("pillar", "pillar_shaft", "ray_material")
  register_s3_method("vctrs", "vec_ptype_abbr", "ray_material")
  
  register_s3_method("rayrender", "print", "ray_transform")
  register_s3_method("rayrender", "format", "ray_transform")
  register_s3_method("pillar", "pillar_shaft", "ray_transform")
  register_s3_method("vctrs", "vec_ptype_abbr", "ray_transform")
  
  register_s3_method("rayrender", "print", "ray_animated_transform")
  register_s3_method("rayrender", "format", "ray_animated_transform")
  register_s3_method("pillar", "pillar_shaft", "ray_animated_transform")
  register_s3_method("vctrs", "vec_ptype_abbr", "ray_animated_transform")
  
  register_s3_method("rayrender", "print", "ray_shape_info")
  register_s3_method("rayrender", "format", "ray_shape_info")
  register_s3_method("pillar", "pillar_shaft", "ray_shape_info")
  register_s3_method("vctrs", "vec_ptype_abbr", "ray_shape_info")
}

ray_environment = new.env(parent = emptyenv())

.onLoad = function(libname, pkgname) {
  assign("keyframes", data.frame(), envir = ray_environment)
  assign("max_material_id", 1, envir = ray_environment)
  assign("init_time", 0, envir = ray_environment)
  assign("prev_time", 0, envir = ray_environment)
}