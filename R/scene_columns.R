# All four descriptor classes store opaque scene records in ordinary lists.
# Shortcut only the exact, unnamed constructor representation. In particular,
# matching class names alone do not establish compatibility when a user adds
# metadata or subclasses a descriptor; leave those cases to vctrs.
#
# These methods deliberately preserve x rather than restoring from to: casting
# must retain the original records, and the shortcut already has identical
# attributes. This avoids constructing empty slices merely to compare types.

#' @keywords internal
#' @export
vec_cast.ray_material.ray_material = function(x, to, ...) {
  attrs = attributes(x)
  if (
    typeof(x) == "list" &&
      typeof(to) == "list" &&
      length(attrs) == 1L &&
      identical(attrs, attributes(to)) &&
      identical(class(x)[-1L], c("vctrs_vctr", "list"))
  ) {
    return(x)
  }
  vctrs::vec_default_cast(x, to, ...)
}

#' @keywords internal
#' @export
vec_cast.ray_shape_info.ray_shape_info = vec_cast.ray_material.ray_material

#' @keywords internal
#' @export
vec_cast.ray_transform.ray_transform = vec_cast.ray_material.ray_material

#' @keywords internal
#' @export
vec_cast.ray_animated_transform.ray_animated_transform = vec_cast.ray_material.ray_material

#' @keywords internal
#' @export
vec_ptype2.ray_material.ray_material = function(x, y, ...) {
  attrs = attributes(x)
  if (
    typeof(x) == "list" &&
      typeof(y) == "list" &&
      length(attrs) == 1L &&
      identical(attrs, attributes(y)) &&
      identical(class(x)[-1L], c("vctrs_vctr", "list"))
  ) {
    return(vctrs::vec_ptype(x))
  }
  vctrs::vec_default_ptype2(x, y, ...)
}

#' @keywords internal
#' @export
vec_ptype2.ray_shape_info.ray_shape_info = vec_ptype2.ray_material.ray_material

#' @keywords internal
#' @export
vec_ptype2.ray_transform.ray_transform = vec_ptype2.ray_material.ray_material

#' @keywords internal
#' @export
vec_ptype2.ray_animated_transform.ray_animated_transform = vec_ptype2.ray_material.ray_material
