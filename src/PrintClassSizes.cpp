#include "core/buildscene.h"
#include "hitables/hitable.h"
#include "hitables/sphere.h"
#include "hitables/hitablelist.h"

#include "math/perlin.h"
#include "materials/texture.h"
#include "hitables/rectangle.h"
#include "hitables/box.h"
#include "materials/constant.h"
#include "hitables/triangle.h"
#include "math/pdf.h"
#include "hitables/trimesh.h"
#include "hitables/disk.h"
#include "hitables/cylinder.h"
#include "hitables/ellipsoid.h"
#include "hitables/curve.h"
#include "hitables/csg.h"
#include "hitables/plymesh.h"
#include "hitables/mesh3d.h"
#include "hitables/raymesh.h"
#include "hitables/instance.h"
#include "math/transform.h"
#include "math/transformcache.h"
#include "hitables/infinite_area_light.h"
#include "core/camera.h"
#include "math/RayMatrix.h"

// [[Rcpp::export]] 
void PrintClassSizes() {
  Rcpp::Rcout << "hitable                  : " << sizeof(hitable) << "\n";
  Rcpp::Rcout << "---sphere                : " << sizeof(sphere) << "\n";
  Rcpp::Rcout << "---xy_rect               : " << sizeof(xy_rect) << "\n";
  Rcpp::Rcout << "---xz_rect               : " << sizeof(xz_rect) << "\n";
  Rcpp::Rcout << "---yz_rect               : " << sizeof(yz_rect) << "\n";
  Rcpp::Rcout << "---box                   : " << sizeof(box) << "\n";
  Rcpp::Rcout << "---AnimatedHitable       : " << sizeof(AnimatedHitable) << "\n";
  Rcpp::Rcout << "---triangle              : " << sizeof(triangle) << "\n";
  Rcpp::Rcout << "---trimesh               : " << sizeof(trimesh) << "\n";
  Rcpp::Rcout << "---disk                  : " << sizeof(disk) << "\n";
  Rcpp::Rcout << "---cylinder              : " << sizeof(cylinder) << "\n";
  Rcpp::Rcout << "---ellipsoid             : " << sizeof(ellipsoid) << "\n";
  Rcpp::Rcout << "---curve                 : " << sizeof(curve) << "\n";
  Rcpp::Rcout << "---csg                   : " << sizeof(csg) << "\n";
  Rcpp::Rcout << "---plymesh               : " << sizeof(plymesh) << "\n";
  Rcpp::Rcout << "---mesh3d                : " << sizeof(mesh3d) << "\n\n";
  Rcpp::Rcout << "---InfiniteAreaLight     : " << sizeof(InfiniteAreaLight) << "\n\n";
  
  Rcpp::Rcout << "Float                    : " << sizeof(Float) << "\n";
  Rcpp::Rcout << "Matrix4x4                : " << sizeof(Matrix4x4) << "\n";
  Rcpp::Rcout << "Transform                : " << sizeof(Transform) << "\n";
  Rcpp::Rcout << "TransformCache           : " << sizeof(TransformCache) << "\n";
  Rcpp::Rcout << "RayMatrix                : " << sizeof(RayMatrix) << "\n";
  Rcpp::Rcout << "vec3f                    : " << sizeof(vec3f) << "\n";
  Rcpp::Rcout << "vec3i                    : " << sizeof(vec3i) << "\n";
  Rcpp::Rcout << "point3f                  : " << sizeof(point3f) << "\n";
  Rcpp::Rcout << "point3f                  : " << sizeof(point3f) << "\n";
  Rcpp::Rcout << "vec2f                    : " << sizeof(vec2f) << "\n";
  Rcpp::Rcout << "vec2i                    : " << sizeof(vec2i) << "\n";
  Rcpp::Rcout << "random_gen               : " << sizeof(random_gen) << "\n";
  Rcpp::Rcout << "aabb                     : " << sizeof(aabb) << "\n\n";
  Rcpp::Rcout << "hitable_list             : " << sizeof(hitable_list) << "\n\n";
  
  Rcpp::Rcout << "RayCamera                : " << sizeof(RayCamera) << "\n";
  Rcpp::Rcout << "---camera                : " << sizeof(camera) << "\n";
  Rcpp::Rcout << "---ortho_camera          : " << sizeof(ortho_camera) << "\n";
  Rcpp::Rcout << "---environment_camera    : " << sizeof(RayCamera) << "\n";
  Rcpp::Rcout << "---RealisticCamera       : " << sizeof(RealisticCamera) << "\n\n";
  Rcpp::Rcout << "material                 : " << sizeof(material) << "\n";
  Rcpp::Rcout << "---lambertian            : " << sizeof(lambertian) << "\n";
  Rcpp::Rcout << "---metal                 : " << sizeof(metal) << "\n";
  Rcpp::Rcout << "---dielectric            : " << sizeof(dielectric) << "\n";
  Rcpp::Rcout << "---diffuse_light         : " << sizeof(diffuse_light) << "\n";
  Rcpp::Rcout << "---spot_light            : " << sizeof(spot_light) << "\n";
  Rcpp::Rcout << "---isotropic             : " << sizeof(isotropic) << "\n";
  Rcpp::Rcout << "---orennayar             : " << sizeof(orennayar) << "\n";
  Rcpp::Rcout << "---MicrofacetReflection  : " << sizeof(MicrofacetReflection) << "\n";
  Rcpp::Rcout << "---MicrofacetTransmission: " << sizeof(MicrofacetTransmission) << "\n";
  Rcpp::Rcout << "---glossy                : " << sizeof(glossy) << "\n";
  Rcpp::Rcout << "---hair                  : " << sizeof(hair) << "\n";
}