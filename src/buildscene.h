#ifndef BUILDSCENEH
#define BUILDSCENEH

#include "hitable.h"
#include "sphere.h"
#include "hitablelist.h"
#include "bvh_node.h"
#include "perlin.h"
#include "texture.h"
#include "rectangle.h"
#include "box.h"
#include "constant.h"
#include "triangle.h"
#include "pdf.h"
#include "trimesh.h"
#include "disk.h"
#include "cylinder.h"
#include "ellipsoid.h"
#include "cone.h"
#include "curve.h"
#include "csg.h"
#include "plymesh.h"
#include "mesh3d.h"
#include <Rcpp.h>
#include <memory>
using namespace Rcpp;


std::shared_ptr<hitable> rotation_order(std::shared_ptr<hitable> entry, NumericVector temprotvec, NumericVector order_rotation);


std::shared_ptr<hitable> build_scene(IntegerVector& type, 
                                     NumericVector& radius, IntegerVector& shape,
                                     List& position_list,
                                     List& properties, List& velocity, LogicalVector& moving,
                                     int n, Float shutteropen, Float shutterclose,
                                     LogicalVector& ischeckered, List& checkercolors, 
                                     List gradient_info,
                                     NumericVector& noise, LogicalVector& isnoise,
                                     NumericVector& noisephase, NumericVector& noiseintensity, List noisecolorlist,
                                     List& angle, 
                                     LogicalVector& isimage, LogicalVector has_alpha,
                                     std::vector<Float* >& alpha_textures, std::vector<int* >& nveca,
                                     std::vector<Float* >& textures, std::vector<int* >& nvec,
                                     LogicalVector has_bump,
                                     std::vector<Float* >& bump_textures, std::vector<int* >& nvecb,
                                     NumericVector& bump_intensity,
                                     NumericVector& lightintensity,
                                     LogicalVector& isflipped,
                                     LogicalVector& isvolume, NumericVector& voldensity,
                                     List& order_rotation_list, 
                                     LogicalVector& isgrouped, List& group_pivot, List& group_translate,
                                     List& group_angle, List& group_order_rotation, List& group_scale,
                                     LogicalVector& tri_normal_bools, LogicalVector& is_tri_color, List& tri_color_vert, 
                                     CharacterVector& fileinfo, CharacterVector& filebasedir,
                                     List& scale_list, NumericVector& sigma,  List &glossyinfo,
                                     IntegerVector& shared_id_mat, LogicalVector& is_shared_mat,
                                     std::vector<std::shared_ptr<material> >* shared_materials, List& image_repeat_list,
                                     List& csg_info, List& mesh_list, int bvh_type,
                                     random_gen& rng);

std::shared_ptr<hitable> build_imp_sample(IntegerVector& type, 
                                          NumericVector& radius, IntegerVector& shape,
                                          List& position_list,
                                          List& properties, List& velocity, 
                                          int n, Float shutteropen, Float shutterclose, 
                                          List& angle, int i, List& order_rotation_list,
                                          LogicalVector& isgrouped, 
                                          List& group_pivot, List& group_translate,
                                          List& group_angle, List& group_order_rotation, List& group_scale,
                                          CharacterVector& fileinfo, CharacterVector& filebasedir,
                                          List& scale_list, List& mesh_list, int bvh_type, LogicalVector& moving,random_gen& rng);

#endif
