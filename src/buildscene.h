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
#include "curve.h"
#include "csg.h"
#include "plymesh.h"
#include "mesh3d.h"
#include "raymesh.h"
#include "transform.h"
#include "transformcache.h"
#include <Rcpp.h>
#include <memory>
using namespace Rcpp;


std::shared_ptr<hitable> rotation_order(std::shared_ptr<hitable> entry, NumericVector temprotvec, NumericVector order_rotation);


std::shared_ptr<hitable> build_scene(List& scene,
                                     IntegerVector& shape,
                                     List& position_list,
                                     Float shutteropen, 
                                     Float shutterclose,
                                     std::vector<Float* >& textures, 
                                     std::vector<int* >& nvec,
                                     std::vector<unsigned char * >& alpha_textures, 
                                     std::vector<int* >& nveca,
                                     std::vector<unsigned char * >& bump_textures, 
                                     std::vector<int* >& nvecb,
                                     std::vector<unsigned char * >& roughness_textures,  
                                     std::vector<int* >& nvecr,
                                     std::vector<std::shared_ptr<material> >* shared_materials, 
                                     int bvh_type,
                                     TransformCache& transformCache, 
                                     hitable_list& imp_sample_objects,
                                     bool verbose,
                                     random_gen& rng); 


#endif
