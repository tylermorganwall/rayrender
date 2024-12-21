#ifndef BUILDSCENEH
#define BUILDSCENEH

#include <Rcpp.h>
#include <memory>
#include "float.h"
using namespace Rcpp;

class TransformCache;
class TextureCache;
class material;
class hitable;
class random_gen;
class hitable_list;
class alpha_texture;
class bump_texture;
class roughness_texture;

std::shared_ptr<hitable> rotation_order(std::shared_ptr<hitable> entry, NumericVector temprotvec, NumericVector order_rotation);

std::shared_ptr<hitable> build_scene(List& scene,
                                     IntegerVector& shape,
                                     Float shutteropen, 
                                     Float shutterclose,
                                     std::vector<Float* >& textures, 
                                     std::vector<unsigned char * >& alpha_textures, 
                                     std::vector<unsigned char * >& bump_textures, 
                                     std::vector<unsigned char * >& roughness_textures,  
                                     std::vector<std::shared_ptr<material> >* shared_materials, 
                                     std::vector<std::shared_ptr<alpha_texture> >& alpha,
                                     std::vector<std::shared_ptr<bump_texture> >& bump,
                                     std::vector<std::shared_ptr<roughness_texture> >& roughness,
                                     int bvh_type,
                                     TransformCache& transformCache, 
                                     TextureCache& texCache,
                                     hitable_list& imp_sample_objects,
                                     std::vector<std::shared_ptr<hitable> >& instanced_objects,
                                     std::vector<std::shared_ptr<hitable_list> >& instance_importance_sampled,
                                     std::vector<int>& texture_idx,
                                     bool verbose,
                                     random_gen& rng); 


#endif
