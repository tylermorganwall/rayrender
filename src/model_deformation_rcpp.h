#ifndef RAYRENDER_MODEL_DEFORMATION_RCPP_H
#define RAYRENDER_MODEL_DEFORMATION_RCPP_H
#include <Rcpp.h>
#include "math/model_animation.h"
void read_model_deformation(const Rcpp::List &, animation::ModelNode &, const animation::ModelAsset &);
#endif
