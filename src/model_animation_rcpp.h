#ifndef RAYRENDER_MODEL_ANIMATION_RCPP_H
#define RAYRENDER_MODEL_ANIMATION_RCPP_H
#include <Rcpp.h>
#include "math/model_animation.h"
std::shared_ptr<const animation::ModelAsset> read_model_asset(SEXP input);
std::shared_ptr<const animation::ModelInstance> read_model_instance(
  const Rcpp::List &player, std::shared_ptr<const animation::ModelAsset> asset);
#endif
