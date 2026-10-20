#ifndef RAYRENDER_ANIMATION_REPORT_H
#define RAYRENDER_ANIMATION_REPORT_H
#include <Rcpp.h>
#include "../math/animatedtransform.h"

inline Rcpp::List AnimationTimingToR(const SceneTimeContext &time) {
  return Rcpp::List::create(
    Rcpp::_["compile_seconds"] = time.compile_seconds,
    Rcpp::_["pose_seconds"] = time.pose_seconds,
    Rcpp::_["refit_seconds"] = time.refit_seconds,
    Rcpp::_["meshes_updated"] = double(time.meshes_updated),
    Rcpp::_["shutter_samples"] = double(time.shutter_samples),
    Rcpp::_["shutter_times"] = time.shutter_times,
    Rcpp::_["shutter_pose_seconds"] = time.shutter_pose_seconds,
    Rcpp::_["shutter_refit_seconds"] = time.shutter_refit_seconds,
    Rcpp::_["geometry_revision"] = double(time.geometry_revision));
}
#endif
