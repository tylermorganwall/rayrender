#ifndef RAYRENDER_ANIMATION_RCPP_H
#define RAYRENDER_ANIMATION_RCPP_H
#include <Rcpp.h>
#include "math/keyframe.h"
#include "core/camera_track.h"
std::shared_ptr<const animation::CameraTrack> read_camera_track(const Rcpp::DataFrame &input,
                                                              double fps, double start);
std::shared_ptr<const animation::Track> read_keyframe_track(const Rcpp::List &input);
std::shared_ptr<const animation::TransformTrack> read_transform_track(const Rcpp::List &input);
#endif
