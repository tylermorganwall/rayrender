#ifndef INTEGRATORH
#define INTEGRATORH

#include "sampler.h"
#include "Rcpp.h"
#include "RcppThread.h"
#include "RProgress.h"
#include "adaptivesampler.h"
#include "camera.h"
#include "hitable.h"
#include "hitablelist.h"
#include "color.h"
#include "mathinline.h"
#include "filter.h"

void pathtracer(size_t numbercores, size_t nx, size_t ny, size_t ns, int debug_channel,
                Float min_variance, size_t min_adaptive_size, 
                Rcpp::NumericMatrix& routput, Rcpp::NumericMatrix& goutput, Rcpp::NumericMatrix& boutput,
                bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                bool verbose, ortho_camera& ocam, camera &cam, environment_camera &ecam, Float fov,
                hitable_list& world, hitable_list& hlist,
                Float clampval, size_t max_depth, size_t roulette_active);

#endif