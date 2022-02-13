#ifndef INTEGRATORH
#define INTEGRATORH

#include "Rcpp.h"
#include "camera.h"
#include "hitablelist.h"
#include "PreviewDisplay.h"


void pathtracer(std::size_t numbercores, std::size_t nx, std::size_t ny, std::size_t ns, int debug_channel,
                Float min_variance, std::size_t min_adaptive_size, 
                Rcpp::NumericMatrix& routput, Rcpp::NumericMatrix& goutput, Rcpp::NumericMatrix& boutput,
                bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                bool verbose, ortho_camera& ocam, camera &cam, environment_camera &ecam, 
                RealisticCamera &rcam,
                Float fov,
                hitable_list& world, hitable_list& hlist,
                Float clampval, std::size_t max_depth, std::size_t roulette_active,
                PreviewDisplay& display);

#endif