#ifndef INTEGRATORH
#define INTEGRATORH

#include "Rcpp.h"
#include "camera.h"
#include "hitablelist.h"
#include "PreviewDisplay.h"
#include "RayMatrix.h"
#include "color.h"


void pathtracer(std::size_t numbercores, std::size_t nx, std::size_t ny, std::size_t ns, int debug_channel,
                Float min_variance, std::size_t min_adaptive_size, 
                RayMatrix& rgb_output, RayMatrix& normalOutput, RayMatrix& albedoOutput, 
                RayMatrix& alpha_output,
                RayMatrix& draw_rgb_output,
                bool progress_bar, int sample_method, Rcpp::NumericVector& stratified_dim,
                bool verbose,RayCamera* cam, 
                Float fov,
                hitable_list& world, hitable_list& hlist,
                Float clampval, std::size_t max_depth, std::size_t roulette_active,
                PreviewDisplay& display, IntegratorType integrator_type);

#endif