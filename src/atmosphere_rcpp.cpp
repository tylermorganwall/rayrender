#include <Rcpp.h>
#include "lights/atmosphere.h"


// Internal reference queries for validation, in the same world coordinates and
// linear RGB units as the renderer. No exposure, compositing, or denoising.
// [[Rcpp::export]]
Rcpp::List query_prague_atmosphere(Rcpp::List description, Rcpp::NumericMatrix positions,
                                  Rcpp::NumericMatrix directions, Rcpp::NumericVector distances,
                                  bool sample = false, bool build_sampler = false) {
  // Each row describes one query. Ordinary queries provide a 3D direction;
  // sampling queries provide two variates for the light's directional sampler.
  int n = positions.nrow();
  if (positions.ncol() != 3 || directions.nrow() != n ||
      directions.ncol() != (sample ? 2 : 3) || distances.size() != n)
    Rcpp::stop("Atmosphere query dimensions do not match.");


  // Sampling requires proposal tables. build_sampler also lets callers evaluate
  // their PDFs at supplied directions; transport-only queries can skip that work.
  PragueInfiniteLight light(description, sample || build_sampler);
  Rcpp::NumericMatrix radiance(n, 3), transmission(n, 3), inscatter(n, 3), wi(n, 3);
  Rcpp::NumericVector pdf(n);


  for (int i = 0; i < n; ++i) {
    // Validate each query while keeping long batches interruptible from R.
    // Positive infinity is a valid distance for transmission to the environment.
    if (i % 64 == 0) Rcpp::checkUserInterrupt();
    for (int c = 0; c < 3; ++c)
      if (!std::isfinite(positions(i, c))) Rcpp::stop("Nonfinite atmosphere query position.");
    for (int c = 0; c < directions.ncol(); ++c)
      if (!std::isfinite(directions(i, c))) Rcpp::stop("Nonfinite atmosphere query direction.");
    if (std::isnan(distances[i]) || distances[i] < 0) Rcpp::stop("Invalid atmosphere query distance.");


    // Resolve the actual direction first, then use its unit vector for every
    // transport and PDF query. The returned direction also exposes sampler output.
    point3f p(positions(i, 0), positions(i, 1), positions(i, 2));
    vec3f w = sample ? light.Sample(p, vec2f(directions(i, 0), directions(i, 1)), 0)
                    : vec3f(directions(i, 0), directions(i, 1), directions(i, 2));
    if (!(w.squared_length() > 0)) Rcpp::stop("Zero atmosphere query direction.");
    w = unit_vector(w);


    // Radiance includes the enabled environment components. In-scattering is
    // the finite segment's haze, or the complete sky-only field at infinity;
    // the unscattered solar disk is not part of that in-scattering term.
    auto le = light.Radiance(p, w, 0);
    auto tr = light.Transmission(p, w, distances[i]);
    auto source = std::isfinite(distances[i]) ? light.Segment(p, w, distances[i]).radiance
                                             : light.SkyRadiance(p, w);
    pdf[i] = light.Pdf(p, w, 0);


    // Preserve the input row order, with RGB and (x,y,z) direction components in
    // columns. These raw values can be compared directly with renderer queries.
    for (int c = 0; c < 3; ++c) {
      wi(i, c) = w[c]; radiance(i, c) = le[c];
      transmission(i, c) = tr[c]; inscatter(i, c) = source[c];
    }
  }


  return Rcpp::List::create(Rcpp::Named("radiance") = radiance,
                            Rcpp::Named("transmission") = transmission,
                            Rcpp::Named("inscatter") = inscatter,
                            Rcpp::Named("direction") = wi, Rcpp::Named("pdf") = pdf);
}
