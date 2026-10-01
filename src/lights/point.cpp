#include "point.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

double PointLight::Falloff(const vec3f &outgoing) const {
  if (!spot)
    return 1;
  const double cosine = dot(direction, outgoing);
  if (cosine < cos_outer)
    return 0;
  if (cosine >= cos_inner || cos_inner == cos_outer)
    return 1;
  // PBRT v4 uses cubic SmoothStep in cos(theta), not the older quartic model.
  const double t = (cosine - cos_outer) / (cos_inner - cos_outer);
  return t * t * (3 - 2 * t);
}
double PointLight::Power() const {
  const double luminance = .2126 * intensity[0] + .7152 * intensity[1] + .0722 * intensity[2];
  return luminance * (spot ? 2 * M_PI * (1 - .5 * (cos_inner + cos_outer)) : 4 * M_PI);
}
PointLightSample PointLight::Sample(const point3f &p) const {
  PointLightSample result;
  result.position = position;
  vec3f difference = position - p;
  const double distance2 = double(difference[0]) * difference[0] + double(difference[1]) * difference[1] +
                           double(difference[2]) * difference[2];
  if (!(distance2 > 0) || !std::isfinite(distance2))
    return result;
  result.wi = difference / Float(std::sqrt(distance2));
  const double scale = Falloff(-result.wi) / distance2;
  for (int c = 0; c < 3; ++c)
    result.radiance[c] = Float(intensity[c] * scale);
  result.pmf = 1;
  return result;
}
PointLightSet::PointLightSet(const Rcpp::List &descriptions) {
  double maximum = 0;
  for (SEXP element : descriptions) {
    Rcpp::List description(element);
    PointLight light;
    std::string type = Rcpp::as<std::string>(description["type"]);
    if (type != "point" && type != "spot")
      Rcpp::stop("Invalid point/spot light type.");
    light.spot = type == "spot";
    Rcpp::NumericVector position = description["position"], color = description["color"];
    double intensity = Rcpp::as<double>(description["intensity"]);
    if (position.size() != 3 || color.size() != 3 || !std::isfinite(intensity) || intensity < 0)
      Rcpp::stop("Invalid point light values.");
    for (int c = 0; c < 3; ++c) {
      light.position[c] = Float(position[c]);
      light.intensity[c] = Float(color[c] * intensity);
      if (!std::isfinite(light.position[c]) || !std::isfinite(light.intensity[c]) || light.intensity[c] < 0)
        Rcpp::stop("Nonfinite or negative point light values.");
    }
    if (light.spot) {
      Rcpp::NumericVector direction = description["direction"];
      double outer = Rcpp::as<double>(description["cone_angle"]);
      double falloff = Rcpp::as<double>(description["falloff_angle"]);
      if (direction.size() != 3 || !std::isfinite(outer) || !std::isfinite(falloff) || outer <= 0 || outer > 180 ||
          falloff < 0 || falloff > outer)
        Rcpp::stop("Invalid spotlight cone.");
      double length = std::hypot(direction[0], direction[1], direction[2]);
      if (!(length > 0) || !std::isfinite(length))
        Rcpp::stop("Invalid spotlight direction.");
      for (int c = 0; c < 3; ++c)
        light.direction[c] = Float(direction[c] / length);
      light.cos_outer = std::cos(outer * M_PI / 180);
      light.cos_inner = std::cos((outer - falloff) * M_PI / 180);
    }
    if (light.Power() > 0) {
      lights.push_back(light);
      maximum = std::max(maximum, light.Power());
    }
  }
  // Normalize before summing so large scenes cannot overflow their total power.
  double sum = 0;
  for (const auto &light : lights) {
    sum += light.Power() / maximum;
    cdf.push_back(sum);
  }
  for (double &value : cdf)
    value /= sum;
  if (!cdf.empty())
    cdf.back() = 1;
}
PointLightSample PointLightSet::Sample(const point3f &p, double u) const {
  if (lights.empty())
    return {};
  size_t index = std::upper_bound(cdf.begin(), cdf.end(), std::clamp(u, 0., std::nextafter(1., 0.))) - cdf.begin();
  index = std::min(index, lights.size() - 1);
  PointLightSample result = lights[index].Sample(p);
  if (result.pmf > 0)
    result.pmf = cdf[index] - (index ? cdf[index - 1] : 0);
  return result;
}
