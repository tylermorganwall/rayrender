#include "point.h"
#include <algorithm>
#include <climits>
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
  if (projection) return projection_power;
  const double luminance = .2126 * intensity[0] + .7152 * intensity[1] + .0722 * intensity[2];
  return luminance * (spot ? 2 * M_PI * (1 - .5 * (cos_inner + cos_outer)) : 4 * M_PI);
}
point3f PointLight::Project(const vec3f &outgoing) const {
  const Float z = dot(direction, outgoing);
  if (!(z > 0)) return point3f(0);
  const Float u = .5f + .5f * dot(right, outgoing) / (z * Float(half_width));
  const Float v = .5f + .5f * dot(up, outgoing) / (z * Float(half_height));
  if (!(u >= 0 && u <= 1 && v >= 0 && v <= 1)) return point3f(0);
  // Texel-center bilinear filtering with clamped edges. The rectangular beam
  // clips before filtering, so opposite image edges never bleed together.
  Float x = u * width - .5f, y = (1 - v) * height - .5f;
  int ix = int(std::floor(x)), iy = int(std::floor(y));
  Float fx = x - ix, fy = y - iy;
  point3f result(0);
  for (int j = 0; j < 2; ++j)
    for (int i = 0; i < 2; ++i)
      result += image[std::clamp(iy + j, 0, height - 1) * width + std::clamp(ix + i, 0, width - 1)] *
                (i ? fx : 1 - fx) * (j ? fy : 1 - fy);
  return result;
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
  const point3f pattern = projection ? Project(-result.wi) : point3f(1);
  for (int c = 0; c < 3; ++c)
    result.radiance[c] = Float(intensity[c] * pattern[c] * scale);
  result.pmf = 1;
  return result;
}
PointLightSet::PointLightSet(const Rcpp::List &descriptions) {
  double maximum = 0;
  for (SEXP element : descriptions) {
    Rcpp::List description(element);
    PointLight light;
    std::string type = Rcpp::as<std::string>(description["type"]);
    if (type != "point" && type != "spot" && type != "projection")
      Rcpp::stop("Invalid point/spot/projection light type.");
    light.spot = type == "spot";
    light.projection = type == "projection";
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
    if (light.projection) {
      Rcpp::NumericVector direction = description["direction"], up = description["up"];
      if (direction.size() != 3 || up.size() != 3)
        Rcpp::stop("Invalid projector orientation.");
      const double direction_length = std::hypot(direction[0], direction[1], direction[2]);
      const double up_length = std::hypot(up[0], up[1], up[2]);
      for (int c = 0; c < 3; ++c) {
        light.direction[c] = direction[c] / direction_length;
        light.up[c] = up[c] / up_length;
      }
      light.right = unit_vector(cross(light.direction, light.up));
      light.up = cross(light.right, light.direction);
      if (description.containsElementNamed("projection_basis")) {
        Rcpp::NumericMatrix basis = description["projection_basis"];
        if (basis.nrow() != 3 || basis.ncol() != 3)
          Rcpp::stop("Invalid projector basis.");
        for (int c = 0; c < 3; ++c) {
          light.right[c] = basis(0, c);
          light.up[c] = basis(1, c);
          light.direction[c] = basis(2, c);
        }
      }
      for (int c = 0; c < 3; ++c)
        if (!std::isfinite(light.right[c]) || !std::isfinite(light.up[c]) ||
            !std::isfinite(light.direction[c])) Rcpp::stop("Invalid projector basis.");
      double fov = Rcpp::as<double>(description["fov"]);
      if (!std::isfinite(fov) || fov <= 0 || fov >= 180)
        Rcpp::stop("Invalid projector field of view.");
      Rcpp::NumericVector pixels = description["image"];
      Rcpp::IntegerVector dims = pixels.attr("dim");
      if (dims.size() != 3 || dims[0] < 1 || dims[1] < 1 || dims[2] != 3 ||
          double(dims[0]) * dims[1] > INT_MAX / 3 || pixels.size() != double(dims[0]) * dims[1] * 3)
        Rcpp::stop("Invalid projector image dimensions.");
      light.height = dims[0]; light.width = dims[1];
      double aspect = double(light.width) / light.height;
      light.half_width = std::tan(fov * M_PI / 360) * std::max(1., aspect);
      light.half_height = std::tan(fov * M_PI / 360) * std::max(1., 1 / aspect);
      light.image.resize(size_t(light.width) * light.height);
      // Integrate texel-center intensity with the perspective solid-angle
      // Jacobian. This estimate guides selection only; radiance is not rescaled.
      double area = 4 * light.half_width * light.half_height / light.image.size();
      for (int y = 0; y < light.height; ++y) {
        double py = (1 - 2 * (y + .5) / light.height) * light.half_height;
        for (int x = 0; x < light.width; ++x) {
          point3f pixel;
          for (int c = 0; c < 3; ++c) {
            double value = pixels[y + size_t(light.height) * (x + size_t(light.width) * c)];
            if (!std::isfinite(value) || value < 0 || value > 1e30)
              Rcpp::stop("Invalid projector image value.");
            pixel[c] = Float(value);
          }
          light.image[y * light.width + x] = pixel;
          double px = (2 * (x + .5) / light.width - 1) * light.half_width;
          double luminance = .2126 * double(pixel[0]) * light.intensity[0] +
                             .7152 * double(pixel[1]) * light.intensity[1] +
                             .0722 * double(pixel[2]) * light.intensity[2];
          light.projection_power += luminance * area / std::pow(1 + px * px + py * py, 1.5);
        }
      }
    }
    if (light.Power() > 0) {
      maximum = std::max(maximum, light.Power());
      lights.push_back(std::move(light));
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
