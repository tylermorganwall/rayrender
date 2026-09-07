#include "medium.h"
#include "../materials/texture.h"
#include "boundary.h"
#include <algorithm>
#include <stdexcept>

namespace {
point3f rgb(const Rcpp::List &d, const char *name) {
  Rcpp::NumericVector v(d[name]);
  if (v.size() != 1 && v.size() != 3 && !(std::string(name) == "emission" && v.hasAttribute("dim")))
    throw std::runtime_error("Medium RGB coefficients require one or three values.");
  return point3f(v[0], v[v.size() == 1 ? 0 : 1], v[v.size() == 1 ? 0 : 2]);
}
SampledField field(SEXP data, int channels = 1) {
  SampledField f;
  if (Rf_isNull(data))
    return f;
  Rcpp::NumericVector v = Rcpp::as<Rcpp::NumericVector>(data);
  f.values.reserve(v.size());
  for (double value : v) {
    if (!std::isfinite(value) || value < 0 || !std::isfinite(Float(value)))
      throw std::runtime_error(
          "Medium fields must be finite, nonnegative, and representable as float values.");
    f.values.push_back(Float(value));
  }
  f.channels = channels;
  if (v.hasAttribute("dim")) {
    Rcpp::IntegerVector dims = v.attr("dim");
    if (dims.size() != (channels == 1 ? 3 : 4) || (channels == 3 && dims[3] != 3))
      throw std::runtime_error("Invalid medium field dimensions.");
    size_t count = channels;
    for (int i = 0; i < 3; ++i) {
      if (dims[i] < 1)
        throw std::runtime_error("Medium grid dimensions must be positive.");
      f.dims[i] = dims[i];
      count *= size_t(dims[i]);
    }
    if (count != f.values.size())
      throw std::runtime_error("Medium array dimensions do not match its data.");
  }
  return f;
}
} // namespace

Float HGPhaseFunction::p(const vec3f &wo, const vec3f &wi) const {
  double cosine = 0, wo2 = 0, wi2 = 0;
  for (int a = 0; a < 3; ++a) {
    cosine += double(wo[a]) * wi[a];
    wo2 += double(wo[a]) * wo[a];
    wi2 += double(wi[a]) * wi[a];
  }
  cosine = std::clamp(cosine / std::sqrt(wo2 * wi2), -1.0, 1.0);
  double gg = g;
  double d = gg >= 0 ? (1 - gg) * (1 - gg) + 2 * gg * (1 + cosine)
                     : (1 + gg) * (1 + gg) - 2 * gg * (1 - cosine);
  return Float((1 - gg * gg) / (4 * M_PI * d * std::sqrt(d)));
}
PhaseFunctionSample HGPhaseFunction::Sample(const vec3f &wo, Float u, Float v) const {
  double gg = g, cos_theta = 1 - 2 * double(u);
  if (std::abs(gg) > 1e-6) {
    double s = (1 - gg * gg) / (1 - gg + 2 * gg * u);
    cos_theta = (1 + gg * gg - s * s) / (2 * gg);
  }
  cos_theta = std::clamp(cos_theta, -1.0, 1.0);
  double sin_theta = std::sqrt(std::max(0.0, 1 - cos_theta * cos_theta));
  onb frame;
  frame.build_from_w(-wo);
  vec3f wi = unit_vector(frame.local(sin_theta * std::cos(2 * M_PI * v),
                                     sin_theta * std::sin(2 * M_PI * v), cos_theta));
  Float value = p(wo, wi);
  return {wi, value, value};
}

RayMajorantIterator::RayMajorantIterator(const Ray &r, double t_max, const point3f &st,
                                         const MajorantGrid *g)
    : grid(g), sigma_t(st), end(t_max) {
  if (!g)
    return;
  for (int a = 0; a < 3; ++a) {
    double o = r.o[a], d = r.d[a];
    if (d == 0) {
      if (o < g->lo[a] || o > g->hi[a])
        finished = true;
      continue;
    }
    double t0 = (g->lo[a] - o) / d, t1 = (g->hi[a] - o) / d;
    if (t0 > t1)
      std::swap(t0, t1);
    t = std::max(t, t0);
    end = std::min(end, t1);
  }
  if (t >= end) {
    finished = true;
    return;
  }
  for (int a = 0; a < 3; ++a) {
    double width = (g->hi[a] - g->lo[a]) / double(g->resolution);
    double pos = (r.o[a] + t * r.d[a] - g->lo[a]) / width;
    cell[a] = std::clamp(int(std::floor(pos)), 0, g->resolution - 1);
    if (r.d[a] == 0) {
      step[a] = 0;
      next[a] = delta[a] = INFINITY;
    } else {
      step[a] = r.d[a] > 0 ? 1 : -1;
      double boundary = g->lo[a] + (cell[a] + (step[a] > 0 ? 1 : 0)) * width;
      next[a] = std::max(t, (boundary - r.o[a]) / r.d[a]);
      delta[a] = std::abs(width / r.d[a]);
    }
  }
}
std::optional<RayMajorantSegment> RayMajorantIterator::Next() {
  if (finished || !(t < end))
    return {};
  if (!grid) {
    finished = true;
    return RayMajorantSegment{t, end, sigma_t};
  }
  while (t < end) {
    double stop = std::min({end, next[0], next[1], next[2]});
    RayMajorantSegment result{t, stop, sigma_t * grid->Get(cell[0], cell[1], cell[2])};
    t = stop;
    for (int a = 0; a < 3; ++a)
      if (next[a] <= stop) {
        cell[a] += step[a];
        next[a] += delta[a];
        if (cell[a] < 0 || cell[a] >= grid->resolution)
          finished = true;
      }
    if (stop >= end)
      finished = true;
    if (result.t_max > result.t_min)
      return result;
    if (finished)
      return {};
  }
  return {};
}

Float SampledField::Lookup(const point3f &p, int channel) const {
  if (values.empty())
    return 0;
  int lower[3], upper[3];
  Float w[3];
  for (int a = 0; a < 3; ++a) {
    Float v = std::clamp(p[a] * dims[a] - 0.5f, Float(0), Float(dims[a] - 1));
    lower[a] = int(std::floor(v));
    upper[a] = std::min(lower[a] + 1, dims[a] - 1);
    w[a] = v - lower[a];
  }
  size_t plane = size_t(dims[0]) * dims[1] * dims[2];
  Float result = 0;
  for (int z = 0; z < 2; ++z)
    for (int y = 0; y < 2; ++y)
      for (int x = 0; x < 2; ++x) {
        size_t i = (x ? upper[0] : lower[0]) +
                   size_t(dims[0]) *
                       ((y ? upper[1] : lower[1]) + size_t(dims[1]) * (z ? upper[2] : lower[2]));
        result += values[i + channel * plane] * (x ? w[0] : 1 - w[0]) * (y ? w[1] : 1 - w[1]) *
                  (z ? w[2] : 1 - w[2]);
      }
  return result;
}

Medium::Medium(const Rcpp::List &d)
    : sigma_a(rgb(d, "sigma_a") * Rcpp::as<Float>(d["density_scale"])),
      sigma_s(rgb(d, "sigma_s") * Rcpp::as<Float>(d["density_scale"])),
      emission(rgb(d, "emission")), g(Rcpp::as<Float>(d["g"])),
      emission_scale(Rcpp::as<Float>(d["emission_scale"])),
      temperature_scale(Rcpp::as<Float>(d["temperature_scale"])),
      temperature_offset(Rcpp::as<Float>(d["temperature_offset"])),
      medium_to_object(Rcpp::as<Rcpp::NumericMatrix>(d["medium_transform"])) {
  ValidateMediumTransform(medium_to_object);
  if (!std::isfinite(emission_scale) || emission_scale < 0 || !std::isfinite(temperature_scale) ||
      temperature_scale < 0 || !std::isfinite(temperature_offset))
    throw std::runtime_error("Medium emission and temperature controls must be finite and "
                             "representable as float values.");
  for (int channel = 0; channel < 3; ++channel)
    if (!std::isfinite(sigma_a[channel] + sigma_s[channel]) || sigma_a[channel] < 0 ||
        sigma_s[channel] < 0)
      throw std::runtime_error("Scaled medium coefficients must be finite and nonnegative.");
  if (!std::isfinite(g) || std::abs(g) >= 1)
    throw std::runtime_error("Medium g must be strictly between -1 and 1.");
  for (double value : Rcpp::NumericVector(d["emission"])) {
    if (value < 0 || !std::isfinite(value * double(emission_scale)) ||
        !std::isfinite(Float(value * double(emission_scale))))
      throw std::runtime_error("Scaled medium emission must be finite and nonnegative.");
    has_rgb_emission |= value > 0;
  }
  has_temperature = !Rf_isNull(d["temperature"]);
  if (has_temperature)
    BlackbodyRGB(6500); // initialize the immutable lookup table before workers start
  if (has_temperature) {
    Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>(d["temperature"]);
    for (double value : values)
      if (!std::isfinite(value) || value < 0 || !std::isfinite(Float(value)) ||
          !std::isfinite(Float((value - temperature_offset) * temperature_scale)))
        throw std::runtime_error(
            "Scaled medium temperatures must be finite and representable as float values.");
    temperature = values[0];
  }
}
MediumProperties Medium::Properties(Float density, const point3f &le, const point3f &p) const {
  MediumProperties out;
  out.sigma_a = sigma_a * density;
  out.sigma_s = sigma_s * density;
  if (legacy_albedo) {
    point3f a = legacy_albedo->value(0, 0, p);
    out.sigma_s = sigma_s * a;
    out.sigma_a = sigma_s * (point3f(1) - a);
  }
  out.Le = le;
  out.phase = HGPhaseFunction(g);
  return out;
}
MediumProperties Medium::SamplePoint(const point3f &p) const {
  return Properties(Density(p), Emission(p), p);
}
RayMajorantIterator Medium::SampleRay(const Ray &r, double t_max) const {
  return RayMajorantIterator(r, t_max, sigma_a + sigma_s, nullptr);
}
point3f Medium::Emission(const point3f &) const {
  return emission_scale *
         (has_temperature ? BlackbodyRGB((temperature - temperature_offset) * temperature_scale)
                          : emission);
}
bool Medium::IsEmissive() const {
  return emission_scale > 0 && (sigma_a[0] > 0 || sigma_a[1] > 0 || sigma_a[2] > 0) &&
         (has_temperature || has_rgb_emission);
}
GridMedium::GridMedium(const Rcpp::List &d) : Medium(d) {
  density = field(d["density"]);
  temperatures = field(d["temperature"]);
  emissions = field(d["emission"], 3);
  Rcpp::NumericMatrix bounds(d["bounds"]);
  for (int a = 0; a < 3; ++a) {
    majorants.lo[a] = bounds(0, a);
    majorants.hi[a] = bounds(1, a);
  }
  for (int a = 0; a < 3; ++a)
    if (!std::isfinite(majorants.lo[a]) || !std::isfinite(majorants.hi[a]) ||
        !(majorants.lo[a] < majorants.hi[a]))
      throw std::runtime_error(
          "Grid bounds must be finite, ordered, and representable as float values.");
  int r = majorants.resolution;
  majorants.density.resize(r * r * r);
  for (int z = 0; z < r; ++z)
    for (int y = 0; y < r; ++y)
      for (int x = 0; x < r; ++x) {
        int cell[3] = {x, y, z}, lo[3], hi[3];
        for (int a = 0; a < 3; ++a) {
          lo[a] = std::clamp(int(std::floor(double(cell[a]) * density.dims[a] / r - 0.5)), 0,
                             density.dims[a] - 1);
          hi[a] = std::clamp(int(std::floor(double(cell[a] + 1) * density.dims[a] / r - 0.5)) + 1,
                             0, density.dims[a] - 1);
        }
        Float m = 0;
        for (int k = lo[2]; k <= hi[2]; ++k)
          for (int j = lo[1]; j <= hi[1]; ++j)
            for (int i = lo[0]; i <= hi[0]; ++i)
              m = std::max(
                  m,
                  density.values[i + size_t(density.dims[0]) * (j + size_t(density.dims[1]) * k)]);
        m = m == 0 ? 0 : std::nextafter(m, INFINITY);
        for (int a = 0; a < 3; ++a)
          if (!std::isfinite(m * (sigma_a[a] + sigma_s[a])))
            throw std::runtime_error(
                "Grid density times extinction exceeds the supported float range.");
        majorants.density[x + r * (y + r * z)] = m;
      }
}
point3f GridMedium::Normalize(const point3f &p) const {
  return point3f((p[0] - majorants.lo[0]) / (majorants.hi[0] - majorants.lo[0]),
                 (p[1] - majorants.lo[1]) / (majorants.hi[1] - majorants.lo[1]),
                 (p[2] - majorants.lo[2]) / (majorants.hi[2] - majorants.lo[2]));
}
Float GridMedium::Density(const point3f &p) const {
  point3f uvw = Normalize(p);
  for (int a = 0; a < 3; ++a)
    if (uvw[a] < 0 || uvw[a] > 1)
      return 0;
  return density.Lookup(uvw);
}
point3f GridMedium::Emission(const point3f &p) const {
  point3f uvw = Normalize(p);
  if (has_temperature)
    return emission_scale *
           BlackbodyRGB((temperatures.Lookup(uvw) - temperature_offset) * temperature_scale);
  return emission_scale *
         point3f(emissions.Lookup(uvw, 0), emissions.Lookup(uvw, 1), emissions.Lookup(uvw, 2));
}
MediumProperties GridMedium::SamplePoint(const point3f &p) const {
  return Properties(Density(p), Emission(p), p);
}
RayMajorantIterator GridMedium::SampleRay(const Ray &r, double t_max) const {
  return RayMajorantIterator(r, t_max, sigma_a + sigma_s, &majorants);
}

std::shared_ptr<const Medium> LoadMedium(const Rcpp::List &d) {
  std::string type = Rcpp::as<std::string>(d["type"]);
  if (type == "homogeneous")
    return std::make_shared<Medium>(d);
  if (type == "grid")
    return std::make_shared<GridMedium>(d);
  if (type == "nanovdb")
    return std::make_shared<NanoVDBMedium>(d);
  throw std::runtime_error("Unknown medium type: " + type);
}

size_t GridMedium::MemoryBytes() const {
  return sizeof(*this) +
         sizeof(Float) * (density.values.capacity() + temperatures.values.capacity() +
                          emissions.values.capacity() + majorants.density.capacity());
}
