#include "atmosphere.h"
#include "../volumes/cie.h"
#include <algorithm>
#include <limits>
#include <mutex>
#include <stdexcept>


// Adapt Prague's spectral model to both an infinite light and finite atmospheric
// segments. Using the same sky field for both lets a path accumulate nearby haze
// and then connect to the distant environment without counting that haze twice.
namespace {
using Model = skymodelr::PragueSkyModel;
using V = Model::Vector3;

// Distances in the model's local frame are meters; angular constants are radians.
// The uniform sampling floor keeps directions outside the tabulated proposals
// reachable, and the native solar radius defines Prague's disk profile.
constexpr double pi = 3.14159265358979323846, earth_radius = 6378000;
constexpr double uniform_fraction = .001, native_sun_radius = .004654793;


// Small vector operations shared by the spherical geometry and sampling code.
double dotv(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
double length(V a) { return std::sqrt(dotv(a, a)); }
V normalized(V a) { return a / length(a); }


// Remove the fitted table's zero-distance bias. Below 100 m, scale optical depth
// from the 100 m value so tiny segments approach unit transmission continuously.
double normalized_transmission(double zero, double value, double distance) {
  value = zero > 0 ? std::clamp(value / zero, 0.0, 1.0) : 0;
  return distance < 100 ? std::pow(value, distance / 100) : value;
}


// Direction helpers use Prague's Z-up frame. The latitude-longitude mappings
// below are inverses; uniform() instead samples solid angle directly.
V crossv(V a, V b) {
  return V(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

V from_uv(vec2f uv) {
  double theta = (1 - uv[1]) * pi, phi = 2 * pi * uv[0];
  return V(std::cos(phi) * std::sin(theta), std::sin(phi) * std::sin(theta),
           std::cos(theta));
}

V uniform(vec2f u) {
  double z = 1 - 2 * u[0], r = std::sqrt(std::max(0.0, 1 - z * z));
  return V(r * std::cos(2 * pi * u[1]), r * std::sin(2 * pi * u[1]), z);
}

vec2f to_uv(V w) {
  double phi = std::atan2(w.y, w.x);
  if (phi < 0) phi += 2 * pi;
  return vec2f(phi / (2 * pi), 1 - std::acos(std::clamp(w.z, -1.0, 1.0)) / pi);
}


std::shared_ptr<const Model> load_model(const std::string &path, double visibility,
                                      bool cache_spectra, bool transmission_table,
                                      double transmission_table_max_mb) {
  // Keep only the most recently used visibility slice between renders. Loading
  // happens before workers start; const queries share the immutable coefficients.
  static std::mutex mutex;
  static std::shared_ptr<const Model> cached;
  static std::string last_path;
  static double last_visibility = -1;
  static bool last_cache_spectra = false, last_transmission_table = false;
  static double last_table_max_mb = -1;


  // Resolve R's registered callable on the main thread, before taking a C++
  // lock. Every subsequent operation on this wrapper stays in native code.
  auto model = std::make_shared<Model>();
  std::lock_guard<std::mutex> guard(mutex);
  if (cached && path == last_path && visibility == last_visibility &&
      cache_spectra == last_cache_spectra && transmission_table == last_transmission_table &&
      transmission_table_max_mb == last_table_max_mb)
    return cached;


  // A ground-only dataset cannot describe elevated surfaces or cloud paths.
  // Publish the new slice to the cache only after checking its altitude coverage.
  model->initialize(path, visibility, cache_spectra, transmission_table, transmission_table_max_mb);
  auto data = model->getAvailableData();
  if (data.altitudeMax < 15000 || data.altitudeMin > 0)
    throw std::runtime_error("Atmospheric sky lights require the full-altitude Prague dataset.");
  cached = model; last_path = path; last_visibility = visibility;
  last_cache_spectra = cache_spectra; last_transmission_table = transmission_table;
  last_table_max_mb = transmission_table_max_mb;
  return model;
}


// Validate the R description during construction, before worker threads can
// issue native queries. Vector fields cover the scene origin and RGB gains.
double scalar(const Rcpp::List &x, const char *name, double lo, double hi) {
  double value = Rcpp::as<double>(x[name]);
  if (!std::isfinite(value) || value < lo || value > hi)
    throw std::runtime_error(std::string("Invalid Prague atmosphere parameter: ") + name);
  return value;
}


// Missing switches preserve descriptions saved before lighting and finite haze
// became independently configurable. Reject NA rather than treating it as true.
bool flag(const Rcpp::List &x, const char *name, bool fallback = true) {
  if (!x.containsElementNamed(name)) return fallback;
  SEXP value = x[name];
  if (TYPEOF(value) != LGLSXP || Rf_xlength(value) != 1 || LOGICAL(value)[0] == NA_LOGICAL)
    throw std::runtime_error(std::string(name) + " must be TRUE or FALSE.");
  return LOGICAL(value)[0] != 0;
}


std::array<double, 3> triple(const Rcpp::List &x, const char *name, bool positive = false) {
  Rcpp::NumericVector values = x[name];
  if (values.size() != 3) throw std::runtime_error(std::string(name) + " must have three values.");
  std::array<double, 3> result;
  for (int c = 0; c < 3; ++c) {
    result[c] = values[c];
    if (!std::isfinite(result[c]) || (positive && result[c] <= 0))
      throw std::runtime_error(std::string("Invalid Prague atmosphere ") + name);
  }
  return result;
}
} // namespace


// Keep one anchor while a ray crosses null events or invisible boundaries.
// Restart when direction changes or a glass/excluded-volume boundary toggles
// transport. Ordinary subdivisions of one air interval keep their common fit.
void AtmosphereRay::Start(const point3f &p, const vec3f &w, bool integrate) {
  bool enabled = model && integrate;
  if (!initialized || active != enabled || w[0] != direction[0] ||
      w[1] != direction[1] || w[2] != direction[2]) {
    initialized = true; active = enabled; origin = p; direction = w;
    cumulative = AtmosphereSegment(); distance = 0;
    cache.ready = cache.near_ready = false;
    pending = false; reservoir_mass = 0;
  }
}


AtmosphereSegment AtmosphereRay::Advance(const point3f &p) {
  AtmosphereSegment result;
  if (!active) return result;


  // Measure progress along the original ray, rather than summing the lengths
  // of subsegments whose origins may have been offset at geometric boundaries.
  double next = 0;
  for (int c = 0; c < 3; ++c) next += (double(p[c]) - origin[c]) * direction[c];
  // Geometric ray offsets can move a boundary point back by its error bound.
  if (next <= distance) return result;


  // Factor out transport already applied by the integrator. For endpoints a,b:
  // T_ab = T_0b / T_0a and L_ab = (L_0b - L_0a) / T_0a. These increments compose
  // back to the same cumulative result regardless of the number of boundaries.
  auto value = model->Segment(origin, direction, next, cache);
  for (int c = 0; c < 3; ++c) {
    if (cumulative.transmission[c] > 0) {
      result.transmission[c] = value.transmission[c] / cumulative.transmission[c];
      result.radiance[c] = (value.radiance[c] - cumulative.radiance[c]) / cumulative.transmission[c];
    } else {
      result.transmission[c] = result.radiance[c] = 0;
    }
  }
  cumulative = value; distance = next;
  return result;
}


// Factor atmospheric transmission out of the per-event throughput. If S(d) is
// cumulative in-scattering at our fixed anchor, a span contributes
// sum_i c_i [S(d_i) - S(d_{i-1})]. Summation by parts leaves an endpoint term
// plus corrections (c_i - c_{i+1}) at intermediate events. Constant weights
// telescope exactly, including arbitrarily many neutral primary-ray null events.
void AtmosphereRay::Accumulate(const point3f &p, const std::array<double, 3> &weight,
                               double uniform) {
  if (!active) return;
  double next = 0;
  for (int c = 0; c < 3; ++c) next += (double(p[c]) - origin[c]) * direction[c];
  if (next <= (pending ? pending_distance : distance)) return;
  std::array<double, 3> coefficient{}, change{};
  double mass = 0;
  for (int c = 0; c < 3; ++c) {
    coefficient[c] = cumulative.transmission[c] > 0 ? weight[c] / cumulative.transmission[c] : 0;
    change[c] = pending_weight[c] - coefficient[c];
    mass = std::max(mass, std::abs(change[c]));
  }


  // Select a correction using only its RGB weight, without evaluating S there.
  // Dividing by its selection probability at Flush preserves the conditional
  // expectation, including signed corrections from colored media or roulette.
  if (pending && mass > 0) {
    reservoir_mass += mass;
    if (reservoir_mass == mass || uniform < mass / reservoir_mass) {
      selected_weight = change;
      selected_distance = pending_distance;
      selected_mass = mass;
    }
  }
  pending = true;
  pending_distance = next;
  pending_weight = coefficient;
}


// Apply deferred extinction at the endpoint and estimate its source integral.
// Subtract the last flushed S from both terms to reduce cancellation on long
// rays. Keep the original fit anchor across flushes and invisible boundaries.
WeightedAtmosphereSegment AtmosphereRay::Flush(double endpoint_uniform, double correction_uniform) {
  WeightedAtmosphereSegment result;
  if (!pending) return result;
  auto endpoint = model->SampledSegment(origin, direction, pending_distance, cache, endpoint_uniform);
  for (int c = 0; c < 3; ++c) {
    result.transmission[c] = cumulative.transmission[c] > 0
        ? endpoint.transmission[c] / cumulative.transmission[c] : 0;
    result.radiance[c] = pending_weight[c] * (double(endpoint.radiance[c]) - cumulative.radiance[c]);
  }
  if (reservoir_mass > 0) {
    auto selected = model->SampledSegment(origin, direction, selected_distance, cache, correction_uniform);
    for (int c = 0; c < 3; ++c)
      result.radiance[c] += selected_weight[c] * (reservoir_mass / selected_mass) *
          (double(selected.radiance[c]) - cumulative.radiance[c]);
  }
  cumulative = endpoint; distance = pending_distance;
  pending = false; reservoir_mass = 0;
  return result;
}


// The environment is evaluated at the fixed anchor and already includes haze.
// Remove the accumulated part so adding the remaining background reconstructs
// that environment exactly, instead of applying the nearby atmosphere twice.
point3f AtmosphereRay::Remaining(const point3f &radiance) const {
  if (!active) return radiance;
  point3f result;
  for (int c = 0; c < 3; ++c)
    result[c] = cumulative.transmission[c] > 0
                  ? (radiance[c] - cumulative.radiance[c]) / cumulative.transmission[c] : 0;
  return result;
}


// Prepare physical settings, spectral conversion, and optional sampling tables
// once. Rendering queries below only read this state and the shared model.
PragueInfiniteLight::PragueInfiniteLight(const Rcpp::List &description, bool build_sampler) {
  // Validate model-domain limits and establish the scene-to-meters conversion.
  meters_per_unit = scalar(description, "meters_per_unit", std::numeric_limits<double>::min(), 1e12);
  altitude = scalar(description, "altitude", 0, 15000);
  visibility = scalar(description, "visibility", 20, 131.8);
  albedo = scalar(description, "albedo", 0, 1);
  elevation = scalar(description, "elevation", -4.2, 90) * pi / 180;
  azimuth = scalar(description, "azimuth", 0, 360) * pi / 180;
  intensity = scalar(description, "intensity", 0, 1e30);
  attenuation = flag(description, "attenuation");
  query_altitude = flag(description, "query_altitude");
  haze_in_volumes = flag(description, "haze_in_volumes");
  deferred_haze = flag(description, "deferred_haze", false);
  if (description.containsElementNamed("haze_correction_probability"))
    haze_correction_probability = scalar(description, "haze_correction_probability", 0, 1);
  if (!(haze_correction_probability > 0))
    throw std::runtime_error("haze_correction_probability must be greater than zero.");
  if (haze_correction_probability < 1 && !deferred_haze)
    throw std::runtime_error("haze_correction_probability < 1 requires deferred_haze = TRUE.");
  if (attenuation && !query_altitude)
    throw std::runtime_error("attenuation = TRUE requires query_altitude = TRUE.");


  // This budget controls optional storage only. Zero skips expansion, and
  // positive infinity requests the complete table without a user memory cap.
  double transmission_table_max_mb = 512;
  if (description.containsElementNamed("transmission_table_max_mb")) {
    SEXP value = description["transmission_table_max_mb"];
    if ((TYPEOF(value) != REALSXP && TYPEOF(value) != INTSXP) || Rf_xlength(value) != 1)
      throw std::runtime_error("transmission_table_max_mb must be a nonnegative number or Inf.");
    transmission_table_max_mb = Rcpp::as<double>(value);
    if (std::isnan(transmission_table_max_mb) || transmission_table_max_mb < 0)
      throw std::runtime_error("transmission_table_max_mb must be a nonnegative number or Inf.");
  }


  // Keep the light's rotation separate from the enclosing environment transform.
  // The origin anchors the local Earth frame; gain calibrates the final RGB.
  double rotation = scalar(description, "rotation", -1e30, 1e30);
  light_to_environment = RotateY(-std::fmod(rotation, 360.0));
  environment_to_light = Inverse(light_to_environment);
  origin = triple(description, "origin"); gain = triple(description, "rgb_gain", true);


  // The enabled components determine both emitted radiance and Sun/sky sampling
  // probabilities. An explicit Sun light can disable this built-in solar disk.
  sun_radius = scalar(description, "angular_diameter", .01, 2) * pi / 360;
  sun_solid_angle = 4 * pi * std::pow(std::sin(sun_radius / 2), 2);
  std::string mode = Rcpp::as<std::string>(description["render_mode"]);
  if (mode != "all" && mode != "sun" && mode != "atmosphere")
    throw std::runtime_error("Unknown Prague sky render mode.");
  include_sky = mode != "sun"; include_sun = mode != "atmosphere";
  if (description.containsElementNamed("include_sun"))
    include_sun = include_sun && Rcpp::as<bool>(description["include_sun"]);
  sun_fraction = include_sun ? (include_sky ? .5 : 1) : 0;


  // Integrate at the dataset's visible channel centres. Interpolate the 1 nm CIE
  // tables, then convert XYZ matching functions to linear RGB quadrature weights.
  model = load_model(Rcpp::as<std::string>(description["filename"]), visibility,
                     flag(description, "cache_spectra"), flag(description, "transmission_table"),
                     transmission_table_max_mb);
  auto data = model->getAvailableData();
  std::array<double, 3> totals{0, 0, 0};
  for (int i = 0; i < data.channels; ++i) {
    double wavelength = data.channelStart + (i + .5) * data.channelWidth;
    if (wavelength < 380 || wavelength > 740) continue;
    wavelengths.push_back(wavelength);
    double index = wavelength - 360;
    int j = int(index); double f = index - j;
    double x = (1 - f) * CIE_X[j] + f * CIE_X[j + 1];
    double y = (1 - f) * CIE_Y[j] + f * CIE_Y[j + 1];
    double z = (1 - f) * CIE_Z[j] + f * CIE_Z[j + 1];
    std::array<double, 3> rgb{
      3.2404542 * x - 1.5371385 * y - .4985314 * z,
      -.9692660 * x + 1.8760108 * y + .0415560 * z,
      .0556434 * x - .2040259 * y + 1.0572252 * z};

    // Radiance keeps signed RGB matching weights. Generic RGB transmission uses
    // separate nonnegative weights so attenuation stays a bounded average.
    std::array<double, 3> transmission;
    for (int c = 0; c < 3; ++c) {
      rgb[c] *= data.channelWidth;
      transmission[c] = std::max(0.0, rgb[c]);
      totals[c] += transmission[c];
    }
    rgb_weights.push_back(rgb); transmission_weights.push_back(transmission);
  }

  // Match the fixed-capacity spectral buffers and normalize each transmission
  // channel: a spectrum with unit transmission must also produce unit RGB.
  if (wavelengths.empty() || wavelengths.size() > 16)
    throw std::runtime_error("Unsupported spectral channels in Prague dataset.");
  for (auto &weight : transmission_weights)
    for (int c = 0; c < 3; ++c) weight[c] /= totals[c];


  // Celestial images carry known reference spectra. Integrate their filtering
  // spectrally before forming RGB ratios, so a textured Sun retains the native
  // solar color at sunset instead of using generic broadband RGB extinction.
  auto solar = Parameters(V(0, 0, altitude), V(0, 0, 1));
  solar.gamma = 0;
  // The Sun uses Prague's unattenuated spectrum; the Moon uses its 5778 K
  // reflected-light reference. Normalization removes absolute source brightness.
  for (int body = 0; body < 2; ++body) {
    std::array<double, 3> total{};
    for (size_t i = 0; i < wavelengths.size(); ++i) {
      double lambda = wavelengths[i];
      double radiance = body == 0 ? model->sunRadiance(solar, lambda, false)
                                 : std::pow(550.0 / lambda, 5) / std::expm1(1.438776877e7 / (lambda * 5778));
      auto weight = rgb_weights[i];
      for (int c = 0; c < 3; ++c) { weight[c] *= radiance; total[c] += weight[c]; }
      celestial_weights[body].push_back(weight);
    }
    for (auto &weight : celestial_weights[body])
      for (int c = 0; c < 3; ++c) weight[c] /= total[c];
  }


  // Transport-only callers do not need a light sampler. Otherwise, build sky
  // proposals at several heights; actual radiance still comes from native queries
  // at the shading point, independently of this table's angular resolution.
  if (!build_sampler) return;
  int height = int(scalar(description, "resolution", 16, 2048)), width = 2 * height;
  proposal_altitudes = {0, 100, 500, 1500, 4000, 8000, 15000};
  // A fixed observer needs only one proposal; per-position lighting brackets
  // the actual height with neighbouring maps without interpolating radiance.
  if (!query_altitude) proposal_altitudes = {altitude};
  for (double h : proposal_altitudes) {
    std::vector<Float> values(size_t(width) * height);
    double total = 0;
    // Weight luminance by sin(theta) to account for the solid angle represented
    // by each latitude-longitude cell. The solar disk is sampled separately.
    for (int y = 0; y < height; ++y) {
      Rcpp::checkUserInterrupt();
      for (int x = 0; x < width; ++x) {
        vec2f uv((x + .5) / width, (y + .5) / height);
        auto rgb = RGB(Spectrum(V(0, 0, h), from_uv(uv), false, include_sky));
        double value = std::max(0.0, .212671 * rgb[0] + .715160 * rgb[1] + .072169 * rgb[2]) *
                       std::sin(pi * uv[1]);
        values[size_t(y) * width + x] = Float(value);
        total += value;
      }
    }

    // Estimate integrated sky brightness for choosing among scene lights, then
    // normalize this height's proposal. Empty maps receive a usable fallback.
    sampling_weight = std::max(sampling_weight, total * 2 * pi * pi / values.size());
    if (total == 0) std::fill(values.begin(), values.end(), Float(1));
    else for (auto &v : values) v /= total;
    proposals.push_back(std::make_unique<Distribution2D>(values.data(), width, height));
  }


  // Add a separate solar contribution to the light-selection weight, using a
  // high-altitude estimate where extinction does not suppress the disk as much.
  V sun(std::cos(azimuth) * std::cos(elevation), std::sin(azimuth) * std::cos(elevation),
        std::sin(elevation));
  auto rgb = RGB(Spectrum(V(0, 0, query_altitude ? 15000 : altitude), sun, include_sun, false));
  sampling_weight += std::max(0.0, .212671 * rgb[0] + .715160 * rgb[1] + .072169 * rgb[2]) *
                     sun_solid_angle;
}


// Convert between the renderer's Y-up frame and Prague's north/east/up axes.
// Positions also use the configured origin, metre scale, and reference altitude;
// directions only rotate. World() reverses the direction mapping for samples.
PragueInfiniteLight::Vector PragueInfiniteLight::Position(const point3f &p) const {
  vec3f offset(p[0] - origin[0], p[1] - origin[1], p[2] - origin[2]);
  vec3f local = environment_to_light(EnvironmentDirection(offset));
  return V(local[2] * meters_per_unit, -local[0] * meters_per_unit,
           local[1] * meters_per_unit + altitude);
}


// Freeze the entire local observer frame in fixed-altitude mode: its horizon,
// solar visibility, radiance, and sampling must all describe the same observer.
PragueInfiniteLight::Vector PragueInfiniteLight::LightingPosition(const point3f &p) const {
  return query_altitude ? Position(p) : V(0, 0, altitude);
}


PragueInfiniteLight::Vector PragueInfiniteLight::Direction(const vec3f &w) const {
  vec3f local = environment_to_light(EnvironmentDirection(w));
  return normalized(V(local[2], -local[0], local[1]));
}


vec3f PragueInfiniteLight::World(const Vector &w) const {
  return WorldDirection(light_to_environment(vec3f(-w.y, w.z, w.x)));
}


// Earth is centred one radius below the local origin. Radial height accounts for
// curvature when a ray travels far from that origin, unlike a flat Z coordinate.
double PragueInfiniteLight::Height(Vector p) const {
  return length(p + V(0, 0, earth_radius)) - earth_radius;
}


// Keep queries within the fitted altitude range by moving along the local radial
// direction. Sun position, visibility, and ground albedo remain fixed for the sky.
PragueInfiniteLight::Model::Parameters PragueInfiniteLight::Parameters(Vector p, Vector w) const {
  V radial = p + V(0, 0, earth_radius);
  double h = Height(p);
  if (h < 0 || h > 15000)
    p = normalized(radial) * (earth_radius + std::clamp(h, 0.0, 15000.0)) - V(0, 0, earth_radius);
  return model->computeParameters(p, w, elevation, azimuth, visibility, albedo);
}


// Test the Earth sphere using the same altitude clamp as the radiance queries.
// This lets elevated observers see below their horizontal plane while still
// blocking light whose path actually enters the planet.
bool PragueInfiniteLight::PlanetOccludes(Vector p, Vector w, double distance) const {
  V radial = p + V(0, 0, earth_radius);
  double h = Height(p);
  if (h < 0 || h > 15000)
    radial = normalized(radial) * (earth_radius + std::clamp(h, 0.0, 15000.0));


  // Solve t^2 + 2*b*t + c = 0 for a unit direction. Outward and tangent rays
  // do not enter the sphere; the rationalized near root avoids cancellation
  // when the observer is very close to its surface.
  double b = dotv(radial, w);
  if (b >= 0) return false;
  double c = dotv(radial, radial) - earth_radius * earth_radius;
  double disc = b * b - c;
  if (disc <= 0) return false;
  double t = c / (-b + std::sqrt(disc));
  return t < distance && t >= -1e-5;
}


// Cheap, conservative support test for an extended Sun/Moon disk. Include the
// depressed horizon at altitude and the disk radius, so a visible upper limb
// can still contribute after the disk's centre has passed below the horizon.
bool PragueInfiniteLight::MaySeeDisk(const point3f &p, const vec3f &w, double radius) const {
  V position = LightingPosition(p);
  V up = normalized(position + V(0, 0, earth_radius));
  double h = std::clamp(Height(position), 0.0, 15000.0);
  double horizon = -std::acos(earth_radius / (earth_radius + h));
  // Conservative support includes Float rounding at a grazing upper limb.
  return dotv(Direction(w), up) >= std::sin(horizon - radius - 2e-6);
}


// Evaluate the enabled radiance components spectrally before any RGB conversion.
// Finite haze segments request only the sky; direct solar emission is separate.
PragueInfiniteLight::SpectrumValues PragueInfiniteLight::Spectrum(Vector p, Vector w,
                                                                 bool sun, bool sky, bool smooth) const {
  SpectrumValues result{};
  auto params = Parameters(p, w);

  // Map the requested apparent disk size onto Prague's native solar profile,
  // and reject each direction individually when Earth hides it.
  if (sun && params.gamma <= sun_radius && !PlanetOccludes(p, w, INFINITY)) {
    auto solar = params;
    solar.gamma *= native_sun_radius / sun_radius;
    for (size_t i = 0; i < wavelengths.size(); ++i)
      result[i] = model->sunRadiance(solar, wavelengths[i]);
  }
  if (!sky) return result;


  // Section 5.4 of Wilkie et al. recommends averaging nearby vertical directions
  // for finite-distance subtraction near the horizon. Use the same smooth field
  // for the dome and segment endpoints, preserving their transport relationship.
  double near_horizon = std::clamp((5 * pi / 180 - std::abs(params.theta - pi / 2)) /
                                   (pi / 180), 0.0, 1.0);
  if (!smooth) near_horizon = 0;
  V up = normalized(p + V(0, 0, earth_radius));
  V vertical = up - w * dotv(up, w);
  for (int j = -1; j <= 1; ++j) {
    double weight = j == 0 ? 1 - .5 * near_horizon : .25 * near_horizon;
    if (weight == 0) continue;
    double angle = j * .2 * pi / 180;
    auto sample = j == 0 ? params : Parameters(p, normalized(w * std::cos(angle) +
                                                    normalized(vertical) * std::sin(angle)));
    SpectrumValues values{};
    model->skyRadianceSpectrum(sample, wavelengths.data(), wavelengths.size(), values.data());
    for (size_t i = 0; i < wavelengths.size(); ++i)
      result[i] += weight * values[i];
  }
  return result;
}


// Use the same spectral-to-RGB conversion for the environment and finite haze.
// Intensity and channel gains affect emitted radiance, not transmission ratios.
point3f PragueInfiniteLight::RGB(const SpectrumValues &values) const {
  double rgb[3] = {0, 0, 0};
  for (size_t i = 0; i < wavelengths.size(); ++i)
    for (int c = 0; c < 3; ++c) rgb[c] += values[i] * rgb_weights[i][c];
  for (int c = 0; c < 3; ++c) {
    rgb[c] *= intensity * gain[c];
    if (!std::isfinite(rgb[c])) throw std::runtime_error("Nonfinite Prague atmosphere radiance.");
  }
  return point3f(rgb[0], rgb[1], rgb[2]);
}


// Environment queries can include the built-in Sun. The sky-only entry point
// also supports retaining atmospheric haze when the background is transparent.
point3f PragueInfiniteLight::Radiance(const point3f &p, const vec3f &w, Float) const {
  return RGB(Spectrum(LightingPosition(p), Direction(w), include_sun, include_sky));
}


point3f PragueInfiniteLight::SkyRadiance(const point3f &p, const vec3f &w) const {
  return RGB(Spectrum(LightingPosition(p), Direction(w), false, include_sky));
}


// Generic RGB attenuation for a path of the requested length. A zero-length path
// is transparent; an infinite path that intersects Earth is completely blocked.
point3f PragueInfiniteLight::Transmission(const point3f &p, const vec3f &w, double distance) const {
  if (distance <= 0 || (!attenuation && std::isfinite(distance))) return point3f(1);
  // Disabling finite haze does not remove extinction from space to the light's
  // observer. Additional environment images still need that infinite filtering.
  V position = std::isfinite(distance) ? Position(p) : LightingPosition(p), direction = Direction(w);
  if (!std::isfinite(distance) && PlanetOccludes(position, direction, INFINITY)) return point3f(0);
  auto params = Parameters(position, direction);
  double d = std::isfinite(distance) ? distance * meters_per_unit : std::numeric_limits<double>::max();


  // The compressed fit is not exactly one at zero distance. Normalize that
  // endpoint and interpolate optical depth over the first 100 m, where raw fit
  // errors can exceed the attenuation itself. Longer paths use the finite-
  // distance Prague table directly. A fixed ray anchor makes this independent
  // of cloud majorants and invisible boundary counts.
  SpectrumValues zero{}, values{};
  model->transmittanceSpectrum(params, wavelengths.data(), wavelengths.size(), 1e-6, zero.data());
  model->transmittanceSpectrum(params, wavelengths.data(), wavelengths.size(), std::max(100.0, d), values.data());


  // Average attenuation with the nonnegative channel weights prepared above.
  // Each RGB channel therefore stays between opaque and fully transmitting.
  double tr[3] = {0, 0, 0};
  for (size_t i = 0; i < wavelengths.size(); ++i) {
    double value = normalized_transmission(zero[i], values[i], d);
    for (int c = 0; c < 3; ++c) tr[c] += value * transmission_weights[i][c];
  }
  return point3f(std::clamp(tr[0], 0.0, 1.0), std::clamp(tr[1], 0.0, 1.0), std::clamp(tr[2], 0.0, 1.0));
}


// Explicit celestial textures arrive without atmospheric filtering. Their known
// reference spectra let us compute source-specific RGB attenuation to infinity,
// preserving solar/lunar color changes that generic RGB weights would miss.
point3f PragueInfiniteLight::CelestialTransmission(const point3f &p, const vec3f &w,
                                                  InfiniteLightSpectrum spectrum) const {
  V position = LightingPosition(p), direction = Direction(w);
  if (PlanetOccludes(position, direction, INFINITY)) return point3f(0);
  SpectrumValues values{};
  model->transmittanceSpectrum(Parameters(position, direction), wavelengths.data(), wavelengths.size(),
                               std::numeric_limits<double>::max(), values.data());


  // These weights were normalized by the unattenuated source's RGB, so the
  // weighted sums are channel ratios to apply directly to its texture radiance.
  const auto &weights = celestial_weights[spectrum == InfiniteLightSpectrum::Sun ? 0 : 1];
  double rgb[3] = {0, 0, 0};
  for (size_t i = 0; i < wavelengths.size(); ++i)
    for (int c = 0; c < 3; ++c) rgb[c] += weights[i][c] * values[i];
  return point3f(rgb[0], rgb[1], rgb[2]);
}


// Standalone segment queries need a temporary cache. AtmosphereRay instead keeps
// one cache for repeated cumulative queries from an unchanged origin/direction.
AtmosphereSegment PragueInfiniteLight::Segment(const point3f &p, const vec3f &w, double distance) const {
  AtmosphereSegmentCache cache;
  return Segment(p, w, distance, cache);
}


AtmosphereSegment PragueInfiniteLight::Segment(const point3f &p, const vec3f &w, double distance,
                                              AtmosphereSegmentCache &cache) const {
  return EvaluateSegment(p, w, distance, cache, -1);
}

AtmosphereSegment PragueInfiniteLight::SampledSegment(const point3f &p, const vec3f &w,
    double distance, AtmosphereSegmentCache &cache, double uniform) const {
  return EvaluateSegment(p, w, distance, cache, uniform);
}

AtmosphereSegment PragueInfiniteLight::EvaluateSegment(const point3f &p, const vec3f &w,
    double distance, AtmosphereSegmentCache &cache, double uniform) const {
  AtmosphereSegment result;
  if (!attenuation || distance <= 0) return result;
  V position = Position(p), direction = Direction(w);
  double d = distance * meters_per_unit;
  if (!std::isfinite(d)) throw std::runtime_error("Finite atmosphere segment length required.");
  auto params = Parameters(position, direction);


  // Cache the anchor's sky spectrum and zero-distance fit correction. Reusing
  // these endpoints makes repeated subdivisions consistent and avoids queries.
  if (!cache.ready) {
    cache.radiance = Spectrum(position, direction, false, include_sky);
    model->transmittanceSpectrum(params, wavelengths.data(), wavelengths.size(), 1e-6,
                                 cache.zero_transmission.data());
    cache.ready = true;
  }

  // Short paths share a 100 m reference and interpolate its optical depth;
  // longer paths obtain their distance-specific transmission directly below.
  if (d < 100 && !cache.near_ready) {
    model->transmittanceSpectrum(params, wavelengths.data(), wavelengths.size(), 100,
                                 cache.near_transmission.data());
    cache.near_ready = true;
  }


  // Recover the light scattered between endpoints a and b from the common sky
  // field: L_ab = L_a - T_ab * L_b. Subtract spectrally before converting to RGB
  // so wavelength-dependent extinction colors the finite haze consistently.
  // Use the native central direction as a cheap control. Sample only the
  // correction for the two neighboring vertical directions used in smoothing.
  // The correction must be formed AFTER the original spectral clamp, otherwise
  // the nonlinear clamp would bias the sampled source estimate.
  V endpoint = position + direction * d;
  bool sample_correction = false;
  if (include_sky && uniform >= 0 && haze_correction_probability < 1) {
    auto endpoint_params = Parameters(endpoint, direction);
    sample_correction = std::abs(endpoint_params.theta - pi / 2) < 5 * pi / 180;
  }
  bool selected = sample_correction && uniform < haze_correction_probability;
  auto b = Spectrum(endpoint, direction, false, include_sky, !sample_correction);
  SpectrumValues full{};
  if (selected) full = Spectrum(endpoint, direction, false, include_sky);
  SpectrumValues source{};
  SpectrumValues values = cache.near_transmission;
  if (d >= 100)
    model->transmittanceSpectrum(params, wavelengths.data(), wavelengths.size(), d, values.data());
  double tr[3] = {0, 0, 0};
  for (size_t i = 0; i < wavelengths.size(); ++i) {
    double value = normalized_transmission(cache.zero_transmission[i], values[i], d);
    // The fitted data can have small negative residuals; radiance cannot be
    // negative spectrally. RGB conversion follows after this spectral subtraction.
    source[i] = std::max(0.0, cache.radiance[i] - value * b[i]);
    if (selected) {
      double target = std::max(0.0, cache.radiance[i] - value * full[i]);
      source[i] += (target - source[i]) / haze_correction_probability;
    }
    for (int c = 0; c < 3; ++c) tr[c] += value * transmission_weights[i][c];
  }


  // Return dimensionless transmission alongside emitted RGB radiance. Only the
  // latter receives the sky's intensity and color calibration through RGB().
  result.transmission = point3f(std::clamp(tr[0], 0.0, 1.0), std::clamp(tr[1], 0.0, 1.0),
                                std::clamp(tr[2], 0.0, 1.0));
  result.radiance = RGB(source);
  return result;
}


// Bracket the observer's height with two precomputed sampling maps. The fraction
// selects a mixture of their distributions, not an interpolation of sky radiance.
std::pair<size_t, double> PragueInfiniteLight::Proposal(Vector p) const {
  if (proposal_altitudes.size() == 1) return {0, 0};
  double h = std::clamp(Height(p), 0.0, 15000.0);
  auto upper = std::upper_bound(proposal_altitudes.begin(), proposal_altitudes.end(), h);
  size_t i = std::min(proposal_altitudes.size() - 2, size_t(upper - proposal_altitudes.begin() - 1));
  return {i, (h - proposal_altitudes[i]) / (proposal_altitudes[i + 1] - proposal_altitudes[i])};
}


// Convert the height-mixture density from UV area to solid angle using the
// latitude-longitude Jacobian, 2*pi^2*sin(theta). The small uniform component
// keeps support everywhere, including directions the coarse maps assign zero.
Float PragueInfiniteLight::SkyPdf(Vector p, Vector w) const {
  if (proposals.empty()) return 1 / (4 * pi);
  auto mix = Proposal(p);
  vec2f uv = to_uv(w);
  double sine = std::sqrt(std::max(0.0, 1 - w.z * w.z));
  double density = proposals[mix.first]->Pdf(uv);
  if (mix.second > 0)
    density = (1 - mix.second) * density + mix.second * proposals[mix.first + 1]->Pdf(uv);
  double pdf = sine > 0 ? density / (2 * pi * pi * sine) : 0;
  return Float(uniform_fraction / (4 * pi) + (1 - uniform_fraction) * pdf);
}


vec3f PragueInfiniteLight::Sample(const point3f &p, vec2f u, Float) const {
  // Keep variates in [0,1) so each mixture branch can remap its conditional
  // interval without landing on an excluded endpoint or dividing by zero.
  u.xy.x = std::clamp(u[0], Float(0), std::nextafter(Float(1), Float(0)));
  u.xy.y = std::clamp(u[1], Float(0), std::nextafter(Float(1), Float(0)));


  // Sample the small solar cone explicitly and uniformly in solid angle.
  // Form delta = 1-cos(theta) through sin^2 to avoid cancellation at small angles.
  if (u[0] < sun_fraction) {
    double delta = (u[0] / sun_fraction) * 2 * std::pow(std::sin(sun_radius / 2), 2);
    V sun(std::cos(azimuth) * std::cos(elevation), std::sin(azimuth) * std::cos(elevation), std::sin(elevation));
    V right = normalized(crossv(sun, std::abs(sun.z) < .99 ? V(0, 0, 1) : V(1, 0, 0)));
    V up = crossv(right, sun);
    double phi = 2 * pi * u[1];
    V radial = right * std::cos(phi) + up * std::sin(phi);
    vec3f result = World(sun * (1 - delta) + radial * std::sqrt(delta * (2 - delta)));

    // Converting a tiny cone to renderer Float coordinates can round a rim
    // sample just outside its support. Move only such rounded samples inward.
    if (dotv(Direction(result), sun) < std::cos(sun_radius)) {
      double angle = std::max(0.0, sun_radius - 4 * std::numeric_limits<Float>::epsilon());
      result = World(sun * std::cos(angle) + radial * std::sin(angle));
    }
    return result;
  }


  // After choosing the sky branch, reserve a small probability for a uniform
  // sphere sample. This also supplies the fallback when no maps were built.
  u.xy.x = (u[0] - sun_fraction) / (1 - sun_fraction);
  if (proposals.empty()) return World(uniform(u));
  if (u[0] < uniform_fraction) {
    u.xy.x /= uniform_fraction;
    return World(uniform(u));
  }


  // Choose one of the two neighbouring altitude maps and remap the same variate
  // into that map's unit interval. SkyPdf() evaluates the full mixture density.
  u.xy.x = (u[0] - uniform_fraction) / (1 - uniform_fraction);
  auto mix = Proposal(LightingPosition(p));
  size_t index = mix.first;
  if (mix.second > 0 && u[0] < mix.second) { ++index; u.xy.x /= mix.second; }
  else u.xy.x = (u[0] - mix.second) / (1 - mix.second);
  Float unused;
  return World(from_uv(proposals[index]->SampleContinuous(u, &unused)));
}


// Match all branches of Sample(). A direction inside the solar disk can also be
// selected by the sky proposal, so both weighted densities contribute there.
Float PragueInfiniteLight::Pdf(const point3f &p, const vec3f &w, Float) const {
  V d = Direction(w);
  V sun(std::cos(azimuth) * std::cos(elevation), std::sin(azimuth) * std::cos(elevation), std::sin(elevation));
  double solar = dotv(d, sun) >= std::cos(sun_radius) ? sun_fraction / sun_solid_angle : 0;
  return Float(solar + (1 - sun_fraction) * SkyPdf(LightingPosition(p), d));
}


// Include coefficient storage and auxiliary sampling tables in the light's
// memory estimate, even though the loaded coefficient slice may be shared.
size_t PragueInfiniteLight::GetSize() const {
  size_t size = sizeof(*this) + model->memoryUsage();
  for (const auto &weights : celestial_weights) size += weights.capacity() * sizeof(weights[0]);
  for (const auto &proposal : proposals) size += proposal->GetSize();
  return size;
}
