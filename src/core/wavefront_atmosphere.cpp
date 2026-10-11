#include "../lights/atmosphere.h"
#include "../materials/material.h"
#include "../volumes/boundary.h"
#include "../volumes/haze.h"
#include "wavefront.h"
#include <RcppThread.h>

namespace {
point3f point(WFVector p) { return {p.x, p.y, p.z}; }
vec3f direction(WFVector p) { return {p.x, p.y, p.z}; }
double dot3(WFVector a, WFVector b) { return double(a.x) * b.x + double(a.y) * b.y + double(a.z) * b.z; }
uint32_t mix_bits(uint32_t x) {
  x ^= x >> 16;
  x *= 0x7feb352du;
  x ^= x >> 15;
  x *= 0x846ca68bu;
  return x ^ (x >> 16);
}
double uniform(uint32_t pixel, const WFParameters &p) {
  return (double(mix_bits(pixel ^ mix_bits(p.sample) ^ mix_bits(p.depth) ^ p.seed)) + .5) * 0x1p-32;
}
double proposal_pdf(const WFLight &l, WFVector w, const WavefrontScene &scene) {
  constexpr double pi = 3.14159265358979323846;
  if (l.flags & 128)
    return (1 - l.forward.w) / (4 * pi) +
           (dot3(w, l.forward) >= 1 - l.right.w ? l.forward.w / (2 * pi * l.right.w) : 0);
  if (l.type == 1)
    return 1 / (2 * pi * l.forward.w);
  if (l.alias_offset == UINT32_MAX)
    return 1 / (4 * pi);
  double x = dot3(w, l.forward), y = dot3(w, l.right), z = dot3(w, l.up);
  double u = .5 - std::atan2(-x, -z) / (2 * pi);
  u -= std::floor(u);
  double v = .5 + std::asin(std::clamp(y, -1., 1.)) / pi;
  uint32_t column = std::min(uint32_t(u * l.alias_width), l.alias_width - 1);
  uint32_t row = std::min(uint32_t(v * l.alias_height), l.alias_height - 1);
  return scene.texels[l.alias_offset + row * l.alias_width + column].z;
}
// Match the GPU's balance heuristic for both source terms and escaped rays.
// Mixing this with the power heuristic would overcount indirect sky lighting.
double mis_weight(double a, double b) { return a / (a + b); }
Transform medium_inverse(const WFMedium &medium) {
  Float matrix[4][4]{};
  matrix[3][3] = 1;
  for (int a = 0; a < 3; ++a) {
    const auto &row = medium.inverse[a];
    matrix[a][0] = row.x;
    matrix[a][1] = row.y;
    matrix[a][2] = row.z;
    matrix[a][3] = row.w;
  }
  return Transform(matrix);
}
} // namespace

struct WavefrontAtmosphere::Data {
  const WavefrontScene &scene;
  RcppThread::ThreadPool pool;
  std::vector<AtmosphereRay> rays;
  std::vector<Transform> inverses;
  explicit Data(const WavefrontScene &s, size_t capacity)
      : scene(s), pool(std::max(size_t(1), s.workers)), rays(capacity, AtmosphereRay(s.transport_atmosphere)) {
    inverses.emplace_back();
    for (size_t i = 1; i < s.media.size(); ++i)
      inverses.push_back(medium_inverse(s.media[i]));
  }
  // The atmosphere's native spectral attenuation is applied to external image,
  // solar and lunar emitters exactly as InfiniteLightMixture does on the CPU.
  point3f Radiance(size_t index, const point3f &p, const vec3f &w) const {
    const auto *source = scene.light_sources[index];
    point3f value = source->Radiance(p, w, 0);
    if (!source->GetAtmosphere() && (value[0] != 0 || value[1] != 0 || value[2] != 0))
      value *= scene.atmosphere->CelestialTransmission(p, w, source->RadianceSpectrum());
    return value;
  }
  point3f ShadowRadiance(const WFShadow &shadow, const WFParameters &parameters) const {
    point3f origin = point(shadow.origin);
    vec3f w = direction(shadow.direction);
    uint32_t light = shadow.padding[1];
    bool infinite = light < scene.environments;
    if (!scene.transport_atmosphere)
      return infinite ? Radiance(light, origin, w) : point3f(1);

    // A sky connection estimates source radiance all along its direction, even
    // before an opaque blocker. Endpoint visibility times sky radiance loses
    // this term, especially inside fog. Evaluate ordered native shadow transport
    // in this hybrid batch; ordinary GPU scenes retain their device shadow walk.
    bool sky_source = infinite && scene.light_sources[light]->GetAtmosphere();
    Float time = parameters.camera.shutter.x + (parameters.camera.shutter.y-parameters.camera.shutter.x) *
                 std::max(0.f,parameters.camera.shutter.w);
    Ray ray(origin, w, time);
    VolumePathState state;
    const auto *volumes = scene.volume_scene;
    if (volumes && volumes->boundary_bvh)
      state = volumes->InitialState(ray, nullptr);
    random_gen rng(
        mix_bits(shadow.pixel ^ parameters.seed ^ mix_bits(parameters.sample) ^ (parameters.depth * 0x9e3779b9u)));
    AtmosphereRay air(scene.transport_atmosphere);
    point3f transmission(1), radiance(0);
    double previous = 0;
    while (previous < shadow.direction.w) {
      state.SetRay(ray);
      hit_record hit;
      bool found =
          scene.world->hit(ray, std::nextafter(Float(previous), Float(INFINITY)), shadow.direction.w, hit, rng);
      bool escaped = !found || hit.infinite_area_hit;
      double end = escaped ? shadow.direction.w : hit.t;
      double distance = end - previous;
      const auto *entry = state.Active();
      const Medium *medium = entry ? entry->boundary->medium.get() : nullptr;
      if (medium && medium->subsurface_diffusion)
        return radiance;
      bool enabled =
          state.glass.empty() && (!medium || (medium->haze && scene.transport_atmosphere->IntegrateInVolumes()));
      Ray span(ray(Float(previous)), w, time);
      Ray local = span;
      if (entry) {
        Transform inverse = Inverse(entry->medium_to_world);
        local = Ray(inverse(span.o), inverse(w), time);
      }
      MediumHazeIterator intervals(medium, local, distance, enabled);
      auto interval = intervals.Next();
      double cursor = 0;
      auto advance = [&](double stop) {
        while (interval && cursor < stop) {
          double next = std::min(stop, interval->t_max);
          air.Start(span(Float(cursor)), w, interval->integrate);
          auto segment = air.Advance(span(Float(next)), rng.unif_rand());
          if (sky_source)
            radiance += transmission * segment.radiance;
          transmission *= segment.transmission;
          if (const auto *glass = state.ActiveDielectric())
            for (int c = 0; c < 3; ++c)
              transmission[c] *= std::exp(-glass->attenuation[c] * (next - cursor));
          cursor = next;
          if (cursor >= interval->t_max)
            interval = intervals.Next();
        }
      };
      air.Start(span.o, w, interval && interval->integrate);
      if (medium) {
        // Scalar-majorant ratio tracking matches the GPU's RGB transmittance
        // estimator. Advance haze before changing the connection's weight, so
        // each source interval receives attenuation from the preceding fog.
        auto majorants = medium->SampleRay(local, distance);
        while (auto cell = majorants.Next()) {
          double rate = std::max({cell->sigma_maj[0], cell->sigma_maj[1], cell->sigma_maj[2]});
          if (!(rate > 0))
            continue;
          double t = cell->t_min;
          while (true) {
            double next = t - std::log1p(-rng.unif_rand()) / rate;
            if (next >= cell->t_max)
              break;
            if (!(next > t))
              throw PathFailure(PathFailureKind::Radiance, "Unresolvable atmosphere shadow flight");
            t = next;
            advance(t);
            auto properties = medium->SamplePoint(local(Float(t)));
            for (int c = 0; c < 3; ++c)
              transmission[c] *= std::max(0., 1 - (properties.sigma_a[c] + properties.sigma_s[c]) / rate);
            Float largest = std::max({transmission[0], transmission[1], transmission[2]});
            if (!(largest > 0))
              return radiance;
            if (largest < .05) {
              if (rng.unif_rand() < .5)
                return radiance;
              transmission *= 2;
            }
          }
        }
      }
      if (std::isfinite(distance))
        advance(distance);
      if (escaped) {
        point3f value = infinite ? Radiance(light, air.Origin(span.o), w) : point3f(1);
        if (infinite)
          value = sky_source ? air.Remaining(value) : point3f(air.Remaining(value) - air.Remaining(point3f(0)));
        return radiance + transmission * value;
      }
      if (!(end > previous))
        return radiance;
      previous = end;
      if (hit.alpha_miss)
        continue;
      if (hit.medium_boundary && !hit.medium_boundary->keep_surface) {
        state.Cross(hit, w);
        continue;
      }
      if (hit.mat_ptr && hit.mat_ptr->is_dielectric()) {
        auto interface = state.Interface(hit, w);
        if (interface.Hidden()) {
          state.CrossDielectric(hit, w);
          continue;
        }
        if (interface.Eta() == 1 && typeid(*hit.mat_ptr) == typeid(dielectric)) {
          if (interface.entering && (!hit.medium_boundary || !hit.medium_boundary->medium->subsurface))
            transmission *= static_cast<const dielectric *>(hit.mat_ptr)->albedo;
          state.CrossDielectric(hit, w);
          continue;
        }
      }
      return radiance;
    }
    return radiance;
  }
};

WavefrontAtmosphere::WavefrontAtmosphere(const WavefrontScene &scene, size_t capacity)
    : data(new Data(scene, capacity)) {}
WavefrontAtmosphere::~WavefrontAtmosphere() = default;
void WavefrontAtmosphere::Begin(const std::vector<uint32_t> &active) {
  if (data->scene.motion_changed)
    for (size_t i=1;i<data->scene.media.size();++i)
      data->inverses[i]=medium_inverse(data->scene.media[i]);
  for (uint32_t pixel : active)
    data->rays[pixel] = AtmosphereRay(data->scene.transport_atmosphere);
}

bool WavefrontAtmosphere::Paths(const WFParameters &p, const uint32_t *queue, size_t count,
                                WFAtmosphereRequest *requests, WFPath *paths, WFPixel *pixels,
                                const std::function<bool()> &cancelled) {
  const auto &scene = data->scene;
  for (size_t first = 0; first < count; first += 2048) {
    data->pool.parallelFor(first, std::min(count, first + 2048), [&](size_t i) {
      uint32_t pixel = queue[i];
      const auto &request = requests[pixel];
      if (!request.event)
        return;
      auto &air = data->rays[pixel];
      auto &path = paths[pixel];
      auto &film = pixels[pixel];
      point3f origin = point(request.origin);
      vec3f w = direction(request.direction);
      uint32_t medium_id = request.body ? scene.boundaries[request.body].medium : 0;
      const Medium *medium = medium_id ? scene.native_media[medium_id] : nullptr;
      bool enabled = scene.transport_atmosphere && !request.glass &&
                     (!request.body || (medium && medium->haze && scene.transport_atmosphere->IntegrateInVolumes()));
      air.Start(origin, w, enabled);
      if (request.event == 2) {
        point3f anchor = air.Origin(origin);
        point3f sum(0);
        for (uint32_t light = 0; light < scene.environments; ++light) {
          point3f value = data->Radiance(light, anchor, w);
          double weight = (request.padding & 1)
                              ? 1
                              : mis_weight(request.weight.w,
                                           proposal_pdf(scene.lights[light], request.direction, scene) / p.lights);
          if (scene.light_sources[light]->GetAtmosphere())
            value = air.Remaining(value);
          else
            value = point3f(air.Remaining(value) - air.Remaining(point3f(0)));
          sum += value * Float(weight);
        }
        for (int c = 0; c < 3; ++c)
          (&film.radiance.x)[c] += (&path.beta.x)[c] * sum[c];
        return;
      }
      Ray world(origin, w, Float(0));
      Ray local = medium ? Ray(data->inverses[medium_id](origin), data->inverses[medium_id](w), Float(0)) : world;
      MediumHazeIterator intervals(medium, local, request.direction.w, enabled);
      point3f transmission(1), radiance(0);
      while (auto interval = intervals.Next()) {
        air.Start(world(Float(interval->t_min)), w, interval->integrate);
        auto segment = air.Advance(world(Float(interval->t_max)), uniform(pixel, p));
        radiance += transmission * segment.radiance;
        transmission *= segment.transmission;
      }
      double weight = 1;
      if (!(request.padding & 1))
        for (uint32_t light = 0; light < scene.environments; ++light)
          if (scene.light_sources[light]->GetAtmosphere())
            weight =
                mis_weight(request.weight.w, proposal_pdf(scene.lights[light], request.direction, scene) / p.lights);
      for (int c = 0; c < 3; ++c) {
        float &value = (&film.radiance.x)[c];
        float previous = (&request.radiance.x)[c];
        value = previous + (value - previous) * transmission[c] + (&request.weight.x)[c] * radiance[c] * weight;
        (&path.beta.x)[c] *= transmission[c];
      }
    });
    data->pool.wait();
    if (cancelled())
      return false;
  }
  return true;
}

bool WavefrontAtmosphere::Shadows(const WFParameters &parameters, WFShadow *shadows, size_t count, WFPixel *pixels,
                                  const std::function<bool()> &cancelled) {
  for (size_t first = 0; first < count; first += 2048) {
    data->pool.parallelFor(first, std::min(count, first + 2048), [&](size_t i) {
      auto &shadow = shadows[i];
      if (shadow.contribution.x == 0 && shadow.contribution.y == 0 && shadow.contribution.z == 0)
        return;
      point3f value;
      try {
        value = data->ShadowRadiance(shadow, parameters);
      } catch (const PathFailure &) {
        shadow.padding[2] = 1;
        return;
      }
      for (int c = 0; c < 3; ++c)
        (&pixels[shadow.pixel].radiance.x)[c] += (&shadow.contribution.x)[c] * value[c];
    });
    data->pool.wait();
    if (cancelled())
      return false;
  }
  return true;
}

// RGB quadrature belongs to the renderer; skymodelr owns only the spectral fit.
// Refresh this small record on interactive rotation, retaining coefficient data.
void PragueInfiniteLight::ExportMetal(WFAtmosphereData &out) const {
  for (int a=0;a<3;++a) {
    vec3f basis(0); basis[a]=1;
    auto v=Direction(basis);
    (&out.transform[0].x)[a]=v.x;
    (&out.transform[1].x)[a]=v.y;
    (&out.transform[2].x)[a]=v.z;
  }
  out.origin={float(origin[0]),float(origin[1]),float(origin[2]),float(altitude)};
  out.gain={float(gain[0]*intensity),float(gain[1]*intensity),float(gain[2]*intensity),float(meters_per_unit)};
  out.settings={float(elevation),float(azimuth),float(visibility),float(albedo)};
  out.options={float(sun_radius),0,0,0};
  out.channels=wavelengths.size(); out.sky=include_sky; out.sun=include_sun;
  out.haze=haze; out.query_altitude=query_altitude; out.haze_filter=haze_filter; out.haze_in_volumes=haze_in_volumes;
  for(size_t i=0;i<wavelengths.size();++i) {
    auto weight=rgb_weights[i],tr=transmission_weights[i];
    out.weights[i]={float(weight[0]),float(weight[1]),float(weight[2]),float(wavelengths[i])};
    out.transmission[i]={float(tr[0]),float(tr[1]),float(tr[2]),0};
    for(int c=0;c<3;++c) {
      (&out.solar[i].x)[c]=celestial_weights[0][i][c];
      (&out.lunar[i].x)[c]=celestial_weights[1][i][c];
    }
  }
  double sum=0;
  for(size_t i=0;i<129;++i) {
    double x=-3.+6.*i/128, angle=x*.5*3.14159265358979323846/180;
    double weight=std::exp(-.5*x*x)*((i==0 || i==128) ? .5 : 1);
    out.filter[i]={float(std::cos(angle)),float(std::sin(angle)),float(weight),0};sum+=weight;
  }
  for(auto &node:out.filter) node.z/=sum;
}
