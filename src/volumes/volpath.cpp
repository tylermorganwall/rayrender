// RGB adaptation of PBRT v4's null-scattering path-space MIS estimator.
// https://pbr-book.org/4ed/Light_Transport_II_Volume_Rendering/Volume_Scattering_Integrators
#include "volpath.h"
#include "boundary.h"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace {
// A reproducible unidirectional reference for diagnosing MIS independently of
// direct-light sampling. Set before loading the package in a fresh R process.
bool use_light_sampling() {
  static const bool enabled = !(std::getenv("RAYRENDER_VOLUME_REFERENCE") &&
                                std::string(std::getenv("RAYRENDER_VOLUME_REFERENCE")) == "true");
  return enabled;
}
struct RGB {
  double v[3];
  RGB(double x = 0) : v{x, x, x} {}
  RGB(const point3f &p) : v{p[0], p[1], p[2]} {}
  double &operator[](int i) { return v[i]; }
  double operator[](int i) const { return v[i]; }
  RGB operator+(const RGB &b) const {
    RGB r;
    for (int i = 0; i < 3; ++i)
      r[i] = v[i] + b[i];
    return r;
  }
  RGB operator*(const RGB &b) const {
    RGB r;
    for (int i = 0; i < 3; ++i)
      r[i] = v[i] * b[i];
    return r;
  }
  RGB operator/(double b) const {
    RGB r;
    for (int i = 0; i < 3; ++i)
      r[i] = v[i] / b;
    return r;
  }
  RGB &operator*=(const RGB &b) {
    for (int i = 0; i < 3; ++i)
      v[i] *= b[i];
    return *this;
  }
  RGB &operator+=(const RGB &b) {
    for (int i = 0; i < 3; ++i)
      v[i] += b[i];
    return *this;
  }
  RGB &operator/=(double b) {
    for (auto &x : v)
      x /= b;
    return *this;
  }
  double Average() const { return (v[0] + v[1] + v[2]) / 3; }
  double Max() const { return std::max({v[0], v[1], v[2]}); }
  point3f FloatRGB() const { return point3f(v[0], v[1], v[2]); }
};
RGB transmittance(const point3f &sigma, double distance) {
  RGB r;
  for (int i = 0; i < 3; ++i)
    r[i] = sigma[i] == 0 ? 1 : std::exp(-double(sigma[i]) * distance);
  return r;
}
bool cancelled(const std::atomic<bool> *c) { return c && c->load(std::memory_order_relaxed); }
void rescale(RGB &beta, RGB &ru, RGB &rl) {
  double scale = std::max({beta.Max(), ru.Max(), rl.Max()});
  if (scale > 1e80 || (scale > 0 && scale < 1e-80)) {
    beta /= scale;
    ru /= scale;
    rl /= scale;
  }
}
random_gen tracking_rng(Sampler *sampler) {
  uint32_t a = uint32_t(double(sampler->Get1D()) * 4294967295.0);
  uint32_t b = uint32_t(double(sampler->Get1D()) * 4294967295.0);
  uint32_t seed = a ^ (b + 0x9e3779b9u + (a << 6) + (a >> 2));
  seed ^= seed >> 16;
  seed *= 0x7feb352du;
  seed ^= seed >> 15;
  seed *= 0x846ca68bu;
  seed ^= seed >> 16;
  return random_gen(seed);
}
RGB glass_transmittance(const VolumePathState &state, double distance) {
  const dielectric *active = nullptr;
  for (const auto *d : state.glass)
    if (!active || d->priority < active->priority)
      active = d;
  return active ? transmittance(active->attenuation, distance) : RGB(1);
}
normal3f geometric_normal(const hit_record &h) {
  return h.geometric_normal.squared_length() > 0 ? h.geometric_normal : h.normal;
}
Ray spawn(const hit_record &h, const vec3f &wi, Float time, VolumePathState &state) {
  point3f o = OffsetMediumOrigin(h, wi);
  Ray r(o, wi, time);
  state.SetRay(r);
  return r;
}
// An ordinary surface may coincide with a medium boundary. Its normal offset
// can cross that boundary before the next intersection. Resolve only such nearby
// crossings; ordinary surface spawns retain the incremental path state.
void reconcile_surface_origin(const VolumeScene *scene, const hit_record &h, Ray &ray,
                              VolumePathState &state, const std::atomic<bool> *cancel) {
  if (!scene || !scene->boundary_bvh)
    return;
  Float tolerance = 8 * std::numeric_limits<Float>::epsilon() *
                        std::max({Float(1), std::abs(h.p[0]), std::abs(h.p[1]), std::abs(h.p[2])}) +
                    h.pError.length();
  vec3f offset = ray.o - h.p;
  vec3f direction = offset.squared_length() > 0 ? unit_vector(offset) : unit_vector(ray.d);
  Ray probe(h.p - direction * (2 * tolerance), direction, ray.time());
  probe.segment_absorption = true;
  hit_record contact;
  random_gen rng(0);
  if (scene->boundary_bvh->hit(probe, 0, 4 * tolerance, contact, rng)) {
    state = scene->InitialState(ray, cancel);
    state.SetRay(ray);
  }
}
void cross_if_transmitted(VolumePathState &state, const hit_record &h, const vec3f &incoming,
                          const vec3f &outgoing) {
  if (h.medium_boundary &&
      dot(incoming, h.geometric_normal) * dot(outgoing, h.geometric_normal) > 0)
    state.Cross(h, outgoing);
}
point3f null_coefficient(const MediumProperties &mp, const point3f &majorant) {
  point3f n;
  for (int i = 0; i < 3; ++i) {
    double st = double(mp.sigma_a[i]) + mp.sigma_s[i];
    if (!std::isfinite(st) || mp.sigma_a[i] < 0 || mp.sigma_s[i] < 0 ||
        st > double(majorant[i]) * (1 + 2e-5) + 1e-12)
      throw std::runtime_error("Medium coefficients exceed their majorant or are invalid.");
    n[i] = std::max(0.0, double(majorant[i]) - st);
  }
  return n;
}
// Callback receives world distance and the majorant transmittance since its last
// invocation. The returned value is the residual transmittance to the endpoint.
template <class Callback>
RGB sample_majorant(const MediumEntry &entry, const Ray &ray, double distance, int hero,
                    random_gen &rng, Callback &&callback, const std::atomic<bool> *cancel,
                    VolumeStatistics *stats = nullptr) {
  const Medium &medium = *entry.boundary->medium;
  Transform inverse = Inverse(entry.medium_to_world);
  Ray local(inverse(ray.o), inverse(ray.d), ray.time());
  auto iterator = medium.SampleRay(local, distance);
  RGB T(1);
  while (auto seg = iterator.Next()) {
    if (stats)
      stats->segments.fetch_add(1, std::memory_order_relaxed);
    double t = seg->t_min;
    double rate = seg->sigma_maj[hero];
    if (rate == 0) {
      T *= transmittance(seg->sigma_maj, seg->t_max - t);
      continue;
    }
    while (t < seg->t_max && !cancelled(cancel)) {
      double u = std::clamp(double(rng.unif_rand()), 0.0, std::nextafter(1.0, 0.0));
      double candidate = t - std::log1p(-u) / rate;
      if (candidate <= t)
        candidate = std::nextafter(t, INFINITY);
      if (candidate >= seg->t_max) {
        T *= transmittance(seg->sigma_maj, seg->t_max - t);
        break;
      }
      T *= transmittance(seg->sigma_maj, candidate - t);
      MediumProperties mp = medium.SamplePoint(local(Float(candidate)));
      if (!callback(candidate, mp, seg->sigma_maj, T))
        return RGB(1);
      T = RGB(1);
      t = candidate;
    }
    if (cancelled(cancel))
      return RGB(1);
  }
  return T;
}
RGB segment_opacity_transmittance(const MediumEntry &entry, const Ray &ray, double distance,
                                  random_gen &rng, const std::atomic<bool> *cancel) {
  const Medium &medium = *entry.boundary->medium;
  if (medium.IsHomogeneous())
    return transmittance(medium.sigma_a + medium.sigma_s, distance);
  Transform inverse = Inverse(entry.medium_to_world);
  Ray local(inverse(ray.o), inverse(ray.d), ray.time());
  auto iterator = medium.SampleRay(local, distance);
  RGB tr(1);
  while (auto seg = iterator.Next()) {
    double rate = std::max({seg->sigma_maj[0], seg->sigma_maj[1], seg->sigma_maj[2]});
    if (rate == 0)
      continue;
    double t = seg->t_min;
    while (!cancelled(cancel)) {
      double next =
          t - std::log1p(-std::min(double(rng.unif_rand()), std::nextafter(1.0, 0.0))) / rate;
      t = next <= t ? std::nextafter(t, INFINITY) : next;
      if (t >= seg->t_max)
        break;
      auto mp = medium.SamplePoint(local(Float(t)));
      null_coefficient(mp, seg->sigma_maj);
      for (int i = 0; i < 3; ++i)
        tr[i] *= std::max(0.0, 1 - (double(mp.sigma_a[i]) + mp.sigma_s[i]) / rate);
      if (tr.Max() == 0)
        return tr;
    }
    if (cancelled(cancel))
      break;
  }
  return tr;
}
Float primary_transparency(const Ray &input, VolumePathState state, hitable *world, random_gen &rng,
                           const std::atomic<bool> *cancel) {
  Ray ray(input.o, unit_vector(input.d), input.time());
  state.SetRay(ray);
  RGB tr(1);
  while (!cancelled(cancel)) {
    hit_record h;
    if (!world->hit(ray, 0, MaxT, h, rng))
      return Float(tr.Average());
    double distance = h.t;
    tr *= glass_transmittance(state, distance);
    if (const auto *entry = state.Active())
      tr *= segment_opacity_transmittance(*entry, ray, distance, rng, cancel);
    if (h.medium_boundary && !h.medium_boundary->keep_surface) {
      state.Cross(h, ray.d);
      ray = spawn(h, ray.d, ray.time(), state);
      continue;
    }
    bool invisible = false;
    if (h.mat_ptr)
      h.mat_ptr->emitted(ray, h, h.u, h.v, h.p, invisible);
    if (h.alpha_miss || invisible) {
      ray = spawn(h, ray.d, ray.time(), state);
      continue;
    }
    return h.infinite_area_hit ? Float(tr.Average()) : 0;
  }
  return 0;
}
Float light_pdf(hitable_list *lights, const point3f &p, const vec3f &wi, random_gen &rng,
                Float time) {
  if (lights->volume_scene && lights->volume_scene->light_sampler)
    return lights->volume_scene->light_sampler->Pdf(p, wi, rng, time);
  return lights->size() ? lights->pdf_value(p, wi, rng, time) : 0;
}
struct LightSample {
  vec3f wi;
  hit_record endpoint;
  point3f radiance{0};
  double distance = 0;
  Float pdf = 0;
  const hitable *emitter = nullptr;
};
RGB direct_light(const Ray &parent, const point3f &p, const hit_record *surface, pdf *bsdf_pdf,
                 const HGPhaseFunction *phase, VolumePathState state, const RGB &beta,
                 const RGB &rp, hitable *world, hitable_list *lights, int hero, random_gen &rng,
                 Sampler *sampler, const std::atomic<bool> *cancel) {
  if (!lights->size() || !use_light_sampling())
    return RGB(0);
  VolumeStatistics *stats = lights->volume_scene && lights->volume_scene->collect_statistics
                                ? &lights->volume_scene->statistics
                                : nullptr;
  LightSample sample;
  sample.wi = lights->volume_scene && lights->volume_scene->light_sampler
                  ? lights->volume_scene->light_sampler->Sample(p, sampler, parent.time())
                  : lights->random(p, sampler, parent.time());
  if (!(sample.wi.squared_length() > 0))
    return RGB(0);
  sample.wi = unit_vector(sample.wi);
  sample.pdf = light_pdf(lights, p, sample.wi, rng, parent.time());
  if (!(sample.pdf > 0) || !std::isfinite(sample.pdf))
    return RGB(0);
  RGB f;
  Float ps;
  if (surface) {
    f = RGB(surface->mat_ptr->f(parent, *surface, sample.wi));
    ps = bsdf_pdf->value(sample.wi, rng, parent.time());
    cross_if_transmitted(state, *surface, parent.d, sample.wi);
  } else {
    ps = phase->p(-parent.d, sample.wi);
    f = RGB(ps);
  }
  if (!(f.Max() > 0))
    return RGB(0);
  Ray ray =
      surface ? spawn(*surface, sample.wi, parent.time(), state) : Ray(p, sample.wi, parent.time());
  state.SetRay(ray);
  if (surface)
    reconcile_surface_origin(lights->volume_scene.get(), *surface, ray, state, cancel);
  RGB tr(1), ru(1), rl(1);
  random_gen tracker = tracking_rng(sampler);
  while (!cancelled(cancel)) {
    hit_record h;
    if (!world->hit(ray, 0, MaxT, h, rng))
      return RGB(0);
    double distance = h.t;
    if (const auto *entry = state.Active()) {
      RGB T = sample_majorant(
          *entry, ray, distance, hero, tracker,
          [&](double, const MediumProperties &mp, const point3f &m, const RGB &T) {
            if (stats)
              stats->shadow_candidates.fetch_add(1, std::memory_order_relaxed);
            double pdf = T[hero] * m[hero];
            if (!(pdf > 0)) {
              tr = RGB(0);
              return false;
            }
            RGB n(null_coefficient(mp, m));
            tr *= T * n / pdf;
            ru *= T * n / pdf;
            rl *= T * RGB(m) / pdf;
            rescale(tr, ru, rl);
            double relative = (ru + rl).Average();
            if (relative > 0 && tr.Max() / relative < .05) {
              if (tracker.unif_rand() < .5) {
                tr = RGB(0);
                return false;
              }
              tr *= RGB(2);
            }
            return tr.Max() > 0;
          },
          cancel, stats);
      if (tr.Max() == 0)
        return RGB(0);
      if (T[hero] > 0) {
        RGB w = T / T[hero];
        tr *= w;
        ru *= w;
        rl *= w;
      }
    }
    tr *= glass_transmittance(state, distance);
    if (h.medium_boundary && !h.medium_boundary->keep_surface) {
      state.Cross(h, ray.d);
      ray = spawn(h, ray.d, ray.time(), state);
      continue;
    }
    if (h.alpha_miss) {
      ray = spawn(h, ray.d, ray.time(), state);
      continue;
    }
    bool invisible = false;
    sample.radiance = h.mat_ptr ? h.mat_ptr->emitted(ray, h, h.u, h.v, h.p, invisible) : point3f(0);
    sample.endpoint = h;
    sample.distance = (h.p - p).length();
    sample.emitter = h.shape;
    double denominator = (rl * rp * RGB(sample.pdf) + ru * rp * RGB(ps)).Average();
    return denominator > 0 ? beta * f * tr * RGB(sample.radiance) / denominator : RGB(0);
  }
  return RGB(0);
}
} // namespace

void color_volume(const Ray &input, hitable *world, hitable_list *lights, size_t max_depth,
                  size_t roulette_depth, random_gen &rng, Sampler *sampler, Float &transparency,
                  point3f &radiance, normal3f &normal, point3f &albedo,
                  const std::atomic<bool> *cancel) {
  radiance = point3f(0);
  normal = normal3f(0);
  albedo = point3f(0);
  transparency = 0;
  const VolumeScene *scene = lights->volume_scene.get();
  VolumeStatistics *stats = scene && scene->collect_statistics ? &scene->statistics : nullptr;
  if (stats)
    stats->paths.fetch_add(1, std::memory_order_relaxed);
  VolumePathState state = scene ? scene->InitialState(input, cancel) : VolumePathState{};
  if (!scene && input.pri_stack)
    state.glass = *input.pri_stack;
  int hero = std::min(2, int(sampler->Get1D() * 3));
  if (scene && scene->has_media) {
    auto alpha_rng = tracking_rng(sampler);
    transparency = primary_transparency(input, state, world, alpha_rng, cancel);
  }
  Ray ray(input.o, unit_vector(input.d), input.time());
  state.SetRay(ray);
  RGB L(0), beta(1), ru(1), rl(1);
  double eta_scale = 1;
  point3f previous_point(0);
  bool specular = true, any_diffuse = false, wrote_feature = false;
  size_t depth = 0;
  while (!cancelled(cancel)) {
    hit_record h;
    if (!world->hit(ray, 0, MaxT, h, rng))
      break;
    double distance = h.t, previous_distance = 0;
    bool scattered = false, terminated = false;
    if (const auto *entry = state.Active()) {
      auto tracker = tracking_rng(sampler);
      RGB T = sample_majorant(
          *entry, ray, distance, hero, tracker,
          [&](double t, const MediumProperties &mp, const point3f &m, const RGB &T) {
            beta *= glass_transmittance(state, t - previous_distance);
            previous_distance = t;
            point3f n = null_coefficient(mp, m);
            double proposal = T[hero] * m[hero];
            if (!(proposal > 0)) {
              terminated = true;
              return false;
            }
            if (depth < max_depth) {
              RGB re = ru * T * RGB(m) / proposal;
              if (re.Average() > 0)
                L += beta * T * RGB(mp.sigma_a) * RGB(mp.Le) / proposal / re.Average();
            }
            double mode = tracker.unif_rand(), pa = mp.sigma_a[hero] / m[hero],
                   ps = mp.sigma_s[hero] / m[hero];
            if (mode < pa) {
              terminated = true;
              return false;
            }
            if (mode < pa + ps) {
              if (depth >= max_depth) {
                terminated = true;
                return false;
              }
              ++depth;
              if (stats)
                stats->scattering_events.fetch_add(1, std::memory_order_relaxed);
              double pdf = T[hero] * mp.sigma_s[hero];
              if (!(pdf > 0)) {
                terminated = true;
                return false;
              }
              RGB w = T * RGB(mp.sigma_s) / pdf;
              beta *= w;
              ru *= w;
              // Reconstructed surface positions can lie off the ray by their
              // error bound. Keep a near-endpoint collision on the segment's
              // interior side when converting its double distance to a Float point.
              double tolerance =
                  8 * std::numeric_limits<Float>::epsilon() *
                      std::max({1.0, std::abs(double(ray.o[0])), std::abs(double(ray.o[1])),
                                std::abs(double(ray.o[2])), distance}) +
                  h.pError.length();
              double position_t = distance - t < tolerance ? std::max(0.0, t - tolerance) : t;
              point3f collision;
              for (int a = 0; a < 3; ++a)
                collision[a] = Float(double(ray.o[a]) + double(ray.d[a]) * position_t);
              MediumInteraction interaction{collision, -ray.d, ray.time(), mp};
              point3f p = interaction.p;
              L += direct_light(ray, p, nullptr, nullptr, &interaction.properties.phase, state,
                                beta, ru, world, lights, hero, rng, sampler, cancel);
              vec2f u = sampler->Get2D();
              auto sample = mp.phase.Sample(-ray.d, u.xy.x, u.xy.y);
              beta *= RGB(sample.p / sample.pdf);
              rl = ru / sample.pdf;
              previous_point = p;
              ray = Ray(p, unit_vector(sample.wi), ray.time());
              state.SetRay(ray);
              scattered = true;
              specular = false;
              any_diffuse = true;
              if (!wrote_feature) {
                normal = normal3f(0);
                albedo = point3f(0);
                wrote_feature = true;
              }
              return false;
            }
            if (stats)
              stats->null_events.fetch_add(1, std::memory_order_relaxed);
            double pdf = T[hero] * n[hero];
            if (!(pdf > 0)) {
              terminated = true;
              return false;
            }
            RGB w = T * RGB(n) / pdf;
            beta *= w;
            ru *= w;
            rl *= T * RGB(m) / pdf;
            rescale(beta, ru, rl);
            if (!(beta.Max() > 0) || !(ru.Max() > 0)) {
              terminated = true;
              return false;
            }
            return true;
          },
          cancel, stats);
      if (terminated)
        break;
      if (scattered) {
        if (depth > roulette_depth) {
          double survival = std::min(1.0, (beta * RGB(eta_scale) / ru.Average()).Max());
          if (!(survival > 0) || sampler->Get1D() >= survival)
            break;
          beta /= survival;
        }
        continue;
      }
      if (T[hero] > 0) {
        RGB w = T / T[hero];
        beta *= w;
        ru *= w;
        rl *= w;
      }
    }
    beta *= glass_transmittance(state, distance - previous_distance);
    if (h.medium_boundary && !h.medium_boundary->keep_surface) {
      state.Cross(h, ray.d);
      ray = spawn(h, ray.d, ray.time(), state);
      continue;
    }
    if (h.alpha_miss) {
      ray = spawn(h, ray.d, ray.time(), state);
      continue;
    }
    bool invisible = false;
    point3f le = h.mat_ptr ? h.mat_ptr->emitted(ray, h, h.u, h.v, h.p, invisible) : point3f(0);
    if (invisible && !any_diffuse) {
      ray = spawn(h, ray.d, ray.time(), state);
      continue;
    }
    if (depth == 0 && h.infinite_area_hit && !(scene && scene->has_media))
      transparency = 1;
    bool hide_background =
        scene && scene->transparent_background && depth == 0 && h.infinite_area_hit;
    if (!hide_background) {
      double denom = ru.Average();
      if (!specular && use_light_sampling())
        denom =
            (ru + rl * RGB(light_pdf(lights, previous_point, ray.d, rng, ray.time()))).Average();
      if (denom > 0)
        L += beta * RGB(le) / denom;
    }
    if (depth >= max_depth || !(beta.Max() > 0) || !(ru.Average() > 0) || !h.mat_ptr)
      break;
    scatter_record s;
    // Dielectric membership uses the oriented geometric normal, not a shading normal.
    hit_record shading_hit = h;
    if (h.mat_ptr->is_dielectric())
      shading_hit.normal = geometric_normal(h);
    if (!h.mat_ptr->scatter(ray, shading_hit, s, sampler))
      break;
    if (!wrote_feature) {
      normal = h.normal;
      albedo = h.mat_ptr->get_albedo(h);
      wrote_feature = true;
    }
    if (!s.is_passthrough)
      ++depth;
    if (s.is_specular) {
      beta *= RGB(s.attenuation);
      if (s.is_transmission && !s.is_passthrough) {
        double eta2 = double(s.eta) * s.eta;
        beta /= eta2;
        eta_scale *= eta2;
      }
      cross_if_transmitted(state, h, ray.d, s.specular_ray.d);
      ray = spawn(h, unit_vector(s.specular_ray.d), ray.time(), state);
      if (!s.is_passthrough)
        specular = true;
    } else {
      L += direct_light(ray, h.p, &h, s.pdf_ptr, nullptr, state, beta, ru, world, lights, hero, rng,
                        sampler, cancel);
      vec3f wi = s.pdf_ptr->generate(sampler, any_diffuse, ray.time());
      if (!(wi.squared_length() > 0))
        break;
      wi = unit_vector(wi);
      Float p = s.pdf_ptr->value(wi, rng, ray.time());
      if (!(p > 0) || !std::isfinite(p))
        break;
      beta *= RGB(h.mat_ptr->f(ray, h, wi)) / p;
      rl = ru / p;
      previous_point = h.p;
      specular = false;
      any_diffuse = true;
      cross_if_transmitted(state, h, ray.d, wi);
      ray = spawn(h, wi, ray.time(), state);
    }
    reconcile_surface_origin(scene, h, ray, state, cancel);
    rescale(beta, ru, rl);
    if (!s.is_passthrough && depth > roulette_depth) {
      double survival = std::min(1.0, (beta * RGB(eta_scale) / ru.Average()).Max());
      if (!(survival > 0) || sampler->Get1D() >= survival)
        break;
      beta /= survival;
    }
  }
  radiance = L.FloatRGB();
}
