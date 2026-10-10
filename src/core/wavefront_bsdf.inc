R"MSL(
  // BSDF operations use directions pointing away from the hit. Evaluation returns
  // f * abs(cos), sampling returns its density in solid angle (or discrete mass).
  // Keep the existing CPU models, including their legacy glossy/reflection factors.
  constant uint bsdf_metal = 1, bsdf_microfacet = 2, bsdf_transmission = 3, bsdf_glossy = 4,
                bsdf_translucent = 5, bsdf_hair = 6, bsdf_openpbr = 7;
  struct BsdfValue {
    float3 value;
    float pdf;
  };
  struct BsdfSample {
    float3 direction, weight;
    float pdf, eta_squared;
    bool delta;
  };
  struct Frame {
    float3 x, y, z;
    float3 local(float3 v) const { return float3(dot(x, v), dot(y, v), dot(z, v)); }
    float3 world(float3 v) const { return x * v.x + y * v.y + z * v.z; }
  };
  Frame normal_frame(float3 n) {
    float3 y = normalize(cross(n, abs(n.x) > .9999999f ? float3(0, 1, 0) : float3(1, 0, 0)));
    return {cross(y, n), y, n};
  }
  float3 surface_tangent(Triangle t) {
    float2 a = t.uv01.zw - t.uv01.xy, b = t.uv2.xy - t.uv01.xy;
    float det = a.x * b.y - a.y * b.x;
    return abs(det) > 1e-12f
               ? (b.y * (t.p[1].xyz - t.p[0].xyz) - a.y * (t.p[2].xyz - t.p[0].xyz)) / det
               : normal_frame(t.geometric.xyz).x;
  }
  float3 conductor_fresnel(float cosine, float3 eta, float3 k) {
    float3 q = eta * eta + k * k, twice = 2 * eta * cosine;
    float c2 = cosine * cosine;
    return .5f * ((q * c2 - twice + 1) / (q * c2 + twice + 1) + (q - twice + c2) / (q + twice + c2));
  }
  float pow5(float x) { return x * x * x * x * x; }
  float micro_D(float3 h, float2 alpha, uint distribution) {
    if (h.z == 0)
      return 0;
    float z2 = h.z * h.z, e = dot(h.xy / alpha, h.xy / alpha) / z2;
    float denominator = pi * alpha.x * alpha.y * z2 * z2;
    return distribution ? exp(-e) / denominator : 1 / (denominator * (1 + e) * (1 + e));
  }
  float micro_lambda(float3 v, float2 alpha, uint distribution) {
    if (v.z == 0)
      return 0;
    float t2 = dot(v.xy * alpha, v.xy * alpha) / (v.z * v.z);
    if (!distribution)
      return .5f * (sqrt(1 + t2) - 1);
    float a = rsqrt(t2);
    return a >= 1.6f ? 0 : (1 - 1.259f * a + .396f * a * a) / (3.535f * a + 2.181f * a * a);
  }
  float micro_G(float3 view, float3 wi, float3 h) {
    float denominator = abs(dot(wi, h)); // CPU Pdf passes the sampled direction first.
    return min(1.f, min(2 * abs(h.z * wi.z) / denominator, 2 * abs(h.z * view.z) / denominator));
  }
  float erf_inverse(float a) {
    float t = log(max(fma(a, -a, 1.f), FLT_MIN)), p;
    if (abs(t) > 6.125f) {
      p = 3.03697567e-10f;
      p = fma(p, t, 2.93243101e-8f);
      p = fma(p, t, 1.22150334e-6f);
      p = fma(p, t, 2.84108955e-5f);
      p = fma(p, t, 3.93552968e-4f);
      p = fma(p, t, 3.02698812e-3f);
      p = fma(p, t, 4.83185798e-3f);
      p = fma(p, t, -2.64646143e-1f);
      p = fma(p, t, 8.40016484e-1f);
    } else {
      p = 5.43877832e-9f;
      p = fma(p, t, 1.43286059e-7f);
      p = fma(p, t, 1.22775396e-6f);
      p = fma(p, t, 1.12962631e-7f);
      p = fma(p, t, -5.61531961e-5f);
      p = fma(p, t, -1.47697705e-4f);
      p = fma(p, t, 2.31468701e-3f);
      p = fma(p, t, 1.15392562e-2f);
      p = fma(p, t, -2.32015476e-1f);
      p = fma(p, t, 8.86226892e-1f);
    }
    return a * p;
  }
  float2 beckmann_slopes(float cosine, float2 u) {
    if (cosine > .9999f) {
      float r = sqrt(-log(1 - u.x)), phi = 2 * pi * u.y;
      return r * float2(cos(phi), sin(phi));
    }
    float tangent = sqrt(max(0.f, 1 - cosine * cosine)) / cosine, cotangent = 1 / tangent;
    // Same Abramowitz--Stegun approximation as the CPU's Erf helper; MSL has no erf.
    float q = 1 / (1 + .3275911f * cotangent);
    float c =
        1 - (((((1.061405429f * q - 1.453152027f) * q) + 1.421413741f) * q - .284496736f) * q +
             .254829592f) *
                q * exp(-cotangent * cotangent);
    float a = -1, theta = acos(cosine);
    float fit = 1 + theta * (-.876f + theta * (.4265f - .0594f * theta));
    float x = max(u.x, 1e-6f), b = c - (1 + c) * pow(1 - x, fit);
    float norm = 1 / (1 + c + rsqrt(pi) * tangent * exp(-cotangent * cotangent));
    for (uint j = 0; j < 9; ++j) {
      if (!(b >= a && b <= c))
        b = .5f * (a + c);
      float inv = erf_inverse(b), value = norm * (1 + b + rsqrt(pi) * tangent * exp(-inv * inv)) - x;
      if (abs(value) < 1e-5f)
        break;
      if (value > 0)
        c = b;
      else
        a = b;
      b -= value / (norm * (1 - inv * tangent));
    }
    return float2(erf_inverse(b), erf_inverse(2 * max(u.y, 1e-6f) - 1));
  }
  float3 micro_normal(float3 view, float2 alpha, uint distribution, float2 u) {
    bool flip = view.z < 0;
    float3 v = normalize(float3(alpha * view.xy, view.z)) * (flip ? -1.f : 1.f), h;
    if (distribution) {
      float2 slope = beckmann_slopes(v.z, u);
      float r = length(v.xy), c = r > 0 ? v.x / r : 1, s = r > 0 ? v.y / r : 0;
      slope = float2(c * slope.x - s * slope.y, s * slope.x + c * slope.y) * alpha;
      h = normalize(float3(-slope, 1));
    } else {
      // Heitz's disk construction samples the same GGX visible-normal density
      // as the CPU slope inversion; changing the proposal mapping does not change f/pdf.
      float3 x = v.z < .99999f ? normalize(cross(float3(0, 0, 1), v)) : float3(1, 0, 0);
      float3 y = cross(v, x);
      float r = sqrt(u.x), phi = 2 * pi * u.y, px = r * cos(phi), py = r * sin(phi),
            s = .5f * (1 + v.z);
      py = mix(sqrt(max(0.f, 1 - px * px)), py, s);
      float3 nh = px * x + py * y + sqrt(max(0.f, 1 - px * px - py * py)) * v;
      h = normalize(float3(alpha * nh.xy, max(0.f, nh.z)));
    }
    return flip ? -h : h;
  }
  BsdfValue micro_eval(Material m, float3 view, float3 wi, float3 rho, float distance) {
    bool reflection = view.z * wi.z > 0;
    if (view.z == 0 || wi.z == 0)
      return {0, 0};
    float2 alpha = m.parameters.xy;
    if (m.type == bsdf_transmission) {
      float ratio = view.z > 0 ? 1 / m.eta.x : m.eta.x;
      float3 sum = reflection ? view + wi : ratio * view + wi;
      if (dot(sum, sum) == 0)
        return {0, 0};
      float3 h = normalize(sum);
      if (h.z < 0)
        h = -h;
      float vo = dot(view, h), vi = dot(wi, h);
      if (vo * view.z <= 0 || vi * wi.z <= 0)
        return {0, 0};
      float F = fresnel(abs(vo), 1 / ratio), D = micro_D(h, alpha, m.distribution);
      float G = 1 / (1 + micro_lambda(view, alpha, m.distribution) +
                     micro_lambda(wi, alpha, m.distribution));
      float density = D * abs(vo) / (abs(view.z) * (1 + micro_lambda(view, alpha, m.distribution)));
      if (reflection)
        return {float3(F * D * G / (4 * abs(view.z))), F * density / (4 * abs(vo))};
      float denominator = (vi + ratio * vo) * (vi + ratio * vo);
      float3 tint = rho * (view.z < 0 ? exp(-distance * m.k.xyz) : float3(1));
      return {tint * ((1 - F) * D * G * abs(vo * vi) * ratio * ratio / (abs(view.z) * denominator)),
              (1 - F) * density * abs(vi) / denominator};
    }
    if (!reflection)
      return {0, 0};
    float3 h = normalize(view + wi);
    float D = micro_D(h, alpha, m.distribution), G = micro_G(view, wi, h);
    float pdf = D * G * abs(dot(wi, h)) / (abs(wi.z) * 4 * dot(wi, h));
    if (m.type == bsdf_microfacet)
      return {rho * conductor_fresnel(abs(wi.z), m.eta.xyz, m.k.xyz) * G * D * abs(view.z) /
                  (4 * wi.z * view.z),
              pdf};
    if (dot(h, view) < 0)
      return {0, 0};
    float3 diffuse = (28.f / (23 * pi)) * m.k.xyz * rho * (1 - m.eta.xyz) *
                     (1 - pow5(1 - .5f * abs(view.z))) * (1 - pow5(1 - .5f * abs(wi.z)));
    float3 specular = D / (4 * abs(dot(view, h)) * max(abs(view.z), abs(wi.z))) *
                      (m.eta.xyz + pow5(1 - dot(wi, h)) * (1 - m.eta.xyz));
    return {(diffuse + specular) * dot(h, view), .5f * (abs(view.z) / pi + pdf)};
  }
  // Four-lobe fiber model, matching math/mathinline.h and hair_pdf. Hair uses a
  // longitudinal frame and has no surface-hemisphere rejection or extra cosine.
  struct HairState {
    float3 attenuation[4];
    float4 probabilities, variance;
    float sine, cosine, phi, gamma_o, gamma_t, logistic, alpha;
  };
  float hair_i0(float x) {
    float value = 0, term = 1;
    for (uint j = 0; j < 10; ++j) {
      value += term;
      term *= x * x / (4 * (j + 1) * (j + 1));
    }
    return value;
  }
  float hair_longitudinal(float si, float ci, float so, float co, float variance) {
    float a = ci * co / variance, b = si * so / variance;
    float log_i = a > 12 ? a + .5f * (-log(2 * pi) + log(1 / a) + 1 / (8 * a)) : log(hair_i0(a));
    return variance <= .1f ? exp(log_i - b - 1 / variance + .6931f + log(1 / (2 * variance)))
                           : exp(-b) * hair_i0(a) / (sinh(1 / variance) * 2 * variance);
  }
  float logistic_cdf(float x, float s) { return 1 / (1 + exp(-x / s)); }
  float hair_azimuth(float phi, uint lobe, HairState h) {
    if (lobe == 3)
      return 1 / (2 * pi);
    float d = phi - (2 * lobe * h.gamma_t - 2 * h.gamma_o + lobe * pi);
    d -= 2 * pi * floor((d + pi) / (2 * pi));
    float e = exp(-abs(d) / h.logistic);
    return e / (h.logistic * (1 + e) * (1 + e) *
                (logistic_cdf(pi, h.logistic) - logistic_cdf(-pi, h.logistic)));
  }
  float2 hair_tilt(HairState h, uint lobe) {
    float angle = h.alpha * (lobe == 0 ? -2.f : (lobe == 1 ? 1.f : (lobe == 2 ? 4.f : 0.f)));
    return float2(h.sine * cos(angle) + h.cosine * sin(angle),
                  abs(h.cosine * cos(angle) - h.sine * sin(angle)));
  }
  HairState prepare_hair(Material m, float3 view, float impact) {
    HairState h;
    float eta = m.parameters.x, bm = m.parameters.y, bn = m.parameters.z;
    float v = .726f * bm + .812f * bm * bm + 3.7f * pow(bm, 20.f);
    h.variance = max(float4(v * v, .25f * v * v, 4 * v * v, 4 * v * v), float4(1e-12f));
    h.logistic = max(sqrt(pi / 8) * (.265f * bn + 1.194f * bn * bn + 5.372f * pow(bn, 22.f)), 1e-8f);
    h.alpha = m.parameters.w * pi / 180;
    h.sine = view.x;
    h.cosine = sqrt(max(0.f, 1 - view.x * view.x));
    h.phi = atan2(view.z, view.y);
    h.gamma_o = asin(clamp(impact, -1.f, 1.f));
    float ct = sqrt(max(FLT_MIN, 1 - view.x * view.x / (eta * eta)));
    float gt = impact * h.cosine / sqrt(max(FLT_MIN, eta * eta - view.x * view.x));
    h.gamma_t = asin(clamp(gt, -1.f, 1.f));
    float3 T = exp(-m.eta.xyz * (2 * sqrt(max(0.f, 1 - gt * gt)) / ct));
    float f = fresnel(h.cosine * sqrt(max(0.f, 1 - impact * impact)), eta);
    h.attenuation[0] = f;
    h.attenuation[1] = (1 - f) * (1 - f) * T;
    h.attenuation[2] = h.attenuation[1] * T * f;
    float3 ratio = T * f;
    h.attenuation[3] =
        select(float3(0), h.attenuation[2] * ratio / max(1 - ratio, float3(FLT_MIN)), ratio < 1);
    for (uint j = 0; j < 4; ++j)
      h.probabilities[j] = dot(h.attenuation[j], float3(1.f / 3));
    float total = dot(h.probabilities, float4(1));
    h.probabilities = total > 0 ? h.probabilities / total : float4(1, 0, 0, 0);
    return h;
  }
  BsdfValue hair_eval(HairState h, float3 wi) {
    float ci = sqrt(max(0.f, 1 - wi.x * wi.x)), phi = atan2(wi.z, wi.y) - h.phi;
    BsdfValue out = {0, 0};
    for (uint j = 0; j < 4; ++j) {
      float2 tilt = hair_tilt(h, j);
      float density =
          hair_longitudinal(wi.x, ci, tilt.x, tilt.y, h.variance[j]) * hair_azimuth(phi, j, h);
      out.value += h.attenuation[j] * density;
      out.pdf += h.probabilities[j] * density;
    }
    return out;
  }
  float3 hair_sample(HairState h, float4 u) {
    uint lobe = 0;
    while (lobe < 3 && u.x >= h.probabilities[lobe]) {
      u.x -= h.probabilities[lobe];
      ++lobe;
    }
    float2 tilt = hair_tilt(h, lobe);
    float ct = 1 + h.variance[lobe] *
                       log(max(u.y, 1e-5f) + (1 - max(u.y, 1e-5f)) * exp(-2 / h.variance[lobe]));
    float si =
        clamp(-ct * tilt.x + sqrt(max(0.f, 1 - ct * ct)) * cos(2 * pi * u.z) * tilt.y, -1.f, 1.f);
    float ci = sqrt(max(0.f, 1 - si * si)), delta;
    if (lobe == 3)
      delta = 2 * pi * u.w;
    else {
      float low = logistic_cdf(-pi, h.logistic), mass = logistic_cdf(pi, h.logistic) - low;
      float x = -h.logistic * log(1 / (u.w * mass + low) - 1);
      delta = 2 * lobe * h.gamma_t - 2 * h.gamma_o + lobe * pi + clamp(x, -pi, pi);
    }
    return float3(si, ci * cos(h.phi + delta), ci * sin(h.phi + delta));
  }
  struct SurfaceBsdf {
    Material material;
    Frame frame;
    float3 view, rho, transmission;
    float distance, impact;
  };
  BsdfValue surface_eval(SurfaceBsdf b, float3 direction) {
    float3 wi = b.frame.local(direction);
    Material m = b.material;
    if (m.type == bsdf_metal)
      return {0, 0};
    if (m.type == bsdf_hair)
      return hair_eval(prepare_hair(m, b.view, b.impact), wi);
    if (m.type == bsdf_translucent) {
      float r = dot(b.rho, float3(1)), t = dot(b.transmission, float3(1));
      float probability = (r + t) > 0 ? r / (r + t) : .5f;
      bool reflection = b.view.z * wi.z > 0;
      return {(reflection ? b.rho : b.transmission) * (abs(wi.z) / pi),
              (reflection ? probability : 1 - probability) * abs(wi.z) / pi};
    }
    return micro_eval(m, b.view, wi, b.rho, b.distance);
  }
  BsdfSample surface_sample(SurfaceBsdf b, float4 u) {
    Material m = b.material;
    BsdfSample out = {0, 0, 0, 1, false};
    float3 wi;
    if (m.type == bsdf_metal) {
      float z = 1 - 2 * u.y, r = sqrt(max(0.f, 1 - z * z)), phi = 2 * pi * u.z;
      wi = normalize(reflect(-b.view, float3(0, 0, 1)) +
                     m.parameters.z * float3(r * cos(phi), r * sin(phi), z));
      out.direction = b.frame.world(wi);
      out.weight = b.rho * conductor_fresnel(abs(b.view.z), m.eta.xyz, m.k.xyz);
      out.pdf = 1;
      out.delta = true;
      return out;
    }
    if (m.type == bsdf_hair)
      wi = hair_sample(prepare_hair(m, b.view, b.impact), u);
    else if (m.type == bsdf_translucent) {
      float r = dot(b.rho, float3(1)), t = dot(b.transmission, float3(1));
      float probability = (r + t) > 0 ? r / (r + t) : .5f;
      wi = float3(sqrt(u.y) * cos(2 * pi * u.z), sqrt(u.y) * sin(2 * pi * u.z),
                  sqrt(max(0.f, 1 - u.y)));
      wi.z *= b.view.z < 0 ? -1.f : 1.f;
      if (u.x >= probability)
        wi.z = -wi.z;
    } else if (m.type == bsdf_transmission && m.eta.x == 1) {
      out.direction = b.frame.world(-b.view);
      out.pdf = 1;
      out.delta = true;
      out.weight = b.rho * (b.view.z < 0 ? exp(-b.distance * m.k.xyz) : float3(1));
      return out;
    } else if (m.type == bsdf_glossy && u.x >= .5f) {
      wi = float3(sqrt(u.y) * cos(2 * pi * u.z), sqrt(u.y) * sin(2 * pi * u.z),
                  sqrt(max(0.f, 1 - u.y)));
      if (b.view.z < 0)
        wi.z = -wi.z;
    } else {
      float3 h = micro_normal(b.view, m.parameters.xy, m.distribution, u.yz);
      if (m.type == bsdf_transmission) {
        float eta = b.view.z > 0 ? m.eta.x : 1 / m.eta.x;
        bool transmission = u.x >= fresnel(dot(b.view, h), eta);
        wi = transmission ? refract(-b.view, h, 1 / eta) : reflect(-b.view, h);
        if (dot(wi, wi) == 0 || ((wi.z * b.view.z < 0) != transmission))
          return out;
      } else
        wi = reflect(-b.view, h);
    }
    out.direction = b.frame.world(wi);
    BsdfValue value = surface_eval(b, out.direction);
    out.pdf = value.pdf;
    if (out.pdf > 0)
      out.weight = value.value / out.pdf;
    return out;
  }
)MSL"
