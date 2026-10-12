R"MSL(

  // Specialize away all interior state/traversal for scenes without boundaries.
  // This is a pipeline constant, not a per-ray material branch.
  constant bool scene_has_boundaries [[function_constant(0)]];
  constant uint texture_features [[function_constant(1)]];
  constant bool scene_has_quadrics [[function_constant(2)]];
  constant bool scene_has_instances [[function_constant(3)]];
  constant bool scene_has_atmosphere [[function_constant(4)]];
  constant bool defer_shadow_lighting [[function_constant(5)]];
  constant bool primary_coverage [[function_constant(6)]];
#include <metal_raytracing>
#include <metal_stdlib>
  using namespace metal;
  using namespace raytracing;
  constant float pi = 3.14159265358979323846f;
  struct Triangle {
    float4 p[3], n[3], uv01, uv2, geometric;
    uint material, light, boundary, flags;
  };
  struct Texture {
    float4 a, b, c, mapping;
    uint type, offset, width, height;
  };
  struct Material {
    float4 diffuse, optical, attenuation;
    uint texture, emissive, dielectric, priority;
    float4 eta, k, parameters, extra;
    uint type, distribution, second_texture, data;
  };
  struct Boundary {
    float4 sigma_a, sigma_s;
    uint material, subsurface, diffusion_offset, medium;
    float4 color, radius, lo, hi;
  };
  // Bounded membership storage, not a walk-event limit. Overflow terminates only
  // the affected ray and is reported to the host; it never silently loses a medium.
  struct VolumeState {
    uint ids[16], count, padding[3];
  };
  struct Light {
    float4 forward, right, up, color;
    uint type, texture, primitive, flags;
    uint instance, alias_offset, alias_width, alias_height;
  };
  struct Camera {
    float4 origin, lower_left, horizontal, vertical, right, up, forward;
    uint orthographic, flip_x, padding[2];
    float4 shutter, motion_origin[2], motion_rotation[2], lens;
    uint moving, elements, pupils, simple_weighting;
  };
  struct Params {
    Camera camera;
    uint width, height, active, sample, lights, depth, roulette, sampler, seed, max_depth,
        boundaries, environments;
    float4 boundary_lo, boundary_hi;
  };
  struct Pixel {
    float4 radiance, normal, albedo;
  };
  struct Differential {
    float4 ox, oy, dx, dy;
  };
  struct Path {
    // origin.w stores the last surface primitive plus one as uint bits. A ray
    // leaving a planar triangle cannot physically hit that same triangle again
    // before another interaction. Suppress Metal's occasional positive-distance
    // self hits without changing the free-flight distance or scene-scale bias.
    float4 origin, direction, beta, previous;
    uint dimension, depth, internal, flags;
  };
  struct Hit {
    uint pixel, primitive;
    float u, v;
    float distance, w;
    uint padding[2];
  };
  struct Shadow {
    float4 origin, direction, contribution;
    uint pixel, padding[3];
    VolumeState state;
  };
  // Counters: current ray queue, hit queue, shadow queue, next ray queue, error.
  struct Counts {
    atomic_uint rays, hits, shadows, next, error, discarded, scattering, csg_failed;
  };
  static_assert(sizeof(Triangle) == 160 && sizeof(Texture) == 80 && sizeof(Material) == 144 &&
                    sizeof(Boundary) == 112 && sizeof(VolumeState) == 80 && sizeof(Light) == 96 &&
                    sizeof(Camera) == 240 && sizeof(Params) == 320 && sizeof(Pixel) == 48 &&
                    sizeof(Path) == 80 && sizeof(Hit) == 32 && sizeof(Shadow) == 144 &&
                    sizeof(Counts) == 32,
                "wavefront host/device ABI mismatch");

  uint hash_seed(uint x, uint seed) {
    x = 0x6217c6e1u ^ (x + seed * 0x9e3779b9u);
    x ^= x >> 17;
    x *= 0xed5ad4bbu;
    x ^= x >> 11;
    x *= 0xac4c1b51u;
    x ^= x >> 15;
    x *= 0x31848babu;
    return x ^ (x >> 14);
  }
  uint owen(uint x, uint seed) {
    x = reverse_bits(x);
    seed = hash_seed(seed, 0xa14a177du);
    x ^= x * 0x3d20adeau;
    x += seed;
    x *= (seed >> 16) | 1u;
    x ^= x * 0x05526c56u;
    x ^= x * 0x53a22864u;
    return reverse_bits(x);
  }
  uint sample_bits(uint pixel, thread uint &dimension, constant Params &p,
                   device const uint *directions, device const uint *blue,
                   bool continuous = false) {
    uint dim = dimension++, x = pixel / p.height, y = pixel % p.height;
    if (p.sampler == 0)
      return hash_seed(hash_seed(p.sample, pixel ^ p.seed), dim);
    if (p.sampler == 3 && dim < 8) {
      uint tile = ((x & 127) + (y & 127) * 128) * 8 + dim;
      uint ranked = (p.sample & 255) ^ blue[65536 + tile];
      uint stratum = blue[ranked * 256 + dim] ^ blue[65536 + 131072 + tile];
      // Preserve the blue-noise stratum, restoring its missing low bits only
      // for image-environment proposals. Ordinary callers retain the midpoint.
      uint low =
          continuous ? hash_seed(hash_seed(p.sample, pixel ^ p.seed), dim) & 0xffffffu : 0x800000u;
      return (stratum << 24) | low;
    }
    uint seed = p.sampler == 3 ? x * 0x9e3779b9u ^ y * 0x85ebca6bu : hash_seed(pixel, p.seed);
    seed += dim / 21200;
    dim %= 21200;
    uint index = owen(p.sample, seed), value = 0;
    for (uint bit = 0; index; index >>= 1, ++bit)
      if (index & 1)
        value ^= directions[dim * 32 + bit];
    return owen(value, seed ^ (dim + (seed << 6) + (seed >> 2)));
  }
  float sample1(uint pixel, thread uint &dimension, constant Params &p,
                device const uint *directions, device const uint *blue, bool continuous = false) {
    return min(float(sample_bits(pixel, dimension, p, directions, blue, continuous)) * 0x1p-32f,
               0x1.fffffep-1f);
  }

)MSL"
#include "wavefront_textures.inc"
R"MSL(

  float mis(float a, float b) {
    float ratio = b / a;
    // Native NEE uses the balance heuristic. Keep the same estimator policy
    // across backends, including the legacy glossy material's declared PDF.
    return 1 / (1 + ratio);
  }
  float3 basis_direction(float3 n, float2 uv) {
    float3 tangent = normalize(cross(abs(n.x) > .9f ? float3(0, 1, 0) : float3(1, 0, 0), n));
    float3 bitangent = cross(n, tangent);
    float r = sqrt(uv.x), phi = 2 * pi * uv.y;
    return r * cos(phi) * tangent + r * sin(phi) * bitangent + sqrt(max(0.f, 1 - uv.x)) * n;
  }
  // Same exact EON albedo and multiple-scattering compensation as
  // materials/normalmap.cpp; the coefficient triple is exported from its child.
  float loss(float mu, float b) {
    mu = clamp(mu, 0.f, 1.f);
    float s = sqrt(max(0.f, 1 - mu * mu));
    float g = s * (acos(mu) - s * mu) + (2.f / 3) * (s * mu * (1 + s + s * s) / (1 + s) - s);
    return b * max(0.f, .5f - 2 / (3 * pi) - g / pi);
  }
  float3 brdf(float3 wo, float3 wi, float3 n, float3 rho, Material m) {
    float co = dot(wo, n), ci = dot(wi, n);
    if (co <= 0 || ci <= 0)
      return 0;
    if (m.diffuse.y == 0)
      return rho / pi;
    float t = dot(wo - co * n, wi - ci * n), angular = t > 0 ? t / max(ci, co) : t;
    float3 single = rho * (m.diffuse.x + m.diffuse.y * angular);
    float avg = m.diffuse.z;
    float3 multiple = rho * rho * (1 - avg) / (1 - rho * avg) * (loss(co, m.diffuse.y) / avg) *
                      loss(ci, m.diffuse.y);
    return (single + multiple) / pi;
  }
  float3 offset_origin(float3 position, float3 error, float3 normal) {
    float3 p = position + dot(abs(normal), error) * normal;
    // MSL has no nextafter. Move one IEEE float toward the chosen normal side,
    // including signed zero; the analytical error offset above handles magnitude.
    for (uint k = 0; k < 3; ++k)
      if (normal[k] != 0) {
        uint bits = as_type<uint>(p[k]);
        // Metal flushes subnormal arithmetic; use the smallest normal at zero.
        if (p[k] == 0)
          bits = normal[k] > 0 ? 0x00800000u : 0x80800000u;
        else
          bits += ((p[k] > 0) == (normal[k] > 0)) ? 1u : uint(-1);
        p[k] = as_type<float>(bits);
      }
    return p;
  }
  float2 triangle_uv(Triangle t, float3 b) {
    return b.x * t.uv01.xy + b.y * t.uv01.zw + b.z * t.uv2.xy;
  }
  void surface_point(Triangle t, float2 uv, thread float3 &position, thread float3 &error) {
    float3 b = float3(1 - uv.x - uv.y, uv.x, uv.y);
    position = b.x * t.p[0].xyz + b.y * t.p[1].xyz + b.z * t.p[2].xyz;
    error = 4.172327e-7f * (abs(b.x * t.p[0].xyz) + abs(b.y * t.p[1].xyz) + abs(b.z * t.p[2].xyz));
    float3 extent = abs(t.p[1].xyz - t.p[0].xyz) + abs(t.p[2].xyz - t.p[0].xyz);
    error += 3 * FLT_EPSILON * max(extent.x, max(extent.y, extent.z));
  }

)MSL"
#include "wavefront_intersections.inc"
#include "wavefront_light_sampling.inc"
R"MSL(

  float3 environment(Light l, float3 wi, float3 position, device const Texture *textures,
                     device const float4 *texels, bool evaluated = false) {
    if (scene_has_atmosphere && !evaluated)
      return 1; // Native radiance is multiplied into completed shadow requests.
    if (l.type == 0) {
      float3 local = float3(dot(wi, l.forward.xyz), dot(wi, l.right.xyz), dot(wi, l.up.xyz));
      float2 uv = float2(.5f - atan2(-local.x, -local.z) / (2 * pi),
                         .5f + asin(clamp(local.y, -1.f, 1.f)) / pi);
      return texture_value(l.texture, uv, position, textures, texels);
    }
    if (l.type == 1) {
      float z = dot(wi, l.forward.xyz), delta = l.forward.w;
      if (z <= 0 || z < 1 - delta || (l.flags && dot(wi, l.color.xyz) <= 0))
        return 0;
      float tan_radius = sqrt(delta * (2 - delta)) / (1 - delta);
      float2 uv = .5f + .5f * float2(dot(wi, l.right.xyz), dot(wi, l.up.xyz)) / (z * tan_radius);
      return texture_value(l.texture, uv, position, textures, texels);
    }
    return 0;
  }
  float2 environment_uv(Light l, float3 wi) {
    float3 local = float3(dot(wi, l.forward.xyz), dot(wi, l.right.xyz), dot(wi, l.up.xyz));
    float u = .5f - atan2(-local.x, -local.z) / (2 * pi);
    return float2(fract(u), .5f + asin(clamp(local.y, -1.f, 1.f)) / pi);
  }
  // A tiny solar cone loses several bits in 1-dot(w,axis). Compare its
  // squared chord instead, and use this same predicate for solar radiance.
  // Otherwise a bright rim sample can receive only the uniform-sky density.
  bool solar_disk_contains(float3 wi, float3 axis, float one_minus_cosine) {
    float3 chord = normalize(wi) - normalize(axis);
    return dot(chord, chord) <= 2 * one_minus_cosine;
  }
  float light_pdf(Light l, float3 wi, device const float4 *texels) {
    if (l.flags & 128)
      return (1-l.forward.w)/(4*pi) +
             (solar_disk_contains(wi,l.forward.xyz,l.right.w) ? l.forward.w/(2*pi*l.right.w) : 0);
    if (l.type != 0)
      return 1 / (2 * pi * l.forward.w);
    if (l.alias_offset == UINT_MAX)
      return 1 / (4 * pi);
    float2 uv = environment_uv(l, wi);
    uint x = min(uint(uv.x * l.alias_width), l.alias_width - 1);
    uint y = min(uint(uv.y * l.alias_height), l.alias_height - 1);
    return texels[l.alias_offset + y * l.alias_width + x].z;
  }
  float3 sample_environment(Light l, uint bits, float latitude, device const float4 *texels) {
    if (l.flags & 128) {
      float u = min(float(bits) * 0x1p-32f, 0x1.fffffep-1f);
      float phi = 2*pi*latitude;
      if (u < l.forward.w) {
        float d = (u / l.forward.w) * l.right.w;
        float3 radial = cos(phi)*l.right.xyz + sin(phi)*l.up.xyz;
        float3 direction = normalize((1-d)*l.forward.xyz + sqrt(d*(2-d)) * radial);
        // As on CPU, move only samples rounded outside the cone inward.
        if (!solar_disk_contains(direction,l.forward.xyz,l.right.w)) {
          float angle = max(0.f, 2 * asin(sqrt(.5f*l.right.w)) - 4 * FLT_EPSILON);
          d = 2 * sin(.5f*angle) * sin(.5f*angle);
          direction = normalize((1-d)*l.forward.xyz + sqrt(d*(2-d)) * radial);
        }
        return direction;
      }
      u = (u - l.forward.w) / (1-l.forward.w);
      float y = 1-2*u, r = sqrt(max(0.f, 1-y*y));
      return float3(r*cos(phi), y, r*sin(phi));
    }
    if (l.alias_offset == UINT_MAX) {
      float u = min(float(bits) * 0x1p-32f, 0x1.fffffep-1f);
      float y = 1 - 2 * u, r = sqrt(max(0.f, 1 - y * y)), phi = 2 * pi * latitude;
      return float3(r * cos(phi), y, r * sin(phi));
    }
    uint cells = l.alias_width * l.alias_height;
    uint column_bits = cells > 1 ? 32 - clz(cells - 1) : 0;
    uint coin_bits = min(24u, 32 - column_bits);
    uint cell = column_bits ? bits >> (32 - column_bits) : 0;
    uint coin = (bits >> (32 - column_bits - coin_bits)) & ((1u << coin_bits) - 1);
    float resolution = float(1u << coin_bits);
    float residual = (float(coin) + .5f) / resolution;
    float4 entry = texels[l.alias_offset + cell];
    // Reuse the alias decision's remapped residual for longitude. No rejection
    // or search, and only the two random dimensions needed for a direction.
    if (coin < uint(entry.x * resolution))
      residual /= entry.x;
    else {
      residual = (residual - entry.x) / (1 - entry.x);
      cell = as_type<uint>(entry.y);
    }
    residual = clamp(residual, 0.f, 0x1.fffffep-1f);
    uint x = cell % l.alias_width, y = cell / l.alias_width;
    float phi = -2 * pi * (x + residual) / l.alias_width;
    float c0 = -cos(pi * float(y) / l.alias_height);
    float c1 = -cos(pi * float(y + 1) / l.alias_height);
    float cosine = mix(c0, c1, latitude), radius = sqrt(max(0.f, 1 - cosine * cosine));
    float3 local = float3(radius * sin(phi), cosine, radius * cos(phi));
    return local.x * l.forward.xyz + local.y * l.right.xyz + local.z * l.up.xyz;
  }

  uint active_boundary(thread const VolumeState &state, device const Boundary *boundaries,
                       device const Material *materials) {
    if (!scene_has_boundaries)
      return 0;
    uint winner = 0;
    for (uint j = 0; j < state.count; ++j) {
      uint id = state.ids[j];
      if ((texture_features & 16) && !materials[boundaries[id].material].dielectric)
        continue;
      if ((texture_features & 16) && boundaries[id].subsurface == 3 && boundaries[id].radius.w == 0)
        continue;
      if (!winner || materials[boundaries[id].material].priority <=
                         materials[boundaries[winner].material].priority)
        winner = id;
    }
    return winner;
  }
  uint active_medium(thread const VolumeState &state, device const Boundary *boundaries,
                     device const Material *materials) {
    // Explicit media use nesting, while SSS follows the winning dielectric.
    if (texture_features & 16)
      for (uint j = state.count; j > 0; --j)
        if (boundaries[state.ids[j - 1]].subsurface == 3)
          return state.ids[j - 1];
    uint body = active_boundary(state, boundaries, materials);
    return body && boundaries[body].subsurface ? body : 0;
  }
  bool cross_boundary(thread VolumeState &state, uint id, bool entering) {
    uint found = state.count;
    for (uint j = 0; j < state.count; ++j)
      if (state.ids[j] == id)
        found = j;
    if (entering) {
      if (found < state.count)
        return true; // Membership is a set per placement.
      if (state.count == 16)
        return false;
      state.ids[state.count++] = id;
    } else if (found < state.count) {
      for (uint j = found + 1; j < state.count; ++j)
        state.ids[j - 1] = state.ids[j];
      --state.count;
    }
    return true;
  }
  float boundary_ior(uint id, device const Boundary *boundaries, device const Material *materials) {
    return id ? materials[boundaries[id].material].optical.x : 1;
  }
  float3 extinction(uint id, device const Boundary *boundaries, device const Material *materials) {
    if (!id)
      return 0;
    Boundary b = boundaries[id];
    return b.subsurface ? b.sigma_a.xyz + b.sigma_s.xyz : materials[b.material].attenuation.xyz;
  }
  float fresnel(float cosine, float eta) {
    float sin2 = max(0.f, 1 - cosine * cosine) / (eta * eta);
    if (sin2 >= 1)
      return 1;
    float ct = sqrt(1 - sin2);
    float rs = (cosine - eta * ct) / (cosine + eta * ct);
    float rp = (eta * cosine - ct) / (eta * cosine + ct);
    return .5f * (rs * rs + rp * rp);
  }
  float3 local_to_world(float3 value, float3 n) {
    float3 x = normalize(cross(abs(n.x) > .9f ? float3(0, 1, 0) : float3(1, 0, 0), n));
    return value.x * x + value.y * cross(n, x) + value.z * n;
  }
  float phase_pdf(float3 incoming, float3 wi, float g) {
    float d = (1 - g) * (1 - g) + 2 * g * (1 - clamp(dot(incoming, wi), -1.f, 1.f));
    return (1 - g * g) / (4 * pi * d * sqrt(d));
  }
  float3 phase_sample(float3 incoming, float g, float u, float v) {
    float cosine = 1 - 2 * u;
    if (abs(g) > 1e-6f) {
      float s = (1 - g * g) / (1 - g + 2 * g * u);
      cosine = clamp((1 + g * g - s * s) / (2 * g), -1.f, 1.f);
    }
    float sine = sqrt(max(0.f, 1 - cosine * cosine)), phi = 2 * pi * v;
    return normalize(local_to_world(float3(sine * cos(phi), sine * sin(phi), cosine), incoming));
  }
  float ggx_g1(float cosine, float alpha) {
    return cosine == 0
               ? 0
               : 2 / (1 + sqrt(1 + alpha * alpha * (1 - cosine * cosine) / (cosine * cosine)));
  }
  // Port of SubsurfaceBoundaryBSDF: returns f*abs(cos_i), plus sampling density.
  float2 boundary_eval(float3 wo, float3 wi, float3 n, float eta, float roughness) {
    float co = dot(wo, n), ci = dot(wi, n), alpha = roughness * roughness;
    if (co <= 0 || ci == 0 || alpha == 0 || eta == 1)
      return 0;
    bool reflection = ci > 0;
    float3 sum = wo + (reflection ? 1.f : eta) * wi;
    if (dot(sum, sum) == 0)
      return 0;
    float3 m = normalize(sum);
    if (dot(m, n) < 0)
      m = -m;
    float om = dot(wo, m), im = dot(wi, m), mn = dot(m, n);
    if (om <= 0 || (reflection ? im <= 0 : im >= 0))
      return 0;
    float3 tangent = m - mn * n;
    float q = dot(tangent, tangent) + alpha * alpha * mn * mn;
    float d = alpha * alpha / (pi * q * q), f = fresnel(om, eta);
    float g = ggx_g1(co, alpha) * ggx_g1(abs(ci), alpha);
    if (reflection)
      return float2(f * d * g / (4 * co), d * mn * f / (4 * om));
    float denom = im + om / eta;
    return denom == 0 ? float2(0)
                      : float2((1 - f) * d * g * abs(im * om) / (co * denom * denom * eta * eta),
                               (1 - f) * d * mn * abs(im) / (denom * denom));
  }
  bool boundary_sample(float3 wo, float3 n, float eta, float roughness, float branch, float u,
                       float v, thread float3 &wi, thread float &weight, thread float &pdf) {
    bool delta = roughness == 0 || eta == 1;
    float3 m = n;
    if (!delta) {
      float alpha = roughness * roughness, z = rsqrt(1 + alpha * alpha * u / (1 - u));
      float xy = sqrt(max(0.f, 1 - z * z));
      m = local_to_world(float3(xy * cos(2 * pi * v), xy * sin(2 * pi * v), z), n);
    }
    float om = dot(wo, m);
    if (om <= 0)
      return false;
    float f = fresnel(om, eta);
    bool transmission = branch >= f;
    wi = transmission ? refract(-wo, m, 1 / eta) : reflect(-wo, m);
    if (dot(wi, wi) == 0 || (transmission ? dot(wi, n) >= 0 : dot(wi, n) <= 0))
      return false;
    wi = normalize(wi);
    float2 eval = delta ? float2(transmission ? 1 / (eta * eta) : 1, transmission ? 1 - f : f)
                        : boundary_eval(wo, wi, n, eta, roughness);
    pdf = eval.y;
    weight = delta ? eval.x : (pdf > 0 ? eval.x / pdf : 0);
    return weight > 0 && pdf > 0;
  }

  // Classify the camera origin by replaying all boundary crossings from outside
  // their aggregate bounds. This also works after interactively moving inside a
  // solid. No per-pixel CPU ray query or per-body scan is needed.
  bool initial_membership(float3 origin, thread VolumeState &state, constant Params &p,
                          device const Triangle *triangles, device const Quadric *quadrics,
                          device const Surface *surfaces, device const Boundary *boundaries,
                          GeometryScene scene) {
    state.count = 0;
    if (!scene_has_boundaries || any(origin < p.boundary_lo.xyz) || any(origin > p.boundary_hi.xyz))
      return true;
    float3 direction = normalize(float3(.419f, .613f, .669f));
    if (p.boundary_lo.w != 0) {
      // Diffusion repeatedly classifies points immediately beside interfaces.
      // Query each containing bound independently: closest-hit replay cannot
      // see both members of an exactly coincident pair. The first outward
      // crossing of a closed, oriented body is an exit iff this point is inside.
      ray probe;
      probe.origin = origin;
      probe.direction = direction;
      probe.min_distance = FLT_MIN;
      probe.max_distance = INFINITY;
      float exit_distance[16];
      for (uint id = 1; id <= p.boundaries; ++id) {
        Boundary b = boundaries[id];
        if (any(origin < b.lo.xyz) || any(origin > b.hi.xyz))
          continue;
        uint body = id;
        auto hit = trace_geometry(probe, scene, quadrics, body);
        if (hit.type != intersection_type::none &&
            dot(direction, geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics,
                                            surfaces, hit.instance, scene.instances, probe.origin, probe.direction)
                               .face.geometric.xyz) > 0) {
          if (!cross_boundary(state, id, true))
            return false;
          // Replay farther exits first, matching nested/equal-priority order.
          uint slot = state.count - 1;
          while (slot && exit_distance[slot - 1] < hit.distance) {
            exit_distance[slot] = exit_distance[slot - 1];
            state.ids[slot] = state.ids[slot - 1];
            --slot;
          }
          exit_distance[slot] = hit.distance;
          state.ids[slot] = id;
        }
      }
      return true;
    }
    float radius = length(p.boundary_hi.xyz - p.boundary_lo.xyz);
    float distance = 2 * radius + 64 * FLT_EPSILON * max(1.f, length(origin));
    ray probe;
    probe.origin = origin - distance * direction;
    probe.direction = direction;
    probe.min_distance = FLT_MIN;
    // The classification endpoint is close to a surface by design. Metal's
    // float traversal distance can round across it, especially when the probe
    // started several object diameters away. Query past that uncertainty and
    // decide near-endpoint crossings from the local triangle plane instead.
    float endpoint_guard =
        64 * FLT_EPSILON *
        max(distance, max(1.f, max(abs(origin.x), max(abs(origin.y), abs(origin.z)))));
    probe.max_distance = distance + endpoint_guard;
    while (true) {
      auto hit = trace_geometry(probe, scene, quadrics);
      if (hit.type == intersection_type::none)
        return true;
      GeometrySurface contact = geometry_surface(hit.primitive_id, hit.coordinates, triangles,
                                                 quadrics, surfaces, hit.instance, scene.instances, probe.origin, probe.direction);
      Triangle t = contact.face;
      float facing = dot(direction, t.geometric.xyz);
      bool before_endpoint = hit.distance < distance - endpoint_guard ||
                             dot(origin - t.p[0].xyz, t.geometric.xyz) * facing >= 0;
      if (t.boundary && before_endpoint && !cross_boundary(state, t.boundary, facing < 0))
        return false;
      float next = as_type<float>(as_type<uint>(hit.distance) + 1u);
      if (!(next > probe.min_distance) || !isfinite(next))
        return false;
      probe.min_distance = next;
    }
  }

  void enqueue(uint pixel, Path path, thread const VolumeState &state, constant Params &p,
               device Path *paths, device VolumeState *membership, device uint *next_queue,
               device Counts &counts) {
    paths[pixel] = path;
    if (scene_has_boundaries)
      membership[pixel] = state;
    uint slot = atomic_fetch_add_explicit(&counts.next, 1, memory_order_relaxed);
    if (slot >= p.active)
      atomic_store_explicit(&counts.error, 1, memory_order_relaxed);
    else
      next_queue[slot] = pixel;
  }
  bool roulette(thread Path &path, uint pixel, constant Params &p, device const uint *directions,
                device const uint *blue, bool internal) {
    if (!(internal ? path.internal > 5 : path.depth > p.roulette))
      return true;
    float survival = min(1.f, max(path.beta.x, max(path.beta.y, path.beta.z)) * path.beta.w);
    if (!(survival > 0) || sample1(pixel, path.dimension, p, directions, blue) >= survival)
      return false;
    path.beta.xyz /= survival;
    return all(isfinite(path.beta));
  }

)MSL"
#include "wavefront_bsdf.inc"
#include "wavefront_diffusion.inc"
#include "wavefront_subsurface_sampling.inc"
#include "wavefront_openpbr.inc"
#include "wavefront_bump.inc"
#include "wavefront_media.inc"
#include "wavefront_atmosphere.inc"
#include "wavefront_haze.inc"
#include "wavefront_camera.inc"
#include "wavefront_motion.inc"
    R"MSL(

#define SCENE_ARGUMENTS                                                                            \
  device const Triangle *triangles [[buffer(0)]], device const Material *materials [[buffer(1)]],  \
      device const Texture *textures [[buffer(2)]], device const float4 *texels [[buffer(3)]],     \
      device const Light *lights [[buffer(4)]], device Path *paths [[buffer(5)]],                  \
      device Pixel *pixels [[buffer(6)]], device uint *ray_queue [[buffer(7)]],                    \
      device uint *next_queue [[buffer(8)]], device Hit *hit_queue [[buffer(9)]],                  \
      device Shadow *shadow_queue [[buffer(10)]], device Counts &counts [[buffer(11)]],            \
      constant Params &p [[buffer(12)]], device const uint *active [[buffer(13)]],                 \
      device const uint *directions [[buffer(14)]], device const uint *blue [[buffer(15)]],        \
      primitive_acceleration_structure triangle_scene [[buffer(16)]],                              \
      device const Boundary *boundaries [[buffer(17)]],                                            \
      device VolumeState *membership [[buffer(18)]], device atomic_uint *indirect [[buffer(19)]],  \
      device const float4 *material_data [[buffer(20)]],                                           \
      device const DiffusionExit *diffusion_exits [[buffer(21)]],                                  \
      device const float4 *diffusion_samples [[buffer(22)]],                                       \
      device const uint *diffusion_pairs [[buffer(23)]],                                           \
      device DiffusionEntry *diffusion_entries [[buffer(24)]],                                     \
      intersection_function_table<triangle_data> triangle_functions [[buffer(25)]],                \
      device const Surface *surfaces [[buffer(26)]],                                               \
      device Differential *differentials [[buffer(27)]],                                           \
      device const GeometryStorage *storage [[buffer(28)]],                                        \
      instance_acceleration_structure mixed_scene [[buffer(29)]],                                  \
      intersection_function_table<triangle_data, instancing> mixed_functions [[buffer(30)]],       \
      uint i [[thread_position_in_grid]]

)MSL"
#include "wavefront_atmosphere_paths.inc"
R"MSL(
  kernel void primary(SCENE_ARGUMENTS) {
    device const uchar *bytes = reinterpret_cast<device const uchar *>(storage);
    device const Quadric *quadrics =
        reinterpret_cast<device const Quadric *>(bytes + storage->quadric_offset);
    device const MeshInstance *instances =
        reinterpret_cast<device const MeshInstance *>(bytes + storage->instance_offset);
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions,
                           instances,
                           reinterpret_cast<device const float4 *>(bytes + storage->motion_offset),
                           reinterpret_cast<device const Quadric *>(bytes + storage->motion_offset + 48)};
    if (i == 0) {
      atomic_store_explicit(&counts.rays, p.active, memory_order_relaxed);
      atomic_store_explicit(&counts.hits, 0, memory_order_relaxed);
      atomic_store_explicit(&counts.shadows, 0, memory_order_relaxed);
      atomic_store_explicit(&counts.next, 0, memory_order_relaxed);
      atomic_store_explicit(&counts.error, 0, memory_order_relaxed);
      atomic_store_explicit(&counts.discarded, 0, memory_order_relaxed);
      atomic_store_explicit(&counts.scattering, 0, memory_order_relaxed);
      for (uint j = 0; j < 9; ++j)
        atomic_store_explicit(&indirect[j], 1, memory_order_relaxed);
      atomic_store_explicit(&indirect[0], max(1u, (p.active + 63) / 64), memory_order_relaxed);
    }
    if (i >= p.active)
      return;
    uint pixel = active[i], dim = 0;
    if (texture_features & 256) diffusion_entries[pixel] = {};
    if (scene_has_atmosphere && (p.camera.padding[1] & 2)) {
      device AirRay *air = reinterpret_cast<device AirRay *>(
          reinterpret_cast<device AtmosphereRequest *>(diffusion_entries + p.camera.padding[0]) + p.camera.padding[0]);
      air[pixel] = {};
    }
    float u = (pixel / p.height + sample1(pixel, dim, p, directions, blue)) / p.width;
    float v = (pixel % p.height + sample1(pixel, dim, p, directions, blue)) / p.height;
    if (p.camera.flip_x)
      u = 1 - u;
    float lu = sample1(pixel, dim, p, directions, blue),
          lv = sample1(pixel, dim, p, directions, blue);
    // Match the native thin-lens camera's uniform disk; unused time retains its dimension.
    float time=sample1(pixel, dim, p, directions, blue)*p.camera.shutter.z;
    if(p.camera.shutter.w>=0) time=p.camera.shutter.w;
    float r = sqrt(lu) * p.camera.right.w, phi = 2 * pi * lv;
    float3 lens = r * (cos(phi) * p.camera.right.xyz + sin(phi) * p.camera.up.xyz);
    float3 film = p.camera.lower_left.xyz + u * p.camera.horizontal.xyz + v * p.camera.vertical.xyz;
    Path path;
    path.origin = float4(p.camera.orthographic ? film : p.camera.origin.xyz + lens, 0);
    path.direction =
        float4(normalize(p.camera.orthographic ? p.camera.forward.xyz : film - path.origin.xyz), 0);
    path.beta = float4(1);
    path.previous = path.origin;
    path.dimension = dim;
    path.depth = 0;
    path.internal = 0;
    path.flags = 1; // delta predecessor, no guide yet.
    bool complex_camera=p.camera.orthographic>=2 || p.camera.moving;
    float camera_weight=1;
    device const float4 *camera_data=reinterpret_cast<device const float4 *>(bytes+storage->camera_offset);
    if(complex_camera) {
      float3 origin,direction;
      float weight=camera_ray(p.camera,float2(u,v),float2(lu,lv),time,camera_data,origin,direction);
      camera_weight=weight;
      path.origin=float4(origin,0);path.direction=float4(direction,0);path.beta=float4(weight>0 ? 1.f : 0.f);
      path.previous=path.origin;
    }
    if (texture_features & 4) {
      float scale = p.camera.up.w;
      float3 dx = p.camera.horizontal.xyz * (p.camera.flip_x ? -1.f : 1.f) / p.width * scale;
      float3 dy = p.camera.vertical.xyz / p.height * scale;
      Differential d;
      d.ox = float4(path.origin.xyz + (p.camera.orthographic ? dx : float3(0)), 1);
      d.oy = float4(path.origin.xyz + (p.camera.orthographic ? dy : float3(0)), 0);
      d.dx = float4(
          p.camera.orthographic ? path.direction.xyz : normalize(film + dx - path.origin.xyz), 0);
      d.dy = float4(
          p.camera.orthographic ? path.direction.xyz : normalize(film + dy - path.origin.xyz), 0);
      differentials[pixel] = d;
      if(complex_camera) {
        d.ox.w=0;
        float3 origins[2],directions[2];bool valid=path.beta.x>0;
        for(uint axis=0;axis<2 && valid;++axis) {
          bool found=false;
          for(int sign=1;sign>=-1;sign-=2) {
            float2 shifted=float2(u,v);
            shifted[axis]+=sign*(axis ? 1.f/p.height : (p.camera.flip_x ? -1.f : 1.f)/p.width);
            if(camera_ray(p.camera,shifted,float2(lu,lv),time,camera_data,origins[axis],directions[axis])<=0) continue;
            origins[axis]=path.origin.xyz+sign*scale*(origins[axis]-path.origin.xyz);
            directions[axis]=path.direction.xyz+sign*scale*(directions[axis]-path.direction.xyz);
            found=true;break;
          }
          valid=found;
        }
        if(valid) d={float4(origins[0],1),float4(origins[1],0),float4(directions[0],0),float4(directions[1],0)};
        differentials[pixel]=d;
      }
    }
    if (scene_has_boundaries) {
      VolumeState state = {};
      if (!initial_membership(path.origin.xyz, state, p, triangles, quadrics, surfaces, boundaries,
                              scene)) {
        path.flags |= 4;
      }
      membership[pixel] = state;
    }
    paths[pixel] = path;
    pixels[pixel] = {float4(0), float4(0), float4(0)};
    pixels[pixel].normal.w=camera_weight;
    ray_queue[i] = pixel;
  }

  kernel void intersect_paths(SCENE_ARGUMENTS) {
    device AtmosphereRequest *air_requests = reinterpret_cast<device AtmosphereRequest *>(
        diffusion_entries + p.camera.padding[0]);
    device const uchar *bytes = reinterpret_cast<device const uchar *>(storage);
    device const Quadric *quadrics =
        reinterpret_cast<device const Quadric *>(bytes + storage->quadric_offset);
    device const MeshInstance *instances =
        reinterpret_cast<device const MeshInstance *>(bytes + storage->instance_offset);
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions,
                           instances,
                           reinterpret_cast<device const float4 *>(bytes + storage->motion_offset),
                           reinterpret_cast<device const Quadric *>(bytes + storage->motion_offset + 48)};
    if (i >= atomic_load_explicit(&counts.rays, memory_order_relaxed))
      return;
    uint pixel = ray_queue[i];
    Path path = paths[pixel];
    if (scene_has_atmosphere)
      air_requests[pixel].event = 0;
    if(all(path.beta.xyz==0)) return;
    if (path.flags & 4) {
      atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
      return;
    }
    VolumeState state = {};
    state.count = 0;
    if (scene_has_boundaries)
      state = membership[pixel];
    uint glass = active_boundary(state, boundaries, materials);
    uint body = active_medium(state, boundaries, materials);
    bool scattering = body && boundaries[body].subsurface;
    device const Medium *media = reinterpret_cast<device const Medium *>(bytes + storage->medium_offset);
    device const float *fields = reinterpret_cast<device const float *>(bytes + storage->field_offset);
    Medium medium = {};
    if ((texture_features & 16) && body)
      medium = media[boundaries[body].medium];
    bool heterogeneous = body && medium.controls.w != 0;
    if (body && boundaries[body].subsurface == 2) {
      if (!(path.flags & 8)) {
        atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
        return;
      }
      uint slot = atomic_fetch_add_explicit(&counts.hits, 1, memory_order_relaxed);
      if (slot >= p.active) {
        atomic_store_explicit(&counts.error, 1, memory_order_relaxed);
        return;
      }
      hit_queue[slot] = {pixel, UINT_MAX - 1, 0, 0, 0, 0, {0, 0}};
      atomic_fetch_max_explicit(&indirect[3], (slot + 64) / 64, memory_order_relaxed);
      return;
    }
    float3 rates = scattering ? extinction(body, boundaries, materials) : float3(0);
    bool guided_flight = (texture_features & 256) && (path.flags & 32) && body &&
        boundaries[body].subsurface == 1 && boundaries[body].radius.w != 0 &&
        as_type<uint>(diffusion_entries[pixel].position.w) == body;
    float3 guided_rates = rates, ordinary_weights = float3(1.f/3), guided_weights = 0;
    if (guided_flight)
      subsurface_flight_mixture(boundaries[body].color.xyz,diffusion_entries[pixel].normal.xyz,
                                path.direction.xyz,rates,guided_rates,ordinary_weights,guided_weights);
    float flight = INFINITY;
    float majorant = heterogeneous ? max(rates.x, max(rates.y, rates.z)) * medium.hi.w : 0;
    if (heterogeneous) {
      flight = medium_flight(medium,path.origin.xyz,path.direction.xyz,max(rates.x,max(rates.y,rates.z)),
                             sample1(pixel,path.dimension,p,directions,blue),fields,majorant);
    } else if (scattering) {
      // A fresh RGB-channel mixture at every segment, with the exact mixture
      // density below. Zero-rate channels retain their mass at infinity.
      float selector = sample1(pixel, path.dimension, p, directions, blue);
      uint channel = min(2u, uint(3 * selector));
      float rate = rates[channel];
      if (guided_flight) {
        float cumulative=0;
        rate=guided_rates.z;
        for (uint c=0;c<6;++c) {
          cumulative += c<3 ? ordinary_weights[c] : guided_weights[c-3];
          if (selector < cumulative) { rate=c<3 ? rates[c] : guided_rates[c-3]; break; }
        }
      }
      float u = sample1(pixel, path.dimension, p, directions, blue);
      if (rate > 0)
        flight = -log(1 - u) / rate;
    }
    ray r;
    r.origin = path.origin.xyz;
    r.direction = path.direction.xyz;
    r.min_distance = FLT_MIN;
    float magnitude =
        max(1.f, max(flight, max(abs(r.origin.x), max(abs(r.origin.y), abs(r.origin.z)))));
    // Match SubsurfaceFlightLimit: query beyond the sampled endpoint only to
    // obtain its positional error. This guard never changes the physical flight.
    r.max_distance = isfinite(flight) ? flight + 32 * FLT_EPSILON * magnitude : INFINITY;
    auto hit = trace_geometry(r, scene, quadrics);
    uint departed = as_type<uint>(path.origin.w);
    while (hit.type != intersection_type::none) {
      bool skip = departed && !is_quadric(departed - 1) && hit.primitive_id == departed - 1 &&
                  hit.instance == as_type<uint>(path.previous.w);
      if (!skip && (texture_features & 64) && !(path.flags & 16)) {
        auto contact = geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics,
                                        surfaces, hit.instance, scene.instances, r.origin, r.direction);
        Material emitter = materials[contact.face.material];
        skip = emitter.emissive == 1 && emitter.extra.w != 0;
      }
      if (!skip && (texture_features & 8)) {
        GeometrySurface contact =
            geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces,
                             hit.instance, scene.instances, r.origin, r.direction);
        Triangle face = contact.face;
        uint alpha = surfaces[as_type<uint>(face.uv2.w)].alpha;
        if (alpha != UINT_MAX) {
          float2 uv = hit.coordinates.xy;
          float3 bary = float3(1 - uv.x - uv.y, uv);
          float3 position = contact.position;
          float opacity = texture_value(alpha, surface_context(face, bary, position, surfaces),
                                        textures, texels)
                              .x;
          skip = opacity < sample1(pixel, path.dimension, p, directions, blue);
        }
      }
      if (!skip)
        break;
      float next = as_type<float>(as_type<uint>(hit.distance) + 1u);
      if (!(next > r.min_distance)) {
        atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
        return;
      }
      r.min_distance = next;
      hit = trace_geometry(r, scene, quadrics);
    }
    float surface = hit.type == intersection_type::none ? INFINITY : hit.distance;
    bool collision = flight < surface;
    float distance = collision ? flight : surface;
    // GPU finite segments need atmosphere transport only when haze is enabled.
    // Escaping rays still evaluate the sky from the last real lighting vertex.
    if (scene_has_atmosphere &&
        (!(p.camera.padding[1] & 2) || (p.camera.padding[1] & 1) || !isfinite(distance))) {
      air_requests[pixel].origin = path.origin;
      air_requests[pixel].direction = float4(r.direction, distance);
      air_requests[pixel].weight = path.beta;
      air_requests[pixel].weight.w = path.direction.w;
      air_requests[pixel].radiance = pixels[pixel].radiance;
      air_requests[pixel].body = body;
      air_requests[pixel].glass = glass;
      air_requests[pixel].event = isfinite(distance) ? 1 : 2;
      air_requests[pixel].padding = path.flags;
    }
    bool null_collision = false;
    if (heterogeneous && collision) {
      float3 uvw = medium_coordinates(medium, r.origin + flight * r.direction);
      float density = medium_density(medium, uvw, fields);
      float3 sigma_a = boundaries[body].sigma_a.xyz * density;
      float3 sigma_s = boundaries[body].sigma_s.xyz * density;
      if (path.depth < p.max_depth)
        pixels[pixel].radiance.xyz += path.beta.xyz * sigma_a * medium_emission(medium, uvw, fields) / majorant;
      float accept = clamp(max(rates.x, max(rates.y, rates.z)) * density / majorant, 0.f, 1.f);
      null_collision = sample1(pixel, path.dimension, p, directions, blue) >= accept;
      path.beta.xyz *= null_collision ? (1 - (sigma_a + sigma_s) / majorant) / (1 - accept)
                                      : sigma_s / (majorant * accept);
    } else if (scattering && !heterogeneous) {
      // Remove a common exponential before normalizing. This is the same ratio
      // without underflow in an optically thick RGB segment.
      float low = min(rates.x, min(rates.y, rates.z));
      if (guided_flight) low=min(low,min(guided_rates.x,min(guided_rates.y,guided_rates.z)));
      float3 relative = isfinite(distance) ? exp(-(rates - low) * distance)
                                           : select(float3(0), float3(1), rates == 0);
      float density = dot(relative * (collision ? rates : float3(1)), float3(1.f / 3));
      if (guided_flight) {
        float3 guided_relative = isfinite(distance) ? exp(-(guided_rates-low)*distance)
            : select(float3(0),float3(1),guided_rates == 0);
        density = dot(relative*(collision ? rates : float3(1)),ordinary_weights) +
                  dot(guided_relative*(collision ? guided_rates : float3(1)),guided_weights);
      }
      if (!(density > 0))
        return;
      if ((texture_features & 16) && collision && boundaries[body].medium && path.depth < p.max_depth)
        pixels[pixel].radiance.xyz += path.beta.xyz * relative * boundaries[body].sigma_a.xyz *
                                     medium.emission.xyz / density;
      path.beta.xyz *= relative * (collision ? boundaries[body].sigma_s.xyz : float3(1)) / density;
    }
    if (glass) {
      float3 absorption = materials[boundaries[glass].material].attenuation.xyz;
      path.beta.xyz *= isfinite(distance) ? exp(-absorption * distance)
                                          : select(float3(0), float3(1), absorption == 0);
    }
    if (!all(isfinite(path.beta))) {
      atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
      return;
    }
    if (all(path.beta.xyz == 0))
      return;
    if (collision) {
      float3 position = fma(float3(flight), r.direction, r.origin);
      if (hit.type != intersection_type::none) {
        // Match SubsurfaceCollisionPoint: rounding a near-boundary collision
        // outside its medium corrupts subsequent membership. Keep the sampled
        // distance/PDF and move only an uncertain position onto the incident side.
        GeometrySurface contact =
            geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces,
                             hit.instance, scene.instances, r.origin, r.direction);
        Triangle endpoint = contact.face;
        float3 surface_position = contact.position, error = contact.error;
        error += 3 * FLT_EPSILON * (abs(r.origin) + abs(flight * r.direction));
        float3 g = endpoint.geometric.xyz, delta = position - surface_position;
        float side = dot(r.direction, g), plane = dot(delta, g);
        if (dot(delta, delta) <= dot(error, error)) {
          position = offset_origin(surface_position, error, side > 0 ? -g : g);
        } else if (side != 0 && plane * side >= 0 && abs(plane) <= dot(error, abs(g))) {
          position -= plane * g;
          position = offset_origin(position, error, side > 0 ? -g : g);
        }
      }
      // A sampled tiny step can round to zero. Reject only when even one mean
      // free path is unresolvable; repeated zero steps must not hang a preview.
      float rate = max(rates.x, max(rates.y, rates.z));
      if (all(position == r.origin) &&
          all(fma(float3(1 / rate), r.direction, r.origin) == r.origin)) {
        atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
        return;
      }
      path.origin.xyz = position;
      path.origin.w = 0;
      if (!null_collision) {
        ++path.internal;
        atomic_fetch_add_explicit(&counts.scattering, 1, memory_order_relaxed);
      }
    } else if (hit.type == intersection_type::none) {
      if (scene_has_atmosphere) {
        paths[pixel] = path;
        return;
      }
      if ((p.camera.padding[1] & 4) && path.depth == 0)
        return;
      float3 radiance = 0;
      for (uint j = 0; j < p.environments; ++j) {
        Light l = lights[j];
        if (l.type > 1)
          continue;
        float weight = (path.flags & 1)
                           ? 1
                           : mis(path.direction.w, light_pdf(l, r.direction, texels) / p.lights);
        radiance += weight * environment(l, r.direction, path.previous.xyz, textures, texels);
      }
      pixels[pixel].radiance.xyz += path.beta.xyz * radiance;
      return;
    }
    paths[pixel] = path;
    uint slot = atomic_fetch_add_explicit(&counts.hits, 1, memory_order_relaxed);
    if (slot >= p.active) {
      atomic_store_explicit(&counts.error, 1, memory_order_relaxed);
      return;
    }
    hit_queue[slot] = {pixel,
                       collision ? (null_collision ? UINT_MAX - 2 : UINT_MAX) : hit.primitive_id,
                       collision ? 0 : hit.coordinates.x,
                       collision ? 0 : hit.coordinates.y,
                       distance,
                       collision ? 0 : hit.coordinates.z,
                       {collision ? 0 : hit.instance, 0}};
    atomic_fetch_max_explicit(&indirect[3], (slot + 64) / 64, memory_order_relaxed);
  }

  template <ushort features> kernel void shade(SCENE_ARGUMENTS) {
    device const uchar *bytes = reinterpret_cast<device const uchar *>(storage);
    device const LightNode *light_nodes =
        reinterpret_cast<device const LightNode *>(bytes + storage->light_offset);
    device const Quadric *quadrics =
        reinterpret_cast<device const Quadric *>(bytes + storage->quadric_offset);
    device const MeshInstance *instances =
        reinterpret_cast<device const MeshInstance *>(bytes + storage->instance_offset);
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions,
                           instances,
                           reinterpret_cast<device const float4 *>(bytes + storage->motion_offset),
                           reinterpret_cast<device const Quadric *>(bytes + storage->motion_offset + 48)};
    if (i >= atomic_load_explicit(&counts.hits, memory_order_relaxed))
      return;
    Hit h = hit_queue[i];
    bool is_diffusion = h.primitive == UINT_MAX - 1;
    if (is_diffusion != (features == 4))
      return;
    bool is_pbr =
        h.primitive < UINT_MAX - 2 &&
        materials[geometry_material(h.primitive, triangles, quadrics)].type == bsdf_openpbr;
    if (is_pbr != (features == 2))
      return;
    uint pixel = h.pixel;
    Path path = paths[pixel];
    uint dim = path.dimension;
    VolumeState state = {};
    state.count = 0;
    if (scene_has_boundaries)
      state = membership[pixel];
    if (h.primitive == UINT_MAX - 2) {
      enqueue(pixel, path, state, p, paths, membership, next_queue, counts);
      return;
    }
    VolumeState adjacent = state;
    bool volume = scene_has_boundaries && h.primitive == UINT_MAX, interface = false,
         internal = volume;
    uint body = volume ? active_medium(state, boundaries, materials)
                       : active_boundary(state, boundaries, materials);
    float3 position = path.origin.xyz, error = 0, g = 0, n = 0, rho = 0, wo = -path.direction.xyz;
    float eta = 1, roughness = 0, anisotropy = body ? boundaries[body].sigma_s.w : 0;
    bool guided_scatter = (texture_features & 256) && volume &&
        boundaries[body].subsurface == 1 && boundaries[body].radius.w != 0 &&
        as_type<uint>(diffusion_entries[pixel].position.w) == body;
    bool entering = false;
    Triangle t = {};
    GeometrySurface geometry = {};
    Material m = {};
    SurfaceBsdf bsdf = {};
    SurfaceDerivatives derivatives = {};
    MappedDiffuse mapped;
    float3 base_normal = 0;
    bool mapped_surface = false;
    PbrState pbr;
    DiffusionExit exit_lobe;
    bool extended = false;
    if (features == 4) {
      if (path.depth >= p.max_depth)
        return;
      ++path.depth;
      DiffusionSample exit;
      bool invalid = false;
      bool found = diffusion_surface(body, diffusion_entries[pixel], pixel, dim, p, directions,
                                     blue, triangles, quadrics, surfaces, boundaries, materials,
                                     scene, exit, invalid);
      if (invalid)
        atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
      if (!found)
        return;
      atomic_fetch_add_explicit(&counts.scattering, 1, memory_order_relaxed);
      path.flags &= ~8u;
      if (!(path.flags & 2)) {
        pixels[pixel].normal.xyz = 0;
        pixels[pixel].albedo = float4(boundaries[body].color.xyz, 0);
        path.flags |= 2;
      }
      path.beta.xyz *= exit.weight;
      h = exit.hit;
      geometry = geometry_surface(h.primitive, float3(h.u, h.v, h.w), triangles, quadrics, surfaces,
                                  h.padding[0], scene.instances, path.origin.xyz, path.direction.xyz);
      t = geometry.face;
      position = exit.position;
      error = exit.error;
      g = n = exit.normal;
      state = adjacent = exit.outside;
      uint outside = active_boundary(state, boundaries, materials);
      uint key = outside ? boundaries[outside].material + 1 : 0;
      exit_lobe = diffusion_exits[diffusion_pairs[boundaries[body].diffusion_offset + key]];
      if (outside && !boundaries[outside].subsurface) {
        Material glass = materials[boundaries[outside].material];
        if (glass.type != bsdf_openpbr)
          path.beta.xyz *=
              texture_value(glass.texture, triangle_uv(t, float3(1 - h.u - h.v, h.u, h.v)),
                            position, textures, texels);
      }
    } else if (volume) {
      if (boundaries[body].subsurface == 3) {
        if (path.depth >= p.max_depth)
          return;
        ++path.depth;
        internal = false;
      }
      if (!(path.flags & 2)) {
        float3 st = boundaries[body].sigma_a.xyz + boundaries[body].sigma_s.xyz;
        pixels[pixel].normal.xyz = 0;
        pixels[pixel].albedo = float4(
            select(float3(0), boundaries[body].sigma_s.xyz / max(st, float3(FLT_MIN)), st > 0), 0);
        path.flags |= 2;
      }
    } else {
      geometry = geometry_surface(h.primitive, float3(h.u, h.v, h.w), triangles, quadrics, surfaces,
                                  h.padding[0], scene.instances, path.origin.xyz, path.direction.xyz);
      t = geometry.face;
      m = materials[t.material];
      float3 bary = float3(1 - h.u - h.v, h.u, h.v);
      position = geometry.position;
      error = geometry.error;
      g = t.geometric.xyz;
      n = geometry.normal;
      if ((texture_features & 16) && t.boundary && boundaries[t.boundary].subsurface == 3 &&
          boundaries[t.boundary].radius.w == 0) {
        entering = dot(g, wo) > 0;
        if (!adjacent_membership(t.boundary, entering, position, error, g, adjacent, p,
                                 triangles, quadrics, surfaces, boundaries, scene)) {
          atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
          return;
        }
        path.origin = float4(offset_origin(position, error, entering ? -g : g),
                             as_type<float>(h.primitive + 1));
        path.previous.w = as_type<float>(h.padding[0]);
        enqueue(pixel, path, adjacent, p, paths, membership, next_queue, counts);
        return;
      }
      TextureContext texture_hit = surface_context(t, bary, position, surfaces);
      base_normal = n;
      if (texture_features & 4) {
        Surface mapping = surfaces[as_type<uint>(t.uv2.w)];
        derivatives = surface_derivatives(t, mapping, texture_hit, n, differentials[pixel]);
        n = bump_normal(t, mapping, texture_hit, n, derivatives, textures, texels);
        mapped_surface = mapping.bump != UINT_MAX;
      }
      rho = texture_value(m.texture, texture_hit, textures, texels);
      float3 feature_normal = mapped_surface ? n :
          consistent_normal(t, surfaces[as_type<uint>(t.uv2.w)], bary, n, wo, texels);
      // Physical diffuse/OpenPBR and sheet transmission keep their own frame.
      // Legacy reflective BSDFs and emission use triangle::hit's normal.
      if (!mapped_surface && (m.emissive == 1 ||
          (m.type != 0 && m.type != bsdf_openpbr && m.type != bsdf_translucent)))
        n = feature_normal;
      float mapped_roughness = -1;
      if (m.second_texture != UINT_MAX && m.type != bsdf_translucent) {
        float3 rough = texture_value(m.second_texture, texture_hit, textures, texels);
        if (m.type == bsdf_openpbr)
          mapped_roughness = rough.x;
        else
          m.parameters.xy = rough.xy;
      }
      if (m.emissive == 1) {
        if (dot(n, wo) > 0) {
          float pdf =
              geometry_light_pdf(h.primitive, geometry, path.previous.xyz, quadrics) *
              light_selection_pdf(t.light, path.previous.xyz, p, light_nodes);
          float weight = (path.flags & 1) ? 1 : mis(path.direction.w, pdf);
          pixels[pixel].radiance.xyz += path.beta.xyz * rho * m.diffuse.w * weight;
        }
        return;
      }
      if (features == 2) {
        // Priority still owns the medium. Only a visible interface invokes the
        // complete OpenPBR surface; never route this subclass through glass BSDFs.
        entering = dot(g, wo) > 0;
        uint after = body;
        if (scene_has_boundaries && m.dielectric) {
          if (!adjacent_membership(t.boundary, entering, position, error, g, adjacent, p, triangles,
                                   quadrics, surfaces, boundaries, scene)) {
            atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
            return;
          }
          after = active_boundary(adjacent, boundaries, materials);
          bool hidden = body == after ||
                        (body && after && boundaries[body].material == boundaries[after].material);
          if (hidden) {
            path.origin = float4(offset_origin(position, error, entering ? -g : g),
                                 as_type<float>(h.primitive + 1));
            path.previous.w = as_type<float>(h.padding[0]);
            remember_subsurface_entry(path, pixel, adjacent, body, t.boundary, position, g, boundaries, materials,
                                     diffusion_entries);
            enqueue(pixel, path, adjacent, p, paths, membership, next_queue, counts);
            return;
          }
        }
        float exterior = boundary_ior(entering ? body : after, boundaries, materials);
        pbr = prepare_pbr(m, rho, wo, n, g, geometry.tangent, exterior, mapped_roughness,
                          material_data);
        float weight = 1;
        if (m.emissive && !(path.flags & 1)) {
          float pdf =
              geometry_light_pdf(h.primitive, geometry, path.previous.xyz, quadrics) *
              light_selection_pdf(t.light, path.previous.xyz, p, light_nodes);
          weight = mis(path.direction.w, pdf);
        }
        pixels[pixel].radiance.xyz += path.beta.xyz * pbr.prepared.emission * weight;
        internal = m.dielectric && !entering;
        if (internal)
          ++path.internal;
        else {
          if (path.depth >= p.max_depth)
            return;
          ++path.depth;
        }
        // Mixed lobes set a guide only after a non-discrete sample below.
      } else if (scene_has_boundaries && m.dielectric) {
        entering = dot(g, wo) > 0;
        if (!adjacent_membership(t.boundary, entering, position, error, g, adjacent, p, triangles,
                                 quadrics, surfaces, boundaries, scene)) {
          atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
          return;
        }
        uint after = active_boundary(adjacent, boundaries, materials);
        bool hidden = body == after ||
                      (body && after && boundaries[body].material == boundaries[after].material);
        eta =
            boundary_ior(after, boundaries, materials) / boundary_ior(body, boundaries, materials);
        if (hidden || eta == 1) {
          // Hidden and index-matched crossings retain MIS history and depth.
          if (!hidden && entering && !boundaries[t.boundary].subsurface)
            path.beta.xyz *= rho;
          path.origin = float4(offset_origin(position, error, entering ? -g : g),
                               as_type<float>(h.primitive + 1));
          path.previous.w = as_type<float>(h.padding[0]);
          path.dimension = dim;
          remember_subsurface_entry(path, pixel, adjacent, body, t.boundary, position, g, boundaries, materials,
                                   diffusion_entries);
          enqueue(pixel, path, adjacent, p, paths, membership, next_queue, counts);
          return;
        }
        internal = boundaries[t.boundary].subsurface && !entering;
        if (p.boundary_lo.w != 0)
          internal = body && boundaries[body].subsurface && !entering;
        if (internal)
          ++path.internal;
        else {
          if (path.depth >= p.max_depth)
            return;
          ++path.depth;
        }
        roughness = boundaries[t.boundary].sigma_a.w;
        // Tessellated analytic boundaries retain smooth refraction normals;
        // native triangle meshes use their actual geometric interface normals.
        n = ((t.flags & 1) || mapped_surface) ? (dot(n, g) < 0 ? -n : n) : g;
        if (!entering) {
          g = -g;
          n = -n;
        }
        if (roughness == 0) {
          float3 wi;
          float weight, pdf;
          if (!boundary_sample(wo, n, eta, 0, sample1(pixel, dim, p, directions, blue), 0, 0, wi,
                               weight, pdf))
            return;
          bool transmission = dot(wi, g) < 0;
          if ((dot(wi, n) < 0) != transmission)
            return;
          path.beta.xyz *= weight;
          if (entering && !boundaries[t.boundary].subsurface)
            path.beta.xyz *= rho;
          if (transmission) {
            state = adjacent;
            path.beta.w *= eta * eta;
          }
          path.origin = float4(offset_origin(position, error, transmission ? -g : g),
                               as_type<float>(h.primitive + 1));
          path.previous = float4(position, as_type<float>(h.padding[0]));
          if (texture_features & 4)
            differentials[pixel] =
                propagate_differential(differentials[pixel], derivatives, n, path.direction.xyz,
                                       path.origin.xyz, wi, transmission, eta);
          path.direction = float4(wi, pdf);
          remember_subsurface_entry(path, pixel, state, body, t.boundary, position, g, boundaries, materials,
                                   diffusion_entries);
          path.flags |= 1;
          path.dimension = dim;
          if (roulette(path, pixel, p, directions, blue, internal))
            enqueue(pixel, path, state, p, paths, membership, next_queue, counts);
          return;
        }
        interface = true;
        if (!(path.flags & 2)) {
          pixels[pixel].normal.xyz = n;
          pixels[pixel].albedo = float4(1, 1, 1, 0);
          path.flags |= 2;
        }
      } else {
        if (path.depth >= p.max_depth)
          return;
        ++path.depth;
        rho = select(float3(0), clamp(rho, 0.f, 1.f), isfinite(rho));
        extended = (features == 1 || features == 3) && m.type != 0;
        if (!extended && dot(g, wo) * dot(base_normal, wo) <= 0)
          return;
        if (!extended && dot(g, wo) < 0) {
          g = -g;
          n = -n;
        }
        if ((texture_features & 4) && !extended)
          mapped = mapped_diffuse(dot(base_normal, wo) < 0 ? -base_normal : base_normal, n, wo);
        if (extended) {
          bsdf.material = m;
          bsdf.frame = normal_frame(m.type == bsdf_translucent ? g : n);
          if (m.type == bsdf_hair) {
            float3 tangent = normalize(geometry.tangent);
            float3 normal = n - tangent * dot(n, tangent);
            normal = dot(normal, normal) < 1e-12f ? normal_frame(tangent).y : normalize(normal);
            bsdf.frame = {tangent, cross(normal, tangent), normal};
          }
          bsdf.view = bsdf.frame.local(wo);
          bsdf.rho = rho;
          bsdf.distance = h.distance;
          bsdf.impact = clamp(2 * triangle_uv(t, bary).y - 1, -1.f, 1.f);
          if (m.type == bsdf_translucent) {
            bsdf.transmission =
                clamp(texture_value(m.second_texture, texture_hit, textures, texels), 0.f, 1.f);
            float3 total = max(bsdf.rho + bsdf.transmission, float3(1));
            bsdf.rho /= total;
            bsdf.transmission /= total;
            rho = bsdf.rho + bsdf.transmission;
          }
        }
        if (!(path.flags & 2) && !(m.type == bsdf_metal && m.parameters.z == 0) &&
            !(m.type == bsdf_transmission && m.eta.x == 1)) {
          pixels[pixel].normal.xyz = feature_normal;
          pixels[pixel].albedo = float4(rho, 0);
          path.flags |= 2;
        }
      }
    }
    // One light proposal per vertex. Each emitter has its own competing BSDF
    // technique; evaluating emission individually on a miss avoids mixture bias.
    if (p.lights) {
      uint choice = select_light(sample1(pixel, dim, p, directions, blue),position,p,light_nodes);
      float selection_pdf = light_selection_pdf(choice,position,p,light_nodes);
      Light l = lights[choice];
      bool continuous = l.type == 0 && l.alias_offset != UINT_MAX;
      uint environment_bits = sample_bits(pixel, dim, p, directions, blue, continuous);
      float u = min(float(environment_bits) * 0x1p-32f, 0x1.fffffep-1f);
      float v = sample1(pixel, dim, p, directions, blue, continuous);
      float3 wi, emitted = 0;
      float3 light_target = 0, light_error = 0, light_normal = 0;
      float pdf = 0, maximum = INFINITY;
      if (l.type == 0) {
        wi = sample_environment(l, environment_bits, v, texels);
        emitted = environment(l, wi, position, textures, texels);
        pdf = light_pdf(l, wi, texels) * selection_pdf;
      } else if (l.type == 1) {
        float d = u * l.forward.w, s = sqrt(d * (2 - d)), phi = 2 * pi * v;
        wi = (1 - d) * l.forward.xyz + s * (cos(phi) * l.right.xyz + sin(phi) * l.up.xyz);
        emitted = environment(l, wi, position, textures, texels);
        pdf = light_pdf(l, wi, texels) * selection_pdf;
      } else if (l.type == 2) {
        light_target = l.forward.xyz;
        float3 delta = l.forward.xyz - position;
        float distance2 = dot(delta, delta);
        maximum = sqrt(distance2);
        wi = delta / maximum;
        float falloff = 1;
        if (l.flags) {
          float cosine = dot(-wi, l.right.xyz);
          if (cosine <= l.right.w)
            falloff = 0;
          else if (cosine < l.up.w) {
            float f = (cosine - l.right.w) / (l.up.w - l.right.w);
            falloff = f * f * f * f;
          }
        }
        emitted = l.color.xyz * falloff / distance2;
        pdf = selection_pdf;
      } else {
        GeometrySurface source =
            sample_geometry_light(l.primitive, float2(u, v), position, triangles, quadrics,
                                  surfaces, l.instance, scene.instances, &counts.csg_failed);
        Triangle emitter = source.face;
        Material lm = materials[emitter.material];
        float3 b = float3(1 - source.coordinates.x - source.coordinates.y, source.coordinates.xy);
        float3 target = source.position;
        light_target = target;
        light_error = source.error;
        light_normal = source.face.geometric.xyz;
        float3 delta = target - position;
        maximum = length(delta);
        wi = delta / maximum;
        float3 en = source.normal;
        float cosine = dot(emitter.geometric.xyz, -wi);
        if (dot(en, -wi) > 0 && abs(cosine) > 0) {
          if (features >= 2 && lm.type == bsdf_openpbr) {
            // Emission is view dependent under OpenPBR's coating layers.
            auto emitter_pbr = prepare_pbr(
                lm,
                texture_value(lm.texture, surface_context(emitter, b, target, surfaces), textures,
                              texels),
                -wi, en, emitter.geometric.xyz, source.tangent,
                boundary_ior(body, boundaries, materials),
                lm.second_texture == UINT_MAX
                    ? -1.f
                    : texture_value(lm.second_texture,
                                    surface_context(emitter, b, target, surfaces), textures, texels)
                          .x,
                material_data);
            emitted = emitter_pbr.prepared.emission;
          } else
            emitted = texture_value(lm.texture, surface_context(emitter, b, target, surfaces),
                                    textures, texels) *
                      lm.diffuse.w;
          if (texture_features & 8) {
            uint alpha = surfaces[as_type<uint>(emitter.uv2.w)].alpha;
            if (alpha != UINT_MAX)
              emitted *= clamp(texture_value(alpha, surface_context(emitter, b, target, surfaces),
                                             textures, texels)
                                   .x,
                               0.f, 1.f);
          }
          pdf = geometry_light_pdf(l.primitive, source, position, quadrics) * selection_pdf;
        }
      }
      float cosine = dot(n, wi);
      BsdfValue mapped_value = {0, 0};
      if ((texture_features & 4) && !volume && !interface && !extended && features != 2 &&
          features != 4)
        mapped_value = mapped_eval(mapped, wi, rho, m);
      BsdfValue surface = extended ? surface_eval(bsdf, wi) : BsdfValue{0, 0};
      if (features == 2)
        surface = pbr_eval(pbr, wi);
      if (features == 4)
        surface = diffusion_eval(exit_lobe, n, wi, diffusion_samples);
      float2 boundary_value = interface ? boundary_eval(wo, wi, n, eta, roughness) : float2(0);
      float competing = volume ? phase_pdf(-wo, wi, anisotropy)
                               : (interface ? boundary_value.y
                                  : (extended || features == 2 || features == 4)
                                      ? surface.pdf
                                      : ((texture_features & 4) ? mapped_value.pdf : cosine / pi));
      float3 value = volume ? float3(competing)
                            : (interface ? float3(boundary_value.x)
                               : (extended || features == 2 || features == 4)
                                   ? surface.value
                                   : ((texture_features & 4) ? mapped_value.value
                                                             : brdf(wo, wi, n, rho, m) * cosine));
      if (guided_scatter)
        competing = subsurface_direction_pdf(boundaries[body].color.xyz,
                                             diffusion_entries[pixel].normal.xyz,wi);
      if (pdf > 0 && any(value > 0) &&
          (volume || interface || extended || features == 2 || features == 4 ||
           (((texture_features & 4) ? dot(mapped.g, wi) > 0 : cosine > 0) && dot(g, wi) > 0)) &&
          any(emitted > 0)) {
        float weight = l.type == 2 ? 1 : mis(pdf, competing);
        float3 contribution = path.beta.xyz * value * emitted * (weight / pdf);
        uint slot = atomic_fetch_add_explicit(&counts.shadows, 1, memory_order_relaxed);
        if (slot >= p.active)
          atomic_store_explicit(&counts.error, 1, memory_order_relaxed);
        else {
          float3 origin =
              volume ? position : offset_origin(position, error, dot(g, wi) < 0 ? -g : g);
          VolumeState direct_state = (interface && dot(g, wi) < 0) ||
                                             (features == 2 && m.dielectric && pbr_crossed(pbr, wi))
                                         ? adjacent
                                         : state;
          // Connect the two offset endpoints. Merely subtracting the origin's
          // shift from t_max leaves a parallel ray aimed beside the sampled
          // light point, which can hit the light itself. Large ground spheres
          // expose this because their conservative position errors are larger.
          // BSDF/PDF evaluation above retains the original, unshifted segment.
          float3 shadow_direction = wi;
          float max_t = maximum;
          if (isfinite(maximum)) {
            float3 target =
                offset_origin(light_target, light_error,
                              dot(light_normal, -wi) < 0 ? -light_normal : light_normal);
            float3 segment = target - origin;
            float distance = length(segment);
            shadow_direction = segment / distance;
            max_t = distance * (1 - 8 * FLT_EPSILON);
          }
          shadow_queue[slot] = {float4(origin, volume ? 0.f : as_type<float>(h.primitive + 1)),
                                float4(shadow_direction, max_t),
                                float4(contribution, 0),
                                pixel,
                                {volume ? 0 : h.padding[0], choice, 0},
                                direct_state};
          atomic_fetch_max_explicit(&indirect[6], (slot + 64) / 64, memory_order_relaxed);
        }
      }
    }
    float3 wi;
    float pdf;
    bool delta_sample = false;
    float specular_eta = 1;
    float branch =
        (interface || extended || features == 2) ? sample1(pixel, dim, p, directions, blue) : 0;
    float u = sample1(pixel, dim, p, directions, blue),
          v = sample1(pixel, dim, p, directions, blue);
    if (features == 4) {
      wi = diffusion_direction(exit_lobe, n, float2(u, v), diffusion_samples);
      BsdfValue value = diffusion_eval(exit_lobe, n, wi, diffusion_samples);
      pdf = value.pdf;
      if (!(pdf > 0))
        return;
      path.beta.xyz *= value.value / pdf;
      path.beta.w /= exit_lobe.shape.x * exit_lobe.shape.x;
    } else if (volume) {
      wi = phase_sample(-wo, anisotropy, u, v);
      pdf = phase_pdf(-wo, wi, anisotropy);
      if (guided_scatter) {
        wi = subsurface_direction(boundaries[body].color.xyz,diffusion_entries[pixel].normal.xyz,
                                  sample1(pixel,dim,p,directions,blue),u,v,wi);
        pdf = subsurface_direction_pdf(boundaries[body].color.xyz,
                                       diffusion_entries[pixel].normal.xyz,wi);
        path.beta.xyz *= (1/(4*pi))/pdf;
        path.flags |= 32;
      } else path.flags &= ~32u;
    } else if (interface) {
      float weight;
      if (!boundary_sample(wo, n, eta, roughness, branch, u, v, wi, weight, pdf))
        return;
      bool transmission = dot(wi, g) < 0;
      if ((dot(wi, n) < 0) != transmission)
        return;
      path.beta.xyz *= weight;
      if (entering && !boundaries[t.boundary].subsurface)
        path.beta.xyz *= rho;
      if (transmission) {
        state = adjacent;
        path.beta.w *= eta * eta;
      }
    } else if (features == 2) {
      BsdfSample sample = pbr_sample(pbr, float3(branch, u, v));
      if (!(sample.pdf > 0) || !all(isfinite(sample.weight)))
        return;
      wi = sample.direction;
      pdf = sample.pdf;
      delta_sample = sample.delta;
      specular_eta = sqrt(sample.eta_squared);
      path.beta.xyz *= sample.weight;
      path.beta.w *= sample.eta_squared;
      if (m.dielectric && pbr_crossed(pbr, wi))
        state = adjacent;
      if (!(path.flags & 2) && !sample.delta) {
        pixels[pixel].normal.xyz = pbr.normal;
        pixels[pixel].albedo = float4(rho, 0);
        path.flags |= 2;
      }
    } else if (extended) {
      BsdfSample sample =
          surface_sample(bsdf, float4(branch, u, v, sample1(pixel, dim, p, directions, blue)));
      if (!(sample.pdf > 0) || !all(isfinite(sample.weight)))
        return;
      wi = sample.direction;
      pdf = sample.pdf;
      path.beta.xyz *= sample.weight;
      path.beta.w *= sample.eta_squared;
      delta_sample = sample.delta;
      specular_eta = sqrt(sample.eta_squared);
    } else {
      if (texture_features & 4) {
        float facet = sample1(pixel, dim, p, directions, blue),
              escape = sample1(pixel, dim, p, directions, blue);
        wi = mapped_sample(mapped, float2(u, v), facet, escape);
        BsdfValue value = mapped_eval(mapped, wi, rho, m);
        pdf = value.pdf;
        if (dot(g, wi) <= 0 || pdf <= 0)
          return;
        path.beta.xyz *= value.value / pdf;
      } else {
        wi = basis_direction(n, float2(u, v));
        pdf = dot(wi, n) / pi;
        if (dot(g, wi) <= 0 || pdf <= 0)
          return;
        path.beta.xyz *= brdf(wo, wi, n, rho, m) * pi;
      }
    }
    if (!all(isfinite(path.beta))) {
      atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
      return;
    }
    if (all(path.beta.xyz == 0))
      return;
    path.previous = float4(position, as_type<float>(h.padding[0]));
    path.origin =
        float4(volume ? position : offset_origin(position, error, dot(g, wi) < 0 ? -g : g),
               volume ? 0.f : as_type<float>(h.primitive + 1));
    if (texture_features & 4) {
      Differential next = {};
      if (delta_sample && !volume) {
        bool transmitted = dot(g, wo) * dot(g, wi) < 0;
        next = propagate_differential(differentials[pixel], derivatives, n, path.direction.xyz,
                                      path.origin.xyz, wi, transmitted, specular_eta);
      }
      differentials[pixel] = next;
    }
    path.direction = float4(wi, pdf);
    path.dimension = dim;
    path.flags &= ~1u;
    if (delta_sample)
      path.flags |= 1;
    else if (texture_features & 64)
      path.flags |= 16; // Persistent: later mirrors must not hide diffuse-path emission.
    if (!volume)
      remember_subsurface_entry(path, pixel, state, body, t.boundary, position, g, boundaries, materials,
                               diffusion_entries);
    if (roulette(path, pixel, p, directions, blue, internal))
      enqueue(pixel, path, state, p, paths, membership, next_queue, counts);
  }

  template [[host_name("shade")]] kernel void shade<0>(SCENE_ARGUMENTS);
  template [[host_name("shade_materials")]] kernel void shade<1>(SCENE_ARGUMENTS);
  template [[host_name("shade_pbr")]] kernel void shade<2>(SCENE_ARGUMENTS);
  template [[host_name("shade_pbr_lights")]] kernel void shade<3>(SCENE_ARGUMENTS);
  template [[host_name("shade_diffusion")]] kernel void shade<4>(SCENE_ARGUMENTS);

  kernel void shadows(SCENE_ARGUMENTS) {
    // Without finite haze, visibility and RGB transmittance need one thread per
    // connection. Keep their result in shadow_queue for the subsequent spectral
    // pass. Finite haze retains cooperative transport: even blocked connections
    // may contribute atmospheric in-scattering along the segment.
    bool gpu_air = scene_has_atmosphere && (p.camera.padding[1] & 2) &&
                   !defer_shadow_lighting;
    uint lane=gpu_air ? i%32 : 0;
    if(gpu_air) i/=32;
    bool gpu_haze = gpu_air && (p.camera.padding[1] & 1);
    if (scene_has_atmosphere && (p.camera.padding[1] & 1) && !gpu_air) return;
    device const uchar *bytes = reinterpret_cast<device const uchar *>(storage);
    device const Quadric *quadrics =
        reinterpret_cast<device const Quadric *>(bytes + storage->quadric_offset);
    device const MeshInstance *instances =
        reinterpret_cast<device const MeshInstance *>(bytes + storage->instance_offset);
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions,
                           instances,
                           reinterpret_cast<device const float4 *>(bytes + storage->motion_offset),
                           reinterpret_cast<device const Quadric *>(bytes + storage->motion_offset + 48)};
    if (i >= (primary_coverage ? p.active : atomic_load_explicit(&counts.shadows, memory_order_relaxed)))
      return;
    Shadow s;
    if (primary_coverage) {
      // Reuse the exact generated camera ray, including lens and shutter samples.
      // This independent straight-line estimator writes transparency only; it
      // does not follow the refracted/scattered RGB path or consume its samples.
      uint pixel = active[i];
      Path path = paths[pixel];
      if ((path.flags & 4) || pixels[pixel].normal.w <= 0) return;
      s = {};
      s.pixel = pixel;
      s.origin = path.origin;
      s.direction = float4(path.direction.xyz, INFINITY);
      s.contribution = float4(1);
      s.padding[1] = UINT_MAX;
      if (scene_has_boundaries) s.state = membership[pixel];
    } else
      s = shadow_queue[i];
    device const AtmosphereData &air = *reinterpret_cast<device const AtmosphereData *>(bytes + storage->atmosphere_offset);
    device const float *air_data = reinterpret_cast<device const float *>(bytes + storage->atmosphere_data_offset);
    AirRay air_ray = {};
    bool infinite = s.padding[1] < p.environments;
    bool sky_source = infinite && (lights[s.padding[1]].flags & 128);
    if (scene_has_atmosphere && !primary_coverage)
      if(lane==0) shadow_queue[i].contribution.xyz = 0;
    if (s.direction.w <= 0)
      return;
    ray r;
    r.origin = s.origin.xyz;
    r.direction = s.direction.xyz;
    r.min_distance = FLT_MIN;
    r.max_distance = s.direction.w;
    if (!scene_has_boundaries && !(texture_features & (8 | 64)) && !gpu_haze) {
      if (trace_spectral_shadow(r, scene, quadrics, gpu_air, lane, true).type == intersection_type::none) {
        if (primary_coverage) {
          if(lane==0) pixels[s.pixel].radiance.w = 1;
        } else if (gpu_air) {
          float3 value = infinite ? air_light(lights[s.padding[1]],r.origin,r.direction,air,air_data,textures,texels) : float3(1);
          if(lane==0) pixels[s.pixel].radiance.xyz += s.contribution.xyz * value;
        } else if (scene_has_atmosphere)
          shadow_queue[i].contribution = s.contribution;
        else
          if(lane==0) pixels[s.pixel].radiance.xyz += s.contribution.xyz;
      }
      return;
    }
    float previous = 0;
    float3 transmittance = 1;
    VolumeState state = s.state;
    device const Medium *media = reinterpret_cast<device const Medium *>(bytes + storage->medium_offset);
    device const float *fields = reinterpret_cast<device const float *>(bytes + storage->field_offset);
    uint shadow_dimension = primary_coverage ? 0x713ad591u : 0;
    while (true) {
      auto hit = trace_spectral_shadow(r, scene, quadrics, gpu_air, lane);
      uint departed = as_type<uint>(s.origin.w);
      if (hit.type != intersection_type::none && departed && !is_quadric(departed - 1) &&
          hit.primitive_id == departed - 1 && hit.instance == s.padding[0]) {
        float next = as_type<float>(as_type<uint>(hit.distance) + 1u);
        if (!(next > r.min_distance)) {
          if(lane==0) atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
          return;
        }
        r.min_distance = next;
        continue;
      }
      float endpoint = hit.type == intersection_type::none ? r.max_distance : hit.distance;
      uint before = active_boundary(state, boundaries, materials);
      uint medium_body = active_medium(state, boundaries, materials);
      // Diffusion integrates the interior nonlocally; it has no ballistic segment.
      if (before && boundaries[before].subsurface == 2)
        return;
      float3 sigma = medium_body ? extinction(medium_body, boundaries, materials) : float3(0);
      float distance = endpoint - previous;
      Medium medium = {};
      if ((texture_features & 16) && medium_body)
        medium = media[boundaries[medium_body].medium];
      bool air_enabled = gpu_haze && !before && (!medium_body || (medium.emission.w!=0 && air.haze_in_volumes));
      float air_cursor = previous;
      if (gpu_air && !isfinite(endpoint)) air_start(air_ray,r.origin+previous*r.direction,r.direction,air_enabled);
      if (medium.controls.w != 0 || (air_enabled && medium_body)) {
        // Ratio tracking estimates RGB transmittance with a common Poisson
        // process. Its independent stream never changes the continuation path.
        float majorant = max(sigma.x, max(sigma.y, sigma.z)) * (medium.controls.w!=0 ? medium.hi.w : 1.f);
        float t = previous;
        MediumFlightIterator iterator = {};
        if (medium.controls.w != 0)
          iterator = medium_iterator(medium, r.origin + previous * r.direction, r.direction,
                                     max(sigma.x, max(sigma.y, sigma.z)));
        while (majorant > 0) {
          uint bits = hash_seed(hash_seed(s.pixel, p.sample ^ p.seed),
                                shadow_dimension++ ^ (p.depth * 0x9e3779b9u));
          float u = (float(bits >> 8) + .5f) * 0x1p-24f;
          float next = medium.controls.w != 0
              ? previous + medium_next(iterator, min(u, 0x1.fffffep-1f), fields, majorant)
              : t - log(1-u)/majorant;
          if (next >= endpoint)
            break;
          // A valid small flight may round to the current ray parameter. It
          // still contributes a ratio-tracking event there. Reject only when
          // even one mean free path is unresolvable, as in camera transport.
          if (!(next > t) && !(t + 1 / majorant > t)) {
            if(lane==0) atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
            return;
          }
          t = next;
          if (gpu_haze) {
            auto segment=air_walk(air_ray,r.origin+air_cursor*r.direction,r.direction,t-air_cursor,
                                  air_enabled,medium,air_uniform(s.pixel^shadow_dimension,p),air,air_data,fields);
            if(sky_source) if(lane==0) pixels[s.pixel].radiance.xyz+=s.contribution.xyz*transmittance*segment.radiance;
            transmittance*=segment.transmission;
            air_cursor=t;
          }
          float density = medium.controls.w != 0 ? medium_density(medium, medium_coordinates(medium, r.origin + t*r.direction), fields) : 1.f;
          transmittance *= max(float3(0), 1 - sigma * density / majorant);
          if (all(transmittance == 0))
            return;
        }
      } else
        transmittance *= isfinite(distance) ? exp(-sigma * distance)
                                            : select(float3(0), float3(1), sigma == 0);
      if(gpu_haze && isfinite(endpoint) && endpoint>air_cursor) {
        auto segment=air_walk(air_ray,r.origin+air_cursor*r.direction,r.direction,endpoint-air_cursor,
                              air_enabled,medium,air_uniform(s.pixel^shadow_dimension,p),air,air_data,fields);
        if(sky_source) if(lane==0) pixels[s.pixel].radiance.xyz+=s.contribution.xyz*transmittance*segment.radiance;
        transmittance*=segment.transmission;
      }
      if (before) {
        float3 absorption = materials[boundaries[before].material].attenuation.xyz;
        transmittance *= isfinite(distance) ? exp(-absorption * distance)
                                            : select(float3(0), float3(1), absorption == 0);
      }
      if (all(transmittance == 0))
        return;
      if (hit.type == intersection_type::none) {
        if (primary_coverage) {
          if (gpu_haze && air_ray.direction.w)
            transmittance *= air_remaining(air_ray,
                air_infinite_transmission(air_ray.origin.xyz, r.direction, air, air_data), false);
          if(lane==0) pixels[s.pixel].radiance.w = clamp(dot(transmittance, float3(1.f/3)), 0.f, 1.f);
        } else if (gpu_air) {
          float3 anchor=gpu_haze ? (air_ray.direction.w ? air_ray.origin.xyz : r.origin+previous*r.direction)
                                : s.origin.xyz;
          float3 value=infinite ? air_light(lights[s.padding[1]],anchor,r.direction,air,air_data,textures,texels) : float3(1);
          if(infinite) value=air_remaining(air_ray,value,sky_source);
          if(lane==0) pixels[s.pixel].radiance.xyz+=transmittance*s.contribution.xyz*value;
        } else if (scene_has_atmosphere)
          shadow_queue[i].contribution.xyz = transmittance * s.contribution.xyz;
        else
          if(lane==0) pixels[s.pixel].radiance.xyz += transmittance * s.contribution.xyz;
        return;
      }
      GeometrySurface contact = geometry_surface(hit.primitive_id, hit.coordinates, triangles,
                                                 quadrics, surfaces, hit.instance, scene.instances, r.origin, r.direction);
      Triangle t = contact.face;
      Material surface_material = materials[t.material];
      if ((texture_features & 64) && surface_material.emissive == 1 && surface_material.extra.w != 0) {
        previous = endpoint;
        r.min_distance = as_type<float>(as_type<uint>(hit.distance) + 1u);
        continue;
      }
      if (texture_features & 8) {
        uint alpha = surfaces[as_type<uint>(t.uv2.w)].alpha;
        if (alpha != UINT_MAX) {
          float2 uv = hit.coordinates.xy;
          float opacity =
              clamp(texture_value(alpha,
                                  surface_context(t, float3(1 - uv.x - uv.y, uv),
                                                  r.origin + endpoint * r.direction, surfaces),
                                  textures, texels)
                        .x,
                    0.f, 1.f);
          // Opaque coverage has an exact scalar shadow transmittance. At a
          // dielectric boundary, choose the covered/hole branch instead: they
          // can enter different media and cannot share one averaged stack.
          bool hole = true;
          if (t.boundary) {
            uint bits = hash_seed(hash_seed(s.pixel, p.sample),
                                  hit.primitive_id ^ (hit.instance * 0x85ebca6bu) ^ p.seed ^
                                      (p.depth * 0x9e3779b9u));
            hole = float(bits) * 0x1p-32f >= opacity;
          } else
            transmittance *= 1 - opacity;
          if (all(transmittance == 0))
            return;
          if (hole) {
            float next = as_type<float>(as_type<uint>(hit.distance) + 1u);
            if (!(next > r.min_distance))
              return;
            previous = endpoint;
            r.min_distance = next;
            continue;
          }
        }
      }
      if (!t.boundary)
        return;
      if (boundaries[t.boundary].subsurface == 3 && boundaries[t.boundary].radius.w != 0 &&
          !materials[t.material].dielectric)
        return;
      bool entering = dot(r.direction, t.geometric.xyz) < 0;
      float3 position = contact.position, error = contact.error;
      if (!adjacent_membership(t.boundary, entering, position, error, t.geometric.xyz, state, p,
                               triangles, quadrics, surfaces, boundaries, scene)) {
        if(lane==0) atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
        return;
      }
      uint after = active_boundary(state, boundaries, materials);
      bool hidden = before == after ||
                    (before && after && boundaries[before].material == boundaries[after].material);
      if (!hidden) {
        // Camera coverage sees every exposed finite surface, including glass
        // with matched IOR. Boundary-only media are not visible surfaces.
        if (primary_coverage && !(boundaries[t.boundary].subsurface == 3 &&
                                  boundaries[t.boundary].radius.w == 0))
          return;
        if (boundary_ior(before, boundaries, materials) !=
            boundary_ior(after, boundaries, materials))
          return;
        if ((before && materials[boundaries[before].material].type == bsdf_openpbr) ||
            (after && materials[boundaries[after].material].type == bsdf_openpbr))
          return; // Matching IOR does not remove OpenPBR's other scattering lobes.
        if (entering && !boundaries[t.boundary].subsurface) {
          Material m = materials[t.material];
          float2 uv = hit.coordinates.xy;
          transmittance *=
              texture_value(m.texture, triangle_uv(t, float3(1 - uv.x - uv.y, uv.x, uv.y)),
                            r.origin + endpoint * r.direction, textures, texels);
        }
      }
      // Stay on the exact light segment; changing its direction would invalidate
      // the direct-light PDF. Refractive interfaces above block this estimator.
      float next = as_type<float>(as_type<uint>(hit.distance) + 1u);
      if (!(next > r.min_distance) || !isfinite(next)) {
        if(lane==0) atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
        return;
      }
      previous = endpoint;
      r.min_distance = next;
    }
  }

  // Producers publish indirect group counts with atomic max. One end-of-event
  // transition replaces four separate counting dispatches; command-encoder order
  // makes all queue records visible before the consuming indirect dispatch.
  kernel void advance_queues(device Counts &c [[buffer(0)]], device uint *args [[buffer(1)]]) {
    uint count = atomic_load_explicit(&c.next, memory_order_relaxed);
    atomic_store_explicit(&c.rays, count, memory_order_relaxed);
    atomic_store_explicit(&c.hits, 0, memory_order_relaxed);
    atomic_store_explicit(&c.shadows, 0, memory_order_relaxed);
    atomic_store_explicit(&c.next, 0, memory_order_relaxed);
    args[0] = max(1u, (count + 63) / 64);
    args[1] = 1;
    args[2] = 1;
    for (uint j = 3; j < 9; ++j)
      args[j] = 1;
  }

)MSL"
