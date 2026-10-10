R"MSL(

  // Specialize away all interior state/traversal for scenes without boundaries.
  // This is a pipeline constant, not a per-ray material branch.
  constant bool scene_has_boundaries [[function_constant(0)]];
  constant uint texture_features [[function_constant(1)]];
  constant bool scene_has_quadrics [[function_constant(2)]];
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
    uint material, subsurface, diffusion_offset, padding;
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
  };
  struct Camera {
    float4 origin, lower_left, horizontal, vertical, right, up, forward;
    uint orthographic, flip_x, padding[2];
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
    atomic_uint rays, hits, shadows, next, error, discarded, scattering;
  };
  static_assert(sizeof(Triangle) == 160 && sizeof(Texture) == 80 && sizeof(Material) == 144 &&
                    sizeof(Boundary) == 112 && sizeof(VolumeState) == 80 && sizeof(Light) == 80 &&
                    sizeof(Camera) == 128 && sizeof(Params) == 208 && sizeof(Pixel) == 48 &&
                    sizeof(Path) == 80 && sizeof(Hit) == 32 && sizeof(Shadow) == 144 &&
                    sizeof(Counts) == 28,
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
  float sample1(uint pixel, thread uint &dimension, constant Params &p,
                device const uint *directions, device const uint *blue) {
    uint dim = dimension++, x = pixel / p.height, y = pixel % p.height;
    if (p.sampler == 0)
      return min(float(hash_seed(hash_seed(p.sample, pixel ^ p.seed), dim)) * 0x1p-32f,
                 0x1.fffffep-1f);
    if (p.sampler == 3 && dim < 8) {
      uint tile = ((x & 127) + (y & 127) * 128) * 8 + dim;
      uint ranked = (p.sample & 255) ^ blue[65536 + tile];
      return (0.5f + float(blue[ranked * 256 + dim] ^ blue[65536 + 131072 + tile])) / 256.f;
    }
    uint seed = p.sampler == 3 ? x * 0x9e3779b9u ^ y * 0x85ebca6bu : hash_seed(pixel, p.seed);
    seed += dim / 21200;
    dim %= 21200;
    uint index = owen(p.sample, seed), value = 0;
    for (uint bit = 0; index; index >>= 1, ++bit)
      if (index & 1)
        value ^= directions[dim * 32 + bit];
    return min(float(owen(value, seed ^ (dim + (seed << 6) + (seed >> 2)))) * 0x1p-32f,
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
R"MSL(

  float3 environment(Light l, float3 wi, float3 position, device const Texture *textures,
                     device const float4 *texels) {
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
  float light_pdf(Light l) { return l.type == 0 ? 1 / (4 * pi) : 1 / (2 * pi * l.forward.w); }

  uint active_boundary(thread const VolumeState &state, device const Boundary *boundaries,
                       device const Material *materials) {
    if (!scene_has_boundaries)
      return 0;
    uint winner = 0;
    for (uint j = 0; j < state.count; ++j) {
      uint id = state.ids[j];
      if (!winner || materials[boundaries[id].material].priority <=
                         materials[boundaries[winner].material].priority)
        winner = id;
    }
    return winner;
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
            dot(direction,
                geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces)
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
      GeometrySurface contact =
          geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces);
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
#include "wavefront_openpbr.inc"
#include "wavefront_bump.inc"
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
      device const Quadric *quadrics [[buffer(28)]],                                               \
      instance_acceleration_structure mixed_scene [[buffer(29)]],                                  \
      intersection_function_table<triangle_data, instancing> mixed_functions [[buffer(30)]],       \
      uint i [[thread_position_in_grid]]

  kernel void primary(SCENE_ARGUMENTS) {
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions};
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
    float u = (pixel / p.height + sample1(pixel, dim, p, directions, blue)) / p.width;
    float v = (pixel % p.height + sample1(pixel, dim, p, directions, blue)) / p.height;
    if (p.camera.flip_x)
      u = 1 - u;
    float lu = sample1(pixel, dim, p, directions, blue),
          lv = sample1(pixel, dim, p, directions, blue);
    // Match the native thin-lens camera's uniform disk; unused time retains its dimension.
    sample1(pixel, dim, p, directions, blue);
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
    ray_queue[i] = pixel;
  }

  kernel void intersect_paths(SCENE_ARGUMENTS) {
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions};
    if (i >= atomic_load_explicit(&counts.rays, memory_order_relaxed))
      return;
    uint pixel = ray_queue[i];
    Path path = paths[pixel];
    if (path.flags & 4) {
      atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
      return;
    }
    VolumeState state = {};
    state.count = 0;
    if (scene_has_boundaries)
      state = membership[pixel];
    uint body = active_boundary(state, boundaries, materials);
    bool scattering = body && boundaries[body].subsurface;
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
    float flight = INFINITY;
    if (scattering) {
      // A fresh RGB-channel mixture at every segment, with the exact mixture
      // density below. Zero-rate channels retain their mass at infinity.
      uint channel = min(2u, uint(3 * sample1(pixel, path.dimension, p, directions, blue)));
      float u = sample1(pixel, path.dimension, p, directions, blue);
      if (rates[channel] > 0)
        flight = -log(1 - u) / rates[channel];
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
      bool skip = departed && !is_quadric(departed - 1) && hit.primitive_id == departed - 1;
      if (!skip && (texture_features & 8)) {
        GeometrySurface contact =
            geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces);
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
    if (scattering) {
      // Remove a common exponential before normalizing. This is the same ratio
      // without underflow in an optically thick RGB segment.
      float low = min(rates.x, min(rates.y, rates.z));
      float3 relative = isfinite(distance) ? exp(-(rates - low) * distance)
                                           : select(float3(0), float3(1), rates == 0);
      float density = dot(relative * (collision ? rates : float3(1)), float3(1.f / 3));
      if (!(density > 0))
        return;
      path.beta.xyz *= relative * (collision ? boundaries[body].sigma_s.xyz : float3(1)) / density;
    }
    if (body) {
      float3 absorption = materials[boundaries[body].material].attenuation.xyz;
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
            geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces);
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
      ++path.internal;
      atomic_fetch_add_explicit(&counts.scattering, 1, memory_order_relaxed);
    } else if (hit.type == intersection_type::none) {
      float3 radiance = 0;
      for (uint j = 0; j < p.environments; ++j) {
        Light l = lights[j];
        if (l.type > 1)
          continue;
        float weight = (path.flags & 1) ? 1 : mis(path.direction.w, light_pdf(l) / p.lights);
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
                       collision ? UINT_MAX : hit.primitive_id,
                       collision ? 0 : hit.coordinates.x,
                       collision ? 0 : hit.coordinates.y,
                       distance,
                       collision ? 0 : hit.coordinates.z,
                       {0, 0}};
    atomic_fetch_max_explicit(&indirect[3], (slot + 64) / 64, memory_order_relaxed);
  }

  template <ushort features> kernel void shade(SCENE_ARGUMENTS) {
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions};
    if (i >= atomic_load_explicit(&counts.hits, memory_order_relaxed))
      return;
    Hit h = hit_queue[i];
    bool is_diffusion = h.primitive == UINT_MAX - 1;
    if (is_diffusion != (features == 4))
      return;
    bool is_pbr =
        h.primitive < UINT_MAX - 1 &&
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
    VolumeState adjacent = state;
    bool volume = scene_has_boundaries && h.primitive == UINT_MAX, interface = false,
         internal = volume;
    uint body = active_boundary(state, boundaries, materials);
    float3 position = path.origin.xyz, error = 0, g = 0, n = 0, rho = 0, wo = -path.direction.xyz;
    float eta = 1, roughness = 0, anisotropy = body ? boundaries[body].sigma_s.w : 0;
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
        pixels[pixel].normal = 0;
        pixels[pixel].albedo = float4(boundaries[body].color.xyz, 0);
        path.flags |= 2;
      }
      path.beta.xyz *= exit.weight;
      h = exit.hit;
      geometry =
          geometry_surface(h.primitive, float3(h.u, h.v, h.w), triangles, quadrics, surfaces);
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
      if (!(path.flags & 2)) {
        float3 st = boundaries[body].sigma_a.xyz + boundaries[body].sigma_s.xyz;
        pixels[pixel].normal = 0;
        pixels[pixel].albedo = float4(
            select(float3(0), boundaries[body].sigma_s.xyz / max(st, float3(FLT_MIN)), st > 0), 0);
        path.flags |= 2;
      }
    } else {
      geometry =
          geometry_surface(h.primitive, float3(h.u, h.v, h.w), triangles, quadrics, surfaces);
      t = geometry.face;
      m = materials[t.material];
      float3 bary = float3(1 - h.u - h.v, h.u, h.v);
      position = geometry.position;
      error = geometry.error;
      g = t.geometric.xyz;
      n = geometry.normal;
      TextureContext texture_hit = surface_context(t, bary, position, surfaces);
      base_normal = n;
      if (texture_features & 4) {
        Surface mapping = surfaces[as_type<uint>(t.uv2.w)];
        derivatives = surface_derivatives(t, mapping, texture_hit, n, differentials[pixel]);
        n = bump_normal(t, mapping, texture_hit, n, derivatives, textures, texels);
        mapped_surface = mapping.bump != UINT_MAX;
      }
      rho = texture_value(m.texture, texture_hit, textures, texels);
      float mapped_roughness = -1;
      if (m.second_texture != UINT_MAX && m.type != bsdf_translucent) {
        float3 rough = texture_value(m.second_texture, texture_hit, textures, texels);
        if (m.type == bsdf_openpbr)
          mapped_roughness = rough.x;
        else
          m.parameters.xy = rough.xy;
      }
      pixels[pixel].radiance.w = 1;
      if (m.emissive == 1) {
        if (dot(n, wo) > 0) {
          float pdf =
              geometry_light_pdf(h.primitive, geometry, path.previous.xyz, quadrics) / p.lights;
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
            remember_diffusion_entry(path, pixel, adjacent, position, g, boundaries, materials,
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
              geometry_light_pdf(h.primitive, geometry, path.previous.xyz, quadrics) / p.lights;
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
          path.dimension = dim;
          remember_diffusion_entry(path, pixel, adjacent, position, g, boundaries, materials,
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
          path.previous = float4(position, 0);
          if (texture_features & 4)
            differentials[pixel] =
                propagate_differential(differentials[pixel], derivatives, n, path.direction.xyz,
                                       path.origin.xyz, wi, transmission, eta);
          path.direction = float4(wi, pdf);
          remember_diffusion_entry(path, pixel, state, position, g, boundaries, materials,
                                   diffusion_entries);
          path.flags |= 1;
          path.dimension = dim;
          if (roulette(path, pixel, p, directions, blue, internal))
            enqueue(pixel, path, state, p, paths, membership, next_queue, counts);
          return;
        }
        interface = true;
        if (!(path.flags & 2)) {
          pixels[pixel].normal = float4(n, 0);
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
            if (t.flags & 6) {
              float3 view_plane = wo - tangent * dot(tangent, wo);
              float projected_length = length(view_plane);
              if (projected_length == 0)
                return;
              view_plane /= projected_length;
              float impact =
                  (t.flags & 2)
                      ? clamp(dot(normal, cross(tangent, view_plane)), -1.f, 1.f)
                      : (1 - 2 * triangle_uv(t, bary).y) * (dot(normal, wo) >= 0 ? 1.f : -1.f);
              float width = mix(t.geometric.w, t.uv2.z, triangle_uv(t, bary).x);
              error = width * float3(t.p[0].w, t.p[1].w, t.p[2].w);
              if (t.flags & 2) {
                // PBRT's cylindrical curve rotates its frame by 90*h degrees.
                // Recover h from the projected cross-section, not its static UV v.
                normal = cos(.5f * pi * impact) * view_plane +
                         sin(.5f * pi * impact) * cross(tangent, view_plane);
                float radius = dot(bary, float3(t.n[0].w, t.n[1].w, t.n[2].w));
                position -= wo * (radius * sqrt(max(0.f, 1 - impact * impact)) / projected_length);
              } else if (dot(normal, wo) < 0)
                normal = -normal;
              g = normal;
              bsdf.impact = impact;
            }
            bsdf.frame = {tangent, cross(normal, tangent), normal};
          }
          bsdf.view = bsdf.frame.local(wo);
          bsdf.rho = rho;
          bsdf.distance = h.distance;
          if (!(m.type == bsdf_hair && (t.flags & 6)))
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
          pixels[pixel].normal = float4(n, 0);
          pixels[pixel].albedo = float4(rho, 0);
          path.flags |= 2;
        }
      }
    }
    // One light proposal per vertex. Each emitter has its own competing BSDF
    // technique; evaluating emission individually on a miss avoids mixture bias.
    if (p.lights) {
      uint choice = min(uint(sample1(pixel, dim, p, directions, blue) * p.lights), p.lights - 1);
      float u = sample1(pixel, dim, p, directions, blue),
            v = sample1(pixel, dim, p, directions, blue);
      Light l = lights[choice];
      float3 wi, emitted = 0;
      float3 light_target = 0, light_error = 0, light_normal = 0;
      float pdf = 0, maximum = INFINITY;
      if (l.type == 0) {
        float y = 1 - 2 * u, r = sqrt(max(0.f, 1 - y * y)), phi = 2 * pi * v;
        wi = float3(r * cos(phi), y, r * sin(phi));
        emitted = environment(l, wi, position, textures, texels);
        pdf = light_pdf(l) / p.lights;
      } else if (l.type == 1) {
        float d = u * l.forward.w, s = sqrt(d * (2 - d)), phi = 2 * pi * v;
        wi = (1 - d) * l.forward.xyz + s * (cos(phi) * l.right.xyz + sin(phi) * l.up.xyz);
        emitted = environment(l, wi, position, textures, texels);
        pdf = light_pdf(l) / p.lights;
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
        pdf = 1.f / p.lights;
      } else {
        GeometrySurface source = sample_geometry_light(l.primitive, float2(u, v), position,
                                                       triangles, quadrics, surfaces);
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
          pdf = geometry_light_pdf(l.primitive, source, position, quadrics) / p.lights;
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
                                {0, 0, 0},
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
        pixels[pixel].normal = float4(pbr.normal, 0);
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
    path.previous = float4(position, 0);
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
    if (!volume)
      remember_diffusion_entry(path, pixel, state, position, g, boundaries, materials,
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
    GeometryScene scene = {triangle_scene, triangle_functions, mixed_scene, mixed_functions};
    if (i >= atomic_load_explicit(&counts.shadows, memory_order_relaxed))
      return;
    Shadow s = shadow_queue[i];
    if (s.direction.w <= 0)
      return;
    ray r;
    r.origin = s.origin.xyz;
    r.direction = s.direction.xyz;
    r.min_distance = FLT_MIN;
    r.max_distance = s.direction.w;
    if (!scene_has_boundaries && !(texture_features & 8)) {
      if (trace_geometry(r, scene, quadrics, 0, true).type == intersection_type::none)
        pixels[s.pixel].radiance.xyz += s.contribution.xyz;
      return;
    }
    float previous = 0;
    float3 transmittance = 1;
    VolumeState state = s.state;
    while (true) {
      auto hit = trace_geometry(r, scene, quadrics);
      uint departed = as_type<uint>(s.origin.w);
      if (hit.type != intersection_type::none && departed && !is_quadric(departed - 1) &&
          hit.primitive_id == departed - 1) {
        float next = as_type<float>(as_type<uint>(hit.distance) + 1u);
        if (!(next > r.min_distance)) {
          atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
          return;
        }
        r.min_distance = next;
        continue;
      }
      float endpoint = hit.type == intersection_type::none ? r.max_distance : hit.distance;
      uint before = active_boundary(state, boundaries, materials);
      // Diffusion integrates the interior nonlocally; it has no ballistic segment.
      if (before && boundaries[before].subsurface == 2)
        return;
      float3 sigma = extinction(before, boundaries, materials);
      float distance = endpoint - previous;
      transmittance *=
          isfinite(distance) ? exp(-sigma * distance) : select(float3(0), float3(1), sigma == 0);
      if (all(transmittance == 0))
        return;
      if (hit.type == intersection_type::none) {
        pixels[s.pixel].radiance.xyz += transmittance * s.contribution.xyz;
        return;
      }
      GeometrySurface contact =
          geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces);
      Triangle t = contact.face;
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
                                  hit.primitive_id ^ p.seed ^ (p.depth * 0x9e3779b9u));
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
      bool entering = dot(r.direction, t.geometric.xyz) < 0;
      float3 position = contact.position, error = contact.error;
      if (!adjacent_membership(t.boundary, entering, position, error, t.geometric.xyz, state, p,
                               triangles, quadrics, surfaces, boundaries, scene)) {
        atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
        return;
      }
      uint after = active_boundary(state, boundaries, materials);
      bool hidden = before == after ||
                    (before && after && boundaries[before].material == boundaries[after].material);
      if (!hidden) {
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
        atomic_fetch_add_explicit(&counts.discarded, 1, memory_order_relaxed);
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
