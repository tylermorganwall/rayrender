R"MSL(
#include <metal_stdlib>
#include <metal_raytracing>
using namespace metal;
using namespace raytracing;
constant float pi = 3.14159265358979323846f;
struct Triangle {
  float4 p[3], n[3], uv01, uv2, geometric;
  uint material, light, padding[2];
};
struct Texture {
  float4 a, b, c, mapping;
  uint type, offset, width, height;
};
struct Material {
  float4 diffuse;
  uint texture, emissive, padding[2];
};
struct Light {
  float4 forward, right, up, color;
  uint type, texture, triangle, flags;
};
struct Camera {
  float4 origin, lower_left, horizontal, vertical, right, up, forward;
  uint orthographic, flip_x, padding[2];
};
struct Params {
  Camera camera;
  uint width, height, active, sample, lights, depth, roulette, sampler, seed, max_depth, padding[2];
};
struct Pixel {
  float4 radiance, normal, albedo;
};
struct Path {
  float4 origin, direction, beta, previous;
  uint dimension, padding[3];
};
struct Hit {
  uint pixel, triangle;
  float u, v;
};
struct Shadow {
  float4 origin, direction, contribution;
  uint pixel, padding[3];
};
// Counters: current ray queue, hit queue, shadow queue, next ray queue, error.
struct Counts {
  atomic_uint rays, hits, shadows, next, error;
};
static_assert(sizeof(Triangle) == 160 && sizeof(Texture) == 80 && sizeof(Material) == 32 &&
              sizeof(Light) == 80 && sizeof(Camera) == 128 && sizeof(Params) == 176 &&
              sizeof(Pixel) == 48 && sizeof(Path) == 80 && sizeof(Hit) == 16 &&
              sizeof(Shadow) == 64 && sizeof(Counts) == 20, "wavefront host/device ABI mismatch");

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
float sample1(uint pixel, thread uint &dimension, constant Params &p, device const uint *directions,
              device const uint *blue) {
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
float wrap(float v) {
  float w = v - floor(v);
  return w == 0 && v > 0 ? 1 : w;
}
float3 texture_value(uint id, float2 uv, float3 position, device const Texture *textures,
                     device const float4 *texels) {
  Texture t = textures[id];
  if (t.type == 0)
    return t.a.xyz;
  if (t.type == 1) {
    float3 s = sin(pi * position / t.a.w);
    s = select(s, float3(1), s == 0);
    return s.x * s.y * s.z < 0 ? t.b.xyz : t.a.xyz;
  }
  if (t.type == 2 || t.type == 3)
    return mix(t.a.xyz, t.b.xyz, t.type == 2 ? uv.x : uv.y);
  if (t.type == 4)
    return uv.x * t.a.xyz + uv.y * t.b.xyz + (1 - uv.x - uv.y) * t.c.xyz;
  float2 coord;
  int x0, x1, y0, y1;
  if (t.type == 6) {
    coord = float2(fract(uv.x) * t.width - .5f,
                   clamp((1 - uv.y) * t.height - .5f, 0.f, float(t.height - 1)));
    x0 = (int(floor(coord.x)) + int(t.width)) % t.width;
    x1 = (x0 + 1) % t.width;
  } else {
    uv = uv * t.mapping.xy + t.mapping.zw;
    coord = float2(wrap(uv.x) * (t.width - 1), (1 - wrap(uv.y)) * (t.height - 1));
    x0 = int(floor(coord.x));
    x1 = min(x0 + 1, int(t.width) - 1);
  }
  y0 = int(floor(coord.y));
  y1 = min(y0 + 1, int(t.height) - 1);
  float2 w = fract(coord);
  return mix(
      mix(texels[t.offset + y0 * t.width + x0].xyz, texels[t.offset + y0 * t.width + x1].xyz, w.x),
      mix(texels[t.offset + y1 * t.width + x0].xyz, texels[t.offset + y1 * t.width + x1].xyz, w.x),
      w.y);
}
float mis(float a, float b) {
  float ratio = b / a;
  return 1 / (1 + ratio * ratio);
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

#define SCENE_ARGUMENTS                                                                            \
  device const Triangle *triangles [[buffer(0)]], device const Material *materials [[buffer(1)]],  \
      device const Texture *textures [[buffer(2)]], device const float4 *texels [[buffer(3)]],     \
      device const Light *lights [[buffer(4)]], device Path *paths [[buffer(5)]],                  \
      device Pixel *pixels [[buffer(6)]], device uint *ray_queue [[buffer(7)]],                    \
      device uint *next_queue [[buffer(8)]], device Hit *hit_queue [[buffer(9)]],                  \
      device Shadow *shadow_queue [[buffer(10)]], device Counts &counts [[buffer(11)]],            \
      constant Params &p [[buffer(12)]], device const uint *active [[buffer(13)]],                 \
      device const uint *directions [[buffer(14)]], device const uint *blue [[buffer(15)]],        \
      primitive_acceleration_structure scene [[buffer(16)]], uint i [[thread_position_in_grid]]

kernel void primary(SCENE_ARGUMENTS) {
  if (i == 0) {
    atomic_store_explicit(&counts.rays, p.active, memory_order_relaxed);
    atomic_store_explicit(&counts.hits, 0, memory_order_relaxed);
    atomic_store_explicit(&counts.shadows, 0, memory_order_relaxed);
    atomic_store_explicit(&counts.next, 0, memory_order_relaxed);
    atomic_store_explicit(&counts.error, 0, memory_order_relaxed);
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
  path.beta = float4(1, 1, 1, 0);
  path.previous = path.origin;
  path.dimension = dim;
  paths[pixel] = path;
  pixels[pixel] = {float4(0), float4(0), float4(0)};
  ray_queue[i] = pixel;
}

kernel void intersect_paths(SCENE_ARGUMENTS) {
  if (i >= atomic_load_explicit(&counts.rays, memory_order_relaxed))
    return;
  uint pixel = ray_queue[i];
  Path path = paths[pixel];
  ray r;
  r.origin = path.origin.xyz;
  r.direction = path.direction.xyz;
  r.min_distance = FLT_MIN;
  r.max_distance = INFINITY;
  intersector<triangle_data> query;
  query.assume_geometry_type(geometry_type::triangle);
  query.force_opacity(forced_opacity::opaque);
  auto hit = query.intersect(r, scene);
  if (hit.type == intersection_type::none) {
    float3 radiance = 0;
    for (uint j = 0; j < p.lights; ++j) {
      Light l = lights[j];
      if (l.type > 1)
        continue;
      float3 emitted = environment(l, r.direction, path.previous.xyz, textures, texels);
      float weight = p.depth == 0 ? 1 : mis(path.direction.w, light_pdf(l) / p.lights);
      radiance += weight * emitted;
    }
    pixels[pixel].radiance.xyz += path.beta.xyz * radiance;
    return;
  }
  uint slot = atomic_fetch_add_explicit(&counts.hits, 1, memory_order_relaxed);
  if (slot >= p.active) {
    atomic_store_explicit(&counts.error, 1, memory_order_relaxed);
    return;
  }
  hit_queue[slot] = {pixel, hit.primitive_id, hit.triangle_barycentric_coord.x,
                     hit.triangle_barycentric_coord.y};
}

kernel void shade(SCENE_ARGUMENTS) {
  if (i >= atomic_load_explicit(&counts.hits, memory_order_relaxed))
    return;
  Hit h = hit_queue[i];
  uint pixel = h.pixel;
  Path path = paths[pixel];
  uint dim = path.dimension;
  Triangle t = triangles[h.triangle];
  Material m = materials[t.material];
  float3 bary = float3(1 - h.u - h.v, h.u, h.v);
  float3 position = bary.x * t.p[0].xyz + bary.y * t.p[1].xyz + bary.z * t.p[2].xyz;
  float3 error = 4.172327e-7f *
                 (abs(bary.x * t.p[0].xyz) + abs(bary.y * t.p[1].xyz) + abs(bary.z * t.p[2].xyz));
  // Metal's float intersection has error as well as barycentric reconstruction.
  // Include triangle extent, even for a zero coordinate on a planar surface;
  // a subnormal-only offset cannot protect a large triangle through the origin.
  float3 extent = abs(t.p[1].xyz - t.p[0].xyz) + abs(t.p[2].xyz - t.p[0].xyz);
  error += 3 * FLT_EPSILON * max(extent.x, max(extent.y, extent.z));
  float3 g = t.geometric.xyz;
  float3 interpolated = bary.x * t.n[0].xyz + bary.y * t.n[1].xyz + bary.z * t.n[2].xyz;
  float3 n = dot(interpolated, interpolated) > 0 ? normalize(interpolated) : g;
  float3 wo = -path.direction.xyz;
  float3 rho = texture_value(m.texture, triangle_uv(t, bary), position, textures, texels);
  if (p.depth == 0)
    pixels[pixel].radiance.w = 1;
  if (m.emissive) {
    if (dot(n, wo) > 0) {
      float area = .5f * length(cross(t.p[1].xyz - t.p[0].xyz, t.p[2].xyz - t.p[0].xyz));
      float pdf =
          distance_squared(position, path.previous.xyz) / (area * abs(dot(g, wo)) * p.lights);
      float weight = p.depth == 0 ? 1 : mis(path.direction.w, pdf);
      pixels[pixel].radiance.xyz += path.beta.xyz * rho * m.diffuse.w * weight;
    }
    return;
  }
  if (p.depth >= p.max_depth)
    return; // Terminal vertices still contribute emission.
  rho = select(float3(0), clamp(rho, 0.f, 1.f), isfinite(rho));
  if (dot(g, wo) * dot(n, wo) <= 0)
    return;
  if (dot(g, wo) < 0) {
    g = -g;
    n = -n;
  }
  if (p.depth == 0) {
    pixels[pixel].normal = float4(n, 0);
    pixels[pixel].albedo = float4(rho, 0);
  }
  // One light proposal per vertex. Each emitter has its own competing BSDF
  // technique; evaluating emission individually on a miss avoids mixture bias.
  if (p.lights) {
    uint choice = min(uint(sample1(pixel, dim, p, directions, blue) * p.lights), p.lights - 1);
    float u = sample1(pixel, dim, p, directions, blue),
          v = sample1(pixel, dim, p, directions, blue);
    Light l = lights[choice];
    float3 wi, emitted = 0;
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
      Triangle emitter = triangles[l.triangle];
      Material lm = materials[emitter.material];
      float su = sqrt(u);
      float3 b = float3(1 - su, su * (1 - v), su * v);
      float3 target = b.x * emitter.p[0].xyz + b.y * emitter.p[1].xyz + b.z * emitter.p[2].xyz;
      float3 delta = target - position;
      maximum = length(delta);
      wi = delta / maximum;
      float3 en =
          normalize(b.x * emitter.n[0].xyz + b.y * emitter.n[1].xyz + b.z * emitter.n[2].xyz);
      float cosine = dot(emitter.geometric.xyz, -wi);
      float area = .5f * length(cross(emitter.p[1].xyz - emitter.p[0].xyz,
                                      emitter.p[2].xyz - emitter.p[0].xyz));
      if (dot(en, -wi) > 0 && abs(cosine) > 0) {
        emitted = texture_value(lm.texture, triangle_uv(emitter, b), target, textures, texels) *
                  lm.diffuse.w;
        pdf = maximum * maximum / (area * abs(cosine) * p.lights);
      }
    }
    float cosine = dot(n, wi);
    if (pdf > 0 && cosine > 0 && dot(g, wi) > 0 && any(emitted > 0)) {
      float weight = l.type == 2 ? 1 : mis(pdf, cosine / pi);
      float3 contribution =
          path.beta.xyz * brdf(wo, wi, n, rho, m) * cosine * emitted * (weight / pdf);
      uint slot = atomic_fetch_add_explicit(&counts.shadows, 1, memory_order_relaxed);
      if (slot >= p.active)
        atomic_store_explicit(&counts.error, 1, memory_order_relaxed);
      else {
        float3 origin = offset_origin(position, error, g);
        // Shorten only finite light segments, accounting for the origin shift.
        float max_t = isfinite(maximum)
                          ? max(0.f, maximum - dot(origin - position, wi)) * (1 - 8 * FLT_EPSILON)
                          : maximum;
        shadow_queue[slot] = {
            float4(origin, 0), float4(wi, max_t), float4(contribution, 0), pixel, {0, 0, 0}};
      }
    }
  }
  float u = sample1(pixel, dim, p, directions, blue), v = sample1(pixel, dim, p, directions, blue);
  float3 wi = basis_direction(n, float2(u, v));
  float pdf = dot(wi, n) / pi;
  if (dot(g, wi) <= 0 || pdf <= 0)
    return; // Null event: do not resample clipped directions.
  path.beta.xyz *= brdf(wo, wi, n, rho, m) * pi;
  float continuation = min(1.f, max(path.beta.x, max(path.beta.y, path.beta.z)));
  if (p.depth >= p.roulette) {
    if (sample1(pixel, dim, p, directions, blue) >= continuation)
      return;
    path.beta.xyz /= continuation;
  }
  if (!all(isfinite(path.beta.xyz)) || all(path.beta.xyz == 0))
    return;
  path.previous = float4(position, 0);
  path.origin = float4(offset_origin(position, error, g), 0);
  path.direction = float4(wi, pdf);
  path.dimension = dim;
  paths[pixel] = path;
  uint slot = atomic_fetch_add_explicit(&counts.next, 1, memory_order_relaxed);
  if (slot >= p.active)
    atomic_store_explicit(&counts.error, 1, memory_order_relaxed);
  else
    next_queue[slot] = pixel;
}

kernel void shadows(SCENE_ARGUMENTS) {
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
  intersector<triangle_data> query;
  query.assume_geometry_type(geometry_type::triangle);
  query.force_opacity(forced_opacity::opaque);
  query.accept_any_intersection(true);
  if (query.intersect(r, scene).type == intersection_type::none)
    pixels[s.pixel].radiance.xyz += s.contribution.xyz;
}

// A one-thread queue transition emits MTLDispatchThreadgroupsIndirectArguments.
// No queue size crosses to the CPU between bounces. At least one group permits
// empty dispatches on every supported device; consumers guard against count=0.
kernel void dispatch_size(device Counts &c [[buffer(0)]], device uint *args [[buffer(1)]],
                          constant uint &stage [[buffer(2)]]) {
  uint count;
  if (stage == 3) {
    count = atomic_load_explicit(&c.next, memory_order_relaxed);
    atomic_store_explicit(&c.rays, count, memory_order_relaxed);
    atomic_store_explicit(&c.hits, 0, memory_order_relaxed);
    atomic_store_explicit(&c.shadows, 0, memory_order_relaxed);
    atomic_store_explicit(&c.next, 0, memory_order_relaxed);
  } else if (stage == 0)
    count = atomic_load_explicit(&c.rays, memory_order_relaxed);
  else if (stage == 1)
    count = atomic_load_explicit(&c.hits, memory_order_relaxed);
  else
    count = atomic_load_explicit(&c.shadows, memory_order_relaxed);
  args[0] = max(1u, (count + 63) / 64);
  args[1] = 1;
  args[2] = 1;
}
)MSL"
