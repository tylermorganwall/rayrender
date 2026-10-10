R"MSL(

  struct Surface {
    float4 inverse[3], forward[3], bump_inverse[3], bump_forward[3], axes;
    uint alpha, bump, bump_transformed, padding;
    float4 bump_parameters;
  };
  struct TextureContext {
    float2 uv;
    float3 position, object_position, normal, object_normal;
    float4 footprint; // dudx, dvdx, dudy, dvdy; zero means no differentials.
  };
  static_assert(sizeof(Surface) == 240, "wavefront surface ABI mismatch");
  TextureContext texture_context(float2 uv, float3 position) {
    return {uv, position, position, float3(0, 1, 0), float3(0, 1, 0), float4(0)};
  }
  float3 object_point(Surface s, float3 p) {
    return float3(dot(s.inverse[0], float4(p, 1)), dot(s.inverse[1], float4(p, 1)),
                  dot(s.inverse[2], float4(p, 1)));
  }
  TextureContext surface_context(Triangle t, float3 bary, float3 position,
                                 device const Surface *surfaces) {
    TextureContext h =
        texture_context(bary.x * t.uv01.xy + bary.y * t.uv01.zw + bary.z * t.uv2.xy, position);
    h.normal = t.geometric.xyz;
    if (texture_features & 2) {
      Surface s = surfaces[as_type<uint>(t.uv2.w)];
      h.object_position = object_point(s, position);
      // Normals transform by the transpose of object-to-world here.
      h.object_normal = normalize(s.forward[0].xyz * h.normal.x + s.forward[1].xyz * h.normal.y +
                                  s.forward[2].xyz * h.normal.z);
      // Analytic primitives define direction textures using the smooth surface
      // normal rather than the tessellation's facet normal.
      if (s.axes.w == 1 || s.axes.w == 2) {
        h.object_normal = normalize(h.object_position / (s.axes.xyz * s.axes.xyz));
        float3 chart = s.axes.w == 1 ? normalize(h.object_position) : h.object_normal;
        h.uv = float2(1 - (atan2(chart.z, chart.x) + pi) / (2 * pi),
                      .5f + asin(clamp(chart.y, -1.f, 1.f)) / pi);
        float3 local = h.object_normal;
        h.normal = normalize(s.inverse[0].xyz * local.x + s.inverse[1].xyz * local.y +
                             s.inverse[2].xyz * local.z);
        if (dot(h.normal, t.geometric.xyz) < 0) {
          h.normal = -h.normal;
          h.object_normal = -h.object_normal;
        }
      }
    }
    return h;
  }
  float wrap(float v) {
    float w = v - floor(v);
    return w == 0 && v > 0 ? 1 : w;
  }
  float3 image_centers(Texture t, float2 uv, bool repeat, device const float4 *texels) {
    uv = repeat ? fract(uv) : clamp(uv, 0.f, 1.f);
    float2 pixel = float2(uv.x, 1 - uv.y) * float2(t.width, t.height) - .5f;
    int2 lower = int2(floor(pixel));
    float2 w = fract(pixel);
    float3 sum = 0;
    for (int y = 0; y < 2; ++y)
      for (int x = 0; x < 2; ++x) {
        int2 tap = lower + int2(x, y), size = int2(t.width, t.height);
        tap = repeat ? (tap % size + size) % size : clamp(tap, int2(0), size - 1);
        sum += texels[t.offset + tap.y * t.width + tap.x].xyz * (x ? w.x : 1 - w.x) *
               (y ? w.y : 1 - w.y);
      }
    return sum;
  }
  uint noise_hash(uint x) {
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    return x ^ (x >> 16);
  }
  float graph_noise(float3 p, uint seed) {
    float3 base = floor(p), f = p - base;
    uint3 cell = uint3(fmod(fmod(base, 1048576.f) + 1048576.f, 1048576.f));
    f = f * f * f * (f * (f * 6 - 15) + 10);
    float result = 0;
    for (uint z = 0; z < 2; ++z)
      for (uint y = 0; y < 2; ++y)
        for (uint x = 0; x < 2; ++x) {
          uint3 q = (cell + uint3(x, y, z)) & 0xfffffu;
          uint hash = noise_hash(seed ^ noise_hash(q.x) ^ noise_hash(noise_hash(q.y)) ^
                                 noise_hash(noise_hash(noise_hash(q.z))));
          result += float(hash >> 8) / float(0xffffffu) * (x ? f.x : 1 - f.x) *
                    (y ? f.y : 1 - f.y) * (z ? f.z : 1 - f.z);
        }
    return result;
  }
  float3 graph_value(Texture root, TextureContext h, device const Texture *textures,
                     device const float4 *texels) {
    float3 values[128];
    for (uint i = 0; i < root.width; ++i) {
      Texture t = textures[root.offset + i];
      uint op = (t.type & 255) - 16, space = (t.type >> 8) & 3;
      bool flag = (t.type & 1024) != 0, scalar = (t.type & 2048) != 0;
      uint4 reg = uint4(t.mapping);
      float3 p = space == 0 ? float3(h.uv, 0) : space == 1 ? h.object_position : h.position;
      p *= t.b.xyz;
      p = float3(t.b.w * p.x - t.a.w * p.y, t.a.w * p.x + t.b.w * p.y, p.z) + t.c.xyz;
      float3 result = 0;
      switch (op) {
      case 0:
        result = t.a.xyz;
        break;
      case 1: {
        float w = clamp(values[reg.z].x, 0.f, 1.f);
        result = w == 0   ? values[reg.x]
                 : w == 1 ? values[reg.y]
                          : mix(values[reg.x], values[reg.y], w);
        break;
      }
      case 2:
        result = values[reg.x] * values[reg.y].x;
        break;
      case 3:
        result = float3(pow(max(0.f, values[reg.x].x), t.a.x));
        break;
      case 4: {
        float3 c = values[reg.x];
        result = float3(t.width < 3    ? c[t.width]
                        : t.width == 3 ? dot(c, float3(1.f / 3))
                                       : dot(c, float3(.2126f, .7152f, .0722f)));
        break;
      }
      case 5: {
        float3 n = space == 1 ? h.object_normal : h.normal;
        float weight = dot(n, n) > 0 ? dot(normalize(n), t.a.xyz) : 0;
        result = float3(flag ? abs(weight) : max(0.f, weight));
        break;
      }
      case 6:
        result =
            fmod(floor(p.x) + floor(p.y) + floor(p.z), 2.f) == 0 ? values[reg.x] : values[reg.y];
        break;
      case 7:
        result = float3(clamp(p[t.width], 0.f, 1.f));
        break;
      case 8: {
        float sum = 0, amplitude = 1, total = 0;
        for (uint octave = 0; octave < t.width; ++octave) {
          sum += amplitude * graph_noise(p, t.offset);
          total += amplitude;
          amplitude *= .5f;
          p *= 2;
        }
        result = float3(sum / total);
        break;
      }
      case 9:
        result = image_centers(t, p.xy, flag, texels);
        if (scalar)
          result = float3(dot(result, float3(1.f / 3)));
        break;
      }
      values[reg.w] = result;
    }
    return values[root.height];
  }
  float roughness_alpha(float x) {
    x = log(max(x, .0001550155f));
    float mapped = 1.62142f + .819955f * x + .1734f * x * x + .0171201f * x * x * x +
                   .000640711f * x * x * x * x;
    return mapped * mapped;
  }
  float3 texture_value(uint id, TextureContext h, device const Texture *textures,
                       device const float4 *texels) {
    Texture t = textures[id];
    float2 uv = h.uv;
    float3 position = h.position;
    if ((texture_features & 1) && t.type == 8)
      return graph_value(t, h, textures, texels);
    if ((texture_features & 1) && t.type == 9) {
      float x = clamp(graph_value(textures[t.offset], h, textures, texels).x, 0.f, 1.f);
      float y =
          t.height ? clamp(graph_value(textures[t.width], h, textures, texels).x, 0.f, 1.f) : x;
      return float3(
          max(float2(t.height ? x : roughness_alpha(x * x), t.height ? y : roughness_alpha(y * y)),
              1e-6f),
          0);
    }
    if ((texture_features & 4) && t.type == 10) {
      float extent =
          2 * max(max(abs(t.mapping.x * h.footprint.x), abs(t.mapping.x * h.footprint.z)),
                  max(abs(t.mapping.y * h.footprint.y), abs(t.mapping.y * h.footprint.w)));
      uint lod = uint(
          clamp(floor(float(t.width - 1) + log2(max(extent, 1e-8f))), 0.f, float(t.width - 1)));
      t = textures[t.offset + lod];
    }
    if ((texture_features & 4) && t.type == 11)
      return image_centers(t, uv * t.mapping.xy + t.mapping.zw, true, texels);
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
    if (t.type == 7) {
      uv = uv * t.mapping.xy + t.mapping.zw;
      uint x = min(uint(wrap(uv.x) * t.width), t.width - 1),
           y = min(uint((1 - wrap(uv.y)) * t.height), t.height - 1);
      return texels[t.offset + y * t.width + x].xyz;
    }
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
    return mix(mix(texels[t.offset + y0 * t.width + x0].xyz,
                   texels[t.offset + y0 * t.width + x1].xyz, w.x),
               mix(texels[t.offset + y1 * t.width + x0].xyz,
                   texels[t.offset + y1 * t.width + x1].xyz, w.x),
               w.y);
  }
  float3 texture_value(uint id, float2 uv, float3 position, device const Texture *textures,
                       device const float4 *texels) {
    return texture_value(id, texture_context(uv, position), textures, texels);
  }

)MSL"
