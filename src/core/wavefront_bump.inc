R"MSL(

  struct SurfaceDerivatives {
    float3 du, dv, nu, nv, dx, dy;
    float4 uv;
    bool valid;
  };
  float3 transform_vector(Surface s, float3 v) {
    return float3(dot(s.forward[0].xyz, v), dot(s.forward[1].xyz, v), dot(s.forward[2].xyz, v));
  }
  float3 transform_normal(Surface s, float3 n) {
    return s.inverse[0].xyz * n.x + s.inverse[1].xyz * n.y + s.inverse[2].xyz * n.z;
  }
  SurfaceDerivatives surface_derivatives(Triangle t, Surface s, TextureContext h, float3 base,
                                         Differential ray) {
    SurfaceDerivatives d = {};
    d.valid = false;
    float2 uv0 = t.uv01.xy - t.uv2.xy, uv1 = t.uv01.zw - t.uv2.xy;
    if (t.flags & 16) {
      uv0 = float2(t.p[0].w - t.p[2].w, t.n[0].w - t.n[2].w);
      uv1 = float2(t.p[1].w - t.p[2].w, t.n[1].w - t.n[2].w);
    }
    float3 p0 = t.p[0].xyz - t.p[2].xyz, p1 = t.p[1].xyz - t.p[2].xyz;
    float det = uv0.x * uv1.y - uv0.y * uv1.x;
    if (t.flags & 32) {
      // Native curves carry their CPU-compatible derivatives explicitly;
      // their UV copies are equal so material lookup never interpolates u/v.
      d.du = t.p[1].xyz - t.p[0].xyz;
      d.dv = t.p[2].xyz - t.p[0].xyz;
    } else if (abs(det) >= 1e-8f) {
      d.du = (uv1.y * p0 - uv0.y * p1) / det;
      d.dv = (-uv1.x * p0 + uv0.x * p1) / det;
      d.nu = (uv1.y * (t.n[0].xyz - t.n[2].xyz) - uv0.y * (t.n[1].xyz - t.n[2].xyz)) / det;
      d.nv = (-uv1.x * (t.n[0].xyz - t.n[2].xyz) + uv0.x * (t.n[1].xyz - t.n[2].xyz)) / det;
    } else {
      auto frame = normal_frame(base);
      d.du = frame.x;
      d.dv = frame.y;
    }
    if (s.axes.w == 1 || s.axes.w == 2) {
      float3 squared = s.axes.xyz * s.axes.xyz;
      // UV position derivatives keep their chart orientation on reversed
      // surfaces; only derivatives of the oriented normal change sign.
      float3 n = normalize(h.object_position / squared);
      float horizontal = length(n.xz), c = horizontal > 0 ? n.x / horizontal : 1,
            si = horizontal > 0 ? n.z / horizontal : 0;
      float3 nu = 2 * pi * float3(n.z, 0, -n.x), nv = pi * float3(-n.y * c, horizontal, -n.y * si);
      float3 an = squared * n;
      float denominator = sqrt(dot(n, an));
      d.du = transform_vector(s, (squared * nu - an * (dot(an, nu) / (denominator * denominator))) /
                                     denominator);
      d.dv = transform_vector(s, (squared * nv - an * (dot(an, nv) / (denominator * denominator))) /
                                     denominator);
      d.nu = transform_normal(s, nu);
      d.nv = transform_normal(s, nv);
      if (dot(transform_normal(s, n), base) < 0) {
        d.nu = -d.nu;
        d.nv = -d.nv;
      }
    }
    // Intersect adjacent rays with the geometric tangent plane, then solve the
    // UV Jacobian in the two coordinates least parallel to its normal.
    float3 g = (t.flags & 1) ? base : t.geometric.xyz;
    float cx = dot(g, ray.dx.xyz), cy = dot(g, ray.dy.xyz);
    if (ray.ox.w != 0 && abs(cx) > 1e-10f && abs(cy) > 1e-10f) {
      d.dx = ray.ox.xyz + ray.dx.xyz * (dot(g, h.position - ray.ox.xyz) / cx) - h.position;
      d.dy = ray.oy.xyz + ray.dy.xyz * (dot(g, h.position - ray.oy.xyz) / cy) - h.position;
      float3 a = abs(g);
      uint drop = a.x > a.y ? (a.x > a.z ? 0 : 2) : (a.y > a.z ? 1 : 2);
      uint x = (drop + 1) % 3, y = (drop + 2) % 3;
      float determinant = d.du[x] * d.dv[y] - d.du[y] * d.dv[x];
      float scale = max(abs(d.du[x] * d.dv[y]), abs(d.du[y] * d.dv[x]));
      if (abs(determinant) > 1e-12f * scale) {
        d.uv = float4((d.dx[x] * d.dv[y] - d.dx[y] * d.dv[x]) / determinant,
                      (d.du[x] * d.dx[y] - d.du[y] * d.dx[x]) / determinant,
                      (d.dy[x] * d.dv[y] - d.dy[y] * d.dv[x]) / determinant,
                      (d.du[x] * d.dy[y] - d.du[y] * d.dy[x]) / determinant);
        d.valid = all(isfinite(d.uv)) && all(abs(d.uv) <= 1e8f);
      }
    }
    if (!d.valid) {
      d.uv = 0;
      d.dx = d.dy = 0;
    }
    return d;
  }
  float3 bump_normal(Triangle t, Surface s, thread TextureContext &h, float3 base,
                     thread SurfaceDerivatives &d, device const Texture *textures,
                     device const float4 *texels) {
    h.footprint = d.uv;
    if (s.bump == UINT_MAX)
      return base;
    float3 height_normal = base, tu = d.du, tv = d.dv, nu = d.nu, nv = d.nv;
    TextureContext origin = h;
    if (s.bump_transformed) {
      // Native instances transform a completed hit, including its bumped
      // normal. Evaluate height in that same pre-instance frame so nonuniform
      // scaling changes bump dimensions along with the rest of the object.
      origin.position = float3(dot(s.bump_inverse[0], float4(h.position, 1)),
                               dot(s.bump_inverse[1], float4(h.position, 1)),
                               dot(s.bump_inverse[2], float4(h.position, 1)));
      height_normal = normalize(s.bump_forward[0].xyz * base.x + s.bump_forward[1].xyz * base.y +
                                s.bump_forward[2].xyz * base.z);
      tu = float3(dot(s.bump_inverse[0].xyz, d.du), dot(s.bump_inverse[1].xyz, d.du),
                  dot(s.bump_inverse[2].xyz, d.du));
      tv = float3(dot(s.bump_inverse[0].xyz, d.dv), dot(s.bump_inverse[1].xyz, d.dv),
                  dot(s.bump_inverse[2].xyz, d.dv));
      nu = s.bump_forward[0].xyz * d.nu.x + s.bump_forward[1].xyz * d.nu.y +
           s.bump_forward[2].xyz * d.nu.z;
      nv = s.bump_forward[0].xyz * d.nv.x + s.bump_forward[1].xyz * d.nv.y +
           s.bump_forward[2].xyz * d.nv.z;
    }
    if (t.flags & 8) {
      float3 bitangent = cross(height_normal, tu);
      if (dot(bitangent, bitangent) > 0) {
        tu = cross(bitangent, height_normal);
        tv = bitangent;
      } else {
        auto frame = normal_frame(height_normal);
        tu = frame.x;
        tv = frame.y;
      }
    }
    origin.normal = height_normal;
    float du = .5f * (abs(d.uv.x) + abs(d.uv.z)), dv = .5f * (abs(d.uv.y) + abs(d.uv.w));
    if (du == 0)
      du = .0005f;
    if (dv == 0)
      dv = .0005f;
    TextureContext u = origin, v = origin;
    u.uv.x += du;
    v.uv.y += dv;
    u.position += du * tu;
    v.position += dv * tv;
    u.normal = normalize(height_normal + du * nu);
    v.normal = normalize(height_normal + dv * nv);
    float height = s.bump_parameters.x * texture_value(s.bump, origin, textures, texels).x;
    float bu = (s.bump_parameters.x * texture_value(s.bump, u, textures, texels).x - height) / du;
    float bv = (s.bump_parameters.x * texture_value(s.bump, v, textures, texels).x - height) / dv;
    tu += bu * height_normal + height * nu;
    tv += bv * height_normal + height * nv;
    float3 crossed = cross(tu, tv);
    if (!(dot(crossed, crossed) > 0) || !all(isfinite(crossed)))
      return base;
    float3 bumped = normalize(crossed);
    if (dot(bumped, height_normal) < 0)
      bumped = -bumped;
    if (s.bump_transformed) {
      d.du = float3(dot(s.bump_forward[0].xyz, tu), dot(s.bump_forward[1].xyz, tu),
                    dot(s.bump_forward[2].xyz, tu));
      d.dv = float3(dot(s.bump_forward[0].xyz, tv), dot(s.bump_forward[1].xyz, tv),
                    dot(s.bump_forward[2].xyz, tv));
      bumped = normalize(s.bump_inverse[0].xyz * bumped.x + s.bump_inverse[1].xyz * bumped.y +
                         s.bump_inverse[2].xyz * bumped.z);
    } else {
      d.du = tu;
      d.dv = tv;
    }
    return bumped;
  }
  // Same projected-area microgeometry as normalmap::Model on the CPU. The
  // geometric frame here is the unbumped smooth normal. Actual triangle-side
  // clipping remains at the call site and adds null events, not energy.
  struct MappedDiffuse {
    float3 g, p, t, wo;
    float c, s, lambda;
    bool identity;
  };
  MappedDiffuse mapped_diffuse(float3 base, float3 bumped, float3 wo) {
    MappedDiffuse m;
    m.g = base;
    m.p = bumped;
    m.wo = wo;
    if (dot(base, wo) < 0) {
      m.g = -base;
      m.p = -bumped;
    }
    m.c = min(1.f, dot(m.g, m.p));
    float3 tilt = m.p - m.c * m.g;
    m.s = length(tilt);
    m.identity = m.c <= 0 || m.s <= 1e-7f;
    m.lambda = 1;
    m.t = 0;
    if (m.identity) {
      m.p = m.g;
      m.c = 1;
      m.s = 0;
    } else {
      m.t = -tilt / m.s;
      float ap = max(0.f, dot(m.p, wo)), at = m.s * max(0.f, dot(m.t, wo));
      m.lambda = ap / (ap + at);
    }
    return m;
  }
  float mapped_mask(MappedDiffuse m, float3 w) {
    if (m.identity)
      return dot(m.g, w) > 0 ? 1 : 0;
    float area = max(0.f, dot(m.p, w)) + m.s * max(0.f, dot(m.t, w));
    return area > 0 ? min(1.f, max(0.f, dot(m.g, w)) * m.c / area) : 0;
  }
  BsdfValue mapped_eval(MappedDiffuse b, float3 wi, float3 rho, Material m) {
    float cosine = dot(b.g, wi);
    if (cosine <= 0)
      return {0, 0};
    float cp = max(0.f, dot(b.p, wi));
    if (b.identity)
      return {brdf(b.wo, wi, b.g, rho, m) * cosine, cp / pi};
    float3 value = b.lambda * brdf(b.wo, wi, b.p, rho, m) * cp;
    float pdf = b.lambda * cp / pi * mapped_mask(b, wi);
    if (dot(b.t, wi) > 0) {
      float3 reflected = wi - 2 * dot(wi, b.t) * b.t;
      float weight = b.lambda * max(0.f, dot(b.p, reflected)) * (1 - mapped_mask(b, reflected));
      value += weight * brdf(b.wo, reflected, b.p, rho, m);
      pdf += weight / pi;
    }
    if (dot(b.t, b.wo) > 0) {
      float3 reflected = b.wo - 2 * dot(b.wo, b.t) * b.t;
      value += (1 - b.lambda) * brdf(reflected, wi, b.p, rho, m) * cp;
      pdf += (1 - b.lambda) * cp / pi;
    }
    float area = cp + b.s * max(0.f, dot(b.t, wi));
    return {value * min(1.f, cosine * b.c / area), pdf};
  }
  float3 mapped_sample(MappedDiffuse m, float2 uv, float facet, float escape) {
    float3 wi = basis_direction(m.p, uv);
    if (!m.identity && facet < m.lambda && escape >= mapped_mask(m, wi))
      wi -= 2 * dot(wi, m.t) * m.t;
    return wi;
  }
  // The same three facet paths as mapped_eval, with a directional glossy
  // child instead of the diffuse child's view-independent cosine proposal.
  // Child evaluation already includes the primary-facet cosine.
  BsdfValue mapped_glossy_child(SurfaceBsdf b, float3 view, float3 wi) {
    b.view = b.frame.local(view);
    if (b.view.z <= 0 || b.frame.local(wi).z <= 0)
      return {0, 0};
    return surface_eval(b, wi);
  }
  BsdfValue mapped_glossy_eval(MappedDiffuse m, SurfaceBsdf b, float3 wi) {
    if (dot(m.g, wi) <= 0)
      return {0, 0};
    b.frame = normal_frame(m.p);
    BsdfValue primary = mapped_glossy_child(b, m.wo, wi);
    if (m.identity) return primary;
    BsdfValue result = {m.lambda * primary.value,
                        m.lambda * primary.pdf * mapped_mask(m, wi)};
    if (dot(m.t, wi) > 0) {
      float3 reflected = wi - 2 * dot(m.t, wi) * m.t;
      BsdfValue child = mapped_glossy_child(b, m.wo, reflected);
      float weight = m.lambda * (1 - mapped_mask(m, reflected));
      result.value += weight * child.value;
      result.pdf += weight * child.pdf;
    }
    if (dot(m.t, m.wo) > 0) {
      float3 reflected = m.wo - 2 * dot(m.t, m.wo) * m.t;
      BsdfValue child = mapped_glossy_child(b, reflected, wi);
      result.value += (1 - m.lambda) * child.value;
      result.pdf += (1 - m.lambda) * child.pdf;
    }
    result.value *= mapped_mask(m, wi);
    return result;
  }
  BsdfSample mapped_glossy_sample(MappedDiffuse m, SurfaceBsdf b, float4 u) {
    bool primary = m.identity || u.x < m.lambda;
    if (!m.identity)
      u.x = primary ? u.x / m.lambda : (u.x - m.lambda) / (1 - m.lambda);
    b.frame = normal_frame(m.p);
    b.view = b.frame.local(primary ? m.wo : m.wo - 2 * dot(m.t, m.wo) * m.t);
    BsdfSample sampled = surface_sample(b, u);
    if (!(sampled.pdf > 0) || b.frame.local(sampled.direction).z <= 0)
      return {0, 0, 0, 1, false};
    float3 wi = sampled.direction;
    if (!m.identity && primary && u.w >= mapped_mask(m, wi))
      wi -= 2 * dot(m.t, wi) * m.t;
    BsdfValue evaluated = mapped_glossy_eval(m, b, wi);
    return {wi, evaluated.pdf > 0 ? evaluated.value / evaluated.pdf : float3(0),
             evaluated.pdf, 1, false};
  }
  // Conditional 16-hit scattering closure, matching normalmap::Scattering.
  // Prefix randomness is sampled once, before NEE. Both direct-light evaluation
  // and continuation use this same positive expansion and normalized mixture.
  struct MappedLobe {
    float3 view, weight;
    float mass;
    bool mirror_exit;
  };
  struct MappedScattering {
    MappedLobe lobes[16];
    uint count;
    float total;
  };
  float mapped_random(thread uint& state) {
    state = state * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    word = (word >> 22u) ^ word;
    return float(word >> 8u) * 0x1p-24f;
  }
  BsdfValue mapped_child_value(MappedDiffuse m, SurfaceBsdf b, float3 view, float3 wi) {
    if (dot(m.p, view) <= 0 || dot(m.p, wi) <= 0) return {0, 0};
    if (b.material.type == bsdf_glossy) {
      b.frame = normal_frame(m.p);
      return mapped_glossy_child(b, view, wi);
    }
    float cosine = dot(m.p, wi);
    return {brdf(view, wi, m.p, b.rho, b.material) * cosine, cosine / pi};
  }
  BsdfSample mapped_child_sample(MappedDiffuse m, SurfaceBsdf b, float3 view, float3 u) {
    if (dot(m.p, view) <= 0) return {0, 0, 0, 1, false};
    if (b.material.type == bsdf_glossy) {
      b.frame = normal_frame(m.p);
      b.view = b.frame.local(view);
      return surface_sample(b, float4(u, 0));
    }
    float3 wi = basis_direction(m.p, u.yz);
    BsdfValue value = mapped_child_value(m, b, view, wi);
    return {wi, value.pdf > 0 ? value.value / value.pdf : float3(0), value.pdf, 1, false};
  }
  void prepare_mapped_scattering(thread MappedScattering& out, MappedDiffuse m,
                                SurfaceBsdf b, uint seed) {
    out.count = 0;
    out.total = 0;
    for (uint branch = 0; branch < 2; ++branch) {
      if (branch && m.identity) break;
      float3 view = branch ? m.wo - 2 * dot(m.t, m.wo) * m.t : m.wo;
      float3 weight = branch ? 1 - m.lambda : m.lambda;
      uint remaining = m.identity ? 1 : 16 - branch;
      while (remaining > 0) {
        float mass = max(weight.x, max(weight.y, weight.z));
        if (!(mass > 0) || !isfinite(mass)) break;
        out.lobes[out.count++] = {view, weight, mass, remaining > 1};
        out.total += mass;
        if (remaining <= 2) break;
        float u = mapped_random(seed), v = mapped_random(seed), w = mapped_random(seed);
        BsdfSample child = mapped_child_sample(m, b, view, float3(u, v, w));
        if (!(child.pdf > 0) || dot(m.p, child.direction) <= 0) break;
        float3 reflected = child.direction - 2 * dot(m.t, child.direction) * m.t;
        weight *= child.weight * (1 - mapped_mask(m, child.direction)) * (1 - mapped_mask(m, reflected));
        view = -reflected;
        remaining -= 2;
      }
    }
  }
  BsdfValue mapped_scattering_eval(thread const MappedScattering& state, MappedDiffuse m,
                                   SurfaceBsdf b, float3 wi) {
    if (!(state.total > 0) || dot(m.g, wi) <= 0) return {0, 0};
    float3 reflected = wi - 2 * dot(m.t, wi) * m.t;
    float escape = mapped_mask(m, wi);
    BsdfValue result = {0, 0};
    for (uint i = 0; i < state.count; ++i) {
      MappedLobe lobe = state.lobes[i];
      BsdfValue term = mapped_child_value(m, b, lobe.view, wi);
      if (lobe.mirror_exit) {
        term.pdf *= escape;
        if (dot(m.t, wi) > 0) {
          BsdfValue mirror = mapped_child_value(m, b, lobe.view, reflected);
          float probability = 1 - mapped_mask(m, reflected);
          term.value += mirror.value * probability;
          term.pdf += mirror.pdf * probability;
        }
      }
      result.value += lobe.weight * term.value;
      result.pdf += lobe.mass * term.pdf;
    }
    result.value *= escape;
    result.pdf /= state.total;
    return result;
  }
  BsdfSample mapped_scattering_sample(thread const MappedScattering& state, MappedDiffuse m,
                                       SurfaceBsdf b, float component, float escape, float3 u) {
    if (!(state.total > 0)) return {0, 0, 0, 1, false};
    float target = component * state.total;
    uint selected = 0;
    while (selected + 1 < state.count && target >= state.lobes[selected].mass)
      target -= state.lobes[selected++].mass;
    MappedLobe lobe = state.lobes[selected];
    BsdfSample child = mapped_child_sample(m, b, lobe.view, u);
    if (!(child.pdf > 0) || dot(m.p, child.direction) <= 0) return {0, 0, 0, 1, false};
    float3 wi = child.direction;
    if (lobe.mirror_exit && escape >= mapped_mask(m, wi)) wi -= 2 * dot(m.t, wi) * m.t;
    BsdfValue value = mapped_scattering_eval(state, m, b, wi);
    return {wi, value.pdf > 0 ? value.value / value.pdf : float3(0), value.pdf, 1, false};
  }
  Differential propagate_differential(Differential incoming, SurfaceDerivatives d, float3 n,
                                      float3 old_direction, float3 origin, float3 direction,
                                      bool transmission, float eta) {
    Differential out = {};
    if (incoming.ox.w == 0 || !d.valid)
      return out;
    float3 nx = normalize(n + d.nu * d.uv.x + d.nv * d.uv.y),
           ny = normalize(n + d.nu * d.uv.z + d.nv * d.uv.w);
    float3 dirs[3] = {old_direction, normalize(incoming.dx.xyz), normalize(incoming.dy.xyz)};
    float3 normals[3] = {n, nx, ny};
    for (uint k = 0; k < 3; ++k) {
      float3 normal = dot(dirs[k], normals[k]) > 0 ? -normals[k] : normals[k];
      dirs[k] = transmission ? refract(dirs[k], normal, 1 / eta) : reflect(dirs[k], normal);
      if (dot(dirs[k], dirs[k]) == 0 || !all(isfinite(dirs[k])))
        return out;
    }
    out.ox = float4(origin + d.dx, 1);
    out.oy = float4(origin + d.dy, 0);
    out.dx = float4(direction + dirs[1] - dirs[0], 0);
    out.dy = float4(direction + dirs[2] - dirs[0], 0);
    return out;
  }

)MSL"
