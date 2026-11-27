R"MSL(
struct CsgInstruction {
  float4 a, b, c, d;
  uint operation, padding[3];
};
static_assert(sizeof(CsgInstruction) == 80, "CSG instruction ABI mismatch");
float csg_tube_segment(CsgInstruction n, float3 p) {
  float3 v = p - n.a.xyz;
  if (n.operation == 2) return length(v) - n.a.w;
  float3 ba = n.b.xyz - n.a.xyz;
  float l2 = dot(ba, ba), rr = n.a.w - n.b.w, a2 = l2 - rr * rr;
  float y = dot(v, ba), z = y - l2;
  float x2 = dot(v * l2 - ba * y, v * l2 - ba * y), y2 = y * y * l2, z2 = z * z * l2;
  float k = sign(rr) * rr * rr * x2;
  if (sign(z) * a2 * z2 > k) return sqrt(x2 + z2) / l2 - n.b.w;
  if (sign(y) * a2 * y2 < k) return sqrt(x2 + y2) / l2 - n.a.w;
  return (sqrt(x2 * a2 / l2) + y * rr) / l2 - n.a.w;
}
// A postfix stack machine: primitive distances push, Boolean operations pop
// two values, and coordinate modifiers save/restore the current query point.
// The host validates nesting <=32; these stacks cannot overflow for valid code.
float csg_material_random(thread uint &seed) {
  uint x = (seed += 0x9e3779b9u);
  x = (x ^ (x >> 16)) * 0x7feb352du;
  x = (x ^ (x >> 15)) * 0x846ca68bu;
  return float((x ^ (x >> 16)) >> 8) * (1.f / 16777216.f);
}
// Only the final surface query carries material choices. The distance-only
// specialization eliminates this stack and all random draws during marching.
template<bool choose_material>
float csg_evaluate(device const CsgInstruction *code, uint count, float3 p,
                   uint fallback, thread uint &seed, thread uint &selected) {
  float values[64];
  uint choices[choose_material ? 64 : 1];
  float3 points[32];
  uint top = 0, frame = 0;
  for (uint i = 0; i < count; ++i) {
    CsgInstruction n = code[i];
    uint op = n.operation;

    // Two compile-time query modes share one immutable program. Material
    // queries retain original order; distance queries use the merged curve BVH.
    if (op == 25) {
      if (choose_material) i += n.padding[0] - 1;
      continue;
    }
    if (op == 26) {
      i += n.padding[0] - 1;
      continue;
    }
    if (op == 23) {
      float best = INFINITY;
      uint end = i + n.padding[0];
      for (uint node = i + 1; node < end;) {
        CsgInstruction bounds = code[node];
        float3 delta = max(bounds.a.xyz - p, p - bounds.b.xyz);
        float lower = max(delta.x, max(delta.y, delta.z));
        if (lower > best) node += bounds.padding[0];
        else if (bounds.padding[1]) {
          best = min(best, csg_tube_segment(code[node + 1], p));
          node += 2;
        } else ++node;
      }
      if (choose_material) choices[top] = fallback;
      values[top++] = best;
      i = end - 1;
      continue;
    }
    if (op == 22) {
      if (choose_material) choices[top - 1] = n.padding[0];
      continue;
    }
    if (op >= 200) {
      p = points[--frame];
      if (op == 219)
        values[top - 1] *= n.a.w;
      if (op == 216) {
        float3 q = abs(p - n.a.xyz) - n.b.xyz;
        values[top - 1] += min(max(q.x, max(q.y, q.z)), 0.f);
      }
      continue;
    }
    if (op >= 100) {
      points[frame++] = p;
      if (op == 119)
        p /= n.a.w;
      else if (op == 120) {
        float3 q = p - n.a.xyz;
        p = float3(dot(q, n.b.xyz), dot(q, n.c.xyz), dot(q, n.d.xyz)) + n.a.xyz;
      } else if (op == 121)
        p -= n.a.xyz;
      else if (op == 115)
        p -= clamp(p - n.a.xyz, -n.b.xyz, n.b.xyz);
      else if (op == 116)
        p = sign(p - n.a.xyz) * max(abs(p - n.a.xyz) - n.b.xyz, 0.f) + n.a.xyz;
      continue;
    }
    if (op >= 30) {
      float b = values[--top], a = values[top - 1], k = n.a.w;
      if (choose_material) {
        float weight = (op == 30 || op == 33) ? float(a <= b) :
                       (op == 31 || op == 34) ? float(a >= -b) :
                       op == 32 ? float(a >= b) : float(k <= .5f);
        if (!n.padding[1]) {
          if (op == 33 && k > 0) weight = clamp(.5f + .5f * (b - a) / k, 0.f, 1.f);
          if (op == 34 && k > 0) weight = clamp(.5f + .5f * (a + b) / k, 0.f, 1.f);
          if (op == 35) weight = 1 - k;
        }
        if (weight <= 0 || (weight < 1 && csg_material_random(seed) >= weight))
          choices[top - 1] = choices[top];
      }
      if (op == 30)
        a = min(a, b);
      else if (op == 31)
        a = max(a, -b);
      else if (op == 32)
        a = max(a, b);
      else if (op == 35)
        a = mix(a, b, k);
      else {
        float h = k > 0 ? max(k - abs((op == 33 ? b : -b) - a), 0.f) / k : 0.f;
        a = op == 33 ? min(a, b) - h * h * k * .25f : max(a, -b) + h * h * k * .25f;
      }
      values[top - 1] = a;
      continue;
    }
    if (op == 17 || op == 18) {
      values[top - 1] = (op == 18 ? abs(values[top - 1]) : values[top - 1]) - n.a.w;
      continue;
    }
    float3 v = p - n.a.xyz;
    float d = INFINITY;
    if (op == 2)
      d = length(v) - n.a.w;
    else if (op == 3) {
      d = dot(v, n.c.xyz);
      float x = abs(dot(v, n.b.xyz)) - n.b.w * .5f, z = abs(dot(v, n.d.xyz)) - n.d.w * .5f;
      if (x > 0)
        d = sqrt(x * x + d * d);
      if (z > 0)
        d = sqrt(z * z + d * d);
    } else if (op == 4) {
      float3 q = abs(v) - n.b.xyz * .5f;
      d = length(max(q, 0.f)) + min(max(q.x, max(q.y, q.z)), 0.f) - n.a.w;
    } else if (op == 7)
      d = length(float2(length(v.xz) - n.a.w, v.y)) - n.b.w;
    else if (op == 8) {
      float3 ba = n.b.xyz - n.a.xyz;
      d = length(v - ba * clamp(dot(v, ba) / dot(ba, ba), 0.f, 1.f)) - n.a.w;
    } else if (op == 9) {
      float3 ba = n.b.xyz - n.a.xyz;
      float baba = dot(ba, ba), paba = dot(v, ba);
      float x = length(v * baba - ba * paba) - n.a.w * baba, y = abs(paba - baba * .5f) - baba * .5f;
      float x2 = x * x, y2 = y * y * baba;
      float q = max(x, y) < 0 ? -min(x2, y2) : (x > 0 ? x2 : 0) + (y > 0 ? y2 : 0);
      d = sign(q) * sqrt(abs(q)) / baba - n.b.w;
    } else if (op == 10) {
      float3 inv = 1.f / n.b.xyz;
      float k0 = length(v * inv), k1 = length(v * inv * inv);
      d = k0 < 1 ? (k0 - 1) * min(n.b.x, min(n.b.y, n.b.z)) : k0 * (k0 - 1) / k1;
    } else if (op == 11) {
      d = csg_tube_segment(n, p);
    } else if (op == 12) {
      float3 qv = float3(dot(v, n.b.xyz), dot(v, n.c.xyz) - n.b.w, dot(v, n.d.xyz));
      float2 q = float2(n.a.w, -n.b.w), w = float2(length(qv.xz), qv.y);
      float2 a = w - q * clamp(dot(w, q) / dot(q, q), 0.f, 1.f);
      float2 b = w - q * float2(clamp(w.x / q.x, 0.f, 1.f), 1);
      d = sqrt(min(dot(a, a), dot(b, b))) * sign(max(sign(q.y) * (w.x * q.y - w.y * q.x), sign(q.y) * (w.y - q.y)));
    } else if (op == 13) {
      v *= n.b.xyz;
      v.xz = abs(v.xz);
      if (v.z > v.x)
        v = v.zyx;
      v -= float3(.5f, 0, .5f);
      float h = n.a.w, m2 = h * h + .25f;
      float3 q = float3(v.z, h * v.y - .5f * v.x, h * v.x + .5f * v.y);
      float s = max(-q.x, 0.f), t = clamp((q.y - .5f * v.z) / (m2 + .25f), 0.f, 1.f);
      float a = m2 * (q.x + s) * (q.x + s) + q.y * q.y;
      float b = m2 * (q.x + .5f * t) * (q.x + .5f * t) + (q.y - m2 * t) * (q.y - m2 * t);
      float d2 = min(q.y, -q.x * m2 - q.y * .5f) > 0 ? 0.f : min(a, b);
      d = sqrt((d2 + q.z * q.z) / m2) * sign(max(q.z, -v.y));
    } else if (op == 14) {
      float3 ba = n.b.xyz - n.a.xyz, cb = n.c.xyz - n.b.xyz, ac = n.a.xyz - n.c.xyz;
      float3 pb = p - n.b.xyz, pc = p - n.c.xyz, nor = cross(ba, ac);
      if (sign(dot(cross(ba, nor), v)) + sign(dot(cross(cb, nor), pb)) + sign(dot(cross(ac, nor), pc)) < 2) {
        float3 a = ba * clamp(dot(ba, v) / dot(ba, ba), 0.f, 1.f) - v;
        float3 b = cb * clamp(dot(cb, pb) / dot(cb, cb), 0.f, 1.f) - pb;
        float3 c = ac * clamp(dot(ac, pc) / dot(ac, ac), 0.f, 1.f) - pc;
        d = sqrt(min(dot(a, a), min(dot(b, b), dot(c, c))));
      } else
        d = abs(dot(nor, v)) / length(nor);
    }
    if (choose_material) choices[top] = fallback;
    values[top++] = d;
  }
  if (choose_material) selected = choices[0];
  return values[0];
}
float csg_distance(device const CsgInstruction *code, uint count, float3 p) {
  uint seed = 0, selected = 0;
  return csg_evaluate<false>(code, count, p, 0, seed, selected);
}
uint csg_material(device const CsgInstruction *code, uint count, float3 p,
                  uint fallback, uint seed) {
  uint selected = fallback;
  csg_evaluate<true>(code, count, p, fallback, seed, selected);
  return selected;
}
)MSL"
