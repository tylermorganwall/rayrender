R"MSL(
// PBRT-style ray-facing cubic curves, matching hitables/curve.cpp. The CPU's
// depth-limited recursion becomes a depth-first stack of at most ten deferred
// siblings. This is a geometric refinement bound, not a ray termination cap.
struct Curve {
  float4 cp[4], normal[2], width, fallback;
  uint type, padding[3];
};
static_assert(sizeof(Curve) == 144, "wavefront curve ABI mismatch");
Curve curve_data(device const Quadric *quadrics, Quadric q) {
  return *reinterpret_cast<device const Curve *>(
      reinterpret_cast<device const uchar *>(quadrics) + q.curve_offset);
}
float3 curve_eval(thread const float3 *cp, float u, thread float3 &derivative) {
  float3 a = mix(cp[0], cp[1], u), b = mix(cp[1], cp[2], u),
         c = mix(cp[2], cp[3], u);
  float3 d = mix(a, b, u), e = mix(b, c, u);
  derivative = 3 * (e - d);
  if (dot(derivative, derivative) == 0)
    derivative = cp[3] - cp[0];
  return mix(d, e, u);
}
float3 curve_ribbon_normal(Curve c, float u) {
  float a = sin((1 - u) * c.normal[0].w) * c.normal[1].w;
  float b = sin(u * c.normal[0].w) * c.normal[1].w;
  return c.normal[0].w > 0 && isfinite(a) && isfinite(b)
             ? a * c.normal[0].xyz + b * c.normal[1].xyz
             : c.normal[0].xyz;
}
struct CurveSegment {
  float3 cp[4];
  float u0, u1;
  uint depth;
};
bool curve_segment_overlaps(CurveSegment s, Curve c, float maximum_z) {
  float radius = .5f * max(mix(c.width.x, c.width.y, s.u0),
                           mix(c.width.x, c.width.y, s.u1));
  float3 lo = min(min(s.cp[0], s.cp[1]), min(s.cp[2], s.cp[3])) - radius;
  float3 hi = max(max(s.cp[0], s.cp[1]), max(s.cp[2], s.cp[3])) + radius;
  return lo.y <= 0 && hi.y >= 0 && lo.x <= 0 && hi.x >= 0 && hi.z >= 0 &&
         lo.z <= maximum_z;
}
// Returns global (u,v) and world ray parameter t. No normalization of
// transformed ray directions: that would corrupt clipping and nonuniform
// placement distances.
float3 curve_intersection(Quadric q, device const Quadric *quadrics,
                          float3 origin, float3 direction, float minimum,
                          float maximum) {
  Curve c = curve_data(quadrics, q);
  float3 o = quadric_local(q, origin), d = quadric_direction(q, direction);
  float speed = length(d);
  float3 miss = float3(0, 0, INFINITY);
  if (!(speed > 0) || !(max(c.width.x, c.width.y) > 0) || minimum > maximum)
    return miss;
  float3 forward = d / speed;
  float3 up = cross(forward, c.cp[3].xyz - c.cp[0].xyz);
  if (dot(up, up) == 0)
    up =
        normalize(cross(forward, abs(forward.x) > .9999999f ? float3(0, 1, 0)
                                                            : float3(1, 0, 0)));
  else
    up = normalize(up);
  float3 right = normalize(cross(up, forward));
  up = cross(forward, right);
  CurveSegment current;
  for (uint i = 0; i < 4; ++i) {
    float3 relative = c.cp[i].xyz - o;
    current.cp[i] =
        float3(dot(right, relative), dot(up, relative), dot(forward, relative));
  }
  current.u0 = c.width.z;
  current.u1 = c.width.w;
  float curvature = 0;
  for (uint i = 0; i < 2; ++i) {
    float3 second =
        abs(current.cp[i] - 2 * current.cp[i + 1] + current.cp[i + 2]);
    curvature = max(curvature, max(second.x, max(second.y, second.z)));
  }
  float refinement =
      1.41421356237f * 6 * curvature / (8 * .05f * max(c.width.x, c.width.y));
  uint bits = as_type<uint>(refinement);
  int rounded_log = refinement < 1
                        ? 0
                        : int(bits >> 23) - 127 + ((bits & (1u << 22)) ? 1 : 0);
  current.depth = uint(clamp(rounded_log / 2, 0, 10));
  CurveSegment deferred[10];
  uint count = 0;
  float3 closest = miss;
  while (true) {
    if (curve_segment_overlaps(current, c, speed * maximum)) {
      if (current.depth) {
        CurveSegment left, right_segment;
        float3 a = .5f * (current.cp[0] + current.cp[1]);
        float3 b = .25f * (current.cp[0] + 2 * current.cp[1] + current.cp[2]);
        float3 middle = .125f * (current.cp[0] + 3 * current.cp[1] +
                                 3 * current.cp[2] + current.cp[3]);
        float3 e = .25f * (current.cp[1] + 2 * current.cp[2] + current.cp[3]);
        float3 f = .5f * (current.cp[2] + current.cp[3]);
        left.cp[0] = current.cp[0];
        left.cp[1] = a;
        left.cp[2] = b;
        left.cp[3] = middle;
        right_segment.cp[0] = middle;
        right_segment.cp[1] = e;
        right_segment.cp[2] = f;
        right_segment.cp[3] = current.cp[3];
        left.u0 = current.u0;
        left.u1 = .5f * (current.u0 + current.u1);
        right_segment.u0 = left.u1;
        right_segment.u1 = current.u1;
        left.depth = right_segment.depth = current.depth - 1;
        deferred[count++] = right_segment;
        current = left;
        continue;
      }
      float start_edge =
          dot(current.cp[0].xy - current.cp[1].xy, current.cp[0].xy);
      float end_edge =
          dot(current.cp[3].xy - current.cp[2].xy, current.cp[3].xy);
      float2 delta = current.cp[3].xy - current.cp[0].xy;
      float denominator = dot(delta, delta);
      if (start_edge >= 0 && end_edge >= 0 && denominator > 0) {
        float w = dot(-current.cp[0].xy, delta) / denominator;
        float u = clamp(mix(current.u0, current.u1, w), current.u0, current.u1);
        float width = mix(c.width.x, c.width.y, u);
        if (c.type == 3)
          width *= abs(dot(curve_ribbon_normal(c, u), d)) / speed;
        float3 derivative;
        float3 point = curve_eval(current.cp, clamp(w, 0.f, 1.f), derivative);
        float radial2 = dot(point.xy, point.xy);
        if (width > 0 && radial2 <= .25f * width * width &&
            point.z >= speed * minimum && point.z <= speed * maximum) {
          float edge = derivative.x * -point.y + point.x * derivative.y;
          float v = .5f + (edge > 0 ? 1 : -1) * sqrt(radial2) / width;
          maximum = point.z / speed;
          closest = float3(u, v, maximum);
        }
      }
    }
    if (!count)
      break;
    current = deferred[--count];
  }
  return closest;
}
struct CurveFrame {
  float3 normal, du, dv;
  float width;
};
CurveFrame curve_frame(Curve c, float2 uv, float3 direction) {
  float3 cp[4], tangent;
  for (uint i = 0; i < 4; ++i)
    cp[i] = c.cp[i].xyz;
  curve_eval(cp, (uv.x - c.width.z) / (c.width.w - c.width.z), tangent);
  if (dot(tangent, tangent) == 0)
    tangent = c.fallback.xyz;
  tangent = normalize(tangent);
  float3 d = normalize(direction), n, dv;
  float width = mix(c.width.x, c.width.y, uv.x);
  if (c.type == 3) {
    n = curve_ribbon_normal(c, uv.x);
    width *= abs(dot(n, d));
    dv = normalize(cross(n, tangent)) * width;
    if (dot(n, d) > 0)
      n = -n;
  } else if (c.type == 2) {
    float3 initial = normalize(cross(d, tangent));
    float theta = (uv.y - .5f) * pi;
    dv = width * (cos(theta) * initial + sin(theta) * cross(tangent, initial));
    n = normalize(-cross(tangent, dv));
  } else {
    dv = normalize(cross(tangent, -d));
    n = -d;
  }
  return {n, tangent, dv, width};
}
)MSL"
