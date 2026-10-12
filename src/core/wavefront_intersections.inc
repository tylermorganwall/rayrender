#include "wavefront_csg.inc"
R"MSL(

// IDs below 2^31 are triangle indices; the high bit selects an analytic
// primitive. UINT_MAX and UINT_MAX-1 remain volume/diffusion queue sentinels.
constant uint quadric_bit = 0x80000000u;
struct Quadric {
  float4 inverse[3], forward[3], measure;
  uint material, light, boundary, surface;
  uint csg_offset, csg_count, curve_offset, reserved;
  float4 csg_lo, csg_hi;
};
static_assert(sizeof(Quadric) == 176, "wavefront quadric ABI mismatch");
struct MeshInstance {
  float4 inverse[3], forward[3];
  uint triangle_start, triangle_count, surface_base, boundary;
  uint light_base, mesh;
  float orientation;
  uint padding;
};
struct GeometryStorage {
  uint quadric_offset, instance_offset, medium_offset, field_offset;
  uint atmosphere_offset, atmosphere_data_offset, camera_offset, motion_offset;
  uint light_offset, padding[3];
};
static_assert(sizeof(MeshInstance) == 128, "wavefront instance ABI mismatch");
bool is_quadric(uint primitive) { return scene_has_quadrics && (primitive & quadric_bit); }
float3 quadric_point(Quadric q, float3 p) {
  return float3(dot(q.forward[0], float4(p, 1)), dot(q.forward[1], float4(p, 1)),
                dot(q.forward[2], float4(p, 1)));
}
float3 quadric_local(Quadric q, float3 p) {
  return float3(dot(q.inverse[0], float4(p, 1)), dot(q.inverse[1], float4(p, 1)),
                dot(q.inverse[2], float4(p, 1)));
}
float3 quadric_direction(Quadric q, float3 d) {
  return float3(dot(q.inverse[0].xyz, d), dot(q.inverse[1].xyz, d), dot(q.inverse[2].xyz, d));
}
float3 quadric_normal(Quadric q, float3 p) {
  return q.inverse[0].xyz * p.x + q.inverse[1].xyz * p.y + q.inverse[2].xyz * p.z;
}
)MSL"
#include "wavefront_curves.inc"
R"MSL(
device const CsgInstruction *csg_program(device const Quadric *quadrics, Quadric q) {
  return reinterpret_cast<device const CsgInstruction *>(
      reinterpret_cast<device const uchar *>(quadrics) + q.csg_offset);
}
float csg_intersection(Quadric q, device const Quadric *quadrics, float3 origin, float3 direction,
                       float minimum, float maximum, device atomic_uint *failed = nullptr) {
  float3 o = quadric_local(q, origin), d = quadric_direction(q, direction);
  float speed = length(d), lower = minimum, upper = maximum;
  if (!(speed > 0))
    return INFINITY;
  for (uint a = 0; a < 3; ++a) {
    if (d[a] == 0) {
      if (o[a] < q.csg_lo[a] || o[a] > q.csg_hi[a])
        return INFINITY;
    } else {
      float x = (q.csg_lo[a] - o[a]) / d[a], y = (q.csg_hi[a] - o[a]) / d[a];
      lower = max(lower, min(x, y));
      upper = min(upper, max(x, y));
    }
  }
  if (lower > upper)
    return INFINITY;
  auto code = csg_program(quadrics, q);
  // March relative to the box entry, not the distant world-ray origin. At a
  // large t, a perfectly resolvable distance step can be smaller than one ULP
  // of t and round away. Rebasing preserves those local steps without relaxing
  // the surface tolerance or skipping thin geometry.
  o = fma(float3(lower), d, o);
  upper -= lower;
  float epsilon = q.csg_lo.w, t = 0;
  // A convex cutter in a hard difference excludes its entire interior ray
  // interval. Skip that known-empty interval instead of taking thousands of
  // tiny distance steps parallel to a box face. Shrink by two tolerances so
  // both cut surfaces retain their full acceptance band.
  float cut_lo = INFINITY, cut_hi = -INFINITY;
  if (q.csg_count >= 3 && code[q.csg_count - 1].operation == 31) {
    CsgInstruction cutter = code[q.csg_count - 2];
    if (cutter.operation == 4 && cutter.a.w == 0) {
      cut_lo = -INFINITY;
      cut_hi = INFINITY;
      float3 halfwidth = cutter.b.xyz * .5f - 2 * epsilon;
      for (uint a = 0; a < 3; ++a) {
        float x = o[a] - cutter.a[a];
        if (halfwidth[a] <= 0 || (d[a] == 0 && abs(x) >= halfwidth[a])) {
          cut_lo = INFINITY;
          cut_hi = -INFINITY;
          break;
        }
        if (d[a] != 0) {
          float near = (-halfwidth[a] - x) / d[a], far = (halfwidth[a] - x) / d[a];
          cut_lo = max(cut_lo, min(near, far));
          cut_hi = min(cut_hi, max(near, far));
        }
      }
    } else if (cutter.operation == 2 && cutter.a.w > 2 * epsilon) {
      float3 v = o - cutter.a.xyz;
      float a = dot(d, d), b = dot(v, d), r = cutter.a.w - 2 * epsilon;
      float3 closest = v - b / a * d;
      float disc = r * r - dot(closest, closest);
      if (disc > 0) {
        float root = sqrt(disc / a);
        cut_lo = -b / a - root;
        cut_hi = -b / a + root;
      }
    }
  }
  // Interior rays use |distance|. Suppress only the initial self-contact,
  // never a new entry reached by advancing to the object's bounding box.
  bool departing =
      lower == minimum && abs(csg_distance(code, q.csg_count, fma(float3(t), d, o))) <= epsilon;
  for (uint step = 0; step < 8192 && t <= upper; ++step) {
    if (t >= cut_lo && t < cut_hi) {
      t = cut_hi;
      if (t > upper)
        return INFINITY;
    }
    float distance = abs(csg_distance(code, q.csg_count, fma(float3(t), d, o)));
    if (!isfinite(distance)) {
      if (failed)
        atomic_fetch_add_explicit(failed, 1, memory_order_relaxed);
      return INFINITY;
    }
    if (distance > epsilon)
      departing = false;
    if (distance <= epsilon && !departing)
      return lower + t;
    // Leave the entire convergence band before accepting another hit. A
    // fixed one-step nudge counts the same grazing surface repeatedly in
    // diffusion chords and biases the reservoir toward those duplicates.
    if (departing)
      distance = max(distance, 4 * epsilon);
    float next = t + distance / speed;
    if (!(next > t)) {
      if (failed)
        atomic_fetch_add_explicit(failed, 1, memory_order_relaxed);
      return INFINITY;
    }
    t = next;
  }
  if (t <= upper && failed)
    atomic_fetch_add_explicit(failed, 1, memory_order_relaxed);
  return INFINITY;
}
float2 quadric_roots(Quadric q, float3 origin, float3 direction) {
  float3 o = quadric_local(q, origin), d = quadric_direction(q, direction);
  float a = dot(d, d), b = dot(o, d);
  // Closest approach avoids subtracting two large squared terms in the
  // discriminant. Vieta's product recovers the near root without subtracting
  // nearly equal far distances. The direction keeps world-distance units.
  float3 closest = fma(float3(-b / a), d, o);
  float discriminant = 1 - dot(closest, closest);
  if (!(a > 0) || discriminant < 0)
    return float2(INFINITY);
  float root = sqrt(a * discriminant), stable = -b - copysign(root, b);
  if (stable == 0)
    return float2(-b / a);
  float length_o = length(o);
  float c = (length_o - 1) * (length_o + 1);
  float2 roots = float2(stable / a, c / stable);
  return float2(min(roots.x, roots.y), max(roots.x, roots.y));
}
struct QuadricIntersection {
  bool accept [[accept_intersection]];
  float distance [[distance]];
};
[[intersection(bounding_box, triangle_data, instancing)]]
QuadricIntersection
intersect_quadric(float3 origin [[origin]], float3 direction [[direction]],
                  float minimum [[min_distance]], float maximum [[max_distance]],
                  uint primitive [[primitive_id]], ray_data uint &body [[payload]],
                  device const Quadric *quadrics [[buffer(1)]],
                  device atomic_uint *failed [[buffer(3)]]) {
  Quadric q = quadrics[primitive];
  if (body && q.boundary != body)
    return {false, 0};
  if ((texture_features & 512) && q.curve_offset) {
    float distance = curve_intersection(q, quadrics, origin, direction, minimum, maximum).z;
    return {isfinite(distance), distance};
  }
  if ((texture_features & 64) && q.csg_count) {
    float distance = csg_intersection(q, quadrics, origin, direction, minimum, maximum, failed);
    return {isfinite(distance), distance};
  }
  float2 roots = quadric_roots(q, origin, direction);
  float distance = roots.x >= minimum ? roots.x : roots.y;
  return {isfinite(distance) && distance >= minimum && distance <= maximum, distance};
}
[[intersection(triangle, triangle_data)]]
bool containment_face(uint primitive [[primitive_id]], ray_data uint &body [[payload]],
                      device const Triangle *triangles [[buffer(0)]]) {
  return !body || triangles[primitive].boundary == body;
}
[[intersection(triangle, triangle_data, instancing)]]
bool containment_face_instanced(uint primitive [[primitive_id]], uint instance [[instance_id]],
                                ray_data uint &body [[payload]],
                                device const Triangle *triangles [[buffer(0)]],
                                device const MeshInstance *instances [[buffer(2)]]) {
  if (!body)
    return true;
  MeshInstance placement = instances[instance];
  return (placement.triangle_count ? placement.boundary : triangles[primitive].boundary) == body;
}
// Metal's primitive structures are homogeneous. Mixed scenes traverse one
// top-level structure over triangle and bounding-box leaves. The original
// direct triangle path is retained for scenes without analytic geometry.
struct GeometryScene {
  primitive_acceleration_structure triangles;
  intersection_function_table<triangle_data> triangle_functions;
  instance_acceleration_structure mixed;
  intersection_function_table<triangle_data, instancing> mixed_functions;
  device const MeshInstance *instances;
  device const float4 *motion_inverse;
  device const Quadric *reference_quadrics;
};
struct GeometryIntersection {
  intersection_type type;
  uint primitive_id, instance;
  float distance;
  float3 coordinates; // triangle barycentric pair, or unit-sphere point
};
GeometryIntersection trace_geometry(ray r, GeometryScene scene, device const Quadric *quadrics,
                                    uint body = 0, bool any_hit = false) {
  if (texture_features & 32) {
    // Keep the immutable reference BVH tight. An affine ray transform retains
    // its world-distance parameter provided the direction is not normalized.
    float3 origin = r.origin, direction = r.direction;
    for (uint a = 0; a < 3; ++a) {
      r.origin[a] = dot(scene.motion_inverse[a], float4(origin, 1));
      r.direction[a] = dot(scene.motion_inverse[a].xyz, direction);
    }
    quadrics = scene.reference_quadrics;
  }
  GeometryIntersection result;
  if (scene_has_quadrics || scene_has_instances) {
    intersector<triangle_data, instancing> query;
    query.force_opacity(body ? forced_opacity::non_opaque : forced_opacity::opaque);
    query.accept_any_intersection(any_hit);
    auto hit = query.intersect(r, scene.mixed, 0xff, scene.mixed_functions, body);
    result = {hit.type, hit.primitive_id, 0, hit.distance, float3(0)};
    if (hit.type == intersection_type::triangle) {
      MeshInstance placement = scene.instances[hit.instance_id];
      if (placement.triangle_count) {
        result.primitive_id += placement.triangle_start;
        result.instance = hit.instance_id + 1;
      }
      result.coordinates.xy = hit.triangle_barycentric_coord;
    } else if (hit.type == intersection_type::bounding_box) {
      result.primitive_id |= quadric_bit;
      Quadric q = quadrics[hit.primitive_id];
      if ((texture_features & 512) && q.curve_offset) {
        // Intersection callbacks return only distance. Recover attributes for
        // the winning curve, never for visibility-only any-hit queries.
        if (!any_hit)
          result.coordinates = curve_intersection(q, quadrics, r.origin, r.direction,
                                                   r.min_distance, r.max_distance);
      } else {
        result.coordinates =
          fma(float3(hit.distance), quadric_direction(q, r.direction), quadric_local(q, r.origin));
        if (!((texture_features & 64) && q.csg_count))
          result.coordinates = normalize(result.coordinates);
      }
    }
  } else {
    intersector<triangle_data> query;
    query.assume_geometry_type(geometry_type::triangle);
    query.force_opacity(body ? forced_opacity::non_opaque : forced_opacity::opaque);
    query.accept_any_intersection(any_hit);
    intersection_result<triangle_data> hit;
    if (body)
      hit = query.intersect(r, scene.triangles, scene.triangle_functions, body);
    else
      hit = query.intersect(r, scene.triangles);
    result = {hit.type, hit.primitive_id, 0, hit.distance,
              float3(hit.triangle_barycentric_coord, 0)};
  }
  return result;
}
// A spectral shadow assigns one ray to all 32 lanes. Traverse once, then share
// the hit before their cooperative spectral evaluation. Every lane must enter
// this helper together when cooperative is true; ordinary shadows stay scalar.
GeometryIntersection trace_spectral_shadow(ray r, GeometryScene scene,
                                          device const Quadric *quadrics,
                                          bool cooperative, uint lane, bool any_hit = false) {
  GeometryIntersection hit = {};
  if (!cooperative || lane == 0)
    hit = trace_geometry(r, scene, quadrics, 0, any_hit);
  if (cooperative) {
    hit.type = intersection_type(simd_broadcast(uint(hit.type), 0));
    hit.primitive_id = simd_broadcast(hit.primitive_id, 0);
    hit.instance = simd_broadcast(hit.instance, 0);
    hit.distance = simd_broadcast(hit.distance, 0);
    hit.coordinates = simd_broadcast(hit.coordinates, 0);
  }
  return hit;
}
uint geometry_material(uint primitive, device const Triangle *triangles,
                       device const Quadric *quadrics) {
  return is_quadric(primitive) ? quadrics[primitive & ~quadric_bit].material
                               : triangles[primitive].material;
}
struct GeometrySurface {
  Triangle face; // Native triangle/curve attributes, or analytic material/chart metadata.
  float3 position, error, normal, tangent, coordinates;
};
GeometrySurface geometry_surface(uint primitive, float3 coordinates,
                                 device const Triangle *triangles, device const Quadric *quadrics,
                                 device const Surface *surfaces, uint instance,
                                 device const MeshInstance *instances,
                                 float3 ray_origin, float3 ray_direction) {
  GeometrySurface out = {};
  out.coordinates = coordinates;
  if (is_quadric(primitive)) {
    Quadric q = quadrics[primitive & ~quadric_bit];
    if ((texture_features & 512) && q.curve_offset) {
      CurveFrame frame = curve_frame(curve_data(quadrics, q), coordinates.xy,
                                      quadric_direction(q, ray_direction));
      out.position = fma(float3(coordinates.z), ray_direction, ray_origin);
      out.normal = q.measure.y * normalize(quadric_normal(q, frame.normal));
      float3 du, dv;
      float3 local_position = quadric_local(q, out.position);
      for (uint a = 0; a < 3; ++a) {
        du[a] = dot(q.forward[a].xyz, frame.du);
        dv[a] = dot(q.forward[a].xyz, frame.dv);
        out.error[a] = frame.width * dot(abs(q.forward[a].xyz), float3(1)) +
          8 * FLT_EPSILON * (dot(abs(q.forward[a].xyz), abs(local_position)) + abs(q.forward[a].w));
        out.face.n[a] = float4(out.normal, 0);
      }
      out.tangent = du;
      out.face.p[0] = float4(out.position, 0);
      out.face.p[1] = float4(out.position + du, 0);
      out.face.p[2] = float4(out.position + dv, 0);
      out.face.material = q.material;
      out.face.light = q.light;
      out.face.boundary = 0;
      out.face.flags = 1 | 32;
      out.face.geometric = float4(out.normal, 0);
      out.face.uv01 = coordinates.xyxy;
      out.face.uv2 = float4(coordinates.xy, 0, as_type<float>(q.surface));
      return out;
    }
    if ((texture_features & 64) && q.csg_count) {
      auto code = csg_program(quadrics, q);
      float e = q.csg_lo.w, gradient[3];
      for (uint a = 0; a < 3; ++a) {
        float3 delta = 0;
        delta[a] = e;
        gradient[a] = csg_distance(code, q.csg_count, coordinates + delta) -
                      csg_distance(code, q.csg_count, coordinates - delta);
      }
      float3 n = float3(gradient[0], gradient[1], gradient[2]);
      if (dot(n, n) == 0)
        n = float3(0, 1, 0);
      n = normalize(quadric_normal(q, n));
      out.position = quadric_point(q, coordinates);
      out.normal = q.measure.y * n;
      out.face.material = q.material;
      out.face.light = q.light;
      out.face.boundary = q.boundary;
      out.face.geometric = float4(q.measure.w * n, 0);
      out.face.flags = 1;
      out.face.uv01 = float4(.5f);
      out.face.uv2 = float4(.5f, .5f, 0, as_type<float>(q.surface));
      out.tangent = normalize(cross(abs(n.x) > .9f ? float3(0, 1, 0) : float3(1, 0, 0), n));
      for (uint k = 0; k < 3; ++k) {
        out.face.p[k] = float4(out.position, 0);
        out.face.n[k] = float4(out.normal, 0);
        out.error[k] =
            2 * e * length(q.forward[k].xyz) +
            8 * FLT_EPSILON * (dot(abs(q.forward[k].xyz), abs(coordinates)) + abs(q.forward[k].w));
      }
      return out;
    }
    float3 unit = normalize(coordinates);
    out.position = quadric_point(q, unit);
    out.normal = q.measure.y * normalize(quadric_normal(q, unit));
    out.face.material = q.material;
    out.face.light = q.light;
    out.face.boundary = q.boundary;
    out.face.flags = 1;
    out.face.geometric = float4(q.measure.w * normalize(quadric_normal(q, unit)), 0);
    Surface s = surfaces[q.surface];
    float3 chart = s.axes.w == 1 ? unit : normalize(unit / s.axes.xyz);
    float2 uv = float2(1 - (atan2(chart.z, chart.x) + pi) / (2 * pi),
                       .5f + asin(clamp(chart.y, -1.f, 1.f)) / pi);
    out.face.uv01 = float4(uv, uv);
    out.face.uv2 = float4(uv, 0, as_type<float>(q.surface));
    // Analytic charts replace the triangle Jacobian in surface_derivatives().
    // Populate equal UVs/normals for the shared texture/material code only.
    for (uint k = 0; k < 3; ++k) {
      out.face.n[k] = float4(out.normal, 0);
      out.face.p[k] = float4(out.position, 0);
      out.error[k] =
          8 * FLT_EPSILON * (dot(abs(q.forward[k].xyz), abs(unit)) + abs(q.forward[k].w));
    }
    // Ellipsoid UVs parameterize the object normal. Differentiate that chart,
    // rather than a latitude on the unit sphere, for anisotropic materials.
    float3 squared = s.axes.xyz * s.axes.xyz, an = squared * chart;
    float denominator = sqrt(dot(chart, an));
    float3 nu = 2 * pi * float3(chart.z, 0, -chart.x);
    float3 tangent =
        (squared * nu - an * (dot(an, nu) / (denominator * denominator))) / denominator;
    if (dot(tangent, tangent) == 0)
      tangent = float3(1, 0, 0);
    out.tangent = float3(dot(s.forward[0].xyz, tangent), dot(s.forward[1].xyz, tangent),
                         dot(s.forward[2].xyz, tangent));
    return out;
  }
  out.face = triangles[primitive];
  if (scene_has_instances && instance) {
    MeshInstance placement = instances[instance - 1];
    Triangle local = out.face;
    for (uint k = 0; k < 3; ++k) {
      float3 p = local.p[k].xyz, n = local.n[k].xyz;
      out.face.p[k].xyz =
          float3(dot(placement.forward[0], float4(p, 1)), dot(placement.forward[1], float4(p, 1)),
                 dot(placement.forward[2], float4(p, 1)));
      out.face.n[k].xyz = placement.inverse[0].xyz * n.x + placement.inverse[1].xyz * n.y +
                          placement.inverse[2].xyz * n.z;
    }
    float3 g = local.geometric.xyz;
    out.face.geometric.xyz = placement.orientation * normalize(placement.inverse[0].xyz * g.x +
                                                               placement.inverse[1].xyz * g.y +
                                                               placement.inverse[2].xyz * g.z);
    out.face.boundary = placement.boundary;
    out.face.uv2.w = as_type<float>(placement.surface_base + as_type<uint>(local.uv2.w));
    if (local.light != UINT_MAX)
      out.face.light = placement.light_base + local.light;
  }

  float3 bary = float3(1 - coordinates.x - coordinates.y, coordinates.xy);
  surface_point(out.face, coordinates.xy, out.position, out.error);
  if (scene_has_instances && instance) {
    MeshInstance placement = instances[instance - 1];
    Triangle local = triangles[primitive];
    float3 magnitude =
        bary.x * abs(local.p[0].xyz) + bary.y * abs(local.p[1].xyz) + bary.z * abs(local.p[2].xyz);
    for (uint k = 0; k < 3; ++k)
      out.error[k] += 8 * FLT_EPSILON *
                      (dot(abs(placement.forward[k].xyz), magnitude) + abs(placement.forward[k].w));
  }
  float3 normal =
      bary.x * out.face.n[0].xyz + bary.y * out.face.n[1].xyz + bary.z * out.face.n[2].xyz;
  out.normal = dot(normal, normal) > 0 ? normalize(normal) : out.face.geometric.xyz;
  float2 a = out.face.uv01.zw - out.face.uv01.xy, b = out.face.uv2.xy - out.face.uv01.xy;
  float det = a.x * b.y - a.y * b.x;
  float3 g = out.face.geometric.xyz;
  float3 bitangent = normalize(cross(g, abs(g.x) > .9999999f ? float3(0, 1, 0) : float3(1, 0, 0)));
  out.tangent = abs(det) > 1e-12f ? (b.y * (out.face.p[1].xyz - out.face.p[0].xyz) -
                                     a.y * (out.face.p[2].xyz - out.face.p[0].xyz)) /
                                        det
                                  : cross(bitangent, g);
  return out;
}
float geometry_light_pdf(uint primitive, GeometrySurface surface, float3 reference,
                         device const Quadric *quadrics) {
  float3 delta = reference - surface.position;
  float distance2 = dot(delta, delta),
        cosine = abs(dot(surface.face.geometric.xyz, normalize(delta)));
  if (is_quadric(primitive)) {
    Quadric q = quadrics[primitive & ~quadric_bit];
    if (((texture_features & 64) && q.csg_count) || ((texture_features & 512) && q.curve_offset)) {
      float3 center = quadric_point(q, .5f * (q.csg_lo.xyz + q.csg_hi.xyz));
      float3 halfwidth = .5f * (q.csg_hi.xyz - q.csg_lo.xyz);
      float3 extent =
          float3(dot(abs(q.forward[0].xyz), halfwidth), dot(abs(q.forward[1].xyz), halfwidth),
                 dot(abs(q.forward[2].xyz), halfwidth));
      float ratio = dot(extent, extent) / distance_squared(reference, center);
      return ratio >= 1 ? 1 / (4 * pi) : 1 / (2 * pi * (ratio / (1 + sqrt(1 - ratio))));
    }
    float3 center = float3(q.forward[0].w, q.forward[1].w, q.forward[2].w);
    float center2 = distance_squared(reference, center), radius2 = q.measure.z * q.measure.z;
    if (radius2 > 0 && center2 > radius2) {
      float ratio = min(1.f, radius2 / center2);
      float delta_cosine = ratio / (1 + sqrt(max(0.f, 1 - ratio)));
      return 1 / (2 * pi * delta_cosine);
    }
    float3 unit = normalize(quadric_local(q, surface.position));
    float jacobian = q.measure.x * length(quadric_normal(q, unit));
    return distance2 / (4 * pi * jacobian * cosine);
  }
  Triangle t = surface.face;
  float area = .5f * length(cross(t.p[1].xyz - t.p[0].xyz, t.p[2].xyz - t.p[0].xyz));
  return distance2 / (area * cosine);
}
GeometrySurface
sample_geometry_light(uint primitive, float2 u, float3 reference, device const Triangle *triangles,
                      device const Quadric *quadrics, device const Surface *surfaces, uint instance,
                      device const MeshInstance *instances, device atomic_uint *failed) {
  float3 coordinates;
  if (is_quadric(primitive)) {
    Quadric q = quadrics[primitive & ~quadric_bit];
    if (((texture_features & 64) && q.csg_count) || ((texture_features & 512) && q.curve_offset)) {
      float3 center = quadric_point(q, .5f * (q.csg_lo.xyz + q.csg_hi.xyz));
      float3 halfwidth = .5f * (q.csg_hi.xyz - q.csg_lo.xyz);
      float3 extent =
          float3(dot(abs(q.forward[0].xyz), halfwidth), dot(abs(q.forward[1].xyz), halfwidth),
                 dot(abs(q.forward[2].xyz), halfwidth));
      float ratio = dot(extent, extent) / distance_squared(reference, center);
      float dc = u.x * (ratio >= 1 ? 2.f : ratio / (1 + sqrt(1 - ratio)));
      float3 forward =
          distance_squared(reference, center) > 0 ? normalize(center - reference) : float3(0, 0, 1);
      float3 x =
          normalize(cross(abs(forward.x) > .9f ? float3(0, 1, 0) : float3(1, 0, 0), forward));
      float3 direction =
          (1 - dc) * forward +
          sqrt(dc * (2 - dc)) * (cos(2 * pi * u.y) * x + sin(2 * pi * u.y) * cross(forward, x));
      float3 curve_hit = float3(0, 0, INFINITY);
      if ((texture_features & 512) && q.curve_offset)
        curve_hit = curve_intersection(q, quadrics, reference, direction, 0, INFINITY);
      float distance = q.curve_offset ? curve_hit.z :
          csg_intersection(q, quadrics, reference, direction, 0, INFINITY, failed);
      if (!isfinite(distance)) {
        // Rejection is a zero-contribution sample. Do not resample: that would
        // condition on the unknown projected area and invalidate the PDF.
        GeometrySurface miss = {};
        miss.face.material = q.material;
        miss.position = reference + direction;
        return miss;
      }
      coordinates = q.curve_offset ? curve_hit :
          fma(float3(distance), quadric_direction(q, direction), quadric_local(q, reference));
      return geometry_surface(primitive, coordinates, triangles, quadrics, surfaces, instance,
                              instances, reference, direction);
    }
    float3 center = float3(q.forward[0].w, q.forward[1].w, q.forward[2].w);
    float center2 = distance_squared(reference, center), radius2 = q.measure.z * q.measure.z;
    if (radius2 > 0 && center2 > radius2) {
      float ratio = min(1.f, radius2 / center2);
      float dc = u.x * ratio / (1 + sqrt(max(0.f, 1 - ratio)));
      float3 forward = normalize(center - reference);
      float3 x =
          normalize(cross(abs(forward.x) > .9f ? float3(0, 1, 0) : float3(1, 0, 0), forward));
      float3 direction =
          (1 - dc) * forward +
          sqrt(dc * (2 - dc)) * (cos(2 * pi * u.y) * x + sin(2 * pi * u.y) * cross(forward, x));
      float distance = quadric_roots(q, reference, direction).x;
      coordinates = normalize(
          fma(float3(distance), quadric_direction(q, direction), quadric_local(q, reference)));
    } else {
      float y = 1 - 2 * u.x, radius = sqrt(max(0.f, 1 - y * y)), phi = 2 * pi * u.y;
      coordinates = float3(radius * cos(phi), y, radius * sin(phi));
    }
  } else {
    float root = sqrt(u.x);
    coordinates = float3(root * (1 - u.y), root * u.y, 0);
  }
  return geometry_surface(primitive, coordinates, triangles, quadrics, surfaces, instance,
                          instances, float3(0), float3(0));
}
)MSL"
