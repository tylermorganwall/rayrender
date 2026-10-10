R"MSL(

  // IDs below 2^31 are triangle indices; the high bit selects a unit-sphere
  // primitive. UINT_MAX and UINT_MAX-1 remain volume/diffusion queue sentinels.
  constant uint quadric_bit = 0x80000000u;
  struct Quadric {
    float4 inverse[3], forward[3], measure;
    uint material, light, boundary, surface;
  };
  static_assert(sizeof(Quadric) == 128, "wavefront quadric ABI mismatch");
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
                    device const Quadric *quadrics [[buffer(1)]]) {
    Quadric q = quadrics[primitive];
    if (body && q.boundary != body)
      return {false, 0};
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
  bool containment_face_instanced(uint primitive [[primitive_id]], ray_data uint &body [[payload]],
                                  device const Triangle *triangles [[buffer(0)]]) {
    return !body || triangles[primitive].boundary == body;
  }
  // Metal's primitive structures are homogeneous. Mixed scenes traverse one
  // top-level structure over triangle and bounding-box leaves. The original
  // direct triangle path is retained for scenes without analytic geometry.
  struct GeometryScene {
    primitive_acceleration_structure triangles;
    intersection_function_table<triangle_data> triangle_functions;
    instance_acceleration_structure mixed;
    intersection_function_table<triangle_data, instancing> mixed_functions;
  };
  struct GeometryIntersection {
    intersection_type type;
    uint primitive_id;
    float distance;
    float3 coordinates; // triangle barycentric pair, or unit-sphere point
  };
  GeometryIntersection trace_geometry(ray r, GeometryScene scene, device const Quadric *quadrics,
                                      uint body = 0, bool any_hit = false) {
    GeometryIntersection result;
    if (scene_has_quadrics) {
      intersector<triangle_data, instancing> query;
      query.force_opacity(body ? forced_opacity::non_opaque : forced_opacity::opaque);
      query.accept_any_intersection(any_hit);
      auto hit = query.intersect(r, scene.mixed, 0xff, scene.mixed_functions, body);
      result = {hit.type, hit.primitive_id, hit.distance, float3(0)};
      if (hit.type == intersection_type::triangle)
        result.coordinates.xy = hit.triangle_barycentric_coord;
      else if (hit.type == intersection_type::bounding_box) {
        result.primitive_id |= quadric_bit;
        Quadric q = quadrics[hit.primitive_id];
        result.coordinates = normalize(fma(float3(hit.distance), quadric_direction(q, r.direction),
                                           quadric_local(q, r.origin)));
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
      result = {hit.type, hit.primitive_id, hit.distance,
                float3(hit.triangle_barycentric_coord, 0)};
    }
    return result;
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
                                   device const Surface *surfaces) {
    GeometrySurface out = {};
    out.coordinates = coordinates;
    if (is_quadric(primitive)) {
      Quadric q = quadrics[primitive & ~quadric_bit];
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
    float3 bary = float3(1 - coordinates.x - coordinates.y, coordinates.xy);
    surface_point(out.face, coordinates.xy, out.position, out.error);
    float3 normal =
        bary.x * out.face.n[0].xyz + bary.y * out.face.n[1].xyz + bary.z * out.face.n[2].xyz;
    out.normal = dot(normal, normal) > 0 ? normalize(normal) : out.face.geometric.xyz;
    float2 a = out.face.uv01.zw - out.face.uv01.xy, b = out.face.uv2.xy - out.face.uv01.xy;
    float det = a.x * b.y - a.y * b.x;
    float3 g = out.face.geometric.xyz;
    float3 bitangent =
        normalize(cross(g, abs(g.x) > .9999999f ? float3(0, 1, 0) : float3(1, 0, 0)));
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
  GeometrySurface sample_geometry_light(uint primitive, float2 u, float3 reference,
                                        device const Triangle *triangles,
                                        device const Quadric *quadrics,
                                        device const Surface *surfaces) {
    float3 coordinates;
    if (is_quadric(primitive)) {
      Quadric q = quadrics[primitive & ~quadric_bit];
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
    return geometry_surface(primitive, coordinates, triangles, quadrics, surfaces);
  }

)MSL"
