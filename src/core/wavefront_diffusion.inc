R"MSL(

  // Normalized diffusion retains the CPU's three-axis RGB proposal. A chord can
  // intersect any number of components. Weighted reservoir selection avoids a
  // per-ray candidate array without dropping crossings or truncating the profile.
  struct DiffusionEntry {
    float4 position, normal;
  };
  struct DiffusionExit {
    float4 shape, normalization;
    uint offset, padding[3];
  };
  struct DiffusionSample {
    Hit hit;
    float3 position, error, normal, weight;
    VolumeState outside;
  };
  bool adjacent_membership(uint boundary, bool entering, float3 position, float3 error,
                           float3 geometric, thread VolumeState &adjacent, constant Params &p,
                           device const Triangle *triangles, device const Quadric *quadrics,
                           device const Surface *surfaces, device const Boundary *boundaries,
                           GeometryScene scene) {
    if (!cross_boundary(adjacent, boundary, entering))
      return false;
    if (p.boundary_lo.w != 0) {
      VolumeState resolved;
      if (!initial_membership(offset_origin(position, error, entering ? -geometric : geometric),
                              resolved, p, triangles, quadrics, surfaces, boundaries, scene))
        return false;
      // The reported face can be the losing member of an exact contact. Even
      // when the winner is glass on both sides, retain ALL coincident members
      // so a later glass exit exposes the correct diffusion region.
      adjacent = resolved;
    }
    return true;
  }
  void remember_diffusion_entry(thread Path &path, uint pixel, thread const VolumeState &state,
                                float3 position, float3 geometric,
                                device const Boundary *boundaries, device const Material *materials,
                                device DiffusionEntry *entries) {
    uint active = active_boundary(state, boundaries, materials);
    if (active && boundaries[active].subsurface == 2) {
      float3 normal = dot(path.direction.xyz, geometric) > 0 ? -geometric : geometric;
      entries[pixel] = {float4(position, 0), float4(normal, 0)};
      path.flags |= 8;
    } else
      path.flags &= ~8u;
  }
  float log_add(float a, float b) {
    float largest = max(a, b);
    return isinf(largest) ? largest : largest + log(exp(a - largest) + exp(b - largest));
  }
  float profile_log(float radius, float distance) {
    float d = max(distance, FLT_MIN), x = d / radius;
    return log_add(-x, -x / 3) - log(8 * pi) - log(radius) - log(d);
  }
  float3 box_corner(Boundary b, uint corner) {
    return float3((corner & 1) ? b.hi.x : b.lo.x, (corner & 2) ? b.hi.y : b.lo.y,
                  (corner & 4) ? b.hi.z : b.lo.z);
  }
  bool diffusion_surface(uint body, DiffusionEntry entry, uint pixel, thread uint &dimension,
                         constant Params &p, device const uint *directions, device const uint *blue,
                         device const Triangle *triangles, device const Quadric *quadrics,
                         device const Surface *surfaces, device const Boundary *boundaries,
                         device const Material *materials, GeometryScene scene,
                         thread DiffusionSample &result, thread bool &invalid) {
    Boundary b = boundaries[body];
    if (all(b.color.xyz == 0))
      return false;
    Frame frame = normal_frame(entry.normal.xyz);
    float axes = sample1(pixel, dimension, p, directions, blue);
    uint axis = axes < .5f ? 2 : (axes < .75f ? 0 : 1);
    uint channel = min(2u, uint(3 * sample1(pixel, dimension, p, directions, blue)));
    float component = sample1(pixel, dimension, p, directions, blue);
    float uniform = max(sample1(pixel, dimension, p, directions, blue), FLT_EPSILON);
    float radius = -b.radius[channel] * (component < .25f ? 1.f : 3.f) * log(1 - uniform);
    float phi = 2 * pi * sample1(pixel, dimension, p, directions, blue);
    float3 basis[3] = {frame.x, frame.y, frame.z};
    float farthest = 0, magnitude = 1;
    for (uint j = 0; j < 8; ++j)
      farthest = max(farthest, length(box_corner(b, j) - entry.position.xyz));
    if (radius > farthest * (1 + 8 * FLT_EPSILON))
      return false;
    float3 center = entry.position.xyz +
                    radius * (cos(phi) * basis[(axis + 1) % 3] + sin(phi) * basis[(axis + 2) % 3]);
    float3 direction = basis[axis];
    float lower = INFINITY, upper = -INFINITY;
    for (uint j = 0; j < 8; ++j) {
      float3 corner = box_corner(b, j);
      float projection = dot(corner - center, direction);
      lower = min(lower, projection);
      upper = max(upper, projection);
      magnitude = max(magnitude, max(abs(corner.x), max(abs(corner.y), abs(corner.z))));
    }
    float padding = 64 * FLT_EPSILON * magnitude;
    ray probe;
    probe.origin = center + (lower - padding) * direction;
    probe.direction = direction;
    probe.min_distance = FLT_MIN;
    probe.max_distance = upper - lower + 2 * padding;
    float total = -INFINITY, selected_score = -INFINITY;
    float3 last_position = 0, last_error = 0;
    bool have_candidate = false, last_leaving = false;
    while (true) {
      auto hit = trace_geometry(probe, scene, quadrics);
      if (hit.type == intersection_type::none)
        break;
      GeometrySurface contact =
          geometry_surface(hit.primitive_id, hit.coordinates, triangles, quadrics, surfaces);
      Triangle t = contact.face;
      if (t.boundary) {
        float3 position = contact.position, error = contact.error;
        float3 g = t.geometric.xyz;
        float3 forward = dot(direction, g) > 0 ? g : -g;
        VolumeState before, after;
        // Classify BOTH sides: one closest hit alone cannot represent coincident
        // milk/glass faces. The winning region, not mesh identity, defines an exit.
        bool classified = initial_membership(offset_origin(position, error, -forward), before, p,
                                             triangles, quadrics, surfaces, boundaries, scene) &&
                          initial_membership(offset_origin(position, error, forward), after, p,
                                             triangles, quadrics, surfaces, boundaries, scene);
        if (!classified) {
          invalid = true;
          return false;
        }
        bool leaving = active_boundary(before, boundaries, materials) == body;
        bool entering = active_boundary(after, boundaries, materials) == body;
        if (leaving != entering) {
          bool duplicate = have_candidate && last_leaving == leaving &&
                           length(position - last_position) <=
                               length(error) + length(last_error) + 8 * FLT_EPSILON * magnitude;
          if (!duplicate) {
            float score = -INFINITY, distance = length(position - entry.position.xyz);
            for (uint c = 0; c < 3; ++c)
              if (b.color[c] > 0)
                score = log_add(score, log(b.color[c]) + profile_log(b.radius[c], distance));
            float next_total = log_add(total, score);
            float choose = sample1(pixel, dimension, p, directions, blue);
            if (isfinite(score) &&
                (selected_score == -INFINITY || choose < exp(score - next_total))) {
              result.hit = {pixel,
                            hit.primitive_id,
                            hit.coordinates.x,
                            hit.coordinates.y,
                            hit.distance,
                            hit.coordinates.z,
                            {0, 0}};
              result.position = position;
              result.error = error;
              result.normal = leaving ? forward : -forward;
              result.outside = leaving ? after : before;
              selected_score = score;
            }
            total = next_total;
            last_position = position;
            last_error = error;
            last_leaving = leaving;
            have_candidate = true;
          }
        }
      }
      float next = as_type<float>(as_type<uint>(hit.distance) + 1u);
      if (!(next > probe.min_distance) || !isfinite(next)) {
        invalid = true;
        return false;
      }
      probe.min_distance = next;
    }
    if (!isfinite(selected_score) || !isfinite(total))
      return false;
    float3 coordinates = frame.local(result.position - entry.position.xyz);
    float distance = length(coordinates), density = -INFINITY;
    for (uint a = 0; a < 3; ++a) {
      float projected = length(float2(coordinates[(a + 1) % 3], coordinates[(a + 2) % 3]));
      float cosine = abs(dot(result.normal, basis[a]));
      if (cosine == 0)
        continue;
      for (uint c = 0; c < 3; ++c)
        density = log_add(density, profile_log(b.radius[c], projected) +
                                       log((a == 2 ? .5f : .25f) * cosine / 3));
    }
    if (!isfinite(density))
      return false;
    for (uint c = 0; c < 3; ++c)
      result.weight[c] =
          b.color[c] * exp(profile_log(b.radius[c], distance) - density - selected_score + total);
    return all(isfinite(result.weight));
  }
  BsdfValue diffusion_eval(DiffusionExit table, float3 normal, float3 wi,
                           device const float4 *samples) {
    float cosine = clamp(dot(normal, wi), 0.f, 1.f), eta = table.shape.x;
    if (cosine == 0)
      return {0, 0};
    if (eta == 1)
      return {float3(cosine / pi), cosine / pi};
    float critical = table.shape.y, support = table.shape.z, mass = table.shape.w;
    if (cosine * cosine <= critical)
      return {0, 0};
    float t = sqrt(max(0.f, (cosine * cosine - critical) / support)),
          coordinate = min(t, 1.f) * 1024;
    uint bin = min(uint(coordinate), 1023u);
    float density =
        mix(samples[table.offset + bin].x, samples[table.offset + bin + 1].x, coordinate - bin);
    float pdf = density * cosine / (2 * pi * support * t * mass);
    float value = cosine / pi * (1 - fresnel(cosine, eta)) / table.normalization.x * eta * eta;
    return {float3(value), pdf};
  }
  float3 diffusion_direction(DiffusionExit table, float3 normal, float2 u,
                             device const float4 *samples) {
    float cosine;
    if (table.shape.x == 1)
      cosine = sqrt(1 - u.x);
    else {
      float target = (1 - u.x) * table.shape.w;
      uint lo = 0, hi = 1024;
      while (lo + 1 < hi) {
        uint mid = (lo + hi) / 2;
        if (samples[table.offset + mid].y <= target)
          lo = mid;
        else
          hi = mid;
      }
      float2 a = samples[table.offset + lo].xy, b = samples[table.offset + lo + 1].xy;
      float area = max(0.f, (target - a.y) * 1024), delta = b.x - a.x;
      float denominator = a.x + sqrt(max(0.f, a.x * a.x + 2 * delta * area));
      float fraction = denominator > 0 ? 2 * area / denominator : 0;
      float t = (lo + clamp(fraction, 0.f, 1.f)) / 1024;
      cosine = sqrt(table.shape.y + table.shape.z * t * t);
      float edge = sqrt(table.shape.y);
      cosine = clamp(cosine, as_type<float>(as_type<uint>(edge) + 1u), 1.f);
    }
    float radius = sqrt(max(0.f, 1 - cosine * cosine)), phi = 2 * pi * u.y;
    return normal_frame(normal).world(float3(radius * cos(phi), radius * sin(phi), cosine));
  }

)MSL"
