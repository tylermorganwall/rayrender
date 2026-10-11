R"MSL(
// Trilinear interpolation restricted to a ray is cubic. Return every crossing
// in [0,1], isolating roots on monotone intervals before bisection. Fixed arrays
// hold at most three roots; no stochastic flight spacing enters classification.
float haze_polynomial(float4 c, float u) { return ((c.w * u + c.z) * u + c.y) * u + c.x; }
uint haze_edges(thread const float *values, float3 start, float3 delta, float threshold,
                thread float *edges, thread float4 &polynomial) {
  float low = values[0], high = values[0];
  for (uint i = 1; i < 8; ++i) {
    low = min(low, values[i]);
    high = max(high, values[i]);
  }
  // A trilinear interpolant stays within its corner extrema. Match the native
  // solver's whole-cell rejection before constructing a polynomial.
  if (high < threshold || low >= threshold) {
    polynomial = float4(low - threshold, 0, 0, 0);
    edges[0] = 0;
    edges[1] = 1;
    return 2;
  }
  float4 c = 0;
  for (uint corner = 0; corner < 8; ++corner) {
    float4 term = float4(values[corner], 0, 0, 0);
    for (uint axis = 0; axis < 3; ++axis) {
      bool upper = corner & (1u << axis);
      float a = upper ? start[axis] : 1 - start[axis], b = upper ? delta[axis] : -delta[axis];
      term = float4(term.x * a, term.y * a + term.x * b, term.z * a + term.y * b,
                    term.w * a + term.z * b);
    }
    c += term;
  }
  c.x -= threshold;
  float magnitude = max(max(abs(c.x), abs(c.y)), max(abs(c.z), abs(c.w)));
  if (magnitude > 0)
    c /= magnitude;
  polynomial = c;
  float cuts[4] = {0, 1, 1, 1};
  uint count = 2;
  float a = 3 * c.w, b = 2 * c.z, d = c.y;
  if (abs(a) < 1e-7f) {
    float root = abs(b) >= 1e-7f ? -d / b : -1;
    if (root > 0 && root < 1)
      cuts[count++] = root;
  } else {
    float discriminant = b * b - 4 * a * d;
    if (discriminant >= 0) {
      float q = -.5f * (b + copysign(sqrt(discriminant), b));
      float2 roots = q != 0 ? float2(q / a, d / q) : float2(-b / (2 * a), -1);
      for (uint j = 0; j < 2; ++j)
        if (roots[j] > 0 && roots[j] < 1)
          cuts[count++] = roots[j];
    }
  }
  for (uint j = 1; j < count; ++j)
    for (uint k = j; k > 0 && cuts[k] < cuts[k - 1]; --k) {
      float swap = cuts[k];
      cuts[k] = cuts[k - 1];
      cuts[k - 1] = swap;
    }
  uint n = 1;
  edges[0] = 0;
  for (uint j = 1; j < count; ++j) {
    float lo = cuts[j - 1], hi = cuts[j], left = haze_polynomial(c, lo),
          right = haze_polynomial(c, hi);
    // Tangencies do not switch selection. Exact derivative-root crossings
    // are included once, without manufacturing positive-width intervals.
    if (left == 0 && lo > edges[n - 1] && lo < 1)
      edges[n++] = lo;
    if ((left < 0 && right > 0) || (left > 0 && right < 0)) {
      for (uint k = 0; k < 28; ++k) {
        float mid = lo + (hi - lo) * .5f;
        if (mid == lo || mid == hi)
          break;
        if ((haze_polynomial(c, mid) < 0) == (left < 0))
          lo = mid;
        else
          hi = mid;
      }
      float root = lo + (hi - lo) * .5f;
      if (root > edges[n - 1] && root < 1)
        edges[n++] = root;
    }
  }
  edges[n++] = 1;
  return n;
}
void haze_accumulate(thread AirSegment &sum, AirSegment segment) {
  sum.radiance += sum.transmission * segment.radiance;
  sum.transmission *= segment.transmission;
}
// Merge equal selections across cell boundaries: only actual threshold changes
// query the spectrum or reset the cumulative atmosphere anchor. All SIMD lanes
// traverse identical intervals so spectral reductions remain convergent.
void haze_select(thread AirRay &ray, thread AirSegment &sum, float3 point, float3 direction,
                 bool enabled, float uniform, device const AtmosphereData &air,
                 device const float *data) {
  if (ray.origin.w && bool(ray.direction.w) != enabled)
    haze_accumulate(sum, air_advance(ray, point, uniform, air, data));
  air_start(ray, point, direction, enabled);
}
AirSegment air_walk(thread AirRay &ray, float3 origin, float3 direction, float distance,
                    bool enabled, Medium medium, float uniform, device const AtmosphereData &air,
                    device const float *data, device const float *fields) {
  AirSegment sum = {float3(1), float3(0)};
  if (!enabled || medium.lo.w == 0 || medium.controls.w == 0) {
    enabled = enabled && (medium.lo.w == 0 || medium.controls.w != 0 || 1 < medium.lo.w);
    air_start(ray, origin, direction, enabled);
    return air_advance(ray, origin + distance * direction, uniform, air, data);
  }
  // World-distance parameter survives both medium and sparse-grid transforms.
  float3 local, local_d;
  for (uint a = 0; a < 3; ++a) {
    local[a] = dot(medium.inverse[a], float4(origin, 1));
    local_d[a] = dot(medium.inverse[a].xyz, direction);
  }
  float3 uv = (local - medium.lo.xyz) / (medium.hi.xyz - medium.lo.xyz);
  float3 du = local_d / (medium.hi.xyz - medium.lo.xyz);
  float enter = 0, leave = distance;
  for (uint a = 0; a < 3; ++a) {
    if (du[a] == 0) {
      if (uv[a] < 0 || uv[a] > 1)
        leave = 0;
    } else {
      float2 t = float2(-uv[a], 1 - uv[a]) / du[a];
      enter = max(enter, min(t.x, t.y));
      leave = min(leave, max(t.x, t.y));
    }
  }
  enter = min(enter, distance);
  leave = max(enter, leave);
  if (enter > 0 || leave == enter)
    haze_select(ray, sum, origin, direction, true, uniform, air, data);
  float3 index, index_d;
  Field field = medium.density;
  if (field.nx == UINT_MAX) {
    for (uint a = 0; a < 3; ++a) {
      uint at = field.offset + 4 * a;
      float3 row = float3(fields[at], fields[at + 1], fields[at + 2]);
      index[a] = dot(row, local) + fields[at + 3];
      index_d[a] = dot(row, local_d);
    }
  } else {
    float3 dims = float3(field.nx, field.ny, field.nz);
    index = uv * dims - .5f;
    index_d = du * dims;
  }
  float cursor = enter;
  uint lane = simd_prefix_exclusive_sum(1u);
  while (cursor < leave) {
    // Skip entire majorant cells whose maximum lies below the cutoff.
    float3 q = uv + cursor * du;
    uint n = medium.majorant.nx;
    int3 major = clamp(int3(floor(q * n)), int3(0), int3(n - 1));
    float stop = leave;
    for (uint a = 0; a < 3; ++a)
      if (du[a] != 0) {
        float plane = (major[a] + (du[a] > 0)) / float(n);
        float next = (plane - uv[a]) / du[a];
        if (next <= cursor)
          next = (plane + (du[a] > 0 ? 1.f : -1.f) / n - uv[a]) / du[a];
        if (next > cursor)
          stop = min(stop, next);
      }
    if (fields[medium.majorant.offset + major.x + n * (major.y + n * major.z)] < medium.lo.w) {
      haze_select(ray, sum, origin + cursor * direction, direction, true, uniform, air, data);
      cursor = stop;
      continue;
    }
    for (uint a = 0; a < 3; ++a)
      if (index_d[a] != 0) {
        float x = fma(cursor, index_d[a], index[a]);
        float plane = index_d[a] > 0 ? floor(x) + 1 : ceil(x) - 1;
        float next = (plane - index[a]) / index_d[a];
        if (next <= cursor)
          next = (plane + (index_d[a] > 0 ? 1.f : -1.f) - index[a]) / index_d[a];
        if (next > cursor)
          stop = min(stop, next);
      }
    float3 cell = floor(fma(float3(cursor + (stop - cursor) * .5f), index_d, index));
    float values[8];
    // One SIMD group owns the ray. Fetch each corner once in a different lane,
    // then broadcast it to the polynomial solver instead of walking the sparse
    // hierarchy eight times in every lane. Spectral work reuses the same group.
    float corner_value = 0;
    if (lane < 8) {
      int3 p = int3(cell) + int3(lane & 1, (lane >> 1) & 1, (lane >> 2) & 1);
      if (field.nx == UINT_MAX)
        corner_value = sparse_value(field, p, fields);
      else {
        p = clamp(p, int3(0), int3(field.nx - 1, field.ny - 1, field.nz - 1));
        corner_value = fields[field.offset + p.x + field.nx * (p.y + field.ny * p.z)];
      }
    }
    for (uint corner = 0; corner < 8; ++corner)
      values[corner] = simd_shuffle(corner_value, corner);
    float edges[5];
    float4 polynomial;
    uint count = haze_edges(values, fma(float3(cursor), index_d, index) - cell,
                            (stop - cursor) * index_d, medium.lo.w, edges, polynomial);
    for (uint j = 1; j < count; ++j) {
      float start = mix(cursor, stop, edges[j - 1]), end = mix(cursor, stop, edges[j]);
      if (!(end > start))
        continue;
      // The same polynomial that isolated the roots classifies each interval.
      // Reusing it avoids another world/index transform and eight density reads.
      bool select = haze_polynomial(polynomial, (edges[j - 1] + edges[j]) * .5f) < 0;
      haze_select(ray, sum, origin + start * direction, direction, select, uniform, air, data);
    }
    cursor = stop;
  }
  if (leave < distance)
    haze_select(ray, sum, origin + leave * direction, direction, true, uniform, air, data);
  haze_accumulate(sum, air_advance(ray, origin + distance * direction, uniform, air, data));
  return sum;
}
)MSL"
