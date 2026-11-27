R"MSL(
  struct Field { uint offset, nx, ny, nz; };
  struct Medium {
    float4 inverse[3], lo, hi, emission, controls;
    Field density, temperature, red, green, blue, scale, blackbody, majorant;
  };
  static_assert(sizeof(Medium) == 240, "medium ABI mismatch");
  // Keep world-distance ray parameters while traversing the local majorant grid.
  // Shadow ratio tracking samples many events on one straight segment; retaining
  // this state avoids transforming, clipping and locating its cell at every event.
  Medium motion_medium(Medium m, MotionContext motion, uint id) {
    uint owner=motion_owner(motion,4,id);
    if(owner) {
      auto pose=motion_pose(motion,owner);
      motion_multiply(m.inverse,m.inverse,pose.inverse);
    }
    return m;
  }
  struct MediumFlightIterator {
    float t, end, sigma;
    uint n, offset;
    int3 cell, step;
    float3 next, delta;
  };
  MediumFlightIterator medium_iterator(Medium m, float3 origin, float3 direction, float sigma) {
    MediumFlightIterator it = {};
    it.sigma = sigma;
    if (sigma == 0) return it;
    float3 o, d;
    for (uint a = 0; a < 3; ++a) {
      o[a] = dot(m.inverse[a], float4(origin, 1));
      d[a] = dot(m.inverse[a].xyz, direction);
    }
    o = (o - m.lo.xyz) / (m.hi.xyz - m.lo.xyz);
    // Do not normalize d: affine scaling must not change world-space flight lengths.
    d /= m.hi.xyz - m.lo.xyz;
    it.end = INFINITY;
    for (uint a = 0; a < 3; ++a) {
      if (d[a] == 0) {
        if (o[a] < 0 || o[a] > 1) { it.end = 0; return it; }
      } else {
        float near = -o[a] / d[a], far = (1 - o[a]) / d[a];
        it.t = max(it.t, min(near, far));
        it.end = min(it.end, max(near, far));
      }
    }
    if (!(it.t < it.end)) return it;
    it.n = m.majorant.nx;
    it.offset = m.majorant.offset;
    it.cell = clamp(int3(floor((o + it.t * d) * it.n)), int3(0), int3(it.n - 1));
    it.step = select(int3(-1), int3(1), d > 0);
    for (uint a = 0; a < 3; ++a) {
      it.delta[a] = d[a] == 0 ? INFINITY : abs(1 / (it.n * d[a]));
      it.next[a] = d[a] == 0 ? INFINITY :
          max(it.t, ((it.cell[a] + (it.step[a] > 0)) / float(it.n) - o[a]) / d[a]);
    }
    return it;
  }
  // Invert the piecewise-constant optical-depth CDF. Consume one uniform per
  // flight, carrying its residual optical depth through empty and occupied cells.
  // Returned distances are absolute from the iterator origin, not the last event.
  float medium_next(thread MediumFlightIterator &it, float u, device const float *values,
                    thread float &majorant) {
    majorant = 0;
    float optical = -log(1 - u);
    while (it.t < it.end) {
      float stop = min(it.end, min(it.next.x, min(it.next.y, it.next.z)));
      float rate = it.sigma * values[it.offset + it.cell.x + it.n * (it.cell.y + it.n * it.cell.z)];
      if (rate > 0 && optical < rate * (stop - it.t)) {
        majorant = rate;
        it.t += optical / rate;
        return it.t;
      }
      optical -= rate * (stop - it.t);
      it.t = stop;
      // Advance every tied axis at a grid edge/corner. An exhausted iterator
      // remains exhausted, including subsequent calls with a zero-length flight.
      for (uint a = 0; a < 3; ++a) if (it.next[a] <= stop) {
        it.cell[a] += it.step[a];
        it.next[a] += it.delta[a];
        if (it.cell[a] < 0 || it.cell[a] >= int(it.n)) { it.end = it.t; return INFINITY; }
      }
    }
    return INFINITY;
  }
  float medium_flight(Medium m, float3 origin, float3 direction, float sigma,
                      float u, device const float *values, thread float &majorant) {
    auto iterator = medium_iterator(m, origin, direction, sigma);
    return medium_next(iterator, u, values, majorant);
  }
  bool lex_less(int3 a, int3 b) {
    return a.x != b.x ? a.x < b.x : (a.y != b.y ? a.y < b.y : a.z < b.z);
  }
  float sparse_value(Field f, int3 p, device const float *values) {
    int3 key = p >> 12;
    uint root = as_type<uint>(values[f.offset+12]);
    uint count = as_type<uint>(values[f.offset+13]), low = 0, high = count;
    while (low < high) {
      uint mid = (low + high) / 2, at = root + 4*mid;
      int3 found(as_type<int>(values[at]), as_type<int>(values[at+1]), as_type<int>(values[at+2]));
      if (lex_less(found, key)) low = mid + 1;
      else high = mid;
    }
    if (low == count)
      return values[f.offset+14];
    uint at = root + 4*low;
    int3 found(as_type<int>(values[at]), as_type<int>(values[at+1]), as_type<int>(values[at+2]));
    if (any(found != key))
      return values[f.offset+14];
    float entry = values[at+3];
    // NanoVDB's internal and leaf arrays use z-fastest indexing.
    for (uint level = 0; level < 3; ++level) {
      uint bits = as_type<uint>(entry);
      if (!(bits & 0x80000000u)) return entry;
      uint shift = level == 0 ? 7 : (level == 1 ? 3 : 0);
      uint logdim = level == 0 ? 5 : (level == 1 ? 4 : 3);
      uint3 q = (uint3(p) >> shift) & ((1u << logdim)-1);
      uint index = (q.x << (2*logdim)) + (q.y << logdim) + q.z;
      entry = values[(bits & 0x7fffffffu) + index];
    }
    return entry;
  }
  float field_value(Field f, float3 uvw, device const float *values) {
    if (!f.nx)
      return 0;
    if (f.nx == UINT_MAX) {
      float3 q;
      for (uint a = 0; a < 3; ++a) {
        uint o = f.offset + 4*a;
        q[a] = dot(float3(values[o], values[o+1], values[o+2]), uvw) + values[o+3];
      }
      // Keep conversion defined for rays far beyond a sparse grid.
      if (any(abs(q) > 2147483000.f)) return values[f.offset+14];
      int3 lower = int3(floor(q));
      float3 w = q - float3(lower);
      float result = 0;
      for (int z = 0; z < 2; ++z)
        for (int y = 0; y < 2; ++y)
          for (int x = 0; x < 2; ++x)
            result += sparse_value(f, lower + int3(x,y,z), values) *
                       (x ? w.x : 1-w.x) * (y ? w.y : 1-w.y) * (z ? w.z : 1-w.z);
      return result;
    }
    uint3 dims(f.nx, f.ny, f.nz);
    float3 q = clamp(uvw * float3(dims) - .5f, float3(0), float3(dims - 1));
    uint3 lo = uint3(floor(q)), hi = min(lo + 1, dims - 1);
    float3 w = q - float3(lo);
    float value = 0;
    for (uint z = 0; z < 2; ++z)
      for (uint y = 0; y < 2; ++y)
        for (uint x = 0; x < 2; ++x) {
          uint index = (x ? hi.x : lo.x) + dims.x *
                       ((y ? hi.y : lo.y) + dims.y * (z ? hi.z : lo.z));
          value += values[f.offset + index] * (x ? w.x : 1-w.x) *
                   (y ? w.y : 1-w.y) * (z ? w.z : 1-w.z);
        }
    return value;
  }
  float3 medium_coordinates(Medium m, float3 p) {
    float4 q(p, 1);
    float3 local(dot(m.inverse[0], q), dot(m.inverse[1], q), dot(m.inverse[2], q));
    return m.controls.w == 2 ? local : (local - m.lo.xyz) / (m.hi.xyz - m.lo.xyz);
  }
  float medium_density(Medium m, float3 uvw, device const float *values) {
    return m.controls.w != 2 && (any(uvw < 0) || any(uvw > 1)) ? 0 : field_value(m.density, uvw, values);
  }
  float3 medium_emission(Medium m, float3 uvw, device const float *values) {
    if (m.controls.w == 0)
      return m.emission.xyz;
    float scale = m.controls.x * (m.scale.nx ? field_value(m.scale, uvw, values) : 1);
    if (scale == 0)
      return 0;
    if (m.temperature.nx) {
      float kelvin = (field_value(m.temperature, uvw, values) - m.controls.z) * m.controls.y;
      if (kelvin <= 100)
        return 0;
      float q = clamp((kelvin - 100) / 10, 0.f, 1990.f);
      uint i = min(uint(q), 1989u), a = m.blackbody.offset + 3*i;
      float3 low(values[a], values[a+1], values[a+2]);
      float3 high(values[a+3], values[a+4], values[a+5]);
      return scale * mix(low, high, q - i);
    }
    if (m.controls.w == 2)
      return m.emission.xyz;
    return scale * float3(field_value(m.red, uvw, values), field_value(m.green, uvw, values),
                          field_value(m.blue, uvw, values));
  }
)MSL"
