R"MSL(
  struct Field { uint offset, nx, ny, nz; };
  struct Medium {
    float4 inverse[3], lo, hi, emission, controls;
    Field density, temperature, red, green, blue, scale, blackbody, majorant;
  };
  static_assert(sizeof(Medium) == 240, "medium ABI mismatch");
  // Invert the optical-depth CDF of the piecewise-constant majorant grid.
  // The transformed direction is NOT normalized: t remains world distance.
  float medium_flight(Medium m, float3 origin, float3 direction, float sigma,
                      float u, device const float *values, thread float &majorant) {
    majorant = 0;
    if (sigma == 0) return INFINITY;
    float3 o, d;
    for (uint a=0; a<3; ++a) {
      o[a] = dot(m.inverse[a],float4(origin,1));
      d[a] = dot(m.inverse[a].xyz,direction);
    }
    o = (o-m.lo.xyz)/(m.hi.xyz-m.lo.xyz);
    d /= m.hi.xyz-m.lo.xyz;
    float t=0, end=INFINITY;
    for (uint a=0; a<3; ++a) {
      if (d[a] == 0) {
        if (o[a]<0 || o[a]>1) return INFINITY;
      } else {
        float a0=-o[a]/d[a], a1=(1-o[a])/d[a];
        t=max(t,min(a0,a1)); end=min(end,max(a0,a1));
      }
    }
    if (!(t<end)) return INFINITY;
    uint n=m.majorant.nx;
    float3 q=(o+t*d)*n;
    int3 cell=clamp(int3(floor(q)),int3(0),int3(n-1));
    int3 step=select(int3(-1),int3(1),d>0);
    float3 next,delta;
    for (uint a=0;a<3;++a) {
      delta[a]=d[a] == 0 ? INFINITY : abs(1/(n*d[a]));
      next[a]=d[a] == 0 ? INFINITY : max(t,((cell[a]+(step[a]>0))/float(n)-o[a])/d[a]);
    }
    float optical=-log(1-u);
    while (t<end) {
      float stop=min(end,min(next.x,min(next.y,next.z)));
      float rate=sigma*values[m.majorant.offset+cell.x+n*(cell.y+n*cell.z)];
      if (rate>0 && optical < rate*(stop-t)) {
        majorant=rate;
        return t+optical/rate;
      }
      optical -= rate*(stop-t);
      t=stop;
      for (uint a=0;a<3;++a)
        if (next[a]<=stop) {
          cell[a]+=step[a]; next[a]+=delta[a];
          if (cell[a]<0 || cell[a]>=int(n)) return INFINITY;
        }
    }
    return INFINITY;
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
