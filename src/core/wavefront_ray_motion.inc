R"MSL(
// One immutable transform hierarchy is shared by all paths. A context selects
// the camera ray's absolute shutter time; secondary/visibility rays retain it.
struct MotionHeader {
  uint count, owners[5], nodes, bounds;
  uint version, byte_size, reserved[2];
};
struct Animation {
  float4 endpoints[2][3], translation[2], rotation[2], scale[2][3];
  float4 before[3], reference_inverse[3];
  uint parent;
  float start, end;
  uint track;
};
struct TransformTrack {
  uint channels[3], version;
  float4 pivot;
  float4 clock;
  float4 frame;
};
struct PropertyTrack {
  uint count, interpolation, property, keys;
};
struct TrackKey {
  float4 value, incoming, outgoing;
  float time;
  uint reserved[3];
};
struct MotionContext {
  device const uchar *data;
  float time;
};
struct MotionPose {
  float4 forward[3], inverse[3];
};
static_assert(sizeof(Animation) == 368 && sizeof(MotionHeader) == 48,
              "motion ABI mismatch");
static_assert(sizeof(TransformTrack) == 64 && sizeof(PropertyTrack) == 16 && sizeof(TrackKey) == 64,
              "track ABI mismatch");
float4 motion_track(MotionContext m, uint offset, float4 rest) {
  if (!offset) return rest;
  auto track = reinterpret_cast<device const PropertyTrack *>(m.data + offset);
  auto keys = reinterpret_cast<device const TrackKey *>(m.data + track->keys);
  if (m.time <= keys[0].time) return keys[0].value;
  if (m.time >= keys[track->count - 1].time) return keys[track->count - 1].value;
  uint lo = 0, hi = track->count - 1;
  while (hi - lo > 1) {
    uint mid = lo + (hi - lo) / 2;
    if (keys[mid].time <= m.time) lo = mid; else hi = mid;
  }
  TrackKey a = keys[lo], b = keys[lo + 1];
  if (track->interpolation == 0 || m.time == a.time) return a.value;
  float duration = b.time - a.time, u = (m.time - a.time) / duration;
  float4 value;
  if (track->interpolation == 1) {
    if (track->property == 1) {
      float cosine = clamp(dot(a.value, b.value), -1.f, 1.f);
      if (cosine >= .9995f) value = mix(a.value, b.value, u);
      else {
        float theta = acos(cosine);
        value = (sin((1-u)*theta)*a.value + sin(u*theta)*b.value) / sin(theta);
      }
    } else value = mix(a.value, b.value, u);
  } else {
    float4 m0 = duration * a.outgoing, m1 = duration * b.incoming;
    float4 c0 = 2*a.value - 2*b.value + m0 + m1;
    float4 c1 = -3*a.value + 3*b.value - 2*m0 - m1;
    value = ((c0*u + c1)*u + m0)*u + a.value;
  }
  return track->property == 1 ? normalize(value) : value;
}
uint motion_owner(MotionContext m, uint kind, uint index) {
  if (!(texture_features & 1024))
    return 0;
  auto header = reinterpret_cast<device const MotionHeader *>(m.data);
  return reinterpret_cast<device const uint *>(m.data +
                                               header->owners[kind])[index];
}
float4 motion_row(float4 a, thread const float4 *b) {
  return a.x * b[0] + a.y * b[1] + a.z * b[2] + float4(0, 0, 0, a.w);
}
void motion_multiply(thread float4 *result, thread const float4 *a,
                     thread const float4 *b) {
  float4 rows[3];
  for (uint k = 0; k < 3; ++k)
    rows[k] = motion_row(a[k], b);
  for (uint k = 0; k < 3; ++k)
    result[k] = rows[k];
}
MotionPose motion_pose(MotionContext m, uint owner) {
  MotionPose result;
  result.forward[0] = float4(1, 0, 0, 0);
  result.forward[1] = float4(0, 1, 0, 0);
  result.forward[2] = float4(0, 0, 1, 0);
  if (texture_features & 1024) {
    auto header = reinterpret_cast<device const MotionHeader *>(m.data);
    auto nodes =
        reinterpret_cast<device const Animation *>(m.data + header->nodes);
    // Walk leaf-to-root: parent deltas left-multiply the child's exact native
    // T * quaternion-slerp(R) * S interpolation. No baked-pose approximation.
    while (owner) {
      Animation a = nodes[owner - 1];
      float4 value[3];
      if (a.track) {
        auto track = reinterpret_cast<device const TransformTrack *>(m.data + a.track);
        MotionContext local = m;
        if (track->clock.w != 0) {
          float duration = track->clock.z;
          // Keep the owner's clock local: a parent has its own scene-time
          // placement and must never inherit a child's wrapped clip time.
          float phase = track->clock.x + track->clock.y * (m.time - track->frame.x);
          if (duration == 0) local.time = 0;
          else if (track->clock.w == 1) local.time = clamp(phase, 0.f, duration);
          else {
            float period = track->clock.w == 2 ? duration : 2 * duration;
            float wrapped = fmod(phase, period);
            if (wrapped < 0) wrapped += period;
            if (wrapped >= period) wrapped = 0;
            local.time = track->clock.w == 3 && wrapped > duration ? period - wrapped : wrapped;
          }
        }
        float3 p = motion_track(local, track->channels[0], float4(0)).xyz;
        float4 q = normalize(motion_track(local, track->channels[1], float4(0,0,0,1)));
        float3 s = motion_track(local, track->channels[2], float4(1,1,1,0)).xyz;
        float3 x = q.xyz, pivot = track->pivot.xyz;
        value[0] = float4(float3(1-2*(x.y*x.y+x.z*x.z), 2*(x.x*x.y-q.w*x.z), 2*(x.x*x.z+q.w*x.y))*s, 0);
        value[1] = float4(float3(2*(x.x*x.y+q.w*x.z), 1-2*(x.x*x.x+x.z*x.z), 2*(x.y*x.z-q.w*x.x))*s, 0);
        value[2] = float4(float3(2*(x.x*x.z-q.w*x.y), 2*(x.y*x.z+q.w*x.x), 1-2*(x.x*x.x+x.y*x.y))*s, 0);
        for (uint k = 0; k < 3; ++k) value[k].w = p[k] + pivot[k] - dot(value[k].xyz, pivot);
      } else if (m.time <= a.start || m.time >= a.end) {
        uint endpoint = m.time <= a.start ? 0 : 1;
        for (uint k = 0; k < 3; ++k)
          value[k] = a.endpoints[endpoint][k];
      } else {
        float t = (m.time - a.start) / (a.end - a.start);
        float4 q;
        float cosine = dot(a.rotation[0], a.rotation[1]);
        if (cosine > .9995f)
          q = normalize(mix(a.rotation[0], a.rotation[1], t));
        else {
          float angle = acos(clamp(cosine, -1.f, 1.f));
          q = a.rotation[0] * cos(t * angle) +
              normalize(a.rotation[1] - a.rotation[0] * cosine) *
                  sin(t * angle);
        }
        float3 x = q.xyz;
        float4 rotation[3] = {
            float4(1 - 2 * (x.y * x.y + x.z * x.z), 2 * (x.x * x.y - q.w * x.z),
                   2 * (x.x * x.z + q.w * x.y), 0),
            float4(2 * (x.x * x.y + q.w * x.z), 1 - 2 * (x.x * x.x + x.z * x.z),
                   2 * (x.y * x.z - q.w * x.x), 0),
            float4(2 * (x.x * x.z - q.w * x.y), 2 * (x.y * x.z + q.w * x.x),
                   1 - 2 * (x.x * x.x + x.y * x.y), 0)};
        float4 scale[3];
        for (uint k = 0; k < 3; ++k)
          scale[k] = mix(a.scale[0][k], a.scale[1][k], t);
        motion_multiply(value, rotation, scale);
        float3 translation = mix(a.translation[0].xyz, a.translation[1].xyz, t);
        for (uint k = 0; k < 3; ++k)
          value[k].w = translation[k];
      }
      motion_multiply(value, a.before, value);
      motion_multiply(value, value, a.reference_inverse);
      motion_multiply(result.forward, value, result.forward);
      owner = a.parent;
    }
  }
  float3 x = result.forward[0].xyz, y = result.forward[1].xyz,
         z = result.forward[2].xyz;
  float3 c0 = cross(y, z), c1 = cross(z, x), c2 = cross(x, y);
  float determinant = dot(x, c0);
  float3 translation =
      float3(result.forward[0].w, result.forward[1].w, result.forward[2].w);
  for (uint k = 0; k < 3; ++k) {
    float3 row = float3(c0[k], c1[k], c2[k]) / determinant;
    result.inverse[k] = float4(row, -dot(row, translation));
  }
  return result;
}
void motion_compose(MotionPose p, thread float4 *forward,
                    thread float4 *inverse) {
  motion_multiply(forward, p.forward, forward);
  motion_multiply(inverse, inverse, p.inverse);
}
float3 motion_point(MotionPose p, float3 v) {
  return float3(dot(p.forward[0], float4(v, 1)),
                dot(p.forward[1], float4(v, 1)),
                dot(p.forward[2], float4(v, 1)));
}
float3 motion_normal(MotionPose p, float3 n) {
  return p.inverse[0].xyz * n.x + p.inverse[1].xyz * n.y +
         p.inverse[2].xyz * n.z;
}
Triangle motion_triangle(Triangle t, MotionContext m, uint id) {
  uint owner = motion_owner(m, 0, id);
  if (owner) {
    auto p = motion_pose(m, owner);
    for (uint k = 0; k < 3; ++k) {
      t.p[k].xyz = motion_point(p, t.p[k].xyz);
      t.n[k].xyz = motion_normal(p, t.n[k].xyz);
    }
    t.geometric.xyz = normalize(motion_normal(p, t.geometric.xyz));
  }
  return t;
}
)MSL"
