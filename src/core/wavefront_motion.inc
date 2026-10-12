R"MSL(
// Geometry stays resident. Only one pair of affine matrices per animated
// owner crosses from the host at each coherent shutter sample. Every update
// starts from immutable reference data: rounding never accumulates over frames.
struct MotionPose { float4 forward[3], inverse[3]; };
template<typename T> struct MotionReference { uint index, owner; T value; };
struct MotionBox { packed_float3 lo, hi; };
struct MotionInstance {
  packed_float3 columns[4];
  uint options, mask, table_offset, acceleration_index;
};
static_assert(sizeof(MotionReference<Triangle>) == 176 &&
              sizeof(MotionReference<Quadric>) == 192 &&
              sizeof(MotionReference<Surface>) == 256 &&
              sizeof(MotionReference<MeshInstance>) == 144 &&
              sizeof(MotionInstance) == 64, "motion upload ABI mismatch");
float3 motion_point(MotionPose p, float3 v) {
  return float3(dot(p.forward[0], float4(v,1)), dot(p.forward[1], float4(v,1)),
                dot(p.forward[2], float4(v,1)));
}
float3 motion_normal(MotionPose p, float3 n) {
  return p.inverse[0].xyz * n.x + p.inverse[1].xyz * n.y + p.inverse[2].xyz * n.z;
}
float4 motion_row(float4 a, thread const float4 *b) {
  return a.x*b[0] + a.y*b[1] + a.z*b[2] + float4(0,0,0,a.w);
}
void motion_compose(MotionPose p, thread float4 *forward, thread float4 *inverse) {
  float4 f[3], inv[3];
  for (uint a=0; a<3; ++a) {
    f[a] = motion_row(p.forward[a], forward);
    inv[a] = motion_row(inverse[a], p.inverse);
  }
  for (uint a=0; a<3; ++a) { forward[a]=f[a]; inverse[a]=inv[a]; }
}
kernel void motion_triangles(device Triangle *out [[buffer(0)]],
    device const MotionReference<Triangle> *references [[buffer(1)]],
    device const MotionPose *poses [[buffer(2)]], constant uint &count [[buffer(4)]],
    uint id [[thread_position_in_grid]]) {
  if (id >= count) return;
  auto r=references[id]; auto p=poses[r.owner]; auto t=r.value;
  for (uint a=0; a<3; ++a) {
    t.p[a].xyz=motion_point(p,t.p[a].xyz);
    t.n[a].xyz=motion_normal(p,t.n[a].xyz);
  }
  t.geometric.xyz=motion_normal(p,t.geometric.xyz);
  out[r.index]=t;
}
kernel void motion_surfaces(device Surface *out [[buffer(0)]],
    device const MotionReference<Surface> *references [[buffer(1)]],
    device const MotionPose *poses [[buffer(2)]], constant uint &count [[buffer(4)]],
    uint id [[thread_position_in_grid]]) {
  if (id >= count) return;
  auto r=references[id]; auto s=r.value;
  motion_compose(poses[r.owner],s.forward,s.inverse);
  motion_compose(poses[r.owner],s.bump_forward,s.bump_inverse);
  s.bump_transformed=1; out[r.index]=s;
}
kernel void motion_instances(device MeshInstance *out [[buffer(0)]],
    device const MotionReference<MeshInstance> *references [[buffer(1)]],
    device const MotionPose *poses [[buffer(2)]], device MotionInstance *descriptors [[buffer(3)]],
    constant uint &count [[buffer(4)]], constant bool &refit [[buffer(5)]],
    uint id [[thread_position_in_grid]]) {
  if (id >= count) return;
  auto r=references[id]; auto s=r.value;
  motion_compose(poses[r.owner],s.forward,s.inverse); out[r.index]=s;
  if (refit && s.triangle_count) {
    for (uint a=0;a<4;++a)
      descriptors[r.index].columns[a]=packed_float3(s.forward[0][a],s.forward[1][a],s.forward[2][a]);
  }
}
kernel void motion_quadrics(device Quadric *out [[buffer(0)]],
    device const MotionReference<Quadric> *references [[buffer(1)]],
    device const MotionPose *poses [[buffer(2)]], device MotionBox *boxes [[buffer(3)]],
    constant uint &count [[buffer(4)]], constant bool &refit [[buffer(5)]],
    uint id [[thread_position_in_grid]]) {
  if (id >= count) return;
  auto r=references[id]; auto q=r.value;
  motion_compose(poses[r.owner],q.forward,q.inverse);
  float3 x=float3(q.forward[0].x,q.forward[1].x,q.forward[2].x);
  float3 y=float3(q.forward[0].y,q.forward[1].y,q.forward[2].y);
  float3 z=float3(q.forward[0].z,q.forward[1].z,q.forward[2].z);
  float radius=length(x);
  q.measure.x=abs(dot(x,cross(y,z)));
  q.measure.z=abs(length(y)-radius)<1e-6f*radius && abs(length(z)-radius)<1e-6f*radius &&
    abs(dot(x,y))+abs(dot(x,z))+abs(dot(y,z))<1e-6f*radius*radius ? radius : 0;
  out[r.index]=q;
  if (!refit) return;
  float3 lo,hi;
  for (uint a=0;a<3;++a) {
    float4 row=q.forward[a];
    float lower=row.w,upper=row.w,magnitude=abs(row.w);
    if (q.csg_count || q.curve_offset) {
      float3 v=row.xyz*q.csg_lo.xyz,w=row.xyz*q.csg_hi.xyz;
      lower+=dot(min(v,w),float3(1)); upper+=dot(max(v,w),float3(1));
      magnitude+=dot(max(abs(v),abs(w)),float3(1));
    } else {
      float extent=length(row.xyz); lower-=extent; upper+=extent; magnitude+=extent;
    }
    float pad=64*FLT_EPSILON*max(1.f,magnitude);
    lo[a]=lower-pad; hi[a]=upper+pad;
  }
  boxes[r.index]={packed_float3(lo),packed_float3(hi)};
}
)MSL"
