R"MSL(
float4 camera_slerp(float t, float4 a, float4 b) {
  float cosine = dot(a, b);
  if (cosine < 0) {
    b = -b;
    cosine = -cosine;
  }
  if (cosine > .9995f)
    return normalize(mix(a, b, t));
  float angle = acos(clamp(cosine, -1.f, 1.f));
  return a * cos(t * angle) + normalize(b - a * cosine) * sin(t * angle);
}
float3 camera_rotate(float4 q, float3 v) { return v + 2 * cross(q.xyz, cross(q.xyz, v) + q.w * v); }
// Film-to-scene lens tracing matches RealisticCamera: spherical interfaces,
// aperture rejection, total internal reflection, and pupil-area/cos^4 weights.
// A blocked lens sample is a valid zero contribution, never a discarded ray.
float camera_ray(Camera c, float2 film, float2 lens_sample, float time, device const float4 *data,
                 thread float3 &origin, thread float3 &direction) {
  float3 right = c.right.xyz, up = c.up.xyz, forward = c.forward.xyz, center = c.origin.xyz;
  float focal = 1;
  if (c.moving == 1) {
    float4 q = camera_slerp(time, c.motion_rotation[0], c.motion_rotation[1]);
    right = camera_rotate(q, float3(1, 0, 0));
    up = camera_rotate(q, float3(0, 1, 0));
    forward = camera_rotate(q, float3(0, 0, 1));
    float4 position = mix(c.motion_origin[0], c.motion_origin[1], time);
    center = position.xyz;
    focal = position.w > 0 && isfinite(position.w) ? position.w : 1;
  }
  if (c.moving == 2) {
    float scene_time = mix(c.shutter.x, c.shutter.y, time) - c.shutter.w;
    uint count = c.track_count, offset = c.track_offset;
    uint lo = 0, hi = count;
    while (lo < hi) {
      uint mid = lo + (hi-lo)/2;
      if (data[offset+5*mid].x <= scene_time) lo = mid+1; else hi = mid;
    }
    uint left = lo == 0 ? 0 : lo-1, next = min(left+1, count-1);
    float begin = data[offset+5*left].x, end = data[offset+5*next].x;
    float fraction = left == next ? 0 : clamp((scene_time-begin)/(end-begin),0.f,1.f);
    float left_fov = data[offset+5*left+3].w, next_fov = data[offset+5*next+3].w;
    uint left_kind = left_fov < 0 ? 3 : left_fov == 0 ? 1 : left_fov == 360 ? 2 : 0;
    uint next_kind = next_fov < 0 ? 3 : next_fov == 0 ? 1 : next_fov == 360 ? 2 : 0;
    if (left_kind != next_kind) fraction = 0;
    c.orthographic = left_kind;
    float4 position = mix(data[offset+5*left+1],data[offset+5*next+1],fraction);
    float4 target = mix(data[offset+5*left+2],data[offset+5*next+2],fraction);
    float4 camera_up = mix(data[offset+5*left+3],data[offset+5*next+3],fraction);
    float2 size = mix(data[offset+5*left+4].xy,data[offset+5*next+4].xy,fraction);
    center = position.xyz;
    forward = normalize(target.xyz-center);
    right = normalize(cross(camera_up.xyz,forward));
    up = normalize(cross(forward,right));
    if (c.orthographic == 3) {
      // Preserve the physical camera's historical LookAt basis, including its
      // transverse scale when the authored up is not perpendicular to view.
      right = cross(normalize(camera_up.xyz),forward);
      up = cross(forward,right);
    }
    focal = target.w;
    c.right.w = position.w*.5f;
    if (c.orthographic == 1) c.lens.xy = size*.5f;
    else if (c.orthographic == 0) {
      float aspect = c.lens.z;
      c.lens.y = tan(camera_up.w*pi/360.f);
      c.lens.x = aspect*c.lens.y;
    }
  }
  origin = center;
  direction = forward;
  if (c.orthographic == 2) {
    float theta = pi * (1 - film.y), phi = pi + 2 * pi * film.x;
    direction = sin(theta) * sin(phi) * right + cos(theta) * up + sin(theta) * cos(phi) * forward;
    return 1;
  }
  if (c.orthographic < 2) {
    float3 corner = c.lower_left.xyz, horizontal = c.horizontal.xyz, vertical = c.vertical.xyz;
    if (c.moving) {
      float scale = c.orthographic ? 1 : focal;
      corner = center - c.lens.x * scale * right - c.lens.y * scale * up +
               (c.orthographic ? float3(0) : focal * forward);
      horizontal = 2 * c.lens.x * scale * right;
      vertical = 2 * c.lens.y * scale * up;
    }
    float3 point = corner + film.x * horizontal + film.y * vertical;
    float radius = sqrt(lens_sample.x) * c.right.w, phi = 2 * pi * lens_sample.y;
    origin = c.orthographic ? point : center + radius * (cos(phi) * right + sin(phi) * up);
    direction = normalize(c.orthographic ? forward : point - origin);
    return 1;
  }
  // The native integrator passes (1-u,1-v) to physical lens cameras.
  film = 1 - film;
  float3 point = float3((1 - 2 * film.x) * c.lens.x, (2 * film.y - 1) * c.lens.y, 0);
  float radius = length(point.xy);
  uint pupil = min(c.pupils - 1, uint(radius / (c.lens.z * .5f) * c.pupils));
  float4 bounds = data[c.elements + pupil];
  float2 lens = mix(bounds.xy, bounds.zw, lens_sample);
  float2 angle = radius > 0 ? point.xy / radius : float2(1, 0);
  float3 rear =
      float3(angle.x * lens.x - angle.y * lens.y, angle.y * lens.x + angle.x * lens.y, c.lens.w);
  float3 d = normalize(rear - point);
  float cosine = d.z;
  // Camera and lens coordinates differ by their Z orientation.
  float3 o = point;
  d.z = -d.z;
  float z = 0;
  for (int i = int(c.elements) - 1; i >= 0; --i) {
    float4 element = data[i];
    z -= element.y;
    float t;
    float3 normal = 0;
    if (element.x == 0) {
      if (d.z >= 0)
        return 0;
      t = (z - o.z) / d.z;
    } else {
      float3 offset = o - float3(0, 0, z + element.x);
      float a = dot(d, d), b = 2 * dot(d, offset), cc = dot(offset, offset) - element.x * element.x;
      float discriminant = b * b - 4 * a * cc;
      if (discriminant < 0)
        return 0;
      float q = -.5f * (b + copysign(sqrt(discriminant), b));
      float2 roots = q == 0 ? float2(-b / (2 * a)) : float2(q / a, cc / q);
      bool near = (d.z > 0) != (element.x < 0);
      t = near ? min(roots.x, roots.y) : max(roots.x, roots.y);
      if (t < 0)
        return 0;
      normal = normalize(offset + t * d);
      if (dot(normal, -d) < 0)
        normal = -normal;
    }
    o += t * d;
    if (dot(o.xy, o.xy) > element.w * element.w)
      return 0;
    if (element.x != 0) {
      float eta = (i > 0 && data[i - 1].z != 0 ? data[i - 1].z : 1) / element.z;
      float cos_i = dot(normal, -d);
      if (cos_i < 0) {
        eta = 1 / eta;
        cos_i = -cos_i;
        normal = -normal;
      }
      float sin2 = max(0.f, 1 - cos_i * cos_i) / (eta * eta);
      if (sin2 >= 1)
        return 0;
      d = d / eta + (cos_i / eta - sqrt(1 - sin2)) * normal;
    }
  }
  o.z = -o.z;
  d.z = -d.z;
  origin = center + right * o.x + up * o.y + forward * o.z;
  direction = normalize(right * d.x + up * d.y + forward * d.z);
  float area = (bounds.z - bounds.x) * (bounds.w - bounds.y);
  float4 first = data[c.elements];
  return pow(cosine, 4.f) * area *
         (c.simple_weighting ? 1 / ((first.z - first.x) * (first.w - first.y))
                             : c.exposure_weight / (c.lens.w * c.lens.w));
}
)MSL"
