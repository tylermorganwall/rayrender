// Host construction for the finite-light proposal. This is independent of the
// geometry acceleration structure: it bounds emission, not intersections.
// A conservative sphere and approximate flux guide selection. They need not
// predict visibility or texture detail: the GPU retains 5% uniform support.
void RefitWavefrontLightBVH(WavefrontScene &scene) {
  const size_t count = scene.lights.size() - scene.environments;
  std::vector<double> powers(count);
  double maximum = 0;
  for (size_t i = 0; i < count; ++i) {
    const auto &light = scene.lights[scene.environments + i];
    auto &node = scene.light_nodes[i];
    vec3<double> lo, hi;
    double power = 0;
    if (light.type == 2) {
      lo = hi = vec3<double>(light.forward.x, light.forward.y, light.forward.z);
      power = (double(light.color.x) + light.color.y + light.color.z) / 3;
      if (light.flags & 2) power = light.projection_direction.w / (4 * M_PI);
    } else {
      uint32_t material;
      double area;
      if (light.primitive & 0x80000000u) {
        size_t id = light.primitive & 0x7fffffffu;
        const auto &a = scene.quadric_bounds[2 * id];
        const auto &b = scene.quadric_bounds[2 * id + 1];
        lo = vec3<double>(a.x, a.y, a.z);
        hi = vec3<double>(b.x, b.y, b.z);
        vec3<double> extent = hi - lo;
        // Ellipsoid/CSG area estimate; used only for importance, never a PDF.
        area = M_PI * (extent[0] * extent[1] + extent[0] * extent[2] + extent[1] * extent[2]) / 3;
        material = scene.quadrics[id].material;
      } else {
        const auto &triangle = scene.triangles[light.primitive];
        vec3<double> vertices[3];
        for (int j = 0; j < 3; ++j) {
          const auto &v = triangle.p[j];
          vertices[j] = vec3<double>(v.x, v.y, v.z);
          if (light.instance) {
            const auto &transform = scene.instances[light.instance - 1].forward;
            for (int a = 0; a < 3; ++a)
              vertices[j][a] = double(transform[a].x) * v.x + double(transform[a].y) * v.y +
                               double(transform[a].z) * v.z + transform[a].w;
          }
        }
        lo = hi = vertices[0];
        for (int j = 1; j < 3; ++j)
          for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], vertices[j][a]);
            hi[a] = std::max(hi[a], vertices[j][a]);
          }
        area = .5 * cross(vertices[1] - vertices[0], vertices[2] - vertices[0]).length();
        material = triangle.material;
      }
      power = area * scene.emission_estimates[material];
      if (light.primitive & 0x80000000u) {
        const auto &q = scene.quadrics[light.primitive & 0x7fffffffu];
        if (q.reserved) power = area * q.csg_hi.w;
      }
    }
    vec3<double> center = .5 * (lo + hi), radius = .5 * (hi - lo);
    node.center_radius2 = {float(center[0]), float(center[1]), float(center[2]),
                           float(dot(radius, radius))};
    powers[i] = std::isfinite(power) ? std::max(0., power) : 0;
    maximum = std::max(maximum, powers[i]);
  }
  // A common normalization cancels in every branch ratio and prevents flux
  // overflow on bright emitters. Zero-power trees fall back to equal branches.
  for (size_t i = 0; i < count; ++i)
    scene.light_nodes[i].power = maximum > 0 ? powers[i] / maximum : 0;
  for (size_t i = count; i < scene.light_nodes.size(); ++i) {
    auto &node = scene.light_nodes[i];
    const auto &left = scene.light_nodes[node.left];
    const auto &right = scene.light_nodes[node.right];
    const auto &a = left.center_radius2;
    const auto &b = right.center_radius2;
    vec3<double> ca(a.x, a.y, a.z), cb(b.x, b.y, b.z), center = .5 * (ca + cb);
    double radius = .5 * (cb - ca).length() + std::max(std::sqrt(a.w), std::sqrt(b.w));
    node.center_radius2 = {float(center[0]), float(center[1]), float(center[2]),
                           float(radius * radius)};
    node.power = left.power + right.power;
  }
}

void BuildWavefrontLightBVH(WavefrontScene &scene) {
  size_t count = scene.lights.size() - scene.environments;
  // Uniform selection is cheap and adequate for the common one/few-light case.
  if (count < 8)
    return;
  if (count > UINT32_MAX / 2)
    throw std::runtime_error("Metal light hierarchy capacity exceeded");
  scene.light_nodes.resize(count);
  for (size_t i = 0; i < count; ++i)
    scene.light_nodes[i].left = scene.environments + i;
  RefitWavefrontLightBVH(scene);
  std::vector<uint32_t> order(count);
  std::iota(order.begin(), order.end(), 0);
  std::function<uint32_t(size_t, size_t)> build = [&](size_t begin, size_t end) -> uint32_t {
    if (end - begin == 1)
      return order[begin];
    float lo[3] = {INFINITY, INFINITY, INFINITY}, hi[3] = {-INFINITY, -INFINITY, -INFINITY};
    for (size_t i = begin; i < end; ++i)
      for (int a = 0; a < 3; ++a) {
        float v = (&scene.light_nodes[order[i]].center_radius2.x)[a];
        lo[a] = std::min(lo[a], v);
        hi[a] = std::max(hi[a], v);
      }
    int axis = 0;
    for (int a = 1; a < 3; ++a)
      if (hi[a] - lo[a] > hi[axis] - lo[axis])
        axis = a;
    size_t middle = (begin + end) / 2;
    std::nth_element(order.begin() + begin, order.begin() + middle, order.begin() + end,
                     [&](uint32_t a, uint32_t b) {
                       float x = (&scene.light_nodes[a].center_radius2.x)[axis];
                       float y = (&scene.light_nodes[b].center_radius2.x)[axis];
                       return x == y ? a < b : x < y;
                     });
    uint32_t left = build(begin, middle), right = build(middle, end), id = scene.light_nodes.size();
    WFLightNode node;
    node.left = left;
    node.right = right;
    scene.light_nodes.push_back(node);
    scene.light_nodes[left].parent = scene.light_nodes[right].parent = id;
    return id;
  };
  build(0, count);
  RefitWavefrontLightBVH(scene);
  scene.texture_features |= 128;
}
