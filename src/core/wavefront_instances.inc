// Native instances share mesh objects. Export their immutable vertices once,
// then keep only transforms and placement-dependent shading data per copy.
struct MeshPrototype {
  std::vector<WFTriangle> triangles;
  std::vector<WFSurface> surfaces;
  std::vector<uint32_t> emitters;
  WFVector lo{INFINITY, INFINITY, INFINITY, 0}, hi{-INFINITY, -INFINITY, -INFINITY, 0};
  uint32_t material = 0;
  bool mixed_materials = false, dielectric = false;
};
std::unordered_map<const hitable *, uint32_t> mesh_prototypes;
std::vector<MeshPrototype> prototypes;
std::vector<WFTriangle> *prototype_triangles = nullptr;

void Mesh(const hitable &object, const hitable &faces, size_t count, const Transform &placement) {
  if (!current_placement) {
    ReserveTriangles(count);
    Geometry(faces, placement);
    return;
  }
  auto found = mesh_prototypes.find(&object);
  if (found == mesh_prototypes.end()) {
    MeshPrototype prototype;
    prototype.triangles.reserve(count);
    prototype_triangles = &prototype.triangles;
    // A fresh placement key prevents the preceding mesh's mapping cache from
    // being reused while compiling this prototype in its native coordinate frame.
    surface_mesh = nullptr;
    Geometry(faces, Transform());
    prototype_triangles = nullptr;
    surface_mesh = nullptr;
    if (prototype.triangles.empty())
      return;
    std::unordered_map<uint32_t, uint32_t> mappings;
    prototype.material = prototype.triangles.front().material;
    for (size_t i = 0; i < prototype.triangles.size(); ++i) {
      auto &triangle = prototype.triangles[i];
      uint32_t surface;
      std::memcpy(&surface, &triangle.uv2.w, sizeof(surface));
      auto mapping = mappings.emplace(surface, prototype.surfaces.size());
      if (mapping.second)
        prototype.surfaces.push_back(scene.surfaces[surface]);
      SetSurface(triangle, mapping.first->second);
      prototype.mixed_materials |= triangle.material != prototype.material;
      prototype.dielectric |= scene.materials[triangle.material].dielectric != 0;
      if (scene.materials[triangle.material].emissive) {
        triangle.light = prototype.emitters.size();
        prototype.emitters.push_back(i);
      }
      for (const auto &p : triangle.p) {
        prototype.lo.x = std::min(prototype.lo.x, p.x);
        prototype.lo.y = std::min(prototype.lo.y, p.y);
        prototype.lo.z = std::min(prototype.lo.z, p.z);
        prototype.hi.x = std::max(prototype.hi.x, p.x);
        prototype.hi.y = std::max(prototype.hi.y, p.y);
        prototype.hi.z = std::max(prototype.hi.z, p.z);
      }
    }
    found = mesh_prototypes.emplace(&object, prototypes.size()).first;
    prototypes.push_back(std::move(prototype));
  }
  const auto &prototype = prototypes[found->second];
  if (current_boundary) {
    if (prototype.mixed_materials ||
        prototype.material != scene.boundaries[current_boundary].material)
      throw std::runtime_error("mixed materials on a dielectric/SSS boundary");
  } else if (prototype.dielectric)
    throw std::runtime_error("open dielectric surface (a closed boundary is required on Metal)");
  WFMeshInstance out;
  out.mesh = found->second;
  out.boundary = current_boundary;
  out.orientation = boundary_orientation;
  const auto &forward = placement.GetMatrix(), &inverse = placement.GetInverseMatrix();
  for (int k = 0; k < 3; ++k) {
    out.forward[k] = {float(forward.m[k][0]), float(forward.m[k][1]), float(forward.m[k][2]),
                      float(forward.m[k][3])};
    out.inverse[k] = {float(inverse.m[k][0]), float(inverse.m[k][1]), float(inverse.m[k][2]),
                      float(inverse.m[k][3])};
  }
  out.surface_base = scene.surfaces.size();
  for (auto surface : prototype.surfaces) {
    Float matrix[4][4] = {{0}};
    matrix[3][3] = 1;
    for (int k = 0; k < 3; ++k) {
      const auto &row = surface.forward[k];
      matrix[k][0] = row.x;
      matrix[k][1] = row.y;
      matrix[k][2] = row.z;
      matrix[k][3] = row.w;
    }
    Transform composed = placement * Transform(matrix);
    const auto &f = composed.GetMatrix(), &v = composed.GetInverseMatrix();
    for (int k = 0; k < 3; ++k) {
      surface.forward[k] = {float(f.m[k][0]), float(f.m[k][1]), float(f.m[k][2]), float(f.m[k][3])};
      surface.inverse[k] = {float(v.m[k][0]), float(v.m[k][1]), float(v.m[k][2]), float(v.m[k][3])};
      surface.bump_forward[k] = out.forward[k];
      surface.bump_inverse[k] = out.inverse[k];
    }
    surface.bump_transformed = !placement.IsIdentity();
    scene.surfaces.push_back(surface);
  }
  if (current_boundary) {
    auto &body = scene.boundaries[current_boundary];
    for (int corner = 0; corner < 8; ++corner) {
      point3f p = placement(point3f(corner & 1 ? prototype.hi.x : prototype.lo.x,
                                    corner & 2 ? prototype.hi.y : prototype.lo.y,
                                    corner & 4 ? prototype.hi.z : prototype.lo.z));
      for (auto pair : {std::make_pair(&body.lo, &body.hi),
                        std::make_pair(&scene.boundary_lo, &scene.boundary_hi)}) {
        pair.first->x = std::min(pair.first->x, float(p[0]));
        pair.first->y = std::min(pair.first->y, float(p[1]));
        pair.first->z = std::min(pair.first->z, float(p[2]));
        pair.second->x = std::max(pair.second->x, float(p[0]));
        pair.second->y = std::max(pair.second->y, float(p[1]));
        pair.second->z = std::max(pair.second->z, float(p[2]));
      }
    }
  }
  out.light_base = scene.lights.size();
  for (uint32_t face : prototype.emitters) {
    WFLight light;
    light.type = 3;
    light.primitive = face; // Resolve the shared triangle range in FinishMeshes().
    light.instance = scene.instances.size() + 1;
    scene.lights.push_back(light);
    scene.light_sources.push_back(nullptr);
  }
  scene.instances.push_back(out);
}
