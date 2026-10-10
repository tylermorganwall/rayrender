// Cold graph compilation: postorder instructions reuse dead child registers.
// This bounds per-thread GPU storage by expression depth, not graph size. Shared
// image pixels are uploaded once; repeated subexpressions remain immutable.
std::unordered_map<const TextureNode *, uint32_t> graph_roots;
std::unordered_map<const DecodedTextureImage *, uint32_t> graph_images;
std::unordered_map<const alpha_texture *, uint32_t> alpha_maps;
std::unordered_map<const bump_texture *, uint32_t> bump_maps;

const TriangleMesh *surface_mesh = nullptr;
const Transform *surface_object_transform = nullptr;
const alpha_texture *surface_alpha = nullptr;
const bump_texture *surface_bump = nullptr;
uint64_t current_placement = 0, next_placement = 0, surface_placement = 0;
uint32_t mesh_surface = 0;

uint32_t MeshSurface(const triangle &face, int material, const Transform &placement) {
  const auto *alpha = face.mesh->alpha_textures[material].get();
  const auto *bump = face.mesh->bump_textures[material].get();
  // A mesh's faces share their object transform and maps. Geometry assigns a
  // scoped ID to each instance placement: comparing matrices per face is costly,
  // and nested instances can reuse the same temporary Transform address.
  // Key maps rather than material IDs so per-face vertex colors also share it.
  if (surface_mesh != face.mesh || surface_object_transform != face.ObjectToWorld ||
      surface_alpha != alpha || surface_bump != bump || surface_placement != current_placement) {
    mesh_surface = Surface(placement * *face.ObjectToWorld, alpha, bump, {}, placement);
    surface_mesh = face.mesh;
    surface_object_transform = face.ObjectToWorld;
    surface_alpha = alpha;
    surface_bump = bump;
    surface_placement = current_placement;
  }
  return mesh_surface;
}

uint32_t Graph(const TextureNode &root) {
  auto found = graph_roots.find(&root);
  if (found != graph_roots.end())
    return found->second;
  std::vector<WFTexture> program;
  uint32_t next = 0, peak = 0;
  std::function<uint32_t(const TextureNode &)> emit = [&](const TextureNode &node) {
    const auto description = node.Describe();
    WFTexture instruction;
    uint32_t first = next, children[3]{};
    for (int k = 0; k < 3; ++k)
      if (description.children[k])
        children[k] = emit(*description.children[k]);
    next = first + 1;
    peak = std::max(peak, next);
    if (peak > 128 || program.size() >= 65536)
      throw std::runtime_error(
          "texture graph exceeds Metal's 128 live values or 65536 instructions");
    instruction.type = 16 + description.operation + (description.space << 8) +
                       (uint32_t(description.flag) << 10) + (uint32_t(node.scalar) << 11);
    instruction.a = Pack(description.value);
    instruction.a.w = description.sine;
    instruction.b = Pack(description.scale);
    instruction.b.w = description.cosine;
    instruction.c = Pack(description.offset);
    instruction.mapping = {float(children[0]), float(children[1]), float(children[2]),
                           float(first)};
    instruction.width = description.argument;
    instruction.offset = description.seed;
    if (description.image) {
      const auto &image = *description.image;
      auto uploaded = graph_images.find(&image);
      if (uploaded == graph_images.end()) {
        instruction.offset = scene.texels.size();
        graph_images[&image] = instruction.offset;
        for (const auto &pixel : image.pixels)
          scene.texels.push_back(Pack(pixel));
      } else
        instruction.offset = uploaded->second;
      instruction.width = image.width;
      instruction.height = image.height;
    }
    program.push_back(instruction);
    return first;
  };
  uint32_t result = emit(root);
  WFTexture entry;
  entry.type = 8;
  entry.offset = scene.textures.size();
  entry.width = program.size();
  entry.height = result;
  scene.textures.insert(scene.textures.end(), program.begin(), program.end());
  uint32_t index = scene.textures.size();
  scene.textures.push_back(entry);
  scene.texture_features |= 3;
  graph_roots[&root] = index;
  return index;
}

uint32_t Alpha(const alpha_texture *input) {
  if (!input)
    return UINT32_MAX;
  auto found = alpha_maps.find(input);
  if (found != alpha_maps.end())
    return found->second;
  WFTexture out;
  out.a = {float(input->opacity), 0, 0, 0};
  if (input->data) {
    out.type = 7;
    out.mapping = {float(input->repeatu), float(input->repeatv), float(input->offsetu),
                   float(input->offsetv)};
    out.width = input->nx;
    out.height = input->ny;
    out.offset = scene.texels.size();
    for (size_t i = 0; i < size_t(input->nx) * input->ny; ++i)
      scene.texels.push_back(
          {float(input->opacity * input->data[(i + 1) * input->channels - 1] / 255), 0, 0, 0});
  }
  uint32_t index = scene.textures.size();
  scene.textures.push_back(out);
  alpha_maps[input] = index;
  scene.texture_features |= 10;
  return index;
}

uint32_t Bump(const bump_texture *input) {
  if (!input)
    return UINT32_MAX;
  auto found = bump_maps.find(input);
  if (found != bump_maps.end())
    return found->second;
  uint32_t index;
  if (input->height_texture)
    index = Texture(input->height_texture.get());
  else {
    WFTexture out;
    out.type = 10;
    out.offset = scene.textures.size();
    out.width = input->image ? input->image->levels.size() : 1;
    out.mapping = {float(input->repeatu), float(input->repeatv), float(input->offsetu),
                   float(input->offsetv)};
    for (uint32_t level = 0; level < out.width; ++level) {
      WFTexture mip;
      mip.type = 11;
      mip.mapping = out.mapping;
      const HeightLevel *pixels = input->image ? &input->image->levels[level] : nullptr;
      mip.width = pixels ? pixels->width : input->nx;
      mip.height = pixels ? pixels->height : input->ny;
      mip.offset = scene.texels.size();
      for (size_t i = 0; i < size_t(mip.width) * mip.height; ++i)
        scene.texels.push_back(
            {float(pixels ? pixels->pixels[i] : input->data[i * input->channels] / 255.f), 0, 0,
             0});
      scene.textures.push_back(mip);
    }
    index = scene.textures.size();
    scene.textures.push_back(out);
  }
  bump_maps[input] = index;
  scene.texture_features |= 6;
  return index;
}

uint32_t Surface(const Transform &transform, const alpha_texture *alpha, const bump_texture *bump,
                 WFVector axes = {}, const Transform &bump_placement = Transform()) {
  WFSurface out;
  const auto inverse = Inverse(transform).GetMatrix();
  const auto forward = transform.GetMatrix();
  const auto &bump_inverse = bump_placement.GetInverseMatrix();
  const auto &bump_forward = bump_placement.GetMatrix();
  for (int k = 0; k < 3; ++k) {
    out.inverse[k] = {float(inverse.m[k][0]), float(inverse.m[k][1]), float(inverse.m[k][2]),
                      float(inverse.m[k][3])};
    out.forward[k] = {float(forward.m[k][0]), float(forward.m[k][1]), float(forward.m[k][2]),
                      float(forward.m[k][3])};
    out.bump_inverse[k] = {float(bump_inverse.m[k][0]), float(bump_inverse.m[k][1]),
                           float(bump_inverse.m[k][2]), float(bump_inverse.m[k][3])};
    out.bump_forward[k] = {float(bump_forward.m[k][0]), float(bump_forward.m[k][1]),
                           float(bump_forward.m[k][2]), float(bump_forward.m[k][3])};
  }
  out.axes = axes;
  out.bump_transformed = !bump_placement.IsIdentity();
  out.alpha = Alpha(alpha);
  out.bump = Bump(bump);
  if (bump)
    out.bump_parameters.x = bump->intensity;
  // Consecutive faces of a mesh share this record without a per-face hash table.
  if (scene.surfaces.size() > 1 && std::memcmp(&out, &scene.surfaces.back(), sizeof(out)) == 0)
    return scene.surfaces.size() - 1;
  scene.surfaces.push_back(out);
  return scene.surfaces.size() - 1;
}
static void SetSurface(WFTriangle &triangle, uint32_t index) {
  std::memcpy(&triangle.uv2.w, &index, sizeof(index));
}
