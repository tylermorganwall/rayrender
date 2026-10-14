#include "wavefront.h"
#include "wavefront_motion.h"
#include "../hitables/box.h"
#include "../hitables/cylinder.h"
#include "../hitables/csg.h"
#include "../hitables/curve.h"
#include "../hitables/disk.h"
#include "../hitables/ellipsoid.h"
#include "../hitables/infinite_area_light.h"
#include "../hitables/instance.h"
#include "../hitables/mesh3d.h"
#include "../hitables/plymesh.h"
#include "../hitables/raymesh.h"
#include "../hitables/rectangle.h"
#include "../hitables/sphere.h"
#include "../hitables/triangle.h"
#include "../hitables/trimesh.h"
#include "../materials/openpbr.h"
#include "../materials/texturegraph.h"
#include "../materials/ptextexture.h"
#include "../volumes/boundary.h"
#include "../volumes/subsurface.h"
#include "../volumes/lights.h"
#include "../lights/atmosphere.h"
#include "adaptivesampler.h"
#include "bvh.h"
#include "camera.h"
#include <chrono>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

namespace {
template <class T> WFVector Pack(const T &v) { return {float(v[0]), float(v[1]), float(v[2]), 0}; }
using Clock = std::chrono::steady_clock;
#include "wavefront_light_bvh.inc"
} // namespace

// This compiler is the only reader of native implementation details. It runs
// once, before any GPU command is submitted, and rejects the entire snapshot
// if a feature cannot be represented faithfully. The CPU hot path has no hooks.
class WavefrontSceneCompiler {
public:
  explicit WavefrontSceneCompiler(WavefrontScene &scene) : scene(scene) {}
  void Geometry(const hitable &object, const Transform &placement = Transform()) {
    // Triangles dominate large meshes: avoid probing every aggregate class per face.
    if (const auto *tri = dynamic_cast<const triangle *>(&object)) {
      const auto &mesh = *tri->mesh;
      int material_id = mesh.face_material_id[tri->face_number];
      WFTriangle output;
      // TriangleMesh has already applied its object transform. Instances are
      // the only transforms still outstanding; normals use inverse transpose.
      auto a = mesh.p[tri->v[0]], b = mesh.p[tri->v[1]], c = mesh.p[tri->v[2]];
      normal3f geometric = convert_to_normal3(cross(b - a, c - a));
      if (geometric.squared_length() == 0)
        return;
      if (tri->reverseOrientation ^ tri->transformSwapsHandedness)
        geometric = -geometric;
      output.geometric = Pack(unit_vector(placement(geometric)));
      // Native triangles use a flat normal for the whole face when any corner
      // lacks a normal. Do not blend a substituted corner with smooth normals.
      bool smooth = mesh.has_normals && tri->n;
      for (int k = 0; k < 3 && smooth; ++k)
        smooth = tri->n[k] >= 0 && tri->n[k] < mesh.nNormals;
      const bool placed = !placement.IsIdentity();
      const auto &normal_matrix = placement.GetInverseMatrix();
      for (int k = 0; k < 3; ++k) {
        output.p[k] = Pack(placement(mesh.p[tri->v[k]]));
        normal3f n = geometric;
        if (smooth) {
          n = mesh.n[tri->n[k]];
          if (tri->reverseOrientation)
            n = -n;
        }
        if (placed) {
          // CPU instances transform the interpolated normal. Preserve lengths
          // here so nonuniform scaling commutes with interpolation; normalizing
          // each transformed vertex separately changes the shading direction.
          normal3f transformed;
          for (int axis = 0; axis < 3; ++axis)
            transformed[axis] = normal_matrix.m[0][axis] * n[0] +
                                normal_matrix.m[1][axis] * n[1] +
                                normal_matrix.m[2][axis] * n[2];
          output.n[k] = Pack(transformed);
        } else
          output.n[k] = Pack(placement(n));
      }
      point2f uv[3];
      tri->GetUVs(uv);
      if (mesh.has_vertex_colors) {
        // Native vertex colors use barycentrics for lookup, but geometric
        // derivatives still use the mesh's original UV chart.
        output.flags |= 16;
        for (int k = 0; k < 3; ++k) {
          output.p[k].w = uv[k][0];
          output.n[k].w = uv[k][1];
        }
        uv[0] = point2f(1, 0);
        uv[1] = point2f(0, 1);
        uv[2] = point2f(0, 0);
      }
      output.uv01 = {float(uv[0][0]), float(uv[0][1]), float(uv[1][0]), float(uv[1][1])};
      output.uv2 = {float(uv[2][0]), float(uv[2][1]), 0, 0};
      uint32_t face_index=uint32_t(tri->TextureFaceIndex());
      std::memcpy(&output.uv2.z, &face_index, sizeof(face_index));
      output.flags |= 64; // A native mesh has face-local Ptex coordinates.
      if (smooth)
        output.flags |= 8;
      if (smooth && mesh.has_consistent_normals) {
        // Keep the triangle ABI unchanged: the spare geometric lane addresses
        // the three native normal-cone angles in the shared float4 data buffer.
        uint32_t normal_data = scene.texels.size();
        scene.texels.push_back({float(mesh.alpha_v[tri->n[0]]),
                                float(mesh.alpha_v[tri->n[1]]),
                                float(mesh.alpha_v[tri->n[2]]), 0});
        std::memcpy(&output.geometric.w, &normal_data, sizeof(normal_data));
        output.flags |= 128;
      }
      output.material = Material(mesh.mesh_materials[material_id].get());
      SetSurface(output, MeshSurface(*tri, material_id, placement));
      AddTriangle(output);
    } else if (const auto *body = dynamic_cast<const MediumBoundary *>(&object)) {
      const Medium *medium = body->medium.get();
      WFBoundary out;
      out.material = Material(body->geometry->mat_ptr.get());
      if (medium && medium->subsurface && !scene.materials[out.material].dielectric)
        throw std::runtime_error("non-dielectric subsurface boundary");
      if (medium) {
        out.subsurface = medium->subsurface ? (medium->subsurface_diffusion ? 2 : 1) : 3;
        out.sigma_a = Pack(medium->sigma_a);
        out.sigma_s = Pack(medium->sigma_s);
        out.sigma_a.w = medium->subsurface_roughness;
        out.sigma_s.w = medium->g;
        out.color = Pack(medium->diffusion_color);
        out.radius = Pack(medium->diffusion_radius);
        if (medium->subsurface_guided && medium->subsurface_guide_eligible && out.subsurface == 1) {
          out.color = {float(medium->subsurface_pole[0]),float(medium->subsurface_pole[1]),
                        float(medium->subsurface_pole[2]),0};
          out.radius.w = 1;
          scene.texture_features |= 256;
        }
        if (!medium->subsurface) {
          out.radius.w = body->keep_surface;
          out.medium = ExportMedium(*medium, placement * body->medium_to_world);
        }
      }
      uint32_t previous = current_boundary;
      Float previous_orientation = boundary_orientation;
      current_boundary = scene.boundaries.size();
      boundary_orientation = body->orientation_sign;
      scene.boundaries.push_back(out);
      Geometry(*body->geometry, placement);
      current_boundary = previous;
      boundary_orientation = previous_orientation;
    } else if (const auto *list = dynamic_cast<const hitable_list *>(&object)) {
      for (const auto &child : list->objects)
        Geometry(*child, placement);
    } else if (const auto *bvh = dynamic_cast<const BVHAggregate *>(&object)) {
      for (const auto &child : bvh->primitives)
        Geometry(*child, placement);
    } else if (const auto *animated = dynamic_cast<const AnimatedHitable *>(&object)) {
      if (!scene.motion) scene.motion = std::make_shared<WavefrontMotion>();
      auto begin = WavefrontMotion::sizes(scene);
      auto previous_motion = current_motion;
      auto previous_placement = current_placement;
      current_motion = scene.motion->Add(previous_motion, animated->PrimitiveToWorld, placement);
      current_placement = ++next_placement;
      Transform reference; animated->PrimitiveToWorld.Interpolate(0, &reference);
      Geometry(*animated->primitive, placement * reference);
      scene.motion->Claim(begin, WavefrontMotion::sizes(scene), current_motion);
      current_motion = previous_motion; current_placement = previous_placement;
    } else if (const auto *inst = dynamic_cast<const instance *>(&object)) {
      const auto previous_placement = current_placement;
      current_placement = ++next_placement;
      Geometry(*inst->original_scene, placement * *inst->ObjectToWorld);
      current_placement = previous_placement;
    } else if (const auto *mesh = dynamic_cast<const mesh3d *>(&object)) {
      Mesh(object, *mesh->mesh_bvh, mesh->mesh->nTriangles, placement);
    } else if (const auto *mesh = dynamic_cast<const plymesh *>(&object)) {
      Mesh(object, *mesh->ply_mesh_bvh, mesh->mesh->nTriangles, placement);
    } else if (const auto *mesh = dynamic_cast<const raymesh *>(&object)) {
      Mesh(object, *mesh->tri_mesh_bvh, mesh->mesh->nTriangles, placement);
    } else if (const auto *mesh = dynamic_cast<const trimesh *>(&object)) {
      Mesh(object, *mesh->tri_mesh_bvh, mesh->mesh->nTriangles, placement);
    } else if (const auto *shape = dynamic_cast<const box *>(&object)) {
      Geometry(shape->list, placement);
    } else if (const auto *rect = dynamic_cast<const xy_rect *>(&object)) {
      Rectangle(*rect, placement, {rect->x0, rect->y0, rect->k}, {rect->x1, rect->y0, rect->k},
                {rect->x1, rect->y1, rect->k}, {rect->x0, rect->y1, rect->k}, {0, 0, 1});
    } else if (const auto *rect = dynamic_cast<const xz_rect *>(&object)) {
      Rectangle(*rect, placement, {rect->x0, rect->k, rect->z0}, {rect->x1, rect->k, rect->z0},
                {rect->x1, rect->k, rect->z1}, {rect->x0, rect->k, rect->z1}, {0, 1, 0});
    } else if (const auto *rect = dynamic_cast<const yz_rect *>(&object)) {
      Rectangle(*rect, placement, {rect->k, rect->y0, rect->z0}, {rect->k, rect->y1, rect->z0},
                {rect->k, rect->y1, rect->z1}, {rect->k, rect->y0, rect->z1}, {1, 0, 0});
    } else if (const auto *s = dynamic_cast<const sphere *>(&object)) {
      Quadric(object, placement, vec3f(s->radius), s->alpha_mask.get(), s->bump_tex.get(), 1);
    } else if (const auto *s = dynamic_cast<const ellipsoid *>(&object)) {
      Quadric(object, placement, s->axes, s->alpha_mask.get(), s->bump_tex.get(), 2);
    } else if (const auto *s = dynamic_cast<const csg *>(&object)) {
      Quadric(object, placement, vec3f(1), nullptr, nullptr, 0, s);
    } else if (const auto *s = dynamic_cast<const curve *>(&object)) {
      Quadric(object, placement, vec3f(1), nullptr, nullptr, 0, nullptr, s);
    } else if (const auto *light = dynamic_cast<const InfiniteAreaLight *>(&object)) {
      if (light->light) {
        scene.environment_versions.emplace_back(light, light->light.get());
        scene.owned_lights.push_back(light->light);
        Infinite(*light->light);
      } else {
        const auto *emitter = dynamic_cast<const diffuse_light *>(light->mat_ptr.get());
        if (!emitter || emitter->intensity != 1)
          throw std::runtime_error("unrecognized legacy background");
        auto source =
            std::make_shared<ImageInfiniteLight>(emitter->emit, light->width, light->height, 0);
        source->SetEnvironmentTransform(light->ObjectToWorld, light->WorldToObject);
        scene.owned_lights.push_back(source);
        Infinite(*source);
      }
    } else {
      auto triangles = TriangulateWavefrontPrimitive(object, placement);
      if (triangles.empty())
        throw std::runtime_error("geometry '" + object.GetName() + "'");
      uint32_t material = Material(object.mat_ptr.get());
      const alpha_texture *alpha = nullptr;
      const bump_texture *bump = nullptr;
      WFVector axes;
      if (const auto *s = dynamic_cast<const cylinder *>(&object)) {
        alpha = s->alpha_mask.get();
        bump = s->bump_tex.get();
      } else if (const auto *s = dynamic_cast<const disk *>(&object)) {
        alpha = s->alpha_mask.get();
        bump = s->bump_tex.get();
      }
      uint32_t surface = Surface(placement * *object.ObjectToWorld, alpha, bump, axes, placement);
      ReserveTriangles(triangles.size());
      for (auto &tri : triangles) {
        tri.material = material;
        SetSurface(tri, surface);
        AddTriangle(tri);
      }
      ++tessellated_primitives;
    }
  }

  size_t tessellated_primitives = 0;

  void FinishMeshes() {
    FinishPtex();
    for (auto &q : scene.quadrics) if (q.curve_offset) {
      size_t offset = scene.quadrics.size() * sizeof(WFQuadric) +
                      scene.csg_instructions.size() * sizeof(WFCsgInstruction) +
                      size_t(q.curve_offset - 1) * sizeof(WFCurve);
      if (offset > UINT32_MAX)
        throw std::runtime_error("Metal curve data exceeds 32-bit offsets");
      q.curve_offset = uint32_t(offset);
    }
    for (auto &q : scene.quadrics) if (q.csg_count) {
      size_t offset=scene.quadrics.size()*sizeof(WFQuadric)+size_t(q.csg_offset)*sizeof(WFCsgInstruction);
      if(offset>UINT32_MAX) throw std::runtime_error("Metal CSG program exceeds 32-bit offsets");
      q.csg_offset=uint32_t(offset);
    }
    scene.loose_triangles = scene.triangles.size();
    // Identity descriptors precede mesh placements in the TLAS. Their count=0
    // distinguishes flat triangles and baked-world quadrics from shared meshes.
    size_t prefix = (scene.loose_triangles != 0) + !scene.quadrics.empty();
    if (scene.instances.empty())
      return;
    for (auto &prototype : prototypes) {
      WFMeshRange range{uint32_t(scene.triangles.size()), uint32_t(prototype.triangles.size())};
      ReserveTriangles(range.count);
      scene.triangles.insert(scene.triangles.end(), prototype.triangles.begin(),
                             prototype.triangles.end());
      scene.meshes.push_back(range);
    }
    for (auto &instance : scene.instances) {
      instance.triangle_start = scene.meshes[instance.mesh].start;
      instance.triangle_count = scene.meshes[instance.mesh].count;
    }
    for (auto &light : scene.lights)
      if (light.instance) {
        light.primitive += scene.instances[light.instance - 1].triangle_start;
        light.instance += prefix;
      }
    scene.instances.insert(scene.instances.begin(), prefix, WFMeshInstance());
    if (scene.motion) scene.motion->PrefixInstances(prefix);
    prototypes.clear(); // Release staging vertices before allocating GPU resources.
  }

  void DiffusionTables() {

    std::unordered_map<double, uint32_t> tables;
    for (auto &body : scene.boundaries) {
      if (body.subsurface != 2)
        continue;
      scene.boundary_lo.w = 1; // Enable coincident-interface classification on the GPU.
      body.diffusion_offset = scene.diffusion_pairs.size();
      const double inside = scene.materials[body.material].optical.x;
      // Slot zero is vacuum; the remaining slots are keyed by material+1.
      for (size_t outside = 0; outside <= scene.materials.size(); ++outside) {
        const double eta = inside / (outside ? scene.materials[outside - 1].optical.x : 1);
        auto found = tables.find(eta);
        if (found == tables.end()) {
          auto source = ExportDiffusionExitTable(eta);
          WFDiffusionExit table;
          table.shape = {float(eta), float(source.critical2), float(source.support),
                         float(source.mass)};
          table.normalization.x = source.normalization;
          table.offset = scene.diffusion_samples.size();
          for (const auto &v : source.samples)
            scene.diffusion_samples.push_back({v[0], v[1], 0, 0});
          found = tables.emplace(eta, scene.diffusion_exits.size()).first;
          scene.diffusion_exits.push_back(table);
        }
        scene.diffusion_pairs.push_back(found->second);
      }
    }
  }

  void Points(const PointLightSet &set) {
    for (const auto &light : set.lights) {
      WFLight out;
      out.type = 2;
      out.forward = Pack(light.position);
      out.right = Pack(light.direction);
      out.color = Pack(light.intensity);
      out.right.w = light.cos_outer;
      out.up.w = light.cos_inner;
      out.flags = light.spot;
      scene.lights.push_back(out);
      scene.light_sources.push_back(nullptr);
    }
  }

  static void PrepareAtmosphere(WavefrontScene &scene) {
#ifdef RAY_HAS_METAL_BVH
    size_t count = 0;
    for (const auto &light : scene.lights) count += (light.flags & 128) != 0;
    const auto *sky = dynamic_cast<const PragueInfiniteLight *>(scene.atmosphere);
    if (count != 1 || !sky) return;
    scene.atmosphere_snapshot = sky->model->metalSnapshot();
    scene.atmosphere_data.model = scene.atmosphere_snapshot->layout();
    scene.gpu_atmosphere = true;
#endif
  }
  static void RefreshLights(WavefrontScene &scene) {
    for (size_t i = 0; i < scene.lights.size(); ++i) {
      auto &out = scene.lights[i];
      const auto *source = scene.light_sources[i];
      if (const auto *disk = dynamic_cast<const DiskInfiniteLight *>(source)) {
        auto transform = [disk](const vec3<double> &v) {
          return Pack(disk->WorldDirection(vec3f(v[0], v[1], v[2])));
        };
        out.forward = transform(disk->forward);
        out.forward.w = disk->one_minus_cos_radius;
        out.right = transform(disk->right);
        out.up = transform(disk->up);
        // Horizon clipping is expressed in the environment frame.
        out.color = Pack(disk->WorldDirection(vec3f(0, 1, 0)));
        out.color.w = uint32_t(source->RadianceSpectrum());
      } else if (const auto *image = dynamic_cast<const ImageInfiniteLight *>(source)) {
        out.forward = Pack(image->WorldDirection(image->light_to_world(vec3f(1, 0, 0))));
        out.right = Pack(image->WorldDirection(image->light_to_world(vec3f(0, 1, 0))));
        out.up = Pack(image->WorldDirection(image->light_to_world(vec3f(0, 0, 1))));
      } else if (const auto *sky = dynamic_cast<const PragueInfiniteLight *>(source)) {
        sky->ExportMetal(scene.atmosphere_data);
        using V = skymodelr::PragueSkyModel::Vector3;
        vec3f sun = sky->World(V(std::cos(sky->azimuth) * std::cos(sky->elevation),
                                std::sin(sky->azimuth) * std::cos(sky->elevation),
                                std::sin(sky->elevation)));
        vec3f right = unit_vector(cross(std::abs(sun[0]) > .9 ? vec3f(0,1,0) : vec3f(1,0,0), sun));
        out.forward = Pack(sun);
        out.forward.w = sky->sun_fraction;
        out.right = Pack(right);
        out.right.w = 2 * std::pow(std::sin(sky->sun_radius / 2), 2);
        out.up = Pack(cross(sun, right));
      }
    }
  }

private:
#include "wavefront_medium_export.inc"
  WavefrontScene &scene;
  std::unordered_map<const material *, uint32_t> materials;
  std::unordered_map<const texture *, uint32_t> textures;
  uint32_t legacy_noise_offset = UINT32_MAX;
  uint32_t current_boundary = 0;
  uint32_t current_motion = 0;
  Float boundary_orientation = 1;
  std::unordered_map<uint32_t, uint32_t> open_dielectrics;

  // Open interfaces have no containment volume. Share their crossing identity
  // by material, like the CPU dielectric stack, rather than by triangle. Their
  // membership persists across opaque reflections until a matching exit.
  uint32_t OpenDielectric(uint32_t material) {
    if (!scene.materials[material].dielectric)
      return 0;
    auto found = open_dielectrics.find(material);
    if (found != open_dielectrics.end())
      return found->second;
    WFBoundary boundary;
    boundary.material = material;
    boundary.lo.w = 1; // Crossing-only; never use geometric containment queries.
    scene.texture_features |= 2048; // Infer incident glass at unpaired back faces.
    uint32_t id = scene.boundaries.size();
    scene.boundaries.push_back(boundary);
    open_dielectrics.emplace(material, id);
    return id;
  }

#include "wavefront_instances.inc"

  void ReserveTriangles(size_t additional) {
    const size_t required = scene.triangles.size() + additional;
    if (required > UINT32_MAX / 10)
      throw std::runtime_error("Metal vertex index capacity exceeded");
    // Retain geometric growth across many small meshes/instances. Reserving
    // exactly size+count at every mesh would turn exporting them quadratic.
    if (required > scene.triangles.capacity())
      scene.triangles.reserve(std::max(required, 2 * scene.triangles.capacity()));
  }

  uint32_t Texture(const texture *input) {
    if (!input)
      throw std::runtime_error("missing texture");
    auto found = textures.find(input);
    if (found != textures.end())
      return found->second;
    WFTexture out;
    if (const auto *t = dynamic_cast<const graph_texture *>(input)) {
      uint32_t id = Graph(t->Root());
      textures[input] = id;
      return id;
    } else if (const auto *t = dynamic_cast<const constant_texture *>(input)) {
      out.a = Pack(t->color);
    } else if (const auto *t = dynamic_cast<const checker_texture *>(input)) {
      auto *even = dynamic_cast<const constant_texture *>(t->even.get());
      auto *odd = dynamic_cast<const constant_texture *>(t->odd.get());
      if (even && odd) {
        out.type = 1;
        out.a = Pack(even->color);
        out.b = Pack(odd->color);
      } else {
        // Checker evaluation selects one child at the same coordinates. Export
        // children first so the device can follow indices without recursion.
        out.type = 14;
        out.offset = Texture(t->even.get());
        out.width = Texture(t->odd.get());
      }
      out.a.w = t->period;
    } else if (const auto *t = dynamic_cast<const gradient_texture *>(input)) {
      out.type = t->aligned_v ? 2 : 3;
      out.a = Pack(t->gamma_color1);
      out.b = Pack(t->gamma_color2);
      out.height = t->hsv;
    } else if (const auto *t = dynamic_cast<const world_gradient_texture *>(input)) {
      out.type = 13;
      out.a = Pack(t->gamma_color1);
      out.b = Pack(t->gamma_color2);
      out.c = Pack(t->point1);
      out.mapping = Pack(t->dir);
      out.mapping.w = t->inv_trans_length;
      out.height = t->hsv;
    } else if (const auto *t = dynamic_cast<const noise_texture *>(input)) {
      // Copy the native Perlin tables, not a newly seeded approximation or the
      // composable value-noise algorithm. All legacy textures share these tables.
      if (legacy_noise_offset == UINT32_MAX) {
        legacy_noise_offset = scene.texels.size();
        for (int i = 0; i < 256; ++i)
          scene.texels.push_back(Pack(perlin::ranvec[i]));
        for (int i = 0; i < 256; ++i)
          scene.texels.push_back({float(perlin::perm_x[i]), float(perlin::perm_y[i]),
                                  float(perlin::perm_z[i]), 0});
      }
      out.type = 12;
      out.offset = legacy_noise_offset;
      out.a = Pack(t->color);
      out.b = Pack(t->color2);
      out.c = {float(t->scale), float(t->phase), float(t->intensity), 0};
    } else if (const auto *t = dynamic_cast<const triangle_texture *>(input)) {
      out.type = 4;
      out.a = Pack(t->a);
      out.b = Pack(t->b);
      out.c = Pack(t->c);
    } else if (const auto *t = dynamic_cast<const image_texture_float *>(input)) {
      if (t->channels < 3)
        throw std::runtime_error("image with fewer than three channels");
      out.type = dynamic_cast<const latlong_image_texture *>(t) ? 6 : 5;
      out.mapping = {float(t->repeatu), float(t->repeatv), float(t->offsetu), float(t->offsetv)};
      out.width = t->nx;
      out.height = t->ny;
      out.offset = scene.texels.size();
      for (size_t i = 0; i < size_t(t->nx) * t->ny; ++i)
        scene.texels.push_back({float(t->data[t->channels * i] * t->intensity),
                                float(t->data[t->channels * i + 1] * t->intensity),
                                float(t->data[t->channels * i + 2] * t->intensity), 0});
    } else if (const auto *t = dynamic_cast<const image_texture_char *>(input)) {
      if (t->channels < 3)
        throw std::runtime_error("image with fewer than three channels");
      out.type = 5;
      out.mapping = {float(t->repeatu), float(t->repeatv), float(t->offsetu), float(t->offsetv)};
      out.width = t->nx;
      out.height = t->ny;
      out.offset = scene.texels.size();
      for (size_t i = 0; i < size_t(t->nx) * t->ny; ++i) {
        WFVector value;
        float *channels = &value.x;
        for (int k = 0; k < 3; ++k) {
          float c = t->data[t->channels * i + k] * t->intensity / 255;
          channels[k] = c * c;
        }
        scene.texels.push_back(value);
      }
    } else
      throw std::runtime_error("procedural or composable texture");
    if (out.type == 2 || out.type == 3 || out.type == 5)
      scene.texture_features |= 2;
    uint32_t id = scene.textures.size();
    scene.textures.push_back(out);
    textures[input] = id;
    return id;
  }

#include "wavefront_texture_export.inc"

  uint32_t Roughness(const roughness_texture &input, bool alpha, point2f repeat = point2f(1, 1)) {
    WFTexture out;
    out.type = 7; // Native roughness images use nearest-neighbor lookup.
    out.width = input.nx;
    out.height = input.ny;
    out.mapping = {float(repeat[0]), float(repeat[1]), float(input.offsetu), float(input.offsetv)};
    out.offset = scene.texels.size();
    for (size_t i = 0; i < size_t(input.nx) * input.ny; ++i) {
      float x = input.data[input.channels * i] / 255.f;
      float y = input.channels > 1 ? input.data[input.channels * i + 1] / 255.f : x;
      if (alpha) {
        x = roughness_texture::RoughnessToAlpha(x);
        y = roughness_texture::RoughnessToAlpha(y);
        x *= x;
        y *= y;
      }
      scene.texels.push_back({x, y, 0, 0});
    }
    uint32_t index = scene.textures.size();
    scene.textures.push_back(out);
    return index;
  }

  void Distribution(const MicrofacetDistribution *input, WFMaterial &out) {
    const roughness_texture *image = nullptr;
    if (const auto *ggx = dynamic_cast<const TrowbridgeReitzDistribution *>(input)) {
      out.distribution = 0;
      image = ggx->has_roughness ? ggx->roughness.get() : nullptr;
    } else if (const auto *beckmann = dynamic_cast<const BeckmannDistribution *>(input)) {
      out.distribution = 1;
      image = beckmann->has_roughness ? beckmann->roughness.get() : nullptr;
    } else
      throw std::runtime_error("unrecognized microfacet distribution");
    if (image)
      out.second_texture = Roughness(*image, true);
    if (input->roughness_graph) {
      WFTexture map;
      map.type = 9;
      map.offset = Graph(*input->roughness_graph);
      map.width = input->roughness_graph_v ? Graph(*input->roughness_graph_v) : map.offset;
      map.height = input->graph_is_alpha;
      out.second_texture = scene.textures.size();
      scene.textures.push_back(map);
    }
    const auto alphas = input->GetAlphas(0, 0);
    out.parameters.x = alphas[0];
    out.parameters.y = alphas[1];
  }

  uint32_t Material(const material *input) {
    auto found = materials.find(input);
    if (found != materials.end())
      return found->second;
    WFMaterial out;
    if (const auto *diffuse = dynamic_cast<const diffuse_material *>(input)) {
      out.texture = Texture(diffuse->albedo.get());
      out.diffuse = {float(diffuse->child.a), float(diffuse->child.b),
                     float(diffuse->child.average_loss), 0};
    } else if (const auto *light = dynamic_cast<const diffuse_light *>(input)) {
      out.texture = Texture(light->emit.get());
      out.emissive = 1;
      out.diffuse.w = light->intensity;
      // Only emission-only materials use this lane. Keep invisible emitters in
      // the light sampler; visibility depends on the path's scattering history.
      out.extra.w = light->invisible;
      if (light->invisible) scene.texture_features |= 64;
    } else if (const auto *pbr = dynamic_cast<const OpenPBRMaterial *>(input)) {
      // OpenPBR inherits dielectric for priority identity, not its BSDF.
      out.type = 7;
      out.data = scene.material_data.size();
      auto maps = pbr->ExportWavefront(scene.material_data);
      out.texture = Texture(maps.base);
      if (maps.roughness_graph)
        out.second_texture = Graph(*maps.roughness_graph);
      if (maps.roughness)
        out.second_texture = Roughness(*maps.roughness, false, maps.repeat);
      out.dielectric = pbr->is_dielectric();
      out.optical.x = pbr->ref_idx;
      out.priority = pbr->priority;
      out.emissive = pbr->EmissionEstimate() > 0 ? 2 : 0;
    } else if (const auto *conductor = dynamic_cast<const metal *>(input)) {
      out.type = 1;
      out.texture = Texture(conductor->albedo.get());
      out.eta = Pack(conductor->eta);
      out.k = Pack(conductor->k);
      out.parameters.z = conductor->fuzz;
    } else if (const auto *micro = dynamic_cast<const MicrofacetReflection *>(input)) {
      out.type = 2;
      out.texture = Texture(micro->albedo.get());
      out.eta = Pack(micro->eta);
      out.k = Pack(micro->k);
      Distribution(micro->distribution, out);
    } else if (const auto *micro = dynamic_cast<const MicrofacetTransmission *>(input)) {
      out.type = 3;
      out.texture = Texture(micro->albedo.get());
      out.eta.x = micro->eta;
      out.k = Pack(micro->k);
      Distribution(micro->distribution, out);
    } else if (const auto *plastic = dynamic_cast<const glossy *>(input)) {
      out.type = 4;
      out.texture = Texture(plastic->albedo.get());
      out.eta = Pack(plastic->Rs);
      out.k = Pack(plastic->Rd);
      Distribution(plastic->distribution, out);
    } else if (const auto *sheet = dynamic_cast<const translucent_material *>(input)) {
      out.type = 5;
      out.texture = Texture(sheet->reflection.get());
      out.second_texture = Texture(sheet->transmission.get());
    } else if (const auto *fiber = dynamic_cast<const hair *>(input)) {
      out.type = 6;
      out.eta = Pack(fiber->sigma_a);
      out.parameters = {float(fiber->eta), float(fiber->beta_m), float(fiber->beta_n),
                        float(fiber->alpha)};
      WFTexture tint;
      tint.a = Pack(fiber->get_albedo(hit_record()));
      out.texture = scene.textures.size();
      scene.textures.push_back(tint);
    } else if (const auto *glass = dynamic_cast<const dielectric *>(input)) {
      out.dielectric = 1;
      out.optical.x = glass->ref_idx;
      out.attenuation = Pack(glass->attenuation);
      out.priority = glass->priority;
      WFTexture tint;
      tint.a = Pack(glass->albedo);
      out.texture = scene.textures.size();
      scene.textures.push_back(tint);
    } else
      throw std::runtime_error("unsupported material");
    uint32_t id = scene.materials.size();
    scene.materials.push_back(out);
    scene.emission_estimates.push_back(out.emissive ? EstimateLightEmission(input, point3f(0)) : 0);
    materials[input] = id;
    return id;
  }

  void Quadric(const hitable &object, const Transform &placement, vec3f axes,
               const alpha_texture *alpha, const bump_texture *bump, float chart,
               const csg *implicit = nullptr, const curve *fiber = nullptr) {
    if (scene.quadrics.size() >= 0x7ffffffdu)
      throw std::runtime_error("Metal analytic primitive capacity exceeded");
    WFQuadric out;
    out.material = Material(object.mat_ptr.get());
    out.boundary = current_boundary;
    const Transform object_transform = placement * *object.ObjectToWorld;
    const Transform unit_transform = object_transform * Scale(axes[0], axes[1], axes[2]);
    const auto &forward = unit_transform.GetMatrix();
    const auto &inverse = unit_transform.GetInverseMatrix();
    auto mapping = Pack(axes);
    mapping.w = chart;
    out.surface = Surface(object_transform, alpha, bump, mapping, placement);
    out.measure.y = (object.reverseOrientation ^ object.transformSwapsHandedness) ? -1 : 1;
    out.measure.w = out.measure.y * (current_boundary ? boundary_orientation : 1);
    WFVector lo, hi;
    for (int k = 0; k < 3; ++k) {
      out.forward[k] = {float(forward.m[k][0]), float(forward.m[k][1]), float(forward.m[k][2]),
                        float(forward.m[k][3])};
      out.inverse[k] = {float(inverse.m[k][0]), float(inverse.m[k][1]), float(inverse.m[k][2]),
                        float(inverse.m[k][3])};
      // An affine image of the unit ball has axis extent equal to the row norm.
      // Forward and inverse matrices round independently when uploaded. Pad for
      // float arithmetic before rounding outward so grazing roots remain inside.
      const auto &row = out.forward[k];
      double extent =
          std::sqrt(double(row.x) * row.x + double(row.y) * row.y + double(row.z) * row.z);
      double padding = 8 * std::numeric_limits<float>::epsilon() * (std::abs(row.w) + extent);
      (&lo.x)[k] = std::nextafter(float(double(row.w) - extent - padding), -INFINITY);
      (&hi.x)[k] = std::nextafter(float(double(row.w) + extent + padding), INFINITY);
    }
    vec3f a(forward.m[0][0], forward.m[1][0], forward.m[2][0]);
    vec3f b(forward.m[0][1], forward.m[1][1], forward.m[2][1]);
    vec3f c(forward.m[0][2], forward.m[1][2], forward.m[2][2]);
    out.measure.x = std::abs(dot(a, cross(b, c)));
    double radius2 = a.squared_length(), tolerance = 1e-6 * radius2;
    if (std::abs(b.squared_length() - radius2) <= tolerance &&
        std::abs(c.squared_length() - radius2) <= tolerance && std::abs(dot(a, b)) <= tolerance &&
        std::abs(dot(a, c)) <= tolerance && std::abs(dot(b, c)) <= tolerance)
      out.measure.z = std::sqrt(radius2);
    if (!(out.measure.x > 0) || !std::isfinite(out.measure.x))
      throw std::runtime_error("singular or unrepresentable analytic transform");
    if (implicit) {
      scene.texture_features |= 64;
      out.csg_offset=scene.csg_instructions.size();
      ExportWavefrontCSG(*implicit->shapes,scene.csg_instructions);
      out.csg_count=scene.csg_instructions.size()-out.csg_offset;
      out.csg_lo=Pack(implicit->object_bounds.min());
      out.csg_hi=Pack(implicit->object_bounds.max());
      out.csg_lo.w=std::max(1e-7f,1e-6f*float((implicit->object_bounds.max()-implicit->object_bounds.min()).length()));
      // Bounds include the isosurface acceptance band, including flat leaves.
      for(int a=0;a<3;++a) {
        (&out.csg_lo.x)[a]-=2*out.csg_lo.w;
        (&out.csg_hi.x)[a]+=2*out.csg_lo.w;
      }
      out.measure.y=object.reverseOrientation ? -1 : 1;
      out.measure.w=out.measure.y*(current_boundary ? boundary_orientation : 1);
      WavefrontCSGBounds(out,lo,hi);
    }
    if (fiber) {
      if (current_boundary)
        throw std::runtime_error("open curve cannot bound a dielectric/SSS interior on Metal");
      if (scene.materials[out.material].dielectric)
        throw std::runtime_error("dielectric curve interfaces on Metal");
      WFCurve data = ExportWavefrontCurve(*fiber);
      if (!(data.width.z < data.width.w))
        return;
      scene.texture_features |= 512;
      out.curve_offset = uint32_t(scene.curves.size() + 1);
      scene.curves.push_back(data);
      float radius = .5f * std::max(
          (1 - data.width.z) * data.width.x + data.width.z * data.width.y,
          (1 - data.width.w) * data.width.x + data.width.w * data.width.y);
      for (int a = 0; a < 3; ++a) {
        float low = INFINITY, high = -INFINITY;
        for (const auto &p : data.cp) {
          low = std::min(low, (&p.x)[a]);
          high = std::max(high, (&p.x)[a]);
        }
        (&out.csg_lo.x)[a] = low - radius;
        (&out.csg_hi.x)[a] = high + radius;
      }
      out.measure.y = object.reverseOrientation ? -1 : 1;
      out.measure.w = out.measure.y;
      WavefrontCSGBounds(out, lo, hi);
    }
    if (current_boundary) {
      auto &body = scene.boundaries[current_boundary];
      if (body.material != out.material)
        throw std::runtime_error("mixed materials on a dielectric/SSS boundary");
      for (int k = 0; k < 3; ++k) {
        (&scene.boundary_lo.x)[k] = std::min((&scene.boundary_lo.x)[k], (&lo.x)[k]);
        (&scene.boundary_hi.x)[k] = std::max((&scene.boundary_hi.x)[k], (&hi.x)[k]);
        (&body.lo.x)[k] = std::min((&body.lo.x)[k], (&lo.x)[k]);
        (&body.hi.x)[k] = std::max((&body.hi.x)[k], (&hi.x)[k]);
      }
    } else
      out.boundary = OpenDielectric(out.material);
    if (scene.materials[out.material].emissive) {
      WFLight light;
      light.type = 3;
      light.primitive = 0x80000000u | uint32_t(scene.quadrics.size());
      out.light = scene.lights.size();
      scene.lights.push_back(light);
      scene.light_sources.push_back(nullptr);
    }
    scene.quadrics.push_back(out);
    scene.quadric_bounds.push_back(lo);
    scene.quadric_bounds.push_back(hi);
  }

  void AddTriangle(WFTriangle out) {
    if (prototype_triangles) {
      prototype_triangles->push_back(out);
      return;
    }
    // The high primitive bit distinguishes analytic shapes in GPU hit queues.
    if (scene.triangles.size() >= 0x80000000u)
      throw std::runtime_error("Metal triangle capacity exceeded");
    out.boundary = current_boundary;
    if (current_boundary) {
      auto &body = scene.boundaries[current_boundary];
      if (body.material != out.material)
        throw std::runtime_error("mixed materials on a dielectric/SSS boundary");
      out.geometric.x *= boundary_orientation;
      out.geometric.y *= boundary_orientation;
      out.geometric.z *= boundary_orientation;
      for (const auto &p : out.p) {
        auto &lo = scene.boundary_lo;
        auto &hi = scene.boundary_hi;
        lo.x = std::min(lo.x, p.x);
        lo.y = std::min(lo.y, p.y);
        lo.z = std::min(lo.z, p.z);
        hi.x = std::max(hi.x, p.x);
        hi.y = std::max(hi.y, p.y);
        hi.z = std::max(hi.z, p.z);
        body.lo.x = std::min(body.lo.x, p.x);
        body.lo.y = std::min(body.lo.y, p.y);
        body.lo.z = std::min(body.lo.z, p.z);
        body.hi.x = std::max(body.hi.x, p.x);
        body.hi.y = std::max(body.hi.y, p.y);
        body.hi.z = std::max(body.hi.z, p.z);
      }
    } else
      out.boundary = OpenDielectric(out.material);
    if (scene.materials[out.material].emissive) {
      WFLight light;
      light.type = 3;
      light.primitive = scene.triangles.size();
      out.light = scene.lights.size();
      scene.lights.push_back(light);
      scene.light_sources.push_back(nullptr);
    }
    scene.triangles.push_back(out);
  }

  template <class Rect>
  void Rectangle(const Rect &rect, const Transform &placement, point3f a, point3f b, point3f c,
                 point3f d, normal3f normal) {
    Transform transform = placement * *rect.ObjectToWorld;
    uint32_t surface =
        Surface(transform, rect.alpha_mask.get(), rect.bump_tex.get(), {}, placement);
    normal = unit_vector(transform(normal * (rect.reverseOrientation ? -1 : 1)));
    point3f points[4] = {a, b, c, d};
    int indices[6] = {0, 1, 2, 0, 2, 3};
    point2f uv[4] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    if (rect.reverseOrientation)
      for (auto &v : uv)
        v.xy.x = 1 - v[0];
    for (int t = 0; t < 2; ++t) {
      WFTriangle out;
      out.material = Material(rect.mat_ptr.get());
      out.geometric = Pack(normal);
      point2f tex[3];
      for (int k = 0; k < 3; ++k) {
        int j = indices[t * 3 + k];
        out.p[k] = Pack(transform(points[j]));
        out.n[k] = Pack(normal);
        tex[k] = uv[j];
      }
      out.uv01 = {float(tex[0][0]), float(tex[0][1]), float(tex[1][0]), float(tex[1][1])};
      out.uv2 = {float(tex[2][0]), float(tex[2][1]), 0, 0};
      SetSurface(out, surface);
      AddTriangle(out);
    }
  }

  // Alias entries share the immutable texel buffer. Rotation only changes the
  // light basis: both lookup and sampling stay in the image's local frame.
  void EnvironmentAlias(const ImageInfiniteLight &image, WFLight &light) {
    const auto texture = scene.textures[light.texture];
    if (!image.distribution || (texture.type != 5 && texture.type != 6))
      return; // Black/constant environments retain uniform sphere sampling.
    auto found = environment_aliases.find(light.texture);
    if (found != environment_aliases.end()) {
      light.alias_offset = found->second.alias_offset;
      light.alias_width = found->second.alias_width;
      light.alias_height = found->second.alias_height;
      return;
    }
    // Bound table memory while retaining enough bits for each alias decision. Larger maps aggregate
    // all source-cell mass into coarser cells; this changes the proposal only, never the radiance
    // image.
    light.alias_width = std::min(texture.width, 2048u);
    light.alias_height = std::min(texture.height, 1024u);
    size_t cells = size_t(light.alias_width) * light.alias_height;
    size_t count = 1;
    unsigned column_bits = 0;
    while (count < cells) {
      count *= 2;
      ++column_bits;
    }
    // Power-of-two columns split one 32-bit sample exactly into a column and
    // a coin. Padding has zero mass and always aliases to a real image cell.
    // Quantize thresholds to the coin's actual resolution, then derive PDFs
    // from those thresholds. Even very dark cells retain at least one outcome.
    double coin_resolution = double(uint32_t(1) << std::min(24u, 32u - column_bits));
    std::vector<double> mass(count, 0);
    for (uint32_t y = 0; y < texture.height; ++y) {
      if (y % 64 == 0)
        Rcpp::checkUserInterrupt();
      for (uint32_t x = 0; x < texture.width; ++x) {
        double u = (x + .5) / texture.width, v = (y + .5) / texture.height;
        size_t cell =
            size_t(v * light.alias_height) * light.alias_width + size_t(u * light.alias_width);
        mass[cell] += image.distribution->Pdf(vec2f(u, v));
      }
    }
    double total = 0;
    for (double value : mass)
      total += value;
    if (!(total > 0) || !std::isfinite(total))
      return;
    std::vector<double> solid_angle(light.alias_height);
    for (uint32_t y = 0; y < light.alias_height; ++y) {
      solid_angle[y] =
          2 * M_PI / light.alias_width *
          (std::cos(M_PI * y / light.alias_height) - std::cos(M_PI * (y + 1) / light.alias_height));
      for (uint32_t x = 0; x < light.alias_width; ++x) {
        size_t i = size_t(y) * light.alias_width + x;
        // Retain support for bilinear filtering, seams and arbitrarily dark cells.
        mass[i] = .999 * mass[i] / total + .001 * solid_angle[y] / (4 * M_PI);
      }
    }
    std::vector<WFVector> table(count);
    std::vector<double> scaled(count);
    std::vector<uint32_t> small, large;
    for (size_t i = 0; i < count; ++i) {
      scaled[i] = mass[i] * count;
      (scaled[i] < 1 ? small : large).push_back(i);
      table[i].x = 1;
      uint32_t self = i;
      std::memcpy(&table[i].y, &self, sizeof(self));
    }
    // Vose construction: each equally likely column keeps its own mass and
    // borrows the remainder from a larger cell. Sampling requires one lookup.
    while (!small.empty() && !large.empty()) {
      uint32_t a = small.back(), b = large.back();
      small.pop_back();
      large.pop_back();
      table[a].x =
          std::clamp(float(std::ceil(scaled[a] * coin_resolution) / coin_resolution), 0.f, 1.f);
      std::memcpy(&table[a].y, &b, sizeof(b));
      scaled[b] -= 1 - scaled[a];
      (scaled[b] < 1 ? small : large).push_back(b);
    }
    // Store the PDF of the rounded table actually sampled on the GPU, rather
    // than the ideal double-precision weights used to construct it.
    std::fill(mass.begin(), mass.end(), 0);
    for (size_t i = 0; i < count; ++i) {
      uint32_t alias;
      std::memcpy(&alias, &table[i].y, sizeof(alias));
      mass[i] += table[i].x / double(count);
      mass[alias] += (1 - double(table[i].x)) / count;
    }
    for (size_t i = 0; i < cells; ++i)
      table[i].z = mass[i] / solid_angle[i / light.alias_width];
    if (scene.texels.size() + count > UINT32_MAX)
      throw std::runtime_error("Metal environment table capacity exceeded");
    light.alias_offset = scene.texels.size();
    scene.texels.insert(scene.texels.end(), table.begin(), table.end());
    environment_aliases.emplace(light.texture, light);
  }
  std::unordered_map<uint32_t, WFLight> environment_aliases;

  void Infinite(const InfiniteLight &source) {
    if (const auto *mixture = dynamic_cast<const InfiniteLightMixture *>(&source)) {
      for (const auto &light : mixture->lights)
        Infinite(*light);
      return;
    }
    WFLight out;
    if (const auto *disk = dynamic_cast<const DiskInfiniteLight *>(&source)) {
      out.type = 1;
      out.texture = Texture(disk->image.get());
      out.flags = disk->clip_horizon;
    } else if (const auto *image = dynamic_cast<const ImageInfiniteLight *>(&source)) {
      out.type = 0;
      out.texture = Texture(image->image.get());
      EnvironmentAlias(*image, out);
    } else if (const auto *sky = dynamic_cast<const PragueInfiniteLight *>(&source)) {
      out.type = 0;
      out.flags = 128;
      scene.atmosphere = sky;
    } else
      throw std::runtime_error("atmospheric or unrecognized infinite light");
    scene.lights.push_back(out);
    scene.light_sources.push_back(&source);
  }
};

#ifndef RAY_HAS_METAL_BVH
bool MetalWavefrontAvailable() { return false; }
std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &, size_t) {
  throw std::runtime_error("Metal was not enabled in this build");
}
#endif

bool WavefrontLightsChanged(const WavefrontScene &scene) {
  for (const auto &entry : scene.environment_versions)
    if (entry.first->light.get() != entry.second) return true;
  return false;
}

std::unique_ptr<WavefrontSession> PrepareWavefront(const hitable_list &world,
                                                   const hitable_list &lights, RayCamera &camera,
                                                   size_t capacity, int sampler,
                                                   WavefrontScene &scene, WavefrontReport &report,
                                                   size_t workers) {
  auto started = Clock::now();
  report = WavefrontReport();
  report.requested = true;
  try {
    scene.workers = workers;
    scene.volume_scene = lights.volume_scene.get();
    scene.transparent_background = scene.volume_scene && scene.volume_scene->transparent_background;
    scene.world = &world;
    scene.transport_atmosphere = lights.volume_scene ? lights.volume_scene->atmosphere : nullptr;
    if (sampler == 1)
      throw std::runtime_error("stratified sampling");
    if (capacity > UINT32_MAX)
      throw std::runtime_error("more than 2^32 image pixels");
    WavefrontSceneCompiler compiler(scene);
    scene.boundary_lo = {INFINITY, INFINITY, INFINITY, 0};
    scene.boundary_hi = {-INFINITY, -INFINITY, -INFINITY, 0};
    compiler.Geometry(world);
    compiler.FinishMeshes();
    compiler.DiffusionTables();
    if (lights.volume_scene)
      compiler.Points(lights.volume_scene->point_lights);
    // Misses only evaluate infinite emitters. Emissive tessellation can create
    // tens of thousands of finite lights; scanning them on every miss makes
    // otherwise simple scenes device-bound. Keep all direct-light indices valid
    // while grouping the typically tiny infinite-light prefix once on the host.
    auto unordered_lights = std::move(scene.lights);
    auto unordered_sources = std::move(scene.light_sources);
    for (bool infinite : {true, false})
      for (size_t i = 0; i < unordered_lights.size(); ++i)
        if ((unordered_lights[i].type < 2) == infinite) {
          const auto &light = unordered_lights[i];
          if (light.type == 3) {
            if (light.instance) {
              auto &instance = scene.instances[light.instance - 1];
              instance.light_base = scene.lights.size() - scene.triangles[light.primitive].light;
            } else if (light.primitive & 0x80000000u)
              scene.quadrics[light.primitive & 0x7fffffffu].light = scene.lights.size();
            else
              scene.triangles[light.primitive].light = scene.lights.size();
          }
          scene.lights.push_back(light);
          scene.light_sources.push_back(unordered_sources[i]);
          scene.environments += infinite;
        }
    WavefrontSceneCompiler::RefreshLights(scene);
    BuildWavefrontLightBVH(scene);
    WavefrontSceneCompiler::PrepareAtmosphere(scene);
    WFCamera camera_data;
    ExportWavefrontCamera(camera, camera_data, scene.camera_data);
    if (scene.motion) scene.motion->Seal(scene);
    if (scene.triangles.empty() && scene.quadrics.empty())
      throw std::runtime_error("scene without renderable geometry");
    auto session = MakeMetalWavefront(scene, capacity);
    report.used = true;
    report.triangles = scene.triangles.size();
    report.analytic_primitives = scene.quadrics.size();
    report.unique_meshes = scene.meshes.size();
    for (const auto &instance : scene.instances)
      report.mesh_instances += instance.triangle_count != 0;
    report.tessellated_primitives = compiler.tessellated_primitives;
    report.upload_seconds = std::chrono::duration<double>(Clock::now() - started).count();
    return session;
  } catch (const std::exception &error) {
    report.fallback = error.what();
    scene = WavefrontScene();
    Rcpp::warning("Metal wavefront: %s; using the CPU NEE integrator.", report.fallback.c_str());
    return nullptr;
  }
}

void WavefrontContext::BeginFrame(const hitable_list &new_world, const hitable_list &new_lights,
                                 RayCamera &camera, size_t capacity, int sampler,
                                 WavefrontReport &report, size_t workers) {
  bool reusable = session && world == &new_world && lights == &new_lights &&
                  capacity <= pixel_capacity && sampler == sampling_method &&
                  !WavefrontLightsChanged(scene);
  if (reusable) {
    // Lens tables share a packed buffer with lights and atmosphere. A larger
    // camera payload needs a new allocation; camera values and smaller payloads
    // can be uploaded in place by Render(). Projection changes alone are safe.
    try {
      WFCamera descriptor{};
      ExportWavefrontCamera(camera, descriptor, scene.camera_data);
      reusable = scene.camera_data.size() <= camera_capacity;
    } catch (const std::exception &) {
      reusable = false; // PrepareWavefront owns the normal CPU-fallback policy.
    }
  }
  if (reusable) {
    report = setup_report; // Reset per-frame counters, not resident resources.
    report.upload_seconds = 0;
    report.scene_reused = true;
    return;
  }
  // Release the device session before replacing data it references. The same
  // ordering is guaranteed by member destruction at animation exit/cancel.
  session.reset();
  scene = WavefrontScene();
  session = PrepareWavefront(new_world, new_lights, camera, capacity, sampler,
                             scene, report, workers);
  world = &new_world;
  lights = &new_lights;
  pixel_capacity = capacity;
  camera_capacity = scene.camera_data.size();
  sampling_method = sampler;
  setup_report = report;
}

bool RenderWavefrontSample(WavefrontSession &session, WavefrontScene &scene, RayCamera &camera,
                           adaptive_sampler &film, size_t sample, int sampler, uint32_t seed,
                           size_t depth, size_t roulette, float clamp,
                           const std::function<bool()> &cancelled, WavefrontReport &report) {
  auto started = Clock::now();
  WFParameters p;
  p.width = film.nx;
  p.height = film.ny;
  p.sample = sample;
  p.sampler = sampler;
  p.seed = seed;
  p.roulette = roulette;
  p.max_depth = depth;
  p.lights = scene.lights.size();
  p.environments = scene.environments;
  p.boundaries = scene.boundaries.size() - 1;
  p.boundary_lo = scene.boundary_lo;
  p.boundary_hi = scene.boundary_hi;
  ExportWavefrontCamera(camera, p.camera, scene.camera_data);
  p.camera.padding[1] = scene.transparent_background ? 4 : 0;
  // Geometry stays at its reference pose on the host. Its light tree remains
  // a valid positive proposal; device emitter evaluation and PDFs use ray time.
  WavefrontSceneCompiler::RefreshLights(scene);
  std::vector<uint32_t> active;
  for (const auto &block : film.pixel_chunks)
    for (int x = block.startx; x < block.endx; ++x)
      for (int y = block.starty; y < block.endy; ++y)
        active.push_back(y + film.ny * x);
  p.active = active.size();
  p.camera.up.w = std::max(.125f, float(1 / std::sqrt(double(film.ns))));
  const WFPixel *pixels = nullptr;
  size_t previous_discarded = report.discarded_paths;
  size_t previous_csg_failures = report.csg_failed_queries;
  if (!session.Render(p, active, scene.lights, cancelled, pixels, report))
    return false;
  if (!previous_discarded && report.discarded_paths)
    Rcpp::warning("Metal wavefront terminated rays with unresolvable positions, nonfinite weights, "
                  "or more than 16 "
                  "overlapping boundaries; rendering continues. See wavefront$discarded_paths.");
  if (!previous_csg_failures && report.csg_failed_queries)
    Rcpp::warning("Metal skipped CSG intersection queries after nonfinite distances, stalled "
                  "progress, or the 8192-step limit. Rendering continues; "
                  "see wavefront$csg_failed_queries.");
  for (uint32_t index : active) {
    size_t x = index / film.ny, y = index % film.ny;
    const auto &value = pixels[index];
    point3f color(value.radiance.x, value.radiance.y, value.radiance.z);
    Float floor = scene.transport_atmosphere ? -std::numeric_limits<Float>::infinity() : 0;
    color = clamp_point(de_nan(color), floor, clamp) * camera.get_iso();
    color *= value.normal.w;
    film.add_color_main(x, y, color);
    if (sample % 2 == 0)
      film.add_color_sec(x, y, color);
    film.add_albedo(x, y, point3f(value.albedo.x, value.albedo.y, value.albedo.z));
    film.add_normal(x, y, normal3f(value.normal.x, value.normal.y, value.normal.z));
    film.add_alpha_count(x, y, value.radiance.w);
  }
  report.sample_seconds += std::chrono::duration<double>(Clock::now() - started).count();
  ++report.completed_samples;
  return true;
}
