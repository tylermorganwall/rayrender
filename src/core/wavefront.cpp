#include "wavefront.h"
#include "camera.h"
#include "adaptivesampler.h"
#include "bvh.h"
#include "../hitables/triangle.h"
#include "../hitables/mesh3d.h"
#include "../hitables/plymesh.h"
#include "../hitables/trimesh.h"
#include "../hitables/raymesh.h"
#include "../hitables/box.h"
#include "../hitables/rectangle.h"
#include "../hitables/instance.h"
#include "../hitables/infinite_area_light.h"
#include "../volumes/boundary.h"
#include <chrono>
#include <unordered_map>
#include <stdexcept>

namespace {
template <class T> WFVector Pack(const T &v) { return {float(v[0]), float(v[1]), float(v[2]), 0}; }
using Clock = std::chrono::steady_clock;
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
      if (mesh.alpha_textures[material_id] || mesh.bump_textures[material_id])
        throw std::runtime_error("mesh alpha or bump mapping");
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
      for (int k = 0; k < 3; ++k) {
        output.p[k] = Pack(placement(mesh.p[tri->v[k]]));
        normal3f n = geometric;
        if (smooth) {
          n = mesh.n[tri->n[k]];
          if (tri->reverseOrientation)
            n = -n;
        }
        output.n[k] = Pack(placement(n));
      }
      point2f uv[3];
      tri->GetUVs(uv);
      if (mesh.has_vertex_colors) {
        uv[0] = point2f(1, 0);
        uv[1] = point2f(0, 1);
        uv[2] = point2f(0, 0);
      }
      output.uv01 = {float(uv[0][0]), float(uv[0][1]), float(uv[1][0]), float(uv[1][1])};
      output.uv2 = {float(uv[2][0]), float(uv[2][1]), 0, 0};
      output.material = Material(mesh.mesh_materials[material_id].get());
      AddTriangle(output);
    } else if (const auto *list = dynamic_cast<const hitable_list *>(&object)) {
      for (const auto &child : list->objects)
        Geometry(*child, placement);
    } else if (const auto *bvh = dynamic_cast<const BVHAggregate *>(&object)) {
      for (const auto &child : bvh->primitives)
        Geometry(*child, placement);
    } else if (const auto *inst = dynamic_cast<const instance *>(&object)) {
      Geometry(*inst->original_scene, placement * *inst->ObjectToWorld);
    } else if (const auto *mesh = dynamic_cast<const mesh3d *>(&object)) {
      ReserveTriangles(mesh->mesh->nTriangles);
      Geometry(*mesh->mesh_bvh, placement);
    } else if (const auto *mesh = dynamic_cast<const plymesh *>(&object)) {
      ReserveTriangles(mesh->mesh->nTriangles);
      Geometry(*mesh->ply_mesh_bvh, placement);
    } else if (const auto *mesh = dynamic_cast<const raymesh *>(&object)) {
      ReserveTriangles(mesh->mesh->nTriangles);
      Geometry(*mesh->tri_mesh_bvh, placement);
    } else if (const auto *mesh = dynamic_cast<const trimesh *>(&object)) {
      ReserveTriangles(mesh->mesh->nTriangles);
      Geometry(*mesh->tri_mesh_bvh, placement);
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
    } else if (const auto *light = dynamic_cast<const InfiniteAreaLight *>(&object)) {
      if (light->light)
        Infinite(*light->light);
      else {
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
      throw std::runtime_error("geometry '" + object.GetName() + "'");
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
      } else if (const auto *image = dynamic_cast<const ImageInfiniteLight *>(source)) {
        out.forward = Pack(image->WorldDirection(image->light_to_world(vec3f(1, 0, 0))));
        out.right = Pack(image->WorldDirection(image->light_to_world(vec3f(0, 1, 0))));
        out.up = Pack(image->WorldDirection(image->light_to_world(vec3f(0, 0, 1))));
      }
    }
  }

private:
  WavefrontScene &scene;
  std::unordered_map<const material *, uint32_t> materials;
  std::unordered_map<const texture *, uint32_t> textures;

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
    if (const auto *t = dynamic_cast<const constant_texture *>(input)) {
      out.a = Pack(t->color);
    } else if (const auto *t = dynamic_cast<const checker_texture *>(input)) {
      auto *even = dynamic_cast<const constant_texture *>(t->even.get());
      auto *odd = dynamic_cast<const constant_texture *>(t->odd.get());
      if (!even || !odd)
        throw std::runtime_error("nested checker texture");
      out.type = 1;
      out.a = Pack(even->color);
      out.b = Pack(odd->color);
      out.a.w = t->period;
    } else if (const auto *t = dynamic_cast<const gradient_texture *>(input)) {
      if (t->hsv)
        throw std::runtime_error("HSV gradient texture");
      out.type = t->aligned_v ? 2 : 3;
      out.a = Pack(t->gamma_color1);
      out.b = Pack(t->gamma_color2);
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
    uint32_t id = scene.textures.size();
    scene.textures.push_back(out);
    textures[input] = id;
    return id;
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
      if (light->invisible)
        throw std::runtime_error("invisible area light");
      out.texture = Texture(light->emit.get());
      out.emissive = 1;
      out.diffuse.w = light->intensity;
    } else
      throw std::runtime_error("non-diffuse material");
    uint32_t id = scene.materials.size();
    scene.materials.push_back(out);
    materials[input] = id;
    return id;
  }

  void AddTriangle(WFTriangle out) {
    if (scene.triangles.size() >= UINT32_MAX)
      throw std::runtime_error("more than 2^32 triangles");
    if (scene.materials[out.material].emissive) {
      WFLight light;
      light.type = 3;
      light.triangle = scene.triangles.size();
      out.light = scene.lights.size();
      scene.lights.push_back(light);
      scene.light_sources.push_back(nullptr);
    }
    scene.triangles.push_back(out);
  }

  template <class Rect>
  void Rectangle(const Rect &rect, const Transform &placement, point3f a, point3f b, point3f c,
                 point3f d, normal3f normal) {
    if (rect.alpha_mask || rect.bump_tex)
      throw std::runtime_error("rectangle alpha or bump mapping");
    Transform transform = placement * *rect.ObjectToWorld;
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
      AddTriangle(out);
    }
  }

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
    } else
      throw std::runtime_error("atmospheric or unrecognized infinite light");
    scene.lights.push_back(out);
    scene.light_sources.push_back(&source);
  }
};

#ifndef RAY_HAS_METAL_BVH
std::unique_ptr<WavefrontSession> MakeMetalWavefront(const WavefrontScene &, size_t) {
  throw std::runtime_error("Metal was not enabled in this build");
}
#endif

std::unique_ptr<WavefrontSession> PrepareWavefront(const hitable_list &world,
                                                   const hitable_list &lights, RayCamera &camera,
                                                   size_t capacity, int sampler,
                                                   WavefrontScene &scene, WavefrontReport &report) {
  auto started = Clock::now();
  report = WavefrontReport();
  report.requested = true;
  try {
    if (lights.volume_scene && (lights.volume_scene->has_media || lights.volume_scene->atmosphere ||
                                lights.volume_scene->transparent_background))
      throw std::runtime_error("participating media, atmosphere, or transparent background");
    if (camera.get_fov() < 0 || camera.get_fov() == 360 || camera.get_camera_motion_blur())
      throw std::runtime_error("realistic, panoramic, or motion-blurred camera");
    if (sampler == 1)
      throw std::runtime_error("stratified sampling");
    if (capacity > UINT32_MAX)
      throw std::runtime_error("more than 2^32 image pixels");
    WavefrontSceneCompiler compiler(scene);
    compiler.Geometry(world);
    if (lights.volume_scene)
      compiler.Points(lights.volume_scene->point_lights);
    WavefrontSceneCompiler::RefreshLights(scene);
    if (scene.triangles.empty())
      throw std::runtime_error("scene without triangles");
    auto session = MakeMetalWavefront(scene, capacity);
    report.used = true;
    report.triangles = scene.triangles.size();
    report.upload_seconds = std::chrono::duration<double>(Clock::now() - started).count();
    return session;
  } catch (const std::exception &error) {
    report.fallback = error.what();
    scene = WavefrontScene();
    Rcpp::warning("Metal wavefront: %s; using the CPU NEE integrator.", report.fallback.c_str());
    return nullptr;
  }
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
  p.camera.flip_x = camera.get_camera_flip_x();
  if (auto *cam = dynamic_cast<class camera *>(&camera)) {
    p.camera.origin = Pack(cam->origin);
    p.camera.lower_left = Pack(cam->lower_left_corner);
    p.camera.horizontal = Pack(cam->horizontal);
    p.camera.vertical = Pack(cam->vertical);
    p.camera.right = Pack(cam->u);
    p.camera.right.w = cam->lens_radius;
    p.camera.up = Pack(cam->v);
  } else if (auto *cam = dynamic_cast<ortho_camera *>(&camera)) {
    p.camera.orthographic = 1;
    p.camera.origin = Pack(cam->origin);
    p.camera.lower_left = Pack(cam->lower_left_corner);
    p.camera.horizontal = Pack(cam->horizontal);
    p.camera.vertical = Pack(cam->vertical);
    p.camera.forward = Pack(cam->w);
  } else
    throw std::runtime_error("Metal wavefront camera changed to an unsupported type");
  WavefrontSceneCompiler::RefreshLights(scene);
  std::vector<uint32_t> active;
  for (const auto &block : film.pixel_chunks)
    for (int x = block.startx; x < block.endx; ++x)
      for (int y = block.starty; y < block.endy; ++y)
        active.push_back(y + film.ny * x);
  p.active = active.size();
  const WFPixel *pixels = nullptr;
  if (!session.Render(p, active, scene.lights, cancelled, pixels))
    return false;
  for (uint32_t index : active) {
    size_t x = index / film.ny, y = index % film.ny;
    const auto &value = pixels[index];
    point3f color(value.radiance.x, value.radiance.y, value.radiance.z);
    color = clamp_point(de_nan(color), 0, clamp) * camera.get_iso();
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
