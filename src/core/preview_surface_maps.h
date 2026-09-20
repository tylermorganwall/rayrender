#ifndef RAYRENDER_PREVIEW_SURFACE_MAPS_H
#define RAYRENDER_PREVIEW_SURFACE_MAPS_H

#include "preview_textures.h"
#include "../materials/constant.h"
#include "../hitables/instance.h"
#include "../hitables/box.h"
#include "../hitables/trimesh.h"
#include "../hitables/raymesh.h"
#include "../hitables/plymesh.h"
#include "../hitables/mesh3d.h"
#include "../volumes/boundary.h"
#include <set>
#include <vector>
#include "../hitables/sphere.h"
#include "../hitables/ellipsoid.h"
#include "../hitables/rectangle.h"
#include "../hitables/cylinder.h"
#include "../hitables/disk.h"
#include "../hitables/triangle.h"

// Recipes survive disabling a map, so enabling it again restores its pixels.
// They are attached only to editor materials and copied before any edit.
struct PreviewSurfaceMapSettings {
  bool alpha_enabled = false, bump_enabled = false;
  std::string alpha_path, bump_path;
  std::shared_ptr<alpha_texture> alpha;
  std::shared_ptr<bump_texture> bump;
  double bump_intensity = 1;
  std::array<double, 3> bump_repeat{1, 1, 0};

  void Build() {
    if (alpha_enabled && !alpha_path.empty() &&
        (!alpha || alpha_path != alpha->preview_path)) {
      auto owner = std::make_shared<TextureCache>();
      int width, height, channels;
      auto data = owner->LookupChar(
          PreviewTextures::FilePath(alpha_path), width, height, channels, 0);
      // Grayscale/RGB files describe opacity directly; RGBA files use alpha.
      // Store our own one-channel copy so the loader cache can be released.
      auto pixels =
          std::make_shared<std::vector<unsigned char>>(size_t(width) * height);
      for (size_t i = 0; i < pixels->size(); ++i) {
        (*pixels)[i] =
            channels == 2 || channels == 4 ? data[i * channels + channels - 1]
            : channels == 1 ? data[i]
                            : static_cast<unsigned char>(.2126 * data[3 * i] +
                                                         .7152 * data[3 * i + 1] +
                                                         .0722 * data[3 * i + 2]);
      }
      alpha = std::make_shared<alpha_texture>(pixels->data(), width, height, 1);
      alpha->preview_path = alpha_path;
      alpha->preview_owner = pixels;
    }
    if (bump_enabled && !bump_path.empty() &&
        (!bump || bump_path != bump->preview_path)) {
      auto owner = std::make_shared<TextureCache>();
      int width, height, channels;
      auto data = owner->LookupChar(
          PreviewTextures::FilePath(bump_path), width, height, channels, 1);
      if (width < 3 || height < 3) {
        throw std::runtime_error("Bump maps must be at least 3 by 3 pixels.");
      }
      bump = std::make_shared<bump_texture>(data, width, height, 1, bump_intensity);
      bump->preview_path = bump_path;
      bump->preview_owner = owner;
    }
    if (alpha_enabled && !alpha) {
      throw std::runtime_error("Choose an alpha map file before applying.");
    }
    if (bump_enabled && !bump) {
      throw std::runtime_error("Choose a bump map file before applying.");
    }
    if (bump) {
      bump = std::make_shared<bump_texture>(*bump);
      bump->intensity = bump_intensity;
      bump->repeatu = bump_repeat[0];
      bump->repeatv = bump_repeat[1];
    }
  }
};

namespace PreviewTextures {
struct SurfaceTargets {
  std::vector<std::shared_ptr<alpha_texture>*> alpha;
  std::vector<std::shared_ptr<bump_texture>*> bump;
  std::vector<Float*> densities;
  std::set<TriangleMesh*> meshes;
  std::vector<BVHAggregate*> aggregates;

  template <typename Shape> bool Primitive(hitable* root, material* selected) {
    auto p = dynamic_cast<Shape*>(root);
    if (!p) {
      return false;
    }
    if (p->mat_ptr.get() == selected) {
      alpha.push_back(&p->alpha_mask);
      bump.push_back(&p->bump_tex);
    }
    return true;
  }

  void Mesh(TriangleMesh* mesh, material* selected) {
    if (!mesh || !meshes.insert(mesh).second) {
      return;
    }
    for (size_t slot = 0; slot < mesh->mesh_materials.size(); ++slot) {
      if (mesh->mesh_materials[slot].get() == selected) {
        if (slot < mesh->alpha_textures.size()) {
          alpha.push_back(&mesh->alpha_textures[slot]);
        }
        if (slot < mesh->bump_textures.size()) {
          bump.push_back(&mesh->bump_textures[slot]);
        }
      }
    }
  }

  // Locate the geometry-level maps belonging to this slot, including mesh
  // material tables. Repeated references to shared mesh tables are visited once.
  void Gather(hitable* root, material* selected) {
    if (!root) {
      return;
    }
    if (auto p = dynamic_cast<AnimatedHitable*>(root)) {
      Gather(p->primitive.get(), selected);
    } else if (auto p = dynamic_cast<MediumBoundary*>(root)) {
      Gather(p->geometry.get(), selected);
    } else if (auto p = dynamic_cast<instance*>(root)) {
      Gather(p->original_scene, selected);
    } else if (auto p = dynamic_cast<constant_medium*>(root)) {
      if (p->phase_function.get() == selected) {
        densities.push_back(&p->density);
      }
    } else if (auto p = dynamic_cast<trimesh*>(root)) {
      Mesh(p->mesh.get(), selected);
      aggregates.push_back(p->tri_mesh_bvh.get());
    } else if (auto p = dynamic_cast<raymesh*>(root)) {
      Mesh(p->mesh.get(), selected);
      aggregates.push_back(p->tri_mesh_bvh.get());
    } else if (auto p = dynamic_cast<plymesh*>(root)) {
      Mesh(p->mesh.get(), selected);
      aggregates.push_back(p->ply_mesh_bvh.get());
    } else if (auto p = dynamic_cast<mesh3d*>(root)) {
      Mesh(p->mesh.get(), selected);
      aggregates.push_back(p->mesh_bvh.get());
    } else if (auto p = dynamic_cast<triangle*>(root)) {
      Mesh(p->mesh, selected);
    } else if (auto p = dynamic_cast<BVHAggregate*>(root)) {
      for (const auto& child : p->Primitives()) {
        Gather(child.get(), selected);
      }
      aggregates.push_back(p);
    } else if (auto p = dynamic_cast<box*>(root)) {
      Gather(&p->list, selected);
    } else if (auto p = dynamic_cast<hitable_list*>(root)) {
      for (const auto& child : p->objects) {
        Gather(child.get(), selected);
      }
    } else {
      Primitive<sphere>(root, selected) || Primitive<ellipsoid>(root, selected) ||
          Primitive<xy_rect>(root, selected) || Primitive<xz_rect>(root, selected) ||
          Primitive<yz_rect>(root, selected) || Primitive<cylinder>(root, selected) ||
          Primitive<disk>(root, selected);
    }
  }
};
}
#endif
