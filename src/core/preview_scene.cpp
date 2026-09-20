#include "preview_scene.h"
#include "bvh.h"
#include "../hitables/instance.h"
#include "../materials/material.h"
#include <set>
#include <algorithm>
#include <cmath>
#include <limits>

// The GUI uses column-major matrices. Validate before constructing Transform,
// whose inverse must exist for ray intersections and object-space shading.
Transform PreviewMatrix(const std::array<float, 16>& values) {
  Float m[4][4];
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 4; ++col) {
      const double v = values[row + 4 * col];
      if (!std::isfinite(v) || std::abs(v) > 1e12) {
        throw std::runtime_error("Transform values must be finite and bounded.");
      }
      m[row][col] = v;
    }
  }

  if (std::abs(m[3][0]) + std::abs(m[3][1]) + std::abs(m[3][2]) > 1e-6 ||
      std::abs(m[3][3] - 1) > 1e-6) {
    throw std::runtime_error("Object transforms must be affine.");
  }

  const vec3f x(m[0][0], m[1][0], m[2][0]), y(m[0][1], m[1][1], m[2][1]),
      z(m[0][2], m[1][2], m[2][2]);
  // Compare the determinant with the axis lengths to catch nearly parallel axes
  // independently of the overall object size. Reflections remain valid.
  const double volume = x.length() * y.length() * z.length();
  if (volume < 1e-15 || std::abs(dot(x, cross(y, z))) < volume * 1e-6) {
    throw std::runtime_error(
        "Scale must stay nonzero; singular transforms cannot be applied.");
  }

  return Transform(Matrix4x4(m));
}

// Convert the renderer matrix back to the layout expected by the gizmo.
void PreviewMatrix(const Transform& t, std::array<float, 16>& values) {
  const auto& m = t.GetMatrix();
  for (int row = 0; row < 4; ++row) {
    for (int col = 0; col < 4; ++col) {
      values[row + 4 * col] = m.m[row][col];
    }
  }
}

// Numeric controls compose scale, then X/Y/Z rotations in degrees, then translation.
// Round-trip through the GUI representation so typed edits get the same validation.
Transform PreviewNumericTransform(const PreviewObjectState& ui) {
  for (int i = 0; i < 3; ++i) {
    if (!std::isfinite(ui.translation[i]) || !std::isfinite(ui.rotation[i]) ||
        !std::isfinite(ui.scale[i]) || std::abs(ui.scale[i]) < .0001 ||
        std::abs(ui.scale[i]) > 10000) {
      throw std::runtime_error(
          "Use finite transforms and nonzero scales with magnitudes between 0.0001 and 10000.");
    }
  }

  auto t = Translate(vec3f(ui.translation[0], ui.translation[1], ui.translation[2])) *
           RotateZ(ui.rotation[2]) * RotateY(ui.rotation[1]) * RotateX(ui.rotation[0]) *
           Scale(ui.scale[0], ui.scale[1], ui.scale[2]);
  std::array<float, 16> values;
  PreviewMatrix(t, values);
  return PreviewMatrix(values);
}

// Populate the numeric controls from the current gizmo matrix. This extracts
// translation, axis scales and Euler angles; it does not modify the matrix itself.
void PreviewDecompose(PreviewObjectState& ui) {
  auto& m = ui.model;
  for (int i = 0; i < 3; ++i) {
    ui.translation[i] = m[12 + i];
    ui.scale[i] = std::sqrt(m[4 * i] * m[4 * i] + m[4 * i + 1] * m[4 * i + 1] +
                            m[4 * i + 2] * m[4 * i + 2]);
  }

  if (*std::min_element(ui.scale.begin(), ui.scale.end()) < 1e-12) {
    return;
  }

  // A negative determinant represents a reflection. Assign its sign to X so the
  // numeric controls retain a signed scale instead of turning it into a rotation.
  if (dot(vec3f(m[0], m[1], m[2]),
          cross(vec3f(m[4], m[5], m[6]), vec3f(m[8], m[9], m[10]))) < 0) {
    ui.scale[0] = -ui.scale[0];
  }

  const double y = std::asin(std::clamp(-m[2] / ui.scale[0], -1.0, 1.0));
  // At the Euler singularity, fix Z to zero and solve the remaining X rotation.
  const bool singular = std::abs(std::cos(y)) < 1e-6;
  ui.rotation = {
      180 / M_PI *
          (singular ? std::atan2(-m[9] / ui.scale[2], m[5] / ui.scale[1])
                    : std::atan2(m[6] / ui.scale[1], m[10] / ui.scale[2])),
      180 / M_PI * y,
      180 / M_PI * (singular ? 0 : std::atan2(m[1] / ui.scale[0], m[0] / ui.scale[0]))};
}

// R supplies outer-group membership separately from flattened geometry rows.
// Ungrouped rows keep zero here and receive their own row/placement identity.
PreviewScene::PreviewScene(const Rcpp::List& scene) {
  Rcpp::NumericVector x = scene["x"];
  groups.resize(x.size(), 0);
  group_paths.resize(x.size());
  if (scene.containsElementNamed("shape")) {
    shapes = Rcpp::as<std::vector<int>>(scene["shape"]);
  }
  if (scene.containsElementNamed("preview_groups")) {
    groups = Rcpp::as<std::vector<int>>(scene["preview_groups"]);
  }
  if (scene.containsElementNamed("preview_paths")) {
    Rcpp::List paths = scene["preview_paths"];
    for (size_t i = 0; i < group_paths.size(); ++i) {
      group_paths[i] = Rcpp::as<std::vector<int>>(paths[i]);
    }
  } else {
    for (size_t i = 0; i < groups.size(); ++i) {
      if (groups[i]) {
        group_paths[i].push_back(groups[i]);
      }
    }
  }
}

namespace {
constexpr uint64_t GroupBit = UINT64_C(1) << 63;
constexpr uint64_t InstanceSetBit = UINT64_C(1) << 62;
constexpr uint64_t NestedBit = UINT64_C(1) << 61;
uint64_t ObjectId(const PreviewObjectKey& key) {
  return ((uint64_t(key.row) + 1) << 32) | (key.placement + 1);
}

// Source-row order is stable when a moved child changes BVH traversal order.
// Keep mesh slots intact and deduplicate genuinely shared instance materials.
std::vector<material*> RootMaterials(const PreviewObjectRoot& root) {
  if (!root.contents) {
    return PreviewMaterials(root.object.get());
  }
  std::vector<material*> result;
  std::set<material*> seen;
  for (const auto& child : root.contents->roots) {
    for (auto* value : RootMaterials(child)) {
      if (seen.insert(value).second) {
        result.push_back(value);
      }
    }
  }
  return result;
}

// Undo the actual placement wrappers at the ray's shutter time. In particular,
// do not use the gizmo's midpoint frame for animated-instance visibility tests.
Ray ChildRay(const PreviewObjectRoot& root, Ray ray) {
  hitable* object = root.object.get();
  if (auto animated = dynamic_cast<AnimatedHitable*>(object)) {
    Transform transform;
    animated->PrimitiveToWorld.Interpolate(ray.time(), &transform);
    ray = Inverse(transform)(ray);
    object = animated->primitive.get();
  }
  auto placement = dynamic_cast<instance*>(object);
  if (!placement) {
    throw std::runtime_error("Instance editor placement is missing.");
  }
  return (*placement->WorldToObject)(ray);
}
}

bool PreviewSceneSettings::Empty() const {
  if (!transforms.empty() || !materials.empty()) {
    return false;
  }
  return std::all_of(children.begin(), children.end(), [](const auto& child) {
    return child.second.Empty();
  });
}

uint64_t PreviewScene::NestedId(uint64_t placement, uint64_t child) const {
  const auto key = std::make_pair(placement, child);
  auto found = nested_ids.find(key);
  if (found != nested_ids.end()) {
    return found->second;
  }
  const uint64_t id = NestedBit | (nested_nodes.size() + 1);
  nested_ids.emplace(key, id);
  nested_nodes.push_back(key);
  return id;
}

// Resolve UI identity into a source scene plus the enclosing placement path.
// Geometry remains local to its shared source; editor transforms are world-space.
PreviewScene::SelectionContext PreviewScene::Resolve(uint64_t id) const {
  SelectionContext context{this, id, Transform(), {}};
  if (!(id & NestedBit)) {
    return context;
  }
  Hierarchy(); // Rebuilt child snapshots assign their local IDs in scene order.
  const uint64_t index = id & ~NestedBit;
  if (!index || index > nested_nodes.size()) {
    throw std::runtime_error("Selected instance child no longer exists.");
  }
  const auto address = nested_nodes[index - 1];
  for (const auto& root : roots) {
    if (ObjectId(root.key) == address.first && root.contents) {
      context = root.contents->Resolve(address.second);
      context.to_world = root.frame * context.to_world;
      context.placements.insert(context.placements.begin(), &root);
      return context;
    }
  }
  throw std::runtime_error("Selected instance placement no longer exists.");
}

// Selection identity never depends on rebuilt geometry addresses. Tree leaves
// address individual roots; groups and instance sets address their descendants.
bool PreviewScene::Contains(const PreviewObjectRoot& root, uint64_t id) const {
  if (id & NestedBit) {
    const uint64_t index = id & ~NestedBit;
    return index && index <= nested_nodes.size() &&
           nested_nodes[index - 1].first == ObjectId(root.key);
  }
  if (id & GroupBit) {
    const auto& path = group_paths.at(root.key.row);
    return std::find(path.begin(), path.end(), int(id & ~GroupBit)) != path.end();
  }
  if (id & InstanceSetBit) {
    return root.key.row + 1 == (id & ~InstanceSetBit);
  }
  return ObjectId(root.key) == id;
}

std::vector<PreviewHierarchyNode> PreviewScene::Hierarchy() const {
  std::vector<PreviewHierarchyNode> nodes;
  std::set<uint64_t> added;
  for (const auto& root : roots) {
    uint64_t parent = 0;
    for (int group : group_paths.at(root.key.row)) {
      const uint64_t id = GroupBit | uint64_t(group);
      if (added.insert(id).second) {
        nodes.push_back({id, parent, "Group " + std::to_string(group)});
      }
      parent = id;
    }
    const bool instance = root.key.row < shapes.size() && shapes[root.key.row] == 15;
    if (instance) {
      const uint64_t id = InstanceSetBit | (root.key.row + 1);
      if (added.insert(id).second) {
        nodes.push_back({id, parent, "Instances " + std::to_string(root.key.row + 1)});
      }
      parent = id;
    }
    nodes.push_back(
        {ObjectId(root.key),
         parent,
         instance ? "Instance " + std::to_string(root.key.placement + 1)
                  : root.object->GetName() + " " + std::to_string(root.key.row + 1)});
    if (root.contents) {
      const uint64_t placement = ObjectId(root.key);
      for (const auto& child : root.contents->Hierarchy()) {
        nodes.push_back({NestedId(placement, child.id),
                         child.parent ? NestedId(placement, child.parent) : placement,
                         child.label});
      }
    }
  }
  return nodes;
}

Transform PreviewScene::ObjectTransform(size_t row, size_t placement) const {
  auto found = settings.transforms.find({row, placement});
  return found == settings.transforms.end() ? Transform() : found->second;
}

// Transforms only change a placement wrapper. Material edits require a separate
// child scene so shared geometry in sibling placements keeps its original material.
bool PreviewScene::UniqueInstance(size_t row, size_t placement) const {
  const auto child = InstanceSettings(row, placement);
  return settings.materials.count({row, placement}) != 0 || (child && !child->Empty());
}

const PreviewSceneSettings* PreviewScene::InstanceSettings(size_t row,
                                                           size_t placement) const {
  const auto found = settings.children.find({row, placement});
  return found == settings.children.end() ? nullptr : &found->second;
}

// Ancestor material overrides are applied first, then more specific child edits.
// This preserves an individual child color when its enclosing placement is rebuilt.
void PreviewScene::ReplayMaterials() {
  for (const auto& root : roots) {
    const auto edits = settings.materials.find(root.key);
    if (edits != settings.materials.end()) {
      const auto materials = RootMaterials(root);
      for (const auto& edit : edits->second) {
        if (edit.first >= materials.size()) {
          throw std::runtime_error("Material slot no longer exists.");
        }
        PreviewApplyMaterial(materials[edit.first], edit.second, root.object.get());
      }
    }
    if (root.contents) {
      root.contents->ReplayMaterials();
    }
  }
}

// A later parent edit must also update existing, narrower overrides for the same
// material. Leave unrelated child materials and properties at their current values.
void PreviewScene::ForwardMaterial(const PreviewObjectRoot& root,
                                   PreviewSceneSettings& next, material* target,
                                   const PreviewMaterialEdit& edit) const {
  auto settings = next.children.find(root.key);
  if (!root.contents || settings == next.children.end()) {
    return;
  }
  auto& children = settings->second;
  for (const auto& child : root.contents->roots) {
    const auto materials = RootMaterials(child);
    auto saved = children.materials.find(child.key);
    if (saved != children.materials.end()) {
      for (size_t slot = 0; slot < materials.size(); ++slot) {
        if (materials[slot] == target && saved->second.count(slot)) {
          saved->second[slot] = edit;
        }
      }
    }
    root.contents->ForwardMaterial(child, children, target, edit);
  }
}

// Called during scene construction on each selectable root. Apply saved material
// overrides to the newly built geometry before it becomes part of the live scene.
void PreviewScene::Record(size_t row, size_t placement, std::shared_ptr<hitable> object,
                          const Transform& frame,
                          std::shared_ptr<PreviewScene> contents) {
  PreviewObjectKey key{row, placement};
  auto edits = settings.materials.find(key);
  if (edits != settings.materials.end()) {
    auto materials = RootMaterials({key, 0, object, frame, contents});
    for (const auto& entry : edits->second) {
      if (entry.first >= materials.size()) {
        throw std::runtime_error("Material slot no longer exists.");
      }
      PreviewApplyMaterial(materials[entry.first], entry.second, object.get());
    }
  }

  // Reserve the high bit for groups; all members share that selection ID. Other
  // roots encode row and placement, with one-based components to keep zero unused.
  uint64_t id = groups.at(row) > 0 ? (UINT64_C(1) << 63) | uint64_t(groups[row])
                                   : ((uint64_t(row) + 1) << 32) | (placement + 1);
  if (contents) {
    contents->ReplayMaterials();
  }
  roots.push_back({key, id, std::move(object), frame, std::move(contents)});
  selection_bvh.reset();
  ++revision;
}

// Copy material values into pointer-free UI panels. A hit material chooses the
// initial slot, while the root/slot keys identify where a later edit must apply.
void PreviewScene::DescribeMaterials(uint64_t id, PreviewObjectState& ui,
                                     material* picked) {
  if (id & NestedBit) {
    auto context = Resolve(id);
    auto snapshot = *context.scene;
    snapshot.DescribeMaterials(context.id, ui, picked);
    return;
  }
  ui.materials.clear();
  ui.material_slot = 0;
  for (const auto& root : roots) {
    if (Contains(root, id)) {
      auto materials = RootMaterials(root);
      for (size_t slot = 0; slot < materials.size(); ++slot) {
        PreviewMaterialPanel panel;
        panel.row = root.key.row;
        panel.placement = root.key.placement;
        panel.slot = slot;
        panel.targets.push_back({panel.row, panel.placement, panel.slot});
        panel.type = materials[slot]->GetName();
        panel.label = std::to_string(root.key.row + 1) + ":" +
                      std::to_string(root.key.placement + 1) + " / " + panel.type +
                      " " + std::to_string(slot + 1);
        for (const auto& binding :
             PreviewMaterialFields(materials[slot], root.object.get())) {
          panel.fields.push_back(binding.field);
        }
        // An instance parent edits matching source slots across all placements.
        // Keep different slots separate, even when they have the same type.
        auto existing = ui.materials.end();
        if ((id & InstanceSetBit) && !(id & GroupBit)) {
          existing =
              std::find_if(ui.materials.begin(),
                           ui.materials.end(),
                           [&](const PreviewMaterialPanel& other) {
                             if (other.row != panel.row || other.slot != slot ||
                                 other.type != panel.type ||
                                 other.fields.size() != panel.fields.size()) {
                               return false;
                             }
                             for (size_t i = 0; i < panel.fields.size(); ++i) {
                               if (other.fields[i].name != panel.fields[i].name ||
                                   other.fields[i].count != panel.fields[i].count) {
                                 return false;
                               }
                             }
                             return true;
                           });
        }
        size_t panel_index = ui.materials.size();
        if (existing == ui.materials.end()) {
          ui.materials.push_back(std::move(panel));
        } else {
          panel_index = size_t(existing - ui.materials.begin());
          existing->targets.push_back(panel.targets.front());
          existing->label = existing->type + " " + std::to_string(slot + 1);
          for (size_t i = 0; i < panel.fields.size(); ++i) {
            auto& field = existing->fields[i];
            const auto& value = panel.fields[i];
            field.mixed = field.mixed || field.text != value.text;
            for (unsigned j = 0; j < field.count; ++j) {
              field.mixed = field.mixed || field.values[j] != value.values[j];
            }
            field.minimum = std::min(field.minimum, value.minimum);
            field.maximum = std::max(field.maximum, value.maximum);
            field.color = field.color && value.color;
          }
        }
        if (materials[slot] == picked) {
          ui.material_slot = static_cast<int32_t>(panel_index);
        }
      }
    }
  }
}

// Replace the current selection and discard uncommitted UI edits. ID zero (or
// an ID absent from the rebuilt roots) leaves the editor deselected.
void PreviewScene::Describe(uint64_t id, PreviewObjectState& ui, material* picked) {
  if (id & NestedBit) {
    const auto context = Resolve(id);
    // Child snapshots can be shared across untouched placements. Selection must
    // never change that shared editor's pivot or other placement-local UI state.
    auto snapshot = *context.scene;
    snapshot.Describe(context.id, ui, picked);
    if (ui.selected) {
      ui.id = id;
      selected_id = id;
      selected_model = context.to_world * snapshot.selected_model;
      PreviewMatrix(selected_model, ui.model);
      PreviewDecompose(ui);
    }
    return;
  }
  ui.selected = false;
  ui.id = 0;
  ui.materials.clear();
  ui.error.clear();
  ui.cancel_transform = false;
  ui.transform_pending = ui.transform_active = ui.apply_transform =
      ui.material_pending = ui.apply_material = ui.revert = ui.numeric_transform =
          false;
  const PreviewObjectRoot* first = nullptr;
  size_t count = 0;
  aabb bounds;
  for (const auto& root : roots) {
    if (Contains(root, id)) {
      aabb b;
      if (root.object->bounding_box(0, 1, b)) {
        bounds = count ? surrounding_box(bounds, b) : b;
      }
      if (!first) {
        first = &root;
      }
      ++count;
    }
  }

  selected_id = id;
  if (!first) {
    return;
  }

  ui.selected = true;
  ui.id = id;
  const bool group = (id & (GroupBit | InstanceSetBit)) != 0;
  ui.label = first->object->GetName();
  for (const auto& node : Hierarchy()) {
    if (node.id == id) {
      ui.label = node.label;
      break;
    }
  }
  // Groups rotate/scale around their combined bounds center. A single root uses
  // its recorded placement frame, preserving the object or instance origin.
  selected_model = group ? Translate(convert_to_vec3(bounds.Centroid())) : first->frame;
  PreviewMatrix(selected_model, ui.model);
  PreviewDecompose(ui);
  DescribeMaterials(id, ui, picked);
}

namespace {
// These wrappers are only used by editor visibility queries. Replacing shape in
// their hit records recovers the outer root even inside shared meshes/instances.
class SelectionRoot final : public hitable {
public:
  SelectionRoot(std::shared_ptr<hitable> object, size_t index)
      : index(index), object(std::move(object)) {
    ObjectToWorld = WorldToObject = nullptr;
  }
  const bool hit(const Ray& ray, Float minimum, Float maximum, hit_record& rec,
                 random_gen& rng) const override {
    return object->hit(ray, minimum, maximum, rec, rng) && Accept(rec);
  }
  const bool hit(const Ray& ray, Float minimum, Float maximum, hit_record& rec,
                 Sampler* sampler) const override {
    return object->hit(ray, minimum, maximum, rec, sampler) && Accept(rec);
  }
  bool bounding_box(Float begin, Float end, aabb& bounds) const override {
    return object->bounding_box(begin, end, bounds);
  }
  std::string GetName() const override {
    return "Selection root";
  }
  size_t GetSize() override {
    return sizeof(*this);
  }
  void hitable_info_bounds(Float begin, Float end) const override {
  }
  size_t index;

private:
  std::shared_ptr<hitable> object;
  bool Accept(hit_record& rec) const {
    if (rec.alpha_miss || !std::isfinite(rec.t)) {
      return false;
    }
    rec.shape = this;
    return true;
  }
};
}

void PreviewScene::PrepareSelectionBvh() {
  if (selection_bvh || roots.empty()) {
    return;
  }
  std::vector<std::shared_ptr<hitable>> objects;
  objects.reserve(roots.size());
  for (size_t i = 0; i < roots.size(); ++i) {
    objects.push_back(std::make_shared<SelectionRoot>(roots[i].object, i));
  }
  selection_bvh = std::make_shared<BVHAggregate>(std::move(objects), 0, 1, 4, true);
}

const PreviewObjectRoot* PreviewScene::HitRoot(const Ray& ray, material*& picked,
                                               random_gen& rng) {
  PrepareSelectionBvh();
  hit_record rec;
  if (!selection_bvh ||
      !selection_bvh->hit(
          ray, .001f, std::numeric_limits<Float>::infinity(), rec, rng)) {
    return nullptr;
  }
  picked = rec.mat_ptr;
  const auto* wrapper = static_cast<const SelectionRoot*>(rec.shape);
  return &roots[wrapper->index];
}

// Build the hierarchy path of the frontmost hit. Mesh triangles remain a single
// object; instance contents recurse through their own source scene and placement.
std::vector<uint64_t> PreviewScene::HitPath(const Ray& ray, material*& picked,
                                            random_gen& rng,
                                            const std::function<bool()>& cancel) {
  if (cancel && cancel()) {
    return {};
  }
  Rcpp::checkUserInterrupt();
  const auto* root = HitRoot(ray, picked, rng);
  if (!root) {
    return {};
  }
  std::vector<uint64_t> path;
  for (int group : group_paths.at(root->key.row)) {
    path.push_back(GroupBit | uint64_t(group));
  }
  if (root->key.row < shapes.size() && shapes[root->key.row] == 15) {
    path.push_back(InstanceSetBit | (root->key.row + 1));
  }
  const uint64_t placement = ObjectId(root->key);
  path.push_back(placement);
  if (root->contents) {
    for (uint64_t child :
         root->contents->HitPath(ChildRay(*root, ray), picked, rng, cancel)) {
      path.push_back(NestedId(placement, child));
    }
  }
  return path;
}

bool PreviewScene::Pick(const Ray& ray, PreviewObjectState& ui,
                        const std::function<bool()>& cancel) {
  const auto tree = Hierarchy();
  bool cancelled = false;
  auto check_cancel = [&]() {
    cancelled = cancelled || (cancel && cancel());
    return cancelled;
  };
  material* picked = nullptr;
  random_gen rng(1);
  const auto path = HitPath(ray, picked, rng, check_cancel);
  if (cancelled) {
    return false; // Keep the existing selection when the query is interrupted.
  }
  if (path.empty()) {
    Describe(0, ui);
    return false;
  }
  // First selection still stops at the outer group or individual placement.
  // A collection selected in the tree can then drill down to its clicked copy.
  uint64_t id = path.front();
  if ((id & InstanceSetBit) && path.size() > 1) {
    id = path[1];
  }
  if (ui.selected) {
    std::map<uint64_t, uint64_t> parents;
    for (const auto& node : tree) {
      parents.emplace(node.id, node.parent);
    }
    std::vector<uint64_t> current;
    for (uint64_t node = ui.id; node && parents.count(node); node = parents.at(node)) {
      current.push_back(node);
    }
    std::reverse(current.begin(), current.end());
    size_t common = 0;
    while (common < current.size() && common < path.size() &&
           current[common] == path[common]) {
      ++common;
    }
    if (common) {
      // An ancestor selection descends one level. At the same depth, a click
      // on a sibling selects that sibling without jumping back to the parent.
      id = path[std::min(common, path.size() - 1)];
    }
  }
  if (!ui.selected || ui.id != id) {
    Describe(id, ui, picked);
  }
  return true;
}

// Trace nearest visible surfaces, including unselected foreground geometry.
// Reuse an acceleration structure across selection/camera changes; invalidate it
// only after geometry is rebuilt. Rendering RNGs and image pixels stay untouched.
std::vector<uint8_t>
PreviewScene::SelectionMask(uint64_t id, uint32_t width, uint32_t height,
                            const std::function<bool(float, float, Ray&)>& make_ray,
                            const std::function<bool()>& cancel) {
  if (!width || !height) {
    return {};
  }
  std::vector<uint8_t> mask(size_t(width) * height, 0);
  if (!id || roots.empty()) {
    return mask;
  }
  const auto selected = Resolve(id);
  PrepareSelectionBvh();
  random_gen rng(1);
  for (uint32_t y = 0; y < height; ++y) {
    if (y % 8 == 0) {
      if (cancel && cancel()) {
        return {};
      }
      Rcpp::checkUserInterrupt();
    }
    for (uint32_t x = 0; x < width; ++x) {
      Ray ray;
      if (!make_ray((x + .5f) / width, (y + .5f) / height, ray)) {
        continue;
      }
      material* material = nullptr;
      const PreviewObjectRoot* root = HitRoot(ray, material, rng);
      PreviewScene* scene = this;
      // Every ancestor must be the frontmost hit. This keeps unselected objects,
      // including siblings inside an instance, in the occlusion calculation.
      for (const auto* placement : selected.placements) {
        if (root != placement) {
          root = nullptr;
          break;
        }
        ray = ChildRay(*placement, ray);
        scene = placement->contents.get();
        root = scene->HitRoot(ray, material, rng);
      }
      if (root) {
        mask[x + size_t(width) * y] = selected.scene->Contains(*root, selected.id);
      }
    }
  }
  return mask;
}

// Consume explicit apply/reset requests at a checkpoint where render workers have
// drained. Build a candidate scene and UI snapshot before publishing either one.
// History owns value snapshots, not geometry pointers. Rebuild a candidate first
// so missing texture files or other build failures leave the live scene untouched.
std::function<void()> PreviewScene::PrepareRestore(const PreviewSceneSettings& saved,
                                                   uint64_t selection) {
  if (!prepare_rebuild) {
    throw std::runtime_error("Scene editing is not connected to this renderer.");
  }
  auto next = std::make_shared<PreviewScene>(*this);
  next->settings = saved;
  next->roots.clear();
  auto publish = prepare_rebuild(*next);
  PreviewObjectState inspector;
  next->Describe(selection, inspector);
  return [this, next, publish] {
    publish();
    settings.transforms.swap(next->settings.transforms);
    settings.materials.swap(next->settings.materials);
    settings.children.swap(next->settings.children);
    nested_ids.swap(next->nested_ids);
    nested_nodes.swap(next->nested_nodes);
    roots.swap(next->roots);
    selected_id = next->selected_id;
    selected_model = next->selected_model;
    selection_bvh.reset();
    ++revision;
  };
}

bool PreviewScene::Apply(PreviewObjectState& ui) {
  if (!ui.selected || (!ui.apply_transform && !ui.apply_material && !ui.revert)) {
    return false;
  }

  const bool transform = ui.apply_transform, material_edit = ui.apply_material,
             revert = ui.revert;
  ui.apply_transform = ui.apply_material = ui.revert = false;
  try {
    const auto selection = Resolve(ui.id);
    const auto& source = *selection.scene;
    PreviewScene next = *this;
    next.roots.clear();
    auto* target_settings = &next.settings;
    for (const auto* placement : selection.placements) {
      target_settings = &target_settings->children[placement->key];
    }
    Transform model = selected_model;
    if (transform) {
      model =
          ui.numeric_transform ? PreviewNumericTransform(ui) : PreviewMatrix(ui.model);
      // Convert movement of the selection pivot into a world-space delta. Apply
      // the same delta to every group member on top of its previously saved edits.
      const auto world_delta = model * Inverse(selected_model);
      const auto delta = Inverse(selection.to_world) * world_delta * selection.to_world;
      for (const auto& root : source.roots) {
        if (source.Contains(root, selection.id)) {
          target_settings->transforms[root.key] =
              delta * source.ObjectTransform(root.key.row, root.key.placement);
        }
      }
    }
    if (material_edit) {
      for (const auto& panel : ui.materials) {
        if (panel.changed) {
          const bool has_field_edits = std::any_of(
              panel.fields.begin(), panel.fields.end(), [](const PreviewField& field) {
                return field.changed;
              });
          for (const auto& target : panel.targets) {
            const PreviewObjectKey key{target.row, target.placement};
            auto root = std::find_if(
                source.roots.begin(), source.roots.end(), [&](const auto& value) {
                  return value.key == key && source.Contains(value, selection.id);
                });
            if (root == source.roots.end()) {
              throw std::runtime_error("Material target is no longer selected.");
            }
            auto materials = RootMaterials(*root);
            if (target.slot >= materials.size()) {
              throw std::runtime_error("Material slot no longer exists.");
            }
            // Start with each placement's own values. Editing a shared color
            // must not overwrite a roughness override on just one placement.
            PreviewMaterialEdit edit;
            edit.type = panel.type;
            auto bindings =
                PreviewMaterialFields(materials[target.slot], root->object.get());
            if (bindings.size() != panel.fields.size()) {
              throw std::runtime_error("Material fields no longer match.");
            }
            for (size_t i = 0; i < bindings.size(); ++i) {
              auto field = bindings[i].field;
              if (panel.fields[i].changed || !has_field_edits) {
                field = panel.fields[i];
              }
              field.changed = field.mixed = false;
              edit.fields.push_back(std::move(field));
            }
            source.ForwardMaterial(
                *root, *target_settings, materials[target.slot], edit);
            target_settings->materials[key][target.slot] = std::move(edit);
          }
        }
      }
    }
    if (revert) {
      for (const auto& root : source.roots) {
        if (source.Contains(root, selection.id)) {
          target_settings->transforms.erase(root.key);
          target_settings->materials.erase(root.key);
          target_settings->children.erase(root.key);
        }
      }
    }
    if (!prepare_rebuild) {
      throw std::runtime_error("Scene editing is not connected to this renderer.");
    }
    auto publish =
        prepare_rebuild(next); // Build off to the side; errors preserve the live scene.
    // Refresh material values from rebuilt objects before committing: allocations
    // and material inspection can still fail here without replacing the live scene.
    PreviewObjectState updated = ui;
    updated.hierarchy = next.Hierarchy();
    if (revert) {
      next.Describe(ui.id, updated);
    } else {
      const int slot = updated.material_slot;
      next.DescribeMaterials(ui.id, updated, nullptr);
      updated.material_slot =
          std::min(slot, static_cast<int>(updated.materials.size()) - 1);
      PreviewMatrix(model, updated.model);
      PreviewDecompose(updated);
    }
    publish(); // Workers are drained; everything below is a swap or scalar update.
    settings.transforms.swap(next.settings.transforms);
    settings.materials.swap(next.settings.materials);
    settings.children.swap(next.settings.children);
    nested_ids.swap(next.nested_ids);
    nested_nodes.swap(next.nested_nodes);
    roots.swap(next.roots);
    selection_bvh.reset();
    ++revision;
    selected_model = revert ? next.selected_model : model;
    ui = std::move(updated);
    ui.error.clear();
    ui.transform_pending = ui.material_pending = ui.numeric_transform = false;
    return true;
  } catch (const Rcpp::internal::InterruptedException&) {
    // Let R interrupts unwind through renderer cleanup.
    throw;
  } catch (const std::exception& error) {
    // Ordinary edit failures stay in the panel; restore the gizmo to the last
    // committed transform while leaving the live scene unchanged.
    ui.error = error.what();
    PreviewMatrix(selected_model, ui.model);
    PreviewDecompose(ui);
    ui.transform_pending = ui.numeric_transform = false;
    return false;
  }
}

// Return only committed changes for the image's scene_edits attribute. Root and
// material-slot indices become one-based for R. Transform deltas use the enclosing
// scene's coordinates; a children list descends into one instance's source scene.
Rcpp::List PreviewScene::ExportEdits() const {
  Rcpp::List result;

  for (const auto& root : roots) {
    Rcpp::List children = root.contents ? root.contents->ExportEdits() : Rcpp::List();
    if (settings.transforms.count(root.key) || settings.materials.count(root.key) ||
        children.size()) {
      Rcpp::NumericMatrix matrix(4, 4);
      auto transform = ObjectTransform(root.key.row, root.key.placement);
      const auto& m = transform.GetMatrix();
      for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
          matrix(i, j) = m.m[i][j];
        }
      }

      Rcpp::List materials;
      auto edits = settings.materials.find(root.key);
      if (edits != settings.materials.end()) {
        for (const auto& slot : edits->second) {
          Rcpp::List fields;
          for (const auto& field : slot.second.fields) {
            if (field.text_input) {
              fields[field.name] = field.text;
            } else {
              fields[field.name] = Rcpp::NumericVector(
                  field.values.begin(), field.values.begin() + field.count);
            }
          }
          materials.push_back(Rcpp::List::create(Rcpp::_["slot"] = slot.first + 1,
                                                 Rcpp::_["type"] = slot.second.type,
                                                 Rcpp::_["values"] = fields));
        }
      }

      auto entry = Rcpp::List::create(Rcpp::_["row"] = root.key.row + 1,
                                      Rcpp::_["instance"] = root.key.placement + 1,
                                      Rcpp::_["transform"] = matrix,
                                      Rcpp::_["materials"] = materials);
      if (children.size()) {
        entry["children"] = children;
      }
      result.push_back(entry);
    }
  }

  return result;
}

// Dragging only stages a UI matrix, so cancellation needs no geometry rebuild.
void PreviewScene::CancelTransform(PreviewObjectState& ui) {
  PreviewMatrix(selected_model, ui.model);
  PreviewDecompose(ui);
  ui.cancel_transform = ui.transform_pending = ui.transform_active =
      ui.numeric_transform = ui.apply_transform = false;
}

namespace {
// Validate the saved row/placement path before any scene objects are built.
// Each instance owns a separate child settings tree, so restoring one placement
// cannot recolor or transform the other copies of its source geometry.
PreviewSceneSettings ReadPreviewEdits(const Rcpp::List& scene,
                                      const Rcpp::List& edits) {
  PreviewSceneSettings settings;
  Rcpp::IntegerVector shapes = scene["shape"];
  Rcpp::List info = scene["shape_info"];
  std::set<PreviewObjectKey> seen;
  auto index = [](const Rcpp::List& value, const char* name, size_t limit) {
    if (!value.containsElementNamed(name)) {
      throw std::runtime_error(std::string("Missing editor index: ") + name);
    }
    const double number = Rcpp::as<double>(value[name]);
    if (!std::isfinite(number) || number < 1 || number != std::floor(number) ||
        number > double(limit)) {
      throw std::runtime_error(std::string("Invalid editor index: ") + name);
    }
    return size_t(number - 1);
  };
  for (R_xlen_t i = 0; i < edits.size(); ++i) {
    Rcpp::List entry = edits[i];
    const size_t row = index(entry, "row", shapes.size());
    const bool instance = shapes[row] == 15;
    Rcpp::List source, properties;
    size_t placements = 1;
    if (instance) {
      Rcpp::List shape = info[row];
      properties = Rcpp::as<Rcpp::List>(shape["shape_properties"]);
      placements = Rf_xlength(properties["x_values"]);
      source = Rcpp::as<Rcpp::List>(properties["original_scene"])[0];
    }
    const PreviewObjectKey key{row, index(entry, "instance", placements)};
    if (!seen.insert(key).second) {
      throw std::runtime_error("Duplicate object in saved editor changes.");
    }
    Rcpp::NumericMatrix matrix = entry["transform"];
    if (matrix.nrow() != 4 || matrix.ncol() != 4) {
      throw std::runtime_error("Saved editor transforms must be 4 x 4 matrices.");
    }
    std::array<float, 16> values;
    for (int r = 0; r < 4; ++r) {
      for (int c = 0; c < 4; ++c) {
        values[r + 4 * c] = matrix(r, c);
      }
    }
    settings.transforms[key] = PreviewMatrix(values);
    Rcpp::List materials = entry["materials"];
    for (R_xlen_t j = 0; j < materials.size(); ++j) {
      Rcpp::List material = materials[j];
      const size_t slot = index(material, "slot", std::numeric_limits<int>::max());
      if (settings.materials[key].count(slot)) {
        throw std::runtime_error("Duplicate material slot in saved editor changes.");
      }
      PreviewMaterialEdit edit;
      edit.type = Rcpp::as<std::string>(material["type"]);
      Rcpp::List fields = material["values"];
      Rcpp::CharacterVector names = fields.names();
      if (names.size() != fields.size()) {
        throw std::runtime_error("Saved material options must be named.");
      }
      for (R_xlen_t k = 0; k < fields.size(); ++k) {
        PreviewField field;
        field.name = Rcpp::as<std::string>(names[k]);
        SEXP value = fields[k];
        field.text_input = TYPEOF(value) == STRSXP;
        if (field.text_input) {
          field.text = Rcpp::as<std::string>(value);
        } else {
          Rcpp::NumericVector numbers = Rcpp::as<Rcpp::NumericVector>(value);
          if (numbers.size() < 1 || numbers.size() > 3) {
            throw std::runtime_error(
                "Saved material values need one to three numbers.");
          }
          field.count = numbers.size();
          std::copy(numbers.begin(), numbers.end(), field.values.begin());
        }
        edit.fields.push_back(std::move(field));
      }
      // Record() checks the actual material type, slot and option schema, then
      // uses the same setters and texture loading as an interactive commit.
      settings.materials[key][slot] = std::move(edit);
    }
    if (entry.containsElementNamed("children")) {
      if (!instance) {
        throw std::runtime_error("Only instances can contain saved child changes.");
      }
      settings.children[key] = ReadPreviewEdits(source, entry["children"]);
    }
  }
  return settings;
}
}

void PreviewScene::ImportEdits(const Rcpp::List& scene, const Rcpp::List& edits) {
  // Build and validate the complete settings tree before replacing live state.
  settings = ReadPreviewEdits(scene, edits);
}
