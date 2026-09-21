#ifndef RAYRENDER_PREVIEW_SCENE_H
#define RAYRENDER_PREVIEW_SCENE_H
#include "preview_object_state.h"
#include "../hitables/hitable.h"
#include <functional>
#include <map>

// Stable zero-based input row and instance placement; rebuilt object addresses change.
struct PreviewObjectKey {
  size_t row = 0, placement = 0;
  bool operator<(const PreviewObjectKey& b) const {
    return row < b.row || (row == b.row && placement < b.placement);
  }

  bool operator==(const PreviewObjectKey& b) const {
    return row == b.row && placement == b.placement;
  }
};
struct PreviewMaterialEdit {
  std::string type;
  std::vector<PreviewField> fields;
};
// Committed overrides replayed whenever the renderer reconstructs the scene.
struct PreviewSceneSettings {
  std::map<PreviewObjectKey, Transform> transforms;
  std::map<PreviewObjectKey, std::map<size_t, PreviewMaterialEdit>> materials;
  // Descendant overrides belong to a placement, never to the shared source.
  std::map<PreviewObjectKey, PreviewSceneSettings> children;
  bool Empty() const;
};
class MediumBoundary;
// A volume and a retained surface have separate editable slots. Boundary identity
// keeps shared density data from merging independent objects in the inspector.
struct PreviewMaterialSlot {
  material* surface = nullptr;
  MediumBoundary* volume = nullptr;
  std::string Name() const;
  const void* Identity() const;
  bool operator==(const PreviewMaterialSlot& other) const {
    return surface == other.surface && volume == other.volume;
  }
};
class PreviewScene;
// A selectable outer object, its shared group ID, and its world placement frame.
struct PreviewObjectRoot {
  PreviewObjectKey key;
  uint64_t id = 0;
  std::shared_ptr<hitable> object;
  Transform frame;
  std::shared_ptr<PreviewScene> contents;
};
// Optional consumer-owned build hook; absent from ordinary render paths.
class PreviewScene {
public:
  explicit PreviewScene(const Rcpp::List& scene);
  PreviewSceneSettings settings;
  std::vector<PreviewObjectRoot> roots;
  Transform ObjectTransform(size_t row, size_t placement) const;
  bool UniqueInstance(size_t row, size_t placement) const;
  const PreviewSceneSettings* InstanceSettings(size_t row, size_t placement) const;
  void Record(size_t row, size_t placement, std::shared_ptr<hitable> object,
              const Transform& frame, std::shared_ptr<PreviewScene> contents = nullptr);
  bool Pick(const Ray& ray, PreviewObjectState& ui,
            const std::function<bool()>& cancel = {});
  void Describe(uint64_t id, PreviewObjectState& ui, material* hit_material = nullptr);
  bool Apply(PreviewObjectState& ui);
  std::function<void()> PrepareRestore(const PreviewSceneSettings& saved, uint64_t selection);
  // Boolean visibility in display coordinates; an empty vector means cancelled.
  std::vector<uint8_t>
  SelectionMask(uint64_t id, uint32_t width, uint32_t height,
                const std::function<bool(float, float, Ray&)>& make_ray,
                const std::function<bool()>& cancel = {});
  uint64_t Revision() const {
    return revision;
  }
  void BeginTransform();
  void EndTransform();
  bool CancelTransform(PreviewObjectState& ui);
  // Build replacement geometry into the supplied candidate and return its commit
  // closure. Preparation may throw; the commit must publish without further R work
  // or allocations, after all renderer workers have drained.
  std::function<std::function<void()>(PreviewScene&)> prepare_rebuild;
  Rcpp::List ExportEdits() const;
  void ImportEdits(const Rcpp::List& scene, const Rcpp::List& edits);
  std::vector<PreviewHierarchyNode> Hierarchy() const;
  bool Contains(const PreviewObjectRoot& root, uint64_t id) const;

private:
  struct SelectionContext {
    const PreviewScene* scene;
    uint64_t id;
    Transform to_world;
    std::vector<const PreviewObjectRoot*> placements;
  };
  // Local child IDs are mapped into a stable namespace for each placement.
  // A registry avoids hash collisions and survives candidate-scene rebuilds.
  mutable std::map<std::pair<uint64_t, uint64_t>, uint64_t> nested_ids;
  mutable std::vector<std::pair<uint64_t, uint64_t>> nested_nodes;
  uint64_t NestedId(uint64_t placement, uint64_t child) const;
  SelectionContext Resolve(uint64_t id) const;
  void ReplayMaterials();
  void ForwardMaterial(const PreviewObjectRoot& root, PreviewSceneSettings& next,
                       const PreviewMaterialSlot& target,
                       const PreviewMaterialEdit& edit) const;
  void PrepareSelectionBvh();
  const PreviewObjectRoot* HitRoot(const Ray& ray, material*& picked, random_gen& rng);
  std::vector<uint64_t> HitPath(const Ray& ray, material*& picked, random_gen& rng,
                                const std::function<bool()>& cancel);
  std::vector<int> groups, shapes;
  std::vector<std::vector<int>> group_paths;
  uint64_t selected_id = 0, revision = 0;
  std::shared_ptr<hitable> selection_bvh;
  Transform selected_model;
  std::shared_ptr<PreviewSceneSettings> drag_settings;
  Transform drag_model;
  uint64_t drag_revision = 0;
  void DescribeMaterials(uint64_t id, PreviewObjectState& ui, material* hit_material);
};

// Setters refer to native material storage and are short-lived renderer-side
// bindings. Only the value-only field snapshot is copied into the GUI state.
struct PreviewMaterialBinding {
  PreviewField field;
  std::function<void(const std::array<double, 3>&)> set;
  std::function<void(const std::string&)> set_text;
};
std::vector<material*> PreviewMaterials(hitable* root);
std::vector<PreviewMaterialSlot> PreviewMaterialSlots(hitable* root);
std::vector<PreviewMaterialBinding>
PreviewMaterialFields(const PreviewMaterialSlot& slot, hitable* root = nullptr);
void PreviewApplyMaterial(const PreviewMaterialSlot& slot,
                          const PreviewMaterialEdit& edit, hitable* root = nullptr);
std::vector<PreviewMaterialBinding> PreviewMaterialFields(material* mat,
                                                          hitable* root = nullptr);
void PreviewApplyMaterial(material* mat, const PreviewMaterialEdit& edit,
                          hitable* root = nullptr);
Transform PreviewMatrix(const std::array<float, 16>& matrix);
void PreviewMatrix(const Transform& transform, std::array<float, 16>& matrix);
Transform PreviewNumericTransform(const PreviewObjectState& ui);
void PreviewDecompose(PreviewObjectState& ui);
#endif
