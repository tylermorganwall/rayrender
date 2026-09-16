#ifndef RAYRENDER_PREVIEW_SCENE_H
#define RAYRENDER_PREVIEW_SCENE_H
#include "preview_object_state.h"
#include "../hitables/hitable.h"
#include <functional>
#include <map>

struct PreviewObjectKey {
  size_t row=0,placement=0;
  bool operator<(const PreviewObjectKey& b) const { return row<b.row || (row==b.row && placement<b.placement); }
  bool operator==(const PreviewObjectKey& b) const { return row==b.row && placement==b.placement; }
};
struct PreviewMaterialEdit {
  std::string type;
  std::vector<PreviewField> fields;
};
struct PreviewSceneSettings {
  std::map<PreviewObjectKey,Transform> transforms;
  std::map<PreviewObjectKey,std::map<size_t,PreviewMaterialEdit>> materials;
};
struct PreviewObjectRoot {
  PreviewObjectKey key;
  uint64_t id=0;
  std::shared_ptr<hitable> object;
  Transform frame;
};
// Optional consumer-owned build hook; absent from ordinary render paths.
class PreviewScene {
public:
  explicit PreviewScene(const Rcpp::List& scene);
  PreviewSceneSettings settings;
  std::vector<PreviewObjectRoot> roots;
  Transform ObjectTransform(size_t row,size_t placement) const;
  bool UniqueInstance(size_t row,size_t placement) const;
  void Record(size_t row,size_t placement,std::shared_ptr<hitable> object,const Transform& frame);
  bool Pick(const Ray& ray,PreviewObjectState& ui,const std::function<bool()>& cancel={});
  void Describe(uint64_t id,PreviewObjectState& ui,material* hit_material=nullptr);
  bool Apply(PreviewObjectState& ui);
  void CancelTransform(PreviewObjectState& ui);
  std::function<std::function<void()>(PreviewScene&)> prepare_rebuild;
  Rcpp::List ExportEdits() const;
private:
  std::vector<int> groups;
  uint64_t selected_id=0;
  Transform selected_model;
  void DescribeMaterials(uint64_t id,PreviewObjectState& ui,material* hit_material);
};

struct PreviewMaterialBinding {
  PreviewField field;
  std::function<void(const std::array<double,3>&)> set;
};
std::vector<material*> PreviewMaterials(hitable* root);
std::vector<PreviewMaterialBinding> PreviewMaterialFields(material* mat);
void PreviewApplyMaterial(material* mat,const PreviewMaterialEdit& edit);
Transform PreviewMatrix(const std::array<float,16>& matrix);
void PreviewMatrix(const Transform& transform,std::array<float,16>& matrix);
Transform PreviewNumericTransform(const PreviewObjectState& ui);
void PreviewDecompose(PreviewObjectState& ui);
#endif
