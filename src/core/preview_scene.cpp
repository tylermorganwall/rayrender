#include "preview_scene.h"
#include "../materials/material.h"
#include <algorithm>
#include <cmath>
#include <limits>

Transform PreviewMatrix(const std::array<float,16>& values) {
  Float m[4][4];
  for(int row=0;row<4;++row)for(int col=0;col<4;++col) {
    const double v=values[row+4*col];
    if(!std::isfinite(v) || std::abs(v)>1e12)throw std::runtime_error("Transform values must be finite and bounded.");
    m[row][col]=v;
  }
  if(std::abs(m[3][0])+std::abs(m[3][1])+std::abs(m[3][2])>1e-6 || std::abs(m[3][3]-1)>1e-6)
    throw std::runtime_error("Object transforms must be affine.");
  const vec3f x(m[0][0],m[1][0],m[2][0]),y(m[0][1],m[1][1],m[2][1]),z(m[0][2],m[1][2],m[2][2]);
  const double volume=x.length()*y.length()*z.length();
  if(volume<1e-15 || std::abs(dot(x,cross(y,z)))<volume*1e-6)
    throw std::runtime_error("Scale must stay nonzero; singular transforms cannot be applied.");
  return Transform(Matrix4x4(m));
}
void PreviewMatrix(const Transform& t,std::array<float,16>& values) {
  const auto& m=t.GetMatrix();
  for(int row=0;row<4;++row)for(int col=0;col<4;++col)values[row+4*col]=m.m[row][col];
}
Transform PreviewNumericTransform(const PreviewObjectState& ui) {
  for(int i=0;i<3;++i)if(!std::isfinite(ui.translation[i]) || !std::isfinite(ui.rotation[i]) || !std::isfinite(ui.scale[i]) || std::abs(ui.scale[i])<.0001 || std::abs(ui.scale[i])>10000)
    throw std::runtime_error("Use finite transforms and nonzero scales with magnitudes between 0.0001 and 10000.");
  auto t=Translate(vec3f(ui.translation[0],ui.translation[1],ui.translation[2]))*
    RotateZ(ui.rotation[2])*RotateY(ui.rotation[1])*RotateX(ui.rotation[0])*
    Scale(ui.scale[0],ui.scale[1],ui.scale[2]);
  std::array<float,16> values;PreviewMatrix(t,values);return PreviewMatrix(values);
}
void PreviewDecompose(PreviewObjectState& ui) {
  auto& m=ui.model;
  for(int i=0;i<3;++i) {
    ui.translation[i]=m[12+i];
    ui.scale[i]=std::sqrt(m[4*i]*m[4*i]+m[4*i+1]*m[4*i+1]+m[4*i+2]*m[4*i+2]);
  }
  if(*std::min_element(ui.scale.begin(),ui.scale.end())<1e-12)return;
  if(dot(vec3f(m[0],m[1],m[2]),cross(vec3f(m[4],m[5],m[6]),vec3f(m[8],m[9],m[10])))<0)ui.scale[0]=-ui.scale[0];
  const double y=std::asin(std::clamp(-m[2]/ui.scale[0],-1.0,1.0));
  const bool singular=std::abs(std::cos(y))<1e-6;
  ui.rotation={180/M_PI*(singular?std::atan2(-m[9]/ui.scale[2],m[5]/ui.scale[1]):std::atan2(m[6]/ui.scale[1],m[10]/ui.scale[2])),
               180/M_PI*y,180/M_PI*(singular?0:std::atan2(m[1]/ui.scale[0],m[0]/ui.scale[0]))};
}

PreviewScene::PreviewScene(const Rcpp::List& scene) {
  Rcpp::NumericVector x=scene["x"];groups.resize(x.size(),0);
  if(scene.containsElementNamed("preview_groups"))groups=Rcpp::as<std::vector<int>>(scene["preview_groups"]);
}
Transform PreviewScene::ObjectTransform(size_t row,size_t placement) const {
  auto found=settings.transforms.find({row,placement});return found==settings.transforms.end()?Transform():found->second;
}
bool PreviewScene::UniqueInstance(size_t row,size_t placement) const {
  return settings.materials.count({row,placement})!=0;
}
void PreviewScene::Record(size_t row,size_t placement,std::shared_ptr<hitable> object,const Transform& frame) {
  PreviewObjectKey key{row,placement};
  auto edits=settings.materials.find(key);
  if(edits!=settings.materials.end()) {
    auto materials=PreviewMaterials(object.get());
    for(const auto& entry:edits->second) {
      if(entry.first>=materials.size())throw std::runtime_error("Material slot no longer exists.");
      PreviewApplyMaterial(materials[entry.first],entry.second);
    }
  }
  uint64_t id=groups.at(row)>0 ? (UINT64_C(1)<<63)|uint64_t(groups[row]) : ((uint64_t(row)+1)<<32)|(placement+1);
  roots.push_back({key,id,std::move(object),frame});
}
void PreviewScene::DescribeMaterials(uint64_t id,PreviewObjectState& ui,material* picked) {
  ui.materials.clear();ui.material_slot=0;
  for(const auto& root:roots)if(root.id==id) {
    auto materials=PreviewMaterials(root.object.get());
    for(size_t slot=0;slot<materials.size();++slot) {
      PreviewMaterialPanel panel;panel.row=root.key.row;panel.placement=root.key.placement;panel.slot=slot;
      panel.type=materials[slot]->GetName();
      panel.label=std::to_string(root.key.row+1)+":"+std::to_string(root.key.placement+1)+" / "+panel.type+" "+std::to_string(slot+1);
      for(const auto& binding:PreviewMaterialFields(materials[slot]))panel.fields.push_back(binding.field);
      if(materials[slot]==picked)ui.material_slot=static_cast<int32_t>(ui.materials.size());
      ui.materials.push_back(std::move(panel));
    }
  }
}
void PreviewScene::Describe(uint64_t id,PreviewObjectState& ui,material* picked) {
  ui.selected=false;ui.id=0;ui.materials.clear();ui.error.clear();ui.cancel_transform=false;
  ui.transform_pending=ui.transform_active=ui.apply_transform=ui.material_pending=ui.apply_material=ui.revert=ui.numeric_transform=false;
  const PreviewObjectRoot* first=nullptr;size_t count=0;aabb bounds;
  for(const auto& root:roots)if(root.id==id) {
    aabb b;if(root.object->bounding_box(0,1,b))bounds=count?surrounding_box(bounds,b):b;
    if(!first)first=&root;++count;
  }
  selected_id=id;
  if(!first)return;
  ui.selected=true;ui.id=id;
  const bool group=(id>>63)!=0;
  ui.label=group?"Group "+std::to_string(id&~(UINT64_C(1)<<63))+" ("+std::to_string(count)+" objects)":
    first->object->GetName()+" "+std::to_string(first->key.row+1)+" / "+std::to_string(first->key.placement+1);
  selected_model=group?Translate(convert_to_vec3(bounds.Centroid())):first->frame;
  PreviewMatrix(selected_model,ui.model);PreviewDecompose(ui);
  DescribeMaterials(id,ui,picked);
}
bool PreviewScene::Pick(const Ray& ray,PreviewObjectState& ui,const std::function<bool()>& cancel) {
  Float closest=std::numeric_limits<Float>::infinity();uint64_t id=0;material* picked=nullptr;
  random_gen rng(1);
  for(size_t i=0;i<roots.size();++i) {
    if(i%64==0) { if(cancel && cancel())return false;Rcpp::checkUserInterrupt(); }
    hit_record rec;
    if(roots[i].object->hit(ray,.001f,closest,rec,rng) && !rec.alpha_miss && std::isfinite(rec.t)) {
      closest=rec.t;id=roots[i].id;picked=rec.mat_ptr;
    }
  }
  Describe(id,ui,picked);return id!=0;
}
bool PreviewScene::Apply(PreviewObjectState& ui) {
  if(!ui.selected || (!ui.apply_transform && !ui.apply_material && !ui.revert))return false;
  const bool transform=ui.apply_transform, material_edit=ui.apply_material, revert=ui.revert;
  ui.apply_transform=ui.apply_material=ui.revert=false;
  try {
    PreviewScene next=*this;next.roots.clear();
    Transform model=selected_model;
    if(transform) {
      model=ui.numeric_transform?PreviewNumericTransform(ui):PreviewMatrix(ui.model);
      auto delta=model*Inverse(selected_model);
      for(const auto& root:roots)if(root.id==ui.id)next.settings.transforms[root.key]=delta*ObjectTransform(root.key.row,root.key.placement);
    }
    if(material_edit)for(const auto& panel:ui.materials)if(panel.changed)
      next.settings.materials[{panel.row,panel.placement}][panel.slot]={panel.type,panel.fields};
    if(revert)for(const auto& root:roots)if(root.id==ui.id) {
      next.settings.transforms.erase(root.key);next.settings.materials.erase(root.key);
    }
    if(!prepare_rebuild)throw std::runtime_error("Scene editing is not connected to this renderer.");
    auto publish=prepare_rebuild(next); // Build off to the side; errors preserve the live scene.
    PreviewObjectState updated=ui;
    if(revert)next.Describe(ui.id,updated);
    else {
      const int slot=updated.material_slot;next.DescribeMaterials(ui.id,updated,nullptr);
      updated.material_slot=std::min(slot,static_cast<int>(updated.materials.size())-1);
      PreviewMatrix(model,updated.model);PreviewDecompose(updated);
    }
    publish(); // Workers are drained; everything below is a swap or scalar update.
    settings.transforms.swap(next.settings.transforms);settings.materials.swap(next.settings.materials);roots.swap(next.roots);
    selected_model=revert?next.selected_model:model;
    ui=std::move(updated);ui.error.clear();ui.transform_pending=ui.material_pending=ui.numeric_transform=false;
    return true;
  } catch(const Rcpp::internal::InterruptedException&) { throw; }
  catch(const std::exception& error) {
    ui.error=error.what();PreviewMatrix(selected_model,ui.model);PreviewDecompose(ui);
    ui.transform_pending=ui.numeric_transform=false;return false;
  }
}
Rcpp::List PreviewScene::ExportEdits() const {
  Rcpp::List result;
  for(const auto& root:roots)if(settings.transforms.count(root.key) || settings.materials.count(root.key)) {
    Rcpp::NumericMatrix matrix(4,4);auto transform=ObjectTransform(root.key.row,root.key.placement);const auto& m=transform.GetMatrix();
    for(int i=0;i<4;++i)for(int j=0;j<4;++j)matrix(i,j)=m.m[i][j];
    Rcpp::List materials;
    auto edits=settings.materials.find(root.key);
    if(edits!=settings.materials.end())for(const auto& slot:edits->second) {
      Rcpp::List fields;
      for(const auto& field:slot.second.fields)fields[field.name]=Rcpp::NumericVector(field.values.begin(),field.values.begin()+field.count);
      materials.push_back(Rcpp::List::create(Rcpp::_["slot"]=slot.first+1,Rcpp::_["type"]=slot.second.type,Rcpp::_["values"]=fields));
    }
    result.push_back(Rcpp::List::create(Rcpp::_["row"]=root.key.row+1,Rcpp::_["instance"]=root.key.placement+1,
                                       Rcpp::_["transform"]=matrix,Rcpp::_["materials"]=materials));
  }
  return result;
}

void PreviewScene::CancelTransform(PreviewObjectState& ui) {
  PreviewMatrix(selected_model,ui.model);PreviewDecompose(ui);
  ui.cancel_transform=ui.transform_pending=ui.transform_active=ui.numeric_transform=ui.apply_transform=false;
}
