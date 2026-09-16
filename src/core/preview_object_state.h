#ifndef RAYRENDER_PREVIEW_OBJECT_STATE_H
#define RAYRENDER_PREVIEW_OBJECT_STATE_H
#include <array>
#include <string>
#include <vector>
#include <cstdint>

// UI snapshots only. Renderer pointers and R objects never enter draw callbacks.
struct PreviewField {
  std::string name;
  bool color=false;
  unsigned count=1;
  std::array<double,3> values{};
  double minimum=0,maximum=1,speed=.01;
};
struct PreviewMaterialPanel {
  size_t row=0,placement=0,slot=0;
  bool changed=false;
  std::string label,type;
  std::vector<PreviewField> fields;
};
struct PreviewObjectState {
  bool enabled=false,selected=false,pick_pending=false,clear_pending=false;
  bool transform_pending=false,transform_active=false,apply_transform=false;
  bool material_pending=false,apply_material=false,revert=false;
  bool projection_valid=false;
  float pick_u=0,pick_v=0;
  uint64_t id=0;
  std::string label,error;
  int32_t operation=0,mode=1,material_slot=0;
  std::array<float,16> model{1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
  std::array<float,16> view{},projection{};
  uint32_t projection_kind=1;
  std::array<double,3> translation{},rotation{},scale{1,1,1};
  bool numeric_transform=false,cancel_transform=false;
  std::vector<PreviewMaterialPanel> materials;
};
#endif
