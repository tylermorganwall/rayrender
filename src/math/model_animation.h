#ifndef RAYRENDER_MODEL_ANIMATION_H
#define RAYRENDER_MODEL_ANIMATION_H

#include "keyframe.h"
#include "model_deformation.h"
#include <string>

namespace animation {
// IDs are asset-local identities, independent of labels and scene row order.
struct ModelNode {
  std::string id, name;
  int parent = -1;
  TransformTrack rest;
  ModelSkin skin;
  std::vector<std::shared_ptr<const DeformationMesh>> meshes;
  std::vector<double> weights;
};
struct ModelClip {
  std::string id, name;
  double begin = 0, end = 0;
  std::vector<TransformTrack> nodes;
  std::vector<std::vector<std::shared_ptr<const Track>>> weights;
};
struct ModelAsset {
  std::string id;
  std::vector<ModelNode> nodes;
  std::vector<ModelClip> clips;
  std::vector<size_t> order;
  // Validate once before publication; order visits parents before children.
  void Validate();
  size_t NodeIndex(const std::string &id) const;
  size_t ClipIndex(const std::string &id_or_name) const;
};
struct ModelPose {
  double clip_time = 0;
  Transform placement;
  std::vector<Transform> local, global;
  std::vector<std::vector<double>> weights;
};
struct ModelPlayer {
  std::string clip, root;
  double start = 0, offset = 0, speed = 1;
  Loop loop = Loop::Hold;
  bool in_place = false, remove_yaw = false;
  std::array<bool, 3> remove_translation{};
};
// Each instance owns its clock and transform descriptors. The descriptors share
// immutable property tracks with the asset; no pose cursor or R data is retained.
class ModelInstance {
public:
  ModelInstance(std::shared_ptr<const ModelAsset> asset, const ModelPlayer &player);
  ModelPose At(double scene_time, const Transform &placement = Transform()) const;
  DeformedMesh Deform(const ModelPose &, size_t node, size_t primitive) const;
  const std::shared_ptr<const ModelAsset> asset;
  std::shared_ptr<const ClipClock> clock;
  std::vector<std::shared_ptr<const TransformTrack>> nodes;
  std::vector<std::vector<std::shared_ptr<const Track>>> weights;
};
// Placement and local pose have separate dependencies. A future deformation
// consumer can retain its geometry when only placement differs.
struct ModelChanges { bool placement, pose; };
ModelChanges ComparePoses(const ModelPose &before, const ModelPose &after);
} // namespace animation
#endif
