#include "model_animation.h"
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <unordered_set>

namespace animation {
void ModelAsset::Validate() {
  if (id.empty() || nodes.empty()) throw std::runtime_error("A model requires an asset ID and at least one node.");
  std::unordered_set<std::string> ids;
  for (const auto &node : nodes) {
    if (node.id.empty() || !ids.insert(node.id).second)
      throw std::runtime_error("Model node IDs must be nonempty and unique.");
    if (node.parent < -1 || node.parent >= int(nodes.size()))
      throw std::runtime_error("Model parent target does not exist.");
    if (node.rest.clock) throw std::runtime_error("A model rest transform cannot have a player.");
    for (const auto &channel : {node.rest.position, node.rest.rotation, node.rest.scale})
      if (channel && channel->times.size() != 1)
        throw std::runtime_error("Model rest transforms require constant one-key channels.");
    if (node.skin.joints.size() != node.skin.inverse_bind.size())
      throw std::runtime_error("Inverse bind count differs from skin joints.");
    for (size_t joint : node.skin.joints) if (joint >= nodes.size())
      throw std::runtime_error("Skin joint node does not exist.");
    for (double weight : node.weights) if (!std::isfinite(weight))
      throw std::runtime_error("Nonfinite default morph weight.");
    for (const auto &mesh : node.meshes) {
      mesh->Validate(node.skin.joints.size());
      if (mesh->targets.size() != node.weights.size())
        throw std::runtime_error("Default weight count differs from morph targets.");
      if (!node.skin.joints.empty() && mesh->influences.empty())
        throw std::runtime_error("Skinned primitive is missing JOINTS/WEIGHTS attributes.");
    }
  }
  order.clear();
  std::vector<unsigned char> state(nodes.size());
  // Iterative parent walks avoid stack overflow on deeply nested input graphs.
  for (size_t node = 0; node < nodes.size(); ++node) {
    std::vector<size_t> chain;
    int current = int(node);
    while (current >= 0 && state[current] == 0) {
      state[current] = 1;
      chain.push_back(current);
      current = nodes[current].parent;
    }
    if (current >= 0 && state[current] == 1)
      throw std::runtime_error("Model hierarchy contains a cycle.");
    for (auto i = chain.rbegin(); i != chain.rend(); ++i) {
      state[*i] = 2;
      order.push_back(*i);
    }
  }
  ids.clear();
  for (const auto &clip : clips) {
    if (clip.id.empty() || !ids.insert(clip.id).second)
      throw std::runtime_error("Model clip IDs must be nonempty and unique.");
    if (!std::isfinite(clip.begin) || !std::isfinite(clip.end) || clip.end < clip.begin ||
        clip.nodes.size() != nodes.size())
      throw std::runtime_error("Invalid model clip interval or node count.");
    if (!clip.weights.empty()) {
      if (clip.weights.size() != nodes.size()) throw std::runtime_error("Invalid morph channel node count.");
      for (size_t i = 0; i < nodes.size(); ++i) if (!clip.weights[i].empty()) {
        if (clip.weights[i].size() != nodes[i].weights.size())
          throw std::runtime_error("Morph animation width differs from targets.");
        for (const auto &track : clip.weights[i]) if (!track || track->property != Property::Scalar)
          throw std::runtime_error("Morph channels require scalar tracks.");
      }
    }
  }
}

size_t ModelAsset::NodeIndex(const std::string &id) const {
  for (size_t i = 0; i < nodes.size(); ++i) if (nodes[i].id == id) return i;
  throw std::runtime_error("Unknown model node ID: " + id);
}

size_t ModelAsset::ClipIndex(const std::string &selection) const {
  for (size_t i = 0; i < clips.size(); ++i) if (clips[i].id == selection) return i;
  size_t result = clips.size();
  for (size_t i = 0; i < clips.size(); ++i) if (clips[i].name == selection) {
    if (result != clips.size()) throw std::runtime_error("Ambiguous clip name; select a clip ID.");
    result = i;
  }
  if (result == clips.size()) throw std::runtime_error("Unknown model clip: " + selection);
  return result;
}

ModelInstance::ModelInstance(std::shared_ptr<const ModelAsset> asset, const ModelPlayer &player)
  : asset(std::move(asset)) {
  const ModelClip *clip = player.clip.empty() ? nullptr : &this->asset->clips[this->asset->ClipIndex(player.clip)];
  clock = std::make_shared<const ClipClock>(clip ? clip->begin : 0, clip ? clip->end : 0,
    player.start, player.offset, player.speed, player.loop);
  if (player.in_place && (!clip || player.root.empty()))
    throw std::runtime_error("In-place playback requires a selected clip and an explicit locomotion root ID.");
  size_t root = player.in_place ? this->asset->NodeIndex(player.root) : this->asset->nodes.size();
  if (player.in_place && this->asset->nodes[root].parent != -1)
    throw std::runtime_error("The locomotion root must be a top-level model node; axes are in asset space.");
  for (size_t i = 0; i < this->asset->nodes.size(); ++i) {
    auto track = std::make_shared<TransformTrack>(clip ? clip->nodes[i] : this->asset->nodes[i].rest);
    track->clock = clock;
    if (i == root) {
      auto removal = std::make_shared<RootMotion>();
      removal->translation = player.remove_translation;
      removal->position = track->position ? track->position->At(clip->begin) : Value{};
      removal->yaw = player.remove_yaw;
      if (removal->yaw && track->rotation) track->rotation->ValidateYaw();
      Value q = track->rotation ? track->rotation->At(clip->begin) : Value{0, 0, 0, 1};
      double norm = std::hypot(q[1], q[3]);
      if (removal->yaw && norm < 1e-8) throw std::runtime_error("Locomotion yaw is singular at the reference pose.");
      removal->reference_yaw = removal->yaw ? Value{0, q[1]/norm, 0, q[3]/norm} : Value{0, 0, 0, 1};
      track->root_motion = removal;
    }
    nodes.push_back(track);
    weights.push_back(clip && !clip->weights.empty() ? clip->weights[i] :
                     std::vector<std::shared_ptr<const Track>>{});
  }
}

ModelPose ModelInstance::At(double scene_time, const Transform &placement) const {
  ModelPose pose;
  pose.clip_time = clock->At(scene_time);
  pose.placement = placement;
  pose.local.resize(nodes.size());
  pose.global.resize(nodes.size());
  pose.weights.resize(nodes.size());
  for (size_t i : asset->order) {
    pose.local[i] = nodes[i]->At(scene_time);
    int parent = asset->nodes[i].parent;
    pose.global[i] = parent < 0 ? pose.local[i] : pose.global[parent] * pose.local[i];
    pose.weights[i] = asset->nodes[i].weights;
    for (size_t j = 0; j < weights[i].size(); ++j)
      pose.weights[i][j] = weights[i][j]->At(pose.clip_time)[0];
  }
  return pose;
}

DeformedMesh ModelInstance::Deform(const ModelPose &pose, size_t node, size_t primitive) const {
  const auto &source = asset->nodes.at(node);
  std::vector<Transform> palette;
  palette.reserve(source.skin.joints.size());
  for (size_t j = 0; j < source.skin.joints.size(); ++j)
    palette.push_back(Inverse(pose.global.at(node)) * pose.global.at(source.skin.joints[j]) * source.skin.inverse_bind[j]);
  return DeformMesh(*source.meshes.at(primitive), pose.weights.at(node), palette);
}

ModelChanges ComparePoses(const ModelPose &before, const ModelPose &after) {
  return {!(before.placement == after.placement), before.local != after.local || before.weights != after.weights};
}
} // namespace animation
