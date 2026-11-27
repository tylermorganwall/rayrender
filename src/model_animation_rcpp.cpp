#include "model_animation_rcpp.h"
#include "animation_rcpp.h"
#include "model_deformation_rcpp.h"
#include <unordered_set>

std::shared_ptr<const animation::ModelAsset> read_model_asset(SEXP input) {
  if (TYPEOF(input) != ENVSXP) Rcpp::stop("A model asset must be an owned authoring environment.");
  Rcpp::Environment source(input);
  if (!source.exists("version") || (Rcpp::as<int>(source["version"]) != 1 && Rcpp::as<int>(source["version"]) != 2))
    Rcpp::stop("Unsupported model asset schema.");
  auto asset = std::make_shared<animation::ModelAsset>();
  asset->id = Rcpp::as<std::string>(source["id"]);
  Rcpp::List nodes = source["nodes"];
  for (SEXP entry : nodes) {
    Rcpp::List node(entry);
    animation::ModelNode value;
    value.id = Rcpp::as<std::string>(node["id"]);
    value.name = Rcpp::as<std::string>(node["name"]);
    value.rest = *read_transform_track(Rcpp::as<Rcpp::List>(node["transform"]));
    asset->nodes.push_back(std::move(value));
  }
  for (size_t i = 0; i < asset->nodes.size(); ++i) {
    Rcpp::List node = nodes[i];
    if (!Rf_isNull(node["parent"])) asset->nodes[i].parent = asset->NodeIndex(Rcpp::as<std::string>(node["parent"]));
    read_model_deformation(node, asset->nodes[i], *asset);
  }
  Rcpp::List clips = source["clips"];
  for (SEXP entry : clips) {
    Rcpp::List authored(entry), channels = authored["channels"];
    animation::ModelClip clip;
    clip.id = Rcpp::as<std::string>(authored["id"]);
    clip.name = Rcpp::as<std::string>(authored["name"]);
    for (const auto &node : asset->nodes) clip.nodes.push_back(node.rest);
    if (channels.size() && Rf_isNull(channels.names())) Rcpp::stop("Clip channels require node IDs.");
    Rcpp::CharacterVector targets = channels.size() ? Rcpp::CharacterVector(channels.names()) : Rcpp::CharacterVector(0);
    std::unordered_set<std::string> seen;
    double begin = INFINITY, end = -INFINITY;
    for (R_xlen_t i = 0; i < channels.size(); ++i) {
      std::string target = Rcpp::as<std::string>(targets[i]);
      if (!seen.insert(target).second) Rcpp::stop("Duplicate clip target; combine independent properties in one transform_track().");
      size_t node = asset->NodeIndex(target);
      auto track = read_transform_track(Rcpp::as<Rcpp::List>(channels[i]));
      if (track->clock) Rcpp::stop("Clip channels use the instance clock, not a nested player.");
      if (track->pivot != animation::Value{}) Rcpp::stop("Set local pivots on model_node(), not clip channels.");
      if (track->position) clip.nodes[node].position = track->position;
      if (track->rotation) clip.nodes[node].rotation = track->rotation;
      if (track->scale) clip.nodes[node].scale = track->scale;
      for (const auto &property : {track->position, track->rotation, track->scale}) if (property) {
        begin = std::min(begin, property->times.front());
        end = std::max(end, property->times.back());
      }
    }
    clip.weights.resize(asset->nodes.size());
    if (authored.containsElementNamed("weights")) {
      Rcpp::List weights = authored["weights"];
      Rcpp::CharacterVector targets = weights.size() ? Rcpp::CharacterVector(weights.names()) : Rcpp::CharacterVector(0);
      for (R_xlen_t i = 0; i < weights.size(); ++i) {
        size_t node = asset->NodeIndex(Rcpp::as<std::string>(targets[i]));
        Rcpp::List tracks = weights[i];
        for (SEXP track : tracks) {
          auto property = read_keyframe_track(Rcpp::List(track));
          clip.weights[node].push_back(property);
          begin = std::min(begin, property->times.front());
          end = std::max(end, property->times.back());
        }
      }
    }
    clip.begin = std::isfinite(begin) ? begin : 0;
    clip.end = std::isfinite(end) ? end : 0;
    asset->clips.push_back(std::move(clip));
  }
  asset->Validate();
  return asset;
}

std::shared_ptr<const animation::ModelInstance> read_model_instance(
  const Rcpp::List &input, std::shared_ptr<const animation::ModelAsset> asset) {
  if (Rcpp::as<int>(input["version"]) != 1) Rcpp::stop("Unsupported clip player schema.");
  animation::ModelPlayer player;
  if (!Rf_isNull(input["clip"])) player.clip = Rcpp::as<std::string>(input["clip"]);
  auto scalar = [&](const char *name) {
    SEXP value = input[name];
    if (!Rf_isNumeric(value) || Rf_xlength(value) != 1) Rcpp::stop("Clip clock fields must be finite scalars.");
    return Rcpp::as<double>(value);
  };
  player.start = scalar("start_time"); player.offset = scalar("offset"); player.speed = scalar("speed");
  std::string mode = Rcpp::as<std::string>(input["loop"]);
  if (mode == "hold") player.loop = animation::Loop::Hold;
  else if (mode == "repeat") player.loop = animation::Loop::Repeat;
  else if (mode == "pingpong") player.loop = animation::Loop::PingPong;
  else Rcpp::stop("Unknown clip loop mode.");
  std::string policy = Rcpp::as<std::string>(input["root_motion"]);
  if (policy != "preserve" && policy != "in_place") Rcpp::stop("Unsupported root-motion policy; extract is not implemented.");
  player.in_place = policy == "in_place";
  if (player.in_place) {
    if (!input.containsElementNamed("locomotion_root") || Rf_isNull(input["locomotion_root"]))
      Rcpp::stop("In-place playback requires an explicit locomotion root ID.");
    player.root = Rcpp::as<std::string>(input["locomotion_root"]);
    Rcpp::CharacterVector axes = input["translation_axes"];
    for (auto axis : axes) {
      std::string value = Rcpp::as<std::string>(axis);
      if (value != "x" && value != "y" && value != "z") Rcpp::stop("Unknown locomotion translation axis.");
      player.remove_translation[value == "x" ? 0 : value == "y" ? 1 : 2] = true;
    }
    player.remove_yaw = Rcpp::as<bool>(input["remove_yaw"]);
  }
  return std::make_shared<const animation::ModelInstance>(std::move(asset), player);
}

namespace {
Rcpp::NumericMatrix matrix_result(const Transform &transform) {
  Rcpp::NumericMatrix result(4, 4);
  const auto &m = transform.GetMatrix();
  for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c) result(r,c) = m.m[r][c];
  return result;
}
}

// The headless oracle shares the compiled evaluator used by scene construction.
// No native pointer is serialized; an R asset is compiled at the owning boundary.
// [[Rcpp::export]]
Rcpp::List model_pose_native(SEXP asset, Rcpp::List player, Rcpp::List placement,
                                Rcpp::NumericVector times) {
  auto instance = read_model_instance(player, read_model_asset(asset));
  auto path = read_transform_track(placement);
  Rcpp::List result(times.size());
  Rcpp::CharacterVector ids(instance->nodes.size());
  for (size_t i = 0; i < instance->nodes.size(); ++i) ids[i] = instance->asset->nodes[i].id;
  animation::ModelPose previous;
  for (R_xlen_t frame = 0; frame < times.size(); ++frame) {
    auto pose = instance->At(times[frame], path->At(times[frame]));
    Rcpp::List local(ids.size()), global(ids.size()), world(ids.size());
    for (size_t i = 0; i < instance->nodes.size(); ++i) {
      local[i] = matrix_result(pose.local[i]);
      global[i] = matrix_result(pose.global[i]);
      world[i] = matrix_result(pose.placement * pose.global[i]);
    }
    local.attr("names") = ids; global.attr("names") = ids; world.attr("names") = ids;
    auto dirty = frame ? animation::ComparePoses(previous, pose) : animation::ModelChanges{true, true};
    result[frame] = Rcpp::List::create(Rcpp::_ ["clip_time"] = pose.clip_time,
      Rcpp::_ ["local"] = local, Rcpp::_ ["global"] = global, Rcpp::_ ["world"] = world,
      Rcpp::_ ["weights"] = pose.weights,
      Rcpp::_ ["placement_changed"] = dirty.placement, Rcpp::_ ["pose_changed"] = dirty.pose);
    previous = std::move(pose);
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::List animation_model_surface(SEXP asset, Rcpp::List player, double time) {
  auto instance = read_model_instance(player, read_model_asset(asset));
  auto pose = instance->At(time);
  Rcpp::List result;
  for (size_t node=0;node<instance->asset->nodes.size();++node) {
    Rcpp::List meshes;
    for (size_t primitive=0;primitive<instance->asset->nodes[node].meshes.size();++primitive) {
      auto mesh = instance->Deform(pose,node,primitive);
      Rcpp::NumericMatrix positions(mesh.positions.size(),3), normals(mesh.normals.size(),3);
      for (size_t i=0;i<mesh.positions.size();++i) for (int k=0;k<3;++k) positions(i,k)=mesh.positions[i][k];
      for (size_t i=0;i<mesh.normals.size();++i) for (int k=0;k<3;++k) normals(i,k)=mesh.normals[i][k];
      meshes.push_back(Rcpp::List::create(Rcpp::_["positions"]=positions,Rcpp::_["normals"]=normals));
    }
    if (meshes.size()) result[instance->asset->nodes[node].id]=meshes;
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::DataFrame animation_model_clips(SEXP input) {
  auto asset = read_model_asset(input);
  Rcpp::CharacterVector ids(asset->clips.size()), names(asset->clips.size());
  Rcpp::NumericVector begin(ids.size()), end(ids.size());
  for (size_t i = 0; i < asset->clips.size(); ++i) {
    const auto &clip = asset->clips[i];
    ids[i] = clip.id; names[i] = clip.name; begin[i] = clip.begin; end[i] = clip.end;
  }
  return Rcpp::DataFrame::create(Rcpp::_ ["id"] = ids, Rcpp::_ ["name"] = names,
    Rcpp::_ ["start"] = begin, Rcpp::_ ["end"] = end, Rcpp::_ ["duration"] = end - begin);
}
