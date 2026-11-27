#include "model_deformation_rcpp.h"
#include <set>

namespace {
using namespace Rcpp;
NumericMatrix attribute(const List &attributes, const char *name, int rows, int columns, bool required = false) {
  if (!attributes.containsElementNamed(name) || Rf_isNull(attributes[name])) {
    if (required) stop("Missing deformation attribute: %s", name);
    return NumericMatrix(0, columns);
  }
  NumericMatrix result = as<NumericMatrix>(attributes[name]);
  if (result.ncol() != columns || (rows >= 0 && result.nrow() != rows))
    stop("Invalid dimensions for deformation attribute: %s", name);
  for (double value : result) if (!std::isfinite(value)) stop("Nonfinite deformation attribute: %s", name);
  return result;
}
template<class T> std::vector<T> vectors(const NumericMatrix &matrix) {
  std::vector<T> result; result.reserve(matrix.nrow());
  for (int i=0;i<matrix.nrow();++i) result.emplace_back(matrix(i,0),matrix(i,1),matrix(i,2));
  return result;
}
}

void read_model_deformation(const Rcpp::List &input, animation::ModelNode &node,
                            const animation::ModelAsset &asset) {
  using namespace Rcpp;
  if (!input.containsElementNamed("deformation") || Rf_isNull(input["deformation"])) return;
  List source = input["deformation"];
  if (!Rf_isNull(source["skin"])) {
    List skin = source["skin"];
    CharacterVector joints = skin["joints"];
    std::set<size_t> seen;
    for (auto id : joints) {
      size_t joint = asset.NodeIndex(as<std::string>(id));
      if (!seen.insert(joint).second) stop("Skin contains duplicate joint nodes.");
      node.skin.joints.push_back(joint);
    }
    if (joints.size() == 0) stop("A skin requires joint nodes.");
    if (Rf_isNull(skin["inverse_bind"])) node.skin.inverse_bind.resize(joints.size());
    else {
      NumericMatrix binds = skin["inverse_bind"];
      if (binds.nrow() != joints.size() || binds.ncol() != 16) stop("Invalid inverse bind matrix dimensions.");
      for (int i=0;i<binds.nrow();++i) {
        Float m[4][4];
        for (int r=0;r<4;++r) for (int c=0;c<4;++c) {
          double value = binds(i,4*c+r);
          if (!std::isfinite(value)) stop("Nonfinite inverse bind matrix.");
          m[r][c] = value;
        }
        if (m[3][0]!=0 || m[3][1]!=0 || m[3][2]!=0 || m[3][3]!=1) stop("Inverse bind matrix must be affine.");
        node.skin.inverse_bind.emplace_back(m);
      }
    }
  }
  node.weights = as<std::vector<double>>(source["weights"]);
  List primitives = source["primitives"];
  for (SEXP entry : primitives) {
    List primitive(entry), attributes = primitive["attributes"], targets = primitive["targets"];
    auto mesh = std::make_shared<animation::DeformationMesh>();
    if (source.containsElementNamed("recalculate_normals"))
      mesh->recalculate_normals = as<bool>(source["recalculate_normals"]);
    auto positions = attribute(attributes,"POSITION",-1,3,true);
    const int count = positions.nrow();
    mesh->positions = vectors<point3f>(positions);
    mesh->normals = vectors<normal3f>(attribute(attributes,"NORMAL",count,3));
    auto tangents = attribute(attributes,"TANGENT",count,4);
    mesh->tangents = vectors<vec3f>(tangents);
    for (int i=0;i<tangents.nrow();++i) mesh->handedness.push_back(tangents(i,3));
    IntegerVector indices = primitive["indices"];
    for (int index : indices) {
      if (index == NA_INTEGER || index < 1 || index > count) stop("Invalid deformation triangle index.");
      mesh->indices.push_back(index-1);
    }
    for (SEXP entry : targets) {
      List target(entry);
      animation::MorphTarget value;
      value.positions = vectors<vec3f>(attribute(target,"POSITION",count,3));
      value.normals = vectors<vec3f>(attribute(target,"NORMAL",count,3));
      value.tangents = vectors<vec3f>(attribute(target,"TANGENT",count,3));
      mesh->targets.push_back(std::move(value));
    }
    // Enumerate every influence set and reject orphaned or noncontiguous sets.
    // Never silently truncate at four/eight influences per vertex.
    CharacterVector names = attributes.names();
    size_t influence_attributes = 0;
    for (auto name : names) {
      auto text = as<std::string>(name);
      if (text.rfind("JOINTS_",0)==0 || text.rfind("WEIGHTS_",0)==0) ++influence_attributes;
    }
    if (!node.skin.joints.empty()) {
      mesh->influences.resize(count);
      size_t consumed = 0;
      for (size_t set=0;;++set) {
        std::string joint_name="JOINTS_"+std::to_string(set), weight_name="WEIGHTS_"+std::to_string(set);
        if (!attributes.containsElementNamed(joint_name.c_str()) && !attributes.containsElementNamed(weight_name.c_str())) break;
        auto joints=attribute(attributes,joint_name.c_str(),count,4,true);
        auto weights=attribute(attributes,weight_name.c_str(),count,4,true);
        consumed+=2;
        for (int vertex=0;vertex<count;++vertex) for (int slot=0;slot<4;++slot) {
          double joint=joints(vertex,slot), weight=weights(vertex,slot);
          if (joint<0 || joint!=std::floor(joint) || joint>=node.skin.joints.size() || weight<0)
            stop("Invalid joint index or skin weight.");
          mesh->influences[vertex].push_back({size_t(joint),weight});
        }
      }
      if (consumed!=influence_attributes || !consumed) stop("Skin requires contiguous paired JOINTS_n/WEIGHTS_n sets.");
    }
    node.meshes.push_back(std::move(mesh));
  }
}
