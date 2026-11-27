#include <Rcpp.h>
#define CGLTF_IMPLEMENTATION
#include <cgltf.h>
#include "math/quaternion.h"
#include <climits>
#include <cmath>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace {
using namespace Rcpp;

void require_gltf(bool condition, const char *message) {
  if (!condition) stop("glTF: %s", message);
}

// Subtraction/division bounds checks avoid overflow in offset + stride*count.
// cgltf's own validation follows only after these ranges are known safe.
void check_span(const cgltf_buffer_view *view, size_t offset, size_t count,
                size_t stride, size_t width) {
  require_gltf(view && count && width && stride >= width, "Invalid accessor layout.");
  require_gltf(offset <= view->size && width <= view->size - offset,
               "Accessor exceeds its buffer view.");
  require_gltf(count - 1 <= (view->size - offset - width) / stride,
               "Accessor exceeds its buffer view.");
}

cgltf_accessor sparse_indices(const cgltf_accessor &accessor) {
  cgltf_accessor result{};
  result.type = cgltf_type_scalar;
  result.component_type = accessor.sparse.indices_component_type;
  result.buffer_view = accessor.sparse.indices_buffer_view;
  result.offset = accessor.sparse.indices_byte_offset;
  result.count = accessor.sparse.count;
  result.stride = cgltf_component_size(result.component_type);
  return result;
}

void validate_storage(cgltf_data &data) {
  for (size_t i = 0; i < data.buffer_views_count; ++i) {
    const auto &view = data.buffer_views[i];
    require_gltf(!view.has_meshopt_compression, "Meshopt compression is not supported.");
    require_gltf(view.buffer && view.buffer->data, "A buffer view has no loaded buffer.");
    require_gltf(view.offset <= view.buffer->size && view.size <= view.buffer->size - view.offset,
                 "Buffer view exceeds its buffer.");
  }
  for (size_t i = 0; i < data.accessors_count; ++i) {
    const auto &a = data.accessors[i];
    size_t width = cgltf_calc_size(a.type, a.component_type);
    size_t component = cgltf_component_size(a.component_type);
    require_gltf(a.count > 0 && a.count <= INT_MAX / 16 && width && component,
                 "Invalid or oversized accessor count/type.");
    require_gltf(!a.normalized || (a.component_type != cgltf_component_type_r_32f &&
                                   a.component_type != cgltf_component_type_r_32u),
                 "Only 8-bit and 16-bit integer accessors can be normalized.");
    require_gltf(a.stride >= width && a.stride % component == 0 && a.offset % component == 0,
                 "Invalid accessor stride or alignment.");
    if (a.buffer_view) {
      require_gltf(a.buffer_view->offset % component == 0, "Unaligned accessor buffer view.");
      check_span(a.buffer_view, a.offset, a.count, a.stride, width);
    }
    if (a.is_sparse) {
      auto indices = sparse_indices(a);
      require_gltf(indices.component_type == cgltf_component_type_r_8u ||
                   indices.component_type == cgltf_component_type_r_16u ||
                   indices.component_type == cgltf_component_type_r_32u,
                   "Sparse indices must be unsigned integers.");
      require_gltf(indices.count <= a.count, "Too many sparse replacements.");
      check_span(indices.buffer_view, indices.offset, indices.count, indices.stride, indices.stride);
      check_span(a.sparse.values_buffer_view, a.sparse.values_byte_offset, indices.count, width, width);
      require_gltf(indices.offset % indices.stride == 0 && indices.buffer_view->offset % indices.stride == 0 &&
                   a.sparse.values_byte_offset % component == 0 && a.sparse.values_buffer_view->offset % component == 0,
                   "Unaligned sparse accessor.");
      require_gltf(!indices.buffer_view->stride && !a.sparse.values_buffer_view->stride,
                   "Sparse buffer views must be tightly packed.");
      size_t previous = 0;
      for (size_t j = 0; j < indices.count; ++j) {
        size_t index = cgltf_accessor_read_index(&indices, j);
        require_gltf(index < a.count && (!j || index > previous), "Sparse indices must be ordered, unique and in range.");
        previous = index;
      }
    }
  }
}

NumericMatrix unpack(const cgltf_accessor *source) {
  require_gltf(source, "Missing accessor.");
  const size_t components = cgltf_num_components(source->type);
  NumericMatrix result(source->count, components);
  auto base = *source;
  base.is_sparse = false;
  float value[16];
  auto write = [&](size_t row, const cgltf_accessor &accessor, size_t index) {
    require_gltf(cgltf_accessor_read_float(&accessor, index, value, components), "Cannot decode accessor.");
    for (size_t k = 0; k < components; ++k) {
      require_gltf(std::isfinite(value[k]), "Accessor contains nonfinite values.");
      result(row, k) = value[k];
    }
  };
  for (size_t row = 0; row < source->count; ++row) write(row, base, row);
  if (source->is_sparse) {
    auto indices = sparse_indices(*source);
    auto values = base;
    values.buffer_view = source->sparse.values_buffer_view;
    values.offset = source->sparse.values_byte_offset;
    values.stride = cgltf_calc_size(source->type, source->component_type);
    for (size_t i = 0; i < indices.count; ++i)
      write(cgltf_accessor_read_index(&indices, i), values, i);
  }
  return result;
}

List attributes(const cgltf_attribute *input, size_t count) {
  List result(count);
  CharacterVector names(count);
  std::set<std::string> seen;
  for (size_t i = 0; i < count; ++i) {
    require_gltf(input[i].name && seen.insert(input[i].name).second, "Duplicate or unnamed vertex attribute.");
    const auto &a = *input[i].data;
    if (input[i].type == cgltf_attribute_type_joints)
      require_gltf(a.type == cgltf_type_vec4 && !a.normalized &&
        (a.component_type == cgltf_component_type_r_8u || a.component_type == cgltf_component_type_r_16u),
        "Joint attributes require unnormalized unsigned byte/short VEC4 accessors.");
    if (input[i].type == cgltf_attribute_type_weights)
      require_gltf(a.type == cgltf_type_vec4 &&
        (a.component_type == cgltf_component_type_r_32f ||
         (a.normalized && (a.component_type == cgltf_component_type_r_8u || a.component_type == cgltf_component_type_r_16u))),
        "Weight attributes require float or normalized unsigned byte/short VEC4 accessors.");
    names[i] = input[i].name;
    result[i] = unpack(input[i].data);
  }
  result.attr("names") = names;
  return result;
}

List node_rest(const cgltf_node &node) {
  NumericVector position(3), rotation = NumericVector::create(0., 0., 0., 1.), scale(3, 1.);
  if (!node.has_matrix) {
    for (int k = 0; k < 3; ++k) {
      position[k] = node.translation[k]; scale[k] = node.scale[k];
    }
    for (int k = 0; k < 4; ++k) rotation[k] = node.rotation[k];
  } else {
    require_gltf(!node.has_translation && !node.has_rotation && !node.has_scale,
                 "A node cannot specify both matrix and TRS.");
    Matrix4x4 matrix;
    for (int row = 0; row < 4; ++row) for (int col = 0; col < 4; ++col) {
      float value = node.matrix[4*col+row];
      require_gltf(std::isfinite(value), "Node matrix contains nonfinite values.");
      matrix.m[row][col] = value;
    }
    require_gltf(matrix.m[3][0] == 0 && matrix.m[3][1] == 0 && matrix.m[3][2] == 0 && matrix.m[3][3] == 1,
                 "Node matrix must be affine.");
    for (int col = 0; col < 3; ++col) {
      scale[col] = std::hypot(std::hypot(matrix.m[0][col], matrix.m[1][col]), matrix.m[2][col]);
      require_gltf(scale[col] >= 1e-8, "Node matrix has singular scale.");
      position[col] = matrix.m[col][3];
    }
    const double determinant =
      matrix.m[0][0]*(matrix.m[1][1]*matrix.m[2][2]-matrix.m[1][2]*matrix.m[2][1]) -
      matrix.m[0][1]*(matrix.m[1][0]*matrix.m[2][2]-matrix.m[1][2]*matrix.m[2][0]) +
      matrix.m[0][2]*(matrix.m[1][0]*matrix.m[2][1]-matrix.m[1][1]*matrix.m[2][0]);
    if (determinant < 0) scale[0] = -scale[0];
    Matrix4x4 orientation;
    for (int row = 0; row < 3; ++row) for (int col = 0; col < 3; ++col)
      orientation.m[row][col] = matrix.m[row][col] / scale[col];
    for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) {
      double dot = 0;
      for (int k = 0; k < 3; ++k) dot += orientation.m[k][a]*orientation.m[k][b];
      require_gltf(std::abs(dot - (a == b ? 1 : 0)) < 1e-5, "Node matrix contains unsupported shear.");
    }
    Quaternion q = Normalize(Quaternion(orientation));
    for (int k = 0; k < 3; ++k) rotation[k] = q.v[k];
    rotation[3] = q.w;
  }
  return List::create(_["position"] = position, _["rotation"] = rotation, _["scale"] = scale);
}

std::string node_id(const cgltf_data &data, const cgltf_node *node) {
  return "node-" + std::to_string(node - data.nodes);
}
} // namespace

// Decode owned, serializable source data. No parser pointer survives this call.
// [[Rcpp::export]]
Rcpp::List read_gltf_native(std::string filename) {
  using namespace Rcpp;
  cgltf_options options{};
  cgltf_data *raw = nullptr;
  cgltf_result parsed = cgltf_parse_file(&options, filename.c_str(), &raw);
  std::unique_ptr<cgltf_data, decltype(&cgltf_free)> data(raw, cgltf_free);
  if (parsed != cgltf_result_success) stop("glTF parse failed (code %d).", int(parsed));
  require_gltf(!data->extensions_required_count, "Required extensions are not supported by this importer.");
  // Network URI handling is deliberately absent. Check declared resources before
  // cgltf opens files or allocates embedded buffers. R matrices use signed counts.
  for (size_t i = 0; i < data->buffers_count; ++i) {
    auto &buffer = data->buffers[i];
    require_gltf(buffer.size > 0 && buffer.size <= INT_MAX, "Buffer size exceeds the supported range.");
    if (!buffer.uri) continue;
    std::string uri(buffer.uri);
    if (uri.compare(0, 5, "data:") == 0) {
      size_t comma = uri.find(',');
      require_gltf(comma != std::string::npos && comma >= 7 && uri.substr(comma-7, 7) == ";base64",
                   "Embedded buffers require base64 data URIs.");
      size_t bytes = uri.size() - comma - 1;
      require_gltf(bytes >= ((buffer.size + 2)/3)*4, "Embedded buffer data is shorter than byteLength.");
    } else {
      require_gltf(uri.find(':') == std::string::npos && uri.compare(0, 2, "//") != 0 &&
                   uri.compare(0, 2, "\\\\") != 0, "Only local relative buffer URIs are supported; network loading is disabled.");
    }
  }
  cgltf_result loaded = cgltf_load_buffers(&options, data.get(), filename.c_str());
  if (loaded != cgltf_result_success) stop("glTF buffer load failed (code %d).", int(loaded));
  validate_storage(*data);
  for (size_t i = 0; i < data->meshes_count; ++i) {
    require_gltf(data->meshes[i].primitives_count, "Mesh requires primitives.");
    for (size_t j = 0; j < data->meshes[i].primitives_count; ++j) {
      const auto &primitive = data->meshes[i].primitives[j];
      for (size_t k = 0; k < primitive.attributes_count; ++k)
        require_gltf(primitive.attributes[k].data, "Vertex attribute requires an accessor.");
      for (size_t k = 0; k < primitive.targets_count; ++k)
        for (size_t m = 0; m < primitive.targets[k].attributes_count; ++m)
          require_gltf(primitive.targets[k].attributes[m].data, "Morph attribute requires an accessor.");
      require_gltf(!primitive.indices || !primitive.indices->normalized, "Triangle indices cannot be normalized.");
    }
  }
  for (size_t i = 0; i < data->animations_count; ++i)
    for (size_t j = 0; j < data->animations[i].channels_count; ++j) {
      const auto &channel = data->animations[i].channels[j];
      require_gltf(channel.sampler && channel.sampler->input && channel.sampler->output,
                   "Animation sampler requires input and output accessors.");
    }
  // Reject duplicate parent edges before cgltf's cycle/scene validation. Node
  // names are labels; IDs retain the document's stable zero-based node index.
  std::vector<unsigned> parents(data->nodes_count);
  for (size_t i = 0; i < data->nodes_count; ++i) {
    auto &node = data->nodes[i];
    require_gltf(!node.has_mesh_gpu_instancing, "GPU instancing extension is not supported.");
    for (size_t j = 0; j < node.children_count; ++j)
      require_gltf(++parents[node.children[j] - data->nodes] == 1, "Node has duplicate or multiple parent edges.");
  }
  cgltf_result valid = cgltf_validate(data.get());
  if (valid != cgltf_result_success) stop("glTF validation failed (code %d).", int(valid));
  List meshes(data->meshes_count), nodes(data->nodes_count), clips(data->animations_count), skins(data->skins_count);
  for (size_t i = 0; i < data->meshes_count; ++i) {
    auto &mesh = data->meshes[i];
    List primitives(mesh.primitives_count);
    for (size_t j = 0; j < mesh.primitives_count; ++j) {
      const auto &p = mesh.primitives[j];
      require_gltf(!p.has_draco_mesh_compression, "Draco compression is not supported.");
      require_gltf(p.type == cgltf_primitive_type_triangles, "Only triangle primitives are supported.");
      const auto *positions = cgltf_find_accessor(&p, cgltf_attribute_type_position, 0);
      require_gltf(positions && positions->type == cgltf_type_vec3 && positions->component_type == cgltf_component_type_r_32f,
                   "Triangle positions require a float VEC3 accessor.");
      size_t count = p.indices ? p.indices->count : positions->count;
      require_gltf(count % 3 == 0, "Triangle index count must be divisible by three.");
      IntegerVector indices(count);
      // Integer indices are never routed through float, which loses large IDs.
      auto base = p.indices ? *p.indices : cgltf_accessor{};
      base.is_sparse = false;
      for (size_t k = 0; k < count; ++k) {
        size_t index = p.indices ? cgltf_accessor_read_index(&base, k) : k;
        require_gltf(index < positions->count, "Triangle index is out of range.");
        indices[k] = index + 1;
      }
      if (p.indices && p.indices->is_sparse) {
        auto sparse = sparse_indices(*p.indices), values = base;
        values.buffer_view = p.indices->sparse.values_buffer_view;
        values.offset = p.indices->sparse.values_byte_offset;
        values.stride = cgltf_component_size(values.component_type);
        for (size_t k = 0; k < sparse.count; ++k) {
          size_t index = cgltf_accessor_read_index(&values, k);
          require_gltf(index < positions->count, "Sparse triangle index is out of range.");
          indices[cgltf_accessor_read_index(&sparse, k)] = index + 1;
        }
      }
      List targets(p.targets_count);
      for (size_t k = 0; k < p.targets_count; ++k)
        targets[k] = attributes(p.targets[k].attributes, p.targets[k].attributes_count);
      primitives[j] = List::create(_["attributes"] = attributes(p.attributes, p.attributes_count),
                                  _["indices"] = indices, _["targets"] = targets);
    }
    meshes[i] = List::create(_["primitives"] = primitives,
                             _["weights"] = mesh.weights_count ? NumericVector(mesh.weights, mesh.weights + mesh.weights_count) : NumericVector(0));
  }
  for (size_t i = 0; i < data->nodes_count; ++i) {
    auto &node = data->nodes[i];
    nodes[i] = List::create(_["id"] = node_id(*data, &node),
      _["name"] = node.name ? node.name : node_id(*data, &node),
      _["parent"] = node.parent ? wrap(node_id(*data, node.parent)) : R_NilValue,
      _["mesh"] = node.mesh ? int(node.mesh - data->meshes) + 1 : NA_INTEGER,
      _["skin"] = node.skin ? int(node.skin - data->skins) + 1 : NA_INTEGER,
      _["weights"] = node.weights_count ? NumericVector(node.weights, node.weights + node.weights_count) : NumericVector(0),
      _["rest"] = node_rest(node));
  }
  for (size_t i = 0; i < data->animations_count; ++i) {
    auto &clip = data->animations[i];
    List channels(clip.channels_count);
    std::set<std::pair<std::string, int>> targets;
    for (size_t j = 0; j < clip.channels_count; ++j) {
      auto &channel = clip.channels[j];
      require_gltf(channel.target_node && channel.sampler, "Animation requires a node target and sampler.");
      auto &sampler = *channel.sampler;
      require_gltf(sampler.input && sampler.input->type == cgltf_type_scalar &&
                   sampler.input->component_type == cgltf_component_type_r_32f,
                   "Animation times require a float SCALAR accessor.");
      require_gltf(sampler.output && sampler.output->component_type == cgltf_component_type_r_32f,
                   "Animation outputs require floats.");
      std::string target = node_id(*data, channel.target_node), property;
      cgltf_type type;
      switch (channel.target_path) {
        case cgltf_animation_path_type_translation: property = "position"; type = cgltf_type_vec3; break;
        case cgltf_animation_path_type_rotation: property = "rotation"; type = cgltf_type_vec4; break;
        case cgltf_animation_path_type_scale: property = "scale"; type = cgltf_type_vec3; break;
        case cgltf_animation_path_type_weights: property = "weights"; type = cgltf_type_scalar; break;
        default: stop("glTF: Unsupported animation target property.");
      }
      require_gltf(sampler.output->type == type, "Animation output type does not match its target.");
      require_gltf(property == "weights" || !channel.target_node->has_matrix,
                   "A TRS animation cannot target a node with an authored matrix.");
      require_gltf(targets.emplace(target, int(channel.target_path)).second, "Duplicate animation target property.");
      const char *method = sampler.interpolation == cgltf_interpolation_type_step ? "step" :
        sampler.interpolation == cgltf_interpolation_type_linear ? "linear" :
        sampler.interpolation == cgltf_interpolation_type_cubic_spline ? "cubic" : nullptr;
      require_gltf(method, "Unsupported animation interpolation.");
      NumericMatrix times = unpack(sampler.input);
      for (int k = 0; k < times.nrow(); ++k)
        require_gltf(times(k, 0) >= 0 && (!k || times(k, 0) > times(k-1, 0)), "Animation times must increase and be nonnegative.");
      channels[j] = List::create(_["target"] = target, _["property"] = property,
        _["interpolation"] = method, _["time"] = NumericVector(times.begin(), times.end()), _["value"] = unpack(sampler.output));
    }
    clips[i] = List::create(_["id"] = "clip-" + std::to_string(i),
                            _["name"] = clip.name ? clip.name : "clip-" + std::to_string(i), _["channels"] = channels);
  }
  for (size_t i = 0; i < data->skins_count; ++i) {
    auto &skin = data->skins[i];
    require_gltf(!skin.inverse_bind_matrices ||
      (skin.inverse_bind_matrices->type == cgltf_type_mat4 && skin.inverse_bind_matrices->component_type == cgltf_component_type_r_32f),
      "Inverse bind matrices require float MAT4 accessors.");
    CharacterVector joints(skin.joints_count);
    for (size_t j = 0; j < skin.joints_count; ++j) joints[j] = node_id(*data, skin.joints[j]);
    skins[i] = List::create(_["joints"] = joints,
      _["inverse_bind"] = skin.inverse_bind_matrices ? wrap(unpack(skin.inverse_bind_matrices)) : R_NilValue);
  }
  List scenes(data->scenes_count);
  for (size_t i = 0; i < data->scenes_count; ++i) {
    CharacterVector roots(data->scenes[i].nodes_count);
    for (size_t j = 0; j < data->scenes[i].nodes_count; ++j) roots[j] = node_id(*data, data->scenes[i].nodes[j]);
    scenes[i] = roots;
  }
  return List::create(_["nodes"] = nodes, _["meshes"] = meshes, _["clips"] = clips, _["skins"] = skins,
    _["scenes"] = scenes, _["default_scene"] = data->scene ? int(data->scene - data->scenes) + 1 : NA_INTEGER,
    _["materials"] = double(data->materials_count), _["cameras"] = double(data->cameras_count),
    _["lights"] = double(data->lights_count));
}
