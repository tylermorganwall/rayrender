#include "medium.h"
#include <algorithm>
#include <nanovdb/NanoVDB.h>
#include <nanovdb/math/SampleFromVoxels.h>
#include <openvdbr/openvdbr_volume_api.h>
#include <R_ext/Rdynload.h>
#include <stdexcept>

struct NanoVDBMedium::Data {
  Rcpp::Environment provider = Rcpp::Environment::namespace_env("openvdbr");
  const openvdbr_volume_api_v1 *api = nullptr;
  openvdbr_volume_buffer *density = nullptr, *temperature = nullptr;
  uint64_t density_bytes = 0, temperature_bytes = 0;
  const nanovdb::FloatGrid *d = nullptr, *t = nullptr;
  Data() {
    // Resolve on the main R thread. Rendering only reads immutable grids;
    // it performs no R calls or runtime lookups in the sampling hot path.
    auto get_api = reinterpret_cast<openvdbr_get_volume_api_v1_fn>(
        R_GetCCallable("openvdbr", "openvdbr_get_volume_api_v1"));
    api = get_api(1);
    if (!api || api->abi_version != 1 || api->struct_size < sizeof(openvdbr_volume_api_v1) ||
        api->nanovdb_major != NANOVDB_MAJOR_VERSION_NUMBER)
      throw std::runtime_error("Incompatible openvdbr volume runtime API or NanoVDB format.");
  }
  ~Data() {
    if (density) api->destroy(density);
    if (temperature) api->destroy(temperature);
  }
  const nanovdb::FloatGrid *Read(const std::string &filename, const std::string &name,
                                bool optional, openvdbr_volume_buffer *&buffer, uint64_t &size) {
    openvdbr_error error{};
    if (api->read_float_grid(filename.c_str(), name.c_str(), optional, &buffer, &error))
      throw std::runtime_error(std::string("VDB grid ") + name + ": " + error.message);
    if (!buffer) return nullptr;
    const void *bytes = nullptr;
    if (api->data(buffer, &bytes, &size, &error)) throw std::runtime_error(error.message);
    return static_cast<const nanovdb::FloatGrid *>(bytes);
  }
};
namespace {
void validate_render_grid(const nanovdb::FloatGrid *grid, const std::string &name) {
  if (grid->gridClass() != nanovdb::GridClass::FogVolume &&
      grid->gridClass() != nanovdb::GridClass::Unknown)
    throw std::runtime_error("NanoVDB grid '" + name +
                             "' must be a fog or scalar field, not a level set.");
  const auto &map = grid->mMap;
  for (int i = 0; i < 9; ++i)
    if (!std::isfinite(map.mMatD[i]) || !std::isfinite(map.mInvMatD[i]) ||
        !std::isfinite(map.mMatF[i]) || !std::isfinite(map.mInvMatF[i]))
      throw std::runtime_error("NanoVDB native transform must be finite and invertible.");
  for (int i = 0; i < 3; ++i)
    if (!std::isfinite(map.mVecD[i]) || !std::isfinite(map.mVecF[i]))
      throw std::runtime_error("NanoVDB native transform must be finite.");
  auto o = grid->indexToWorld(nanovdb::Vec3d(0));
  auto x = grid->indexToWorld(nanovdb::Vec3d(1, 0, 0)) - o;
  auto y = grid->indexToWorld(nanovdb::Vec3d(0, 1, 0)) - o;
  auto z = grid->indexToWorld(nanovdb::Vec3d(0, 0, 1)) - o;
  double determinant = x.dot(y.cross(z));
  if (!std::isfinite(determinant) || determinant == 0)
    throw std::runtime_error("NanoVDB native transform must be finite and invertible.");
  for (int a = 0; a < 3; ++a)
    if (!std::isfinite(o[a]) || !std::isfinite(Float(o[a])))
      throw std::runtime_error("NanoVDB native transform exceeds the supported coordinate range.");
}
template <class Node, class F> void visit_values(const Node &node, F &f) {
  for (auto it = node.cbeginValueAll(); it; ++it) {
    if constexpr (Node::LEVEL == 0)
      f(it.getCoord(), 1, *it);
    else
      f(it.getCoord(), Node::ChildNodeType::DIM, *it);
  }
  if constexpr (Node::LEVEL > 0)
    for (auto child = node.cbeginChild(); child; ++child)
      visit_values(*child, f);
}
Float lookup(const nanovdb::FloatGrid *grid, const point3f &p) {
  if (!grid)
    return 0;
  auto ijk = grid->worldToIndex(nanovdb::Vec3d(p[0], p[1], p[2]));
  // A temperature grid may have different bounds and transform from density.
  // Avoid undefined float-to-integer conversion when querying far outside it.
  for (int a = 0; a < 3; ++a)
    if (!std::isfinite(ijk[a]) || ijk[a] < double(INT32_MIN) + 1 || ijk[a] > double(INT32_MAX) - 1)
      return grid->tree().background();
  nanovdb::math::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false> sampler(grid->tree());
  return sampler(ijk);
}
} // namespace
NanoVDBMedium::NanoVDBMedium(const Rcpp::List &d) : Medium(d), data(new Data) {
  std::string filename = Rcpp::as<std::string>(d["filename"]);

  data->d = data->Read(filename, Rcpp::as<std::string>(d["density_grid"]), false,
                       data->density, data->density_bytes);
  validate_render_grid(data->d, Rcpp::as<std::string>(d["density_grid"]));
  if (data->d->tree().background() != 0)
    throw std::runtime_error("NanoVDB density background must be zero.");
  if (d.containsElementNamed("temperature_grid") && !Rf_isNull(d["temperature_grid"])) {
    const std::string temperature_name = Rcpp::as<std::string>(d["temperature_grid"]);
    const bool optional = d.containsElementNamed("temperature_optional") &&
                          Rcpp::as<bool>(d["temperature_optional"]);
    data->t = data->Read(filename, temperature_name, optional, data->temperature, data->temperature_bytes);
    if (data->t) {
      validate_render_grid(data->t, temperature_name);
      auto validate = [&](const nanovdb::Coord &, int, float v) {
        if (!std::isfinite(v) || v < 0 ||
            !std::isfinite(Float((double(v) - temperature_offset) * temperature_scale)))
          throw std::runtime_error("NanoVDB temperatures and their scaled values must be finite.");
      };
      validate(nanovdb::Coord(0), 1, data->t->tree().background());
      visit_values(data->t->tree().root(), validate);
    }
  }
  has_temperature = data->t != nullptr;
  majorants.resolution = 64;
  majorants.lo = point3f(INFINITY);
  majorants.hi = point3f(-INFINITY);
  // Visit stored voxels and tiles, never all voxels in the sparse bounding box.
  auto bounds = [&](const nanovdb::Coord &p, int dim, float v) {
    if (!std::isfinite(v) || v < 0)
      throw std::runtime_error("NanoVDB density must be finite and nonnegative.");
    if (v == 0)
      return;
    for (int corner = 0; corner < 8; ++corner) {
      nanovdb::Vec3d q;
      for (int a = 0; a < 3; ++a)
        q[a] = double(p[a]) + ((corner & (1 << a)) ? dim : -1);
      auto w = data->d->indexToWorld(q);
      for (int a = 0; a < 3; ++a) {
        majorants.lo[a] = std::min(majorants.lo[a], Float(w[a]));
        majorants.hi[a] = std::max(majorants.hi[a], Float(w[a]));
      }
    }
  };
  visit_values(data->d->tree().root(), bounds);
  if (!std::isfinite(majorants.lo[0])) {
    majorants.lo = point3f(-0.5);
    majorants.hi = point3f(0.5);
  }
  for (int a = 0; a < 3; ++a) {
    majorants.lo[a] = std::nextafter(majorants.lo[a], Float(-INFINITY));
    majorants.hi[a] = std::nextafter(majorants.hi[a], Float(INFINITY));
    if (!std::isfinite(majorants.lo[a]) || !std::isfinite(majorants.hi[a]) ||
        !(majorants.lo[a] < majorants.hi[a]))
      throw std::runtime_error("NanoVDB coordinates exceed the supported float range.");
  }
  int r = majorants.resolution;
  majorants.density.assign(r * r * r, 0);
  auto stamp = [&](const nanovdb::Coord &p, int dim, float v) {
    if (v == 0)
      return;
    point3f lo(INFINITY), hi(-INFINITY);
    for (int corner = 0; corner < 8; ++corner) {
      nanovdb::Vec3d q;
      for (int a = 0; a < 3; ++a)
        q[a] = double(p[a]) + ((corner & (1 << a)) ? dim : -1);
      auto w = data->d->indexToWorld(q);
      for (int a = 0; a < 3; ++a) {
        lo[a] = std::min(lo[a], Float(w[a]));
        hi[a] = std::max(hi[a], Float(w[a]));
      }
    }
    int lower[3], upper[3];
    for (int a = 0; a < 3; ++a) {
      double scale = r / double(majorants.hi[a] - majorants.lo[a]);
      lower[a] = std::clamp(int(std::floor((lo[a] - majorants.lo[a]) * scale)), 0, r - 1);
      upper[a] = std::clamp(int(std::floor((hi[a] - majorants.lo[a]) * scale)), 0, r - 1);
    }
    v = std::nextafter(v, INFINITY);
    for (int a = 0; a < 3; ++a)
      if (!std::isfinite(v * (sigma_a[a] + sigma_s[a])))
        throw std::runtime_error(
            "NanoVDB density times extinction exceeds the supported float range.");
    for (int z = lower[2]; z <= upper[2]; ++z)
      for (int y = lower[1]; y <= upper[1]; ++y)
        for (int x = lower[0]; x <= upper[0]; ++x) {
          Float &m = majorants.density[x + r * (y + r * z)];
          m = std::max(m, v);
        }
  };
  visit_values(data->d->tree().root(), stamp);
}
NanoVDBMedium::~NanoVDBMedium() = default;
Float NanoVDBMedium::Density(const point3f &p) const { return lookup(data->d, p); }
DensityIndexRay NanoVDBMedium::DensityRay(const Ray &r) const {
  auto o = data->d->worldToIndex(nanovdb::Vec3d(r.o[0], r.o[1], r.o[2]));
  auto d = data->d->worldToIndexDir(nanovdb::Vec3d(r.d[0], r.d[1], r.d[2]));
  return {{o[0], o[1], o[2]}, {d[0], d[1], d[2]}};
}
std::array<double, 8> NanoVDBMedium::DensityCorners(const std::array<double, 3> &cell) const {
  std::array<double, 8> result{};
  for (int corner = 0; corner < 8; ++corner) {
    nanovdb::Coord index;
    bool valid = true;
    for (int a = 0; a < 3; ++a) {
      double v = cell[a] + ((corner >> a) & 1);
      if (v < INT32_MIN || v > INT32_MAX) { valid = false; break; }
      index[a] = int32_t(v);
    }
    if (valid) result[corner] = data->d->tree().getValue(index);
  }
  return result;
}
point3f NanoVDBMedium::Emission(const point3f &p) const {
  return emission_scale *
         (data->t ? BlackbodyRGB((lookup(data->t, p) - temperature_offset) * temperature_scale)
                  : emission);
}
MediumProperties NanoVDBMedium::SamplePoint(const point3f &p) const {
  return Properties(Density(p), Emission(p), p);
}
RayMajorantIterator NanoVDBMedium::SampleRay(const Ray &r, double t_max) const {
  return RayMajorantIterator(r, t_max, sigma_a + sigma_s, &majorants);
}

size_t NanoVDBMedium::MemoryBytes() const {
  return sizeof(*this) + sizeof(Data) + data->density_bytes +
         data->temperature_bytes + sizeof(Float) * majorants.density.capacity();
}
