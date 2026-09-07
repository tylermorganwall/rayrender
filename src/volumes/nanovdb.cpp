#include "medium.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <nanovdb/GridHandle.h>
#include <nanovdb/math/SampleFromVoxels.h>
#include <nanovdb/tools/GridValidator.h>
#include <stdexcept>

struct NanoVDBMedium::Data {
  nanovdb::GridHandle<> density, temperature;
  const nanovdb::FloatGrid *d = nullptr;
  const nanovdb::FloatGrid *t = nullptr;
};
namespace {
struct GridLocation {
  uint64_t offset, size;
  std::string name;
};
uint64_t checked_offset(uint64_t base, int64_t offset, uint64_t size) {
  if (base > size)
    throw std::runtime_error("Invalid NanoVDB metadata offset.");
  if (offset < 0) {
    uint64_t magnitude = uint64_t(-(offset + 1)) + 1;
    if (magnitude > base)
      throw std::runtime_error("Invalid NanoVDB metadata offset.");
    return base - magnitude;
  }
  if (uint64_t(offset) > size - base)
    throw std::runtime_error("Invalid NanoVDB metadata offset.");
  return base + uint64_t(offset);
}
std::string raw_name(std::ifstream &f, uint64_t start, const nanovdb::GridData &grid) {
  if (!grid.mFlags.isMaskOn(nanovdb::GridFlags::HasLongGridName)) {
    if (!std::memchr(grid.mGridName, 0, sizeof(grid.mGridName)))
      throw std::runtime_error("Invalid NanoVDB grid name.");
    return grid.mGridName;
  }
  uint64_t offset = checked_offset(0, grid.mBlindMetadataOffset, grid.mGridSize);
  if (grid.mBlindMetadataCount > (grid.mGridSize - offset) / sizeof(nanovdb::GridBlindMetaData))
    throw std::runtime_error("Invalid NanoVDB name metadata.");
  for (uint32_t i = 0; i < grid.mBlindMetadataCount; ++i) {
    nanovdb::GridBlindMetaData m;
    uint64_t position = offset + i * sizeof(m);
    f.seekg(start + position);
    if (!f.read(reinterpret_cast<char *>(&m), sizeof(m)))
      throw std::runtime_error("Truncated NanoVDB name metadata.");
    if (m.mDataClass != nanovdb::GridBlindDataClass::GridName)
      continue;
    uint64_t data = checked_offset(position, m.mDataOffset, grid.mGridSize);
    if (m.mValueSize != 1 || m.mValueCount == 0 || m.mValueCount > 1048576 ||
        m.mValueCount > grid.mGridSize - data)
      throw std::runtime_error("Invalid NanoVDB long grid name.");
    std::vector<char> name(m.mValueCount);
    f.seekg(start + data);
    if (!f.read(name.data(), name.size()) || name.back() != 0)
      throw std::runtime_error("Invalid NanoVDB long grid name.");
    return name.data();
  }
  throw std::runtime_error("Missing NanoVDB long grid name.");
}
std::vector<GridLocation> preflight(const std::string &path) {
  std::ifstream f(path, std::ios::binary | std::ios::ate);
  if (!f || f.tellg() < 0)
    throw std::runtime_error("Cannot open NanoVDB file: " + path);
  uint64_t size = uint64_t(f.tellg());
  nanovdb::GridData first{};
  f.seekg(0);
  bool raw = bool(f.read(reinterpret_cast<char *>(&first), sizeof(first))) && first.isValid();
  f.clear();
  f.seekg(0);
  std::vector<GridLocation> grids;
  uint32_t raw_count = 0, raw_remaining = 0;
  while (uint64_t(f.tellg()) < size) {
    if (raw) {
      uint64_t start = uint64_t(f.tellg());
      nanovdb::GridData grid;
      if (!f.read(reinterpret_cast<char *>(&grid), sizeof(grid)) || !grid.isValid() ||
          !grid.mVersion.isCompatible() || grid.mGridCount == 0 ||
          grid.mGridIndex >= grid.mGridCount || grid.mGridSize < sizeof(grid) ||
          grid.mGridSize > size - start)
        throw std::runtime_error("Invalid or truncated raw NanoVDB grid.");
      if (raw_remaining == 0)
        raw_count = raw_remaining = grid.mGridCount;
      if (grid.mGridCount != raw_count || grid.mGridIndex != raw_count - raw_remaining)
        throw std::runtime_error("Inconsistent raw NanoVDB grid indices.");
      --raw_remaining;
      grids.push_back({start, grid.mGridSize, raw_name(f, start, grid)});
      f.seekg(start + grid.mGridSize);
      continue;
    }
    nanovdb::io::FileHeader header;
    if (!f.read(reinterpret_cast<char *>(&header), sizeof(header)))
      throw std::runtime_error("Truncated NanoVDB header.");
    auto magic = nanovdb::toMagic(header.magic);
    if (magic != nanovdb::MagicType::NanoVDB && magic != nanovdb::MagicType::NanoFile)
      throw std::runtime_error("Expected an uncompressed .nvdb file, not OpenVDB.");
    if (header.codec != nanovdb::io::Codec::NONE)
      throw std::runtime_error("Compressed NanoVDB is unsupported; convert to uncompressed .nvdb.");
    if (!header.version.isCompatible() || header.gridCount == 0 ||
        header.gridCount > (size - uint64_t(f.tellg())) / sizeof(nanovdb::io::FileMetaData))
      throw std::runtime_error("Invalid or incompatible NanoVDB header.");
    uint64_t payload = 0;
    size_t first_grid = grids.size();
    for (uint32_t i = 0; i < header.gridCount; ++i) {
      nanovdb::io::FileMetaData m;
      if (!f.read(reinterpret_cast<char *>(&m), sizeof(m)))
        throw std::runtime_error("Truncated NanoVDB metadata.");
      if (m.codec != nanovdb::io::Codec::NONE)
        throw std::runtime_error("Compressed NanoVDB grid is unsupported.");
      if (m.nameSize == 0 || m.nameSize > 1048576 || m.nameSize > size - uint64_t(f.tellg()) ||
          m.fileSize != m.gridSize || m.gridSize < sizeof(nanovdb::GridData) ||
          m.gridSize > size - payload)
        throw std::runtime_error("Invalid NanoVDB grid sizes.");
      std::vector<char> name(m.nameSize);
      if (!f.read(name.data(), m.nameSize) || name.back() != '\0')
        throw std::runtime_error("Invalid NanoVDB grid name.");
      grids.push_back({payload, m.gridSize, name.data()});
      payload += m.fileSize;
    }
    uint64_t start = uint64_t(f.tellg());
    if (payload > size - start)
      throw std::runtime_error("Truncated NanoVDB grid data.");
    for (size_t i = first_grid; i < grids.size(); ++i)
      grids[i].offset += start;
    f.seekg(payload, std::ios::cur);
  }
  if (raw_remaining)
    throw std::runtime_error("Truncated raw NanoVDB grid sequence.");
  if (grids.empty())
    throw std::runtime_error("NanoVDB file contains no grids.");
  return grids;
}
// Validate every pointer before upstream validation dereferences sparse nodes.
// The upstream validator assumes a trustworthy allocation and grid header.
template <class Node>
void check_node_storage(const Node *node, uintptr_t begin, uintptr_t end, size_t &remaining,
                        std::array<uint64_t, 3> &counts) {
  uintptr_t p = reinterpret_cast<uintptr_t>(node);
  if (p < begin || p > end || sizeof(Node) > end - p || p % NANOVDB_DATA_ALIGNMENT || !remaining--)
    throw std::runtime_error("Malformed NanoVDB sparse node offsets.");
  ++counts[Node::LEVEL];
  if constexpr (Node::LEVEL > 0)
    for (auto child = node->cbeginChild(); child; ++child)
      check_node_storage(&*child, begin, end, remaining, counts);
}
void check_storage(const nanovdb::FloatGrid *grid, uint64_t size) {
  using Tree = nanovdb::FloatGrid::TreeType;
  using Root = Tree::RootType;
  uintptr_t begin = reinterpret_cast<uintptr_t>(grid), end = begin + size;
  if (size < sizeof(nanovdb::FloatGrid) + sizeof(Tree) || grid->gridSize() != size)
    throw std::runtime_error("Malformed NanoVDB grid size.");
  const auto *root = static_cast<const Root *>(grid->tree().getRoot());
  uintptr_t p = reinterpret_cast<uintptr_t>(root);
  if (p < begin + sizeof(nanovdb::FloatGrid) + sizeof(Tree) || p > end || sizeof(Root) > end - p ||
      p % NANOVDB_DATA_ALIGNMENT)
    throw std::runtime_error("Malformed NanoVDB root offset.");
  if (root->memUsage() > end - p)
    throw std::runtime_error("Malformed NanoVDB root tile count.");
  const auto &tree = grid->tree();
  const uint64_t node_sizes[3] = {sizeof(Tree::Node0), sizeof(Tree::Node1), sizeof(Tree::Node2)};
  for (int level = 0; level < 3; ++level) {
    int64_t offset = tree.mNodeOffset[level];
    uint64_t count = tree.nodeCount(level);
    if (offset < 0 || uint64_t(offset) > size - sizeof(nanovdb::FloatGrid))
      throw std::runtime_error("Malformed NanoVDB node array offset.");
    uint64_t start = sizeof(nanovdb::FloatGrid) + uint64_t(offset);
    if (count && (start < sizeof(nanovdb::FloatGrid) + sizeof(Tree) || start % 32 ||
                  count > (size - start) / node_sizes[level]))
      throw std::runtime_error("Malformed NanoVDB node array size.");
  }
  size_t remaining = size / 32;
  std::array<uint64_t, 3> counts{};
  for (auto child = root->cbeginChild(); child; ++child)
    check_node_storage(&*child, begin, end, remaining, counts);
  for (int level = 0; level < 3; ++level)
    if (counts[level] != tree.nodeCount(level))
      throw std::runtime_error("NanoVDB sparse node counts do not match its tree.");
}
nanovdb::GridHandle<> read_grid(const std::string &filename, const std::string &name,
                                const std::vector<GridLocation> &locations) {
  auto found = std::find_if(locations.begin(), locations.end(),
                            [&](const GridLocation &g) { return g.name == name; });
  if (found == locations.end())
    throw std::runtime_error("NanoVDB grid '" + name + "' is missing.");
  auto buffer = nanovdb::HostBuffer::create(found->size);
  std::ifstream input(filename, std::ios::binary);
  input.seekg(found->offset);
  if (!input.read(reinterpret_cast<char *>(buffer.data()), found->size))
    throw std::runtime_error("Truncated NanoVDB grid: " + name);
  auto *header = reinterpret_cast<nanovdb::GridData *>(buffer.data());
  if (!header->isValid() || !header->mVersion.isCompatible() || header->mGridCount == 0 ||
      header->mGridIndex >= header->mGridCount || header->mGridSize != found->size)
    throw std::runtime_error("Invalid NanoVDB grid header: " + name);
  if (header->mGridType != nanovdb::GridType::Float)
    throw std::runtime_error("NanoVDB grid '" + name + "' is not a float grid.");
  const auto *grid = reinterpret_cast<const nanovdb::FloatGrid *>(buffer.data());
  check_storage(grid, found->size);
  if (grid->gridClass() != nanovdb::GridClass::FogVolume &&
      grid->gridClass() != nanovdb::GridClass::Unknown)
    throw std::runtime_error("NanoVDB grid '" + name +
                             "' must be a fog or scalar field, not a level set.");
  // The upstream full traversal assumes nonempty leaf/lower node arrays and
  // rejects valid tile-only trees. check_storage validates all sparse pointers
  // against the allocation; retain upstream header and full checksum checks.
  if (!nanovdb::tools::isValid<float>(grid, nanovdb::CheckMode::Half, false) ||
      !nanovdb::tools::validateChecksum(grid, nanovdb::CheckMode::Full))
    throw std::runtime_error("Invalid NanoVDB grid: " + name);
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
  // Normalize a selected grid only after validating its original allocation and
  // checksum. GridHandle's constructor can then safely inspect this one grid.
  nanovdb::tools::updateGridCount(header, 0, 1);
  return nanovdb::GridHandle<>(std::move(buffer));
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
  auto locations = preflight(filename);
  data->density = read_grid(filename, Rcpp::as<std::string>(d["density_grid"]), locations);
  data->d = data->density.grid<float>();
  if (data->d->tree().background() != 0)
    throw std::runtime_error("NanoVDB density background must be zero.");
  if (d.containsElementNamed("temperature_grid") && !Rf_isNull(d["temperature_grid"])) {
    data->temperature =
        read_grid(filename, Rcpp::as<std::string>(d["temperature_grid"]), locations);
    data->t = data->temperature.grid<float>();
    auto validate = [&](const nanovdb::Coord &, int, float v) {
      if (!std::isfinite(v) || v < 0 ||
          !std::isfinite(Float((double(v) - temperature_offset) * temperature_scale)))
        throw std::runtime_error("NanoVDB temperatures and their scaled values must be finite.");
    };
    validate(nanovdb::Coord(0), 1, data->t->tree().background());
    visit_values(data->t->tree().root(), validate);
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
  return sizeof(*this) + sizeof(Data) + data->density.bufferSize() +
         data->temperature.bufferSize() + sizeof(Float) * majorants.density.capacity();
}
