// Regenerates native/ZIP/BLOSC fixtures with upstream APIs, independently of
// openvdbr's reader. Compile against the provider's bundled sources and archives.
#define NANOVDB_USE_OPENVDB
#define NANOVDB_USE_TBB
#define NANOVDB_USE_ZIP
#define NANOVDB_USE_BLOSC
#include <openvdb/openvdb.h>
#include <openvdb/io/File.h>
#include <nanovdb/tools/CreateNanoGrid.h>
#include <nanovdb/io/IO.h>
#include <string>
int main(int argc, char** argv) {
  if (argc != 2) return 1;
  const std::string directory = argv[1];
  openvdb::initialize();
  auto density = openvdb::FloatGrid::create(0);
  density->setName("density");
  density->setGridClass(openvdb::GRID_FOG_VOLUME);
  auto transform = openvdb::math::Transform::createLinearTransform(.12);
  transform->postRotate(.35, openvdb::math::Z_AXIS);
  transform->postTranslate(openvdb::Vec3d(.08, -.05, 0));
  density->setTransform(transform);
  auto temperature = openvdb::FloatGrid::create(0);
  temperature->setName("temperature");
  temperature->setTransform(openvdb::math::Transform::createLinearTransform(.24));
  for (int z = -3; z <= 3; ++z) for (int y = -3; y <= 3; ++y) for (int x = -3; x <= 3; ++x) {
    const double radius = std::sqrt(double(x*x + y*y + z*z));
    if (radius < 3.1) density->tree().setValue(openvdb::Coord(x, y, z), float((4-radius)*.2));
    temperature->tree().setValue(openvdb::Coord(x, y, z), float(1800 + 100*(y+3)));
  }
  density->tree().setValueOff(openvdb::Coord(4, 0, 0), .12f);
  openvdb::io::File file(directory + "/transformed.vdb");
  file.write(openvdb::GridPtrVec{density, temperature});
  std::vector<nanovdb::GridHandle<>> handles;
  handles.push_back(nanovdb::tools::createNanoGrid(*density));
  handles.push_back(nanovdb::tools::createNanoGrid(*temperature));
  nanovdb::io::writeGrids(directory + "/transformed-zip.nvdb", handles, nanovdb::io::Codec::ZIP);
  nanovdb::io::writeGrids(directory + "/transformed-blosc.nvdb", handles, nanovdb::io::Codec::BLOSC);
  return 0;
}
