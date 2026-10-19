#if defined(NOT_CRAN) && defined(RAY_HAS_METAL_BVH)
#include "wavefront_deformation.h"
#include "wavefront_motion.h"
#include "camera.h"
#include "../hitables/mesh3d.h"
#include "../materials/texturecache.h"
#include <testthat.h>
#include <cstring>

context("Resident Metal deformation") {
  test_that("bounds refits preserve membership flags and replace old extents") {
    WavefrontScene scene;
    scene.boundaries.resize(2);
    scene.boundary_lo.w = 1;
    scene.boundary_hi.w = 7;
    WFTriangle face;
    face.boundary = 1;
    face.p[0] = {0,0,0,0}; face.p[1] = {1,0,0,0}; face.p[2] = {0,1,0,0};
    scene.triangles = {face};
    WavefrontMotion motion;
    motion.Claim({}, WavefrontMotion::sizes(scene), 0);
    motion.Seal(scene);
    for (float z : {5.f, -3.f, 0.f, 5.f}) {
      for (auto &p : scene.triangles[0].p) p.z = z;
      motion.RefitBounds(scene, {{0,1}});
      expect_true(scene.boundary_lo.w == 1);
      expect_true(scene.boundary_hi.w == 7);
      expect_true(scene.boundaries[1].lo.z == z);
      expect_true(scene.boundaries[1].hi.z == z);
      expect_true(scene.boundary_lo.z == z);
      expect_true(scene.boundary_hi.z == z);
    }
  }
  test_that("refits preserve sessions and face slots through random seeks and collapsed triangles") {
    if (!MetalWavefrontAvailable()) return;
    using namespace animation;
    for (bool initially_collapsed : {false, true}) for (bool placed : {false, true}) {
      auto source = std::make_shared<DeformationMesh>();
      source->positions = {point3f(0,0,0), point3f(1,0,0), point3f(0,1,0)};
      source->normals = {normal3f(0,0,1),normal3f(0,0,1),normal3f(0,0,1)};
      source->indices = {0,1,2};
      MorphTarget target;
      // Move far beyond the first frame's boxes. The collapsed case starts as
      // a zero-area face and must later become hittable without a new BLAS.
      target.positions = {vec3f(0,0,3), vec3f(0,0,3), vec3f(0,0,3)};
      if (initially_collapsed) {
        source->positions[1] = source->positions[2] = point3f(0,0,0);
        target.positions[1][0] = 1; target.positions[2][1] = 1;
      }
      source->targets = {target};
      auto asset = std::make_shared<ModelAsset>(); asset->id = "resident-mesh";
      ModelNode node; node.id = "surface"; node.meshes = {source}; node.weights = {0};
      asset->nodes = {node};
      ModelClip clip; clip.id = "move"; clip.begin = 0; clip.end = 2; clip.nodes.resize(1);
      clip.weights = {{std::make_shared<Track>(Property::Scalar, Interpolation::Linear,
          std::vector<double>{0,2}, std::vector<Value>{{0,0,0,0},{2,0,0,0}})}};
      asset->clips = {clip}; asset->Validate();
      ModelPlayer player; player.clip = "move";
      auto model = std::make_shared<ModelInstance>(asset, player);
      Rcpp::NumericMatrix vertices(3,3), normals(3,3);
      Rcpp::IntegerMatrix indices(1,3);
      for (int i = 0; i < 3; ++i) {
        indices(0,i) = i;
        for (int j = 0; j < 3; ++j) {
          vertices(i,j) = source->positions[i][j]; normals(i,j) = source->normals[i][j];
        }
      }
      auto info = Rcpp::List::create(Rcpp::_["vertices"] = vertices, Rcpp::_["indices"] = indices,
          Rcpp::_["normals"] = normals, Rcpp::_["texcoords"] = Rcpp::NumericMatrix(0,2),
          Rcpp::_["color_vals"] = Rcpp::NumericMatrix(0,3), Rcpp::_["subdivision_levels"] = 0,
          Rcpp::_["texture"] = "", Rcpp::_["bump_texture"] = "", Rcpp::_["bump_intensity"] = 0);
      Transform identity;
      TextureCache textures; random_gen rng;
      auto material = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
      auto mesh = std::make_shared<mesh3d>(info, material, "", 0, false, textures, false, false,
          0, 0, BVHBuildOptions{}, rng, &identity, &identity, false);
      mesh->BindDeformation(model, 0, 0);
      auto clock = std::make_shared<SceneTimeContext>();
      clock->deformation_updates.push_back([&](double time) { return mesh->PrepareDeformation(time); });
      clock->deformation_bounds.push_back([&](Float a, Float b) { mesh->mesh_bvh->Refit(a,b); });
      camera cam({0,0,10}, {0,0,0}, {0,1,0}, 30, 1, 0, 10, 0, 0, 1);
      cam.scene_camera_time = clock;
      hitable_list world, lights;
      std::shared_ptr<hitable> primitive = mesh;
      if (placed) {
        auto track = std::make_shared<TransformTrack>();
        track->position = std::make_shared<Track>(Property::Position, Interpolation::Linear,
            std::vector<double>{0,2}, std::vector<Value>{{0,0,0,0},{0,0,0,0}});
        AnimatedTransform transform(&identity, 0, &identity, 0, track, 0, clock);
        primitive = std::make_shared<AnimatedHitable>(primitive, transform);
      }
      world.add(primitive);
      WavefrontContext context;
      WavefrontReport report;
      clock->UpdateGeometry(0,0,0,0);
      context.BeginFrame(world, lights, cam, 1, 0, report, 1);
      expect_true(bool(context.session));
      if (!context.session) continue;
      auto *session = context.session.get();
      auto *storage = context.scene.triangles.data();
      expect_true(context.scene.triangles.size() == 1);
      for (double time : {1.5,.25,2.,0.,1.,1.}) {
        const auto revision = clock->geometry_revision;
        clock->UpdateGeometry(time,0,0,0);
        context.BeginFrame(world, lights, cam, 1, 0, report, 1);
        expect_true(context.session.get() == session);
        expect_true(context.scene.triangles.data() == storage);
        expect_true(report.scene_reused);
        expect_true(report.acceleration_builds == 0);
        expect_true(report.acceleration_refits ==
            (revision == clock->geometry_revision ? 0u : (placed ? 2u : 1u)));
        expect_true(report.geometry_bytes_updated ==
            (revision == clock->geometry_revision ? 0u : sizeof(WFTriangle)));
        WFGeometryProbe probe;
        probe.origin = {.05f,.05f,10,.0001f}; probe.direction = {0,0,-1,20};
        auto resident = ProbeMetalSessionGeometry(*session, {probe});
        auto fresh = ProbeMetalGeometry(context.scene, {probe});
        const bool hit = !initially_collapsed || time > 0;
        expect_true(bool(resident[0].normal.w) == hit);
        expect_true(resident[0].normal.w == fresh[0].normal.w);
        if (hit) {
          expect_true(std::abs(resident[0].point.w - (10 - 3 * time)) < 1e-5);
          expect_true(std::abs(resident[0].point.w - fresh[0].point.w) < 1e-5);
          expect_true(resident[0].primitive == fresh[0].primitive);
          expect_true(std::abs(resident[0].normal.z - 1) < 1e-5);
        }
      }
    }
  }
}
#endif
