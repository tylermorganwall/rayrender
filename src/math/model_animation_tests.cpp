#ifdef NOT_CRAN
#include "model_animation.h"
#include <testthat.h>

namespace {
using namespace animation;
std::shared_ptr<const Track> positions(std::vector<double> times, std::vector<Value> values) {
  return std::make_shared<const Track>(Property::Position, Interpolation::Linear, times, values);
}
std::shared_ptr<ModelAsset> bird_asset() {
  auto asset = std::make_shared<ModelAsset>();
  asset->id = "bird";
  ModelNode body; body.id = "body";
  body.rest.position = positions({0}, {{0, 2, 0, 0}});
  ModelNode wing; wing.id = "wing"; wing.parent = 0;
  wing.rest.position = positions({0}, {{1, 0, 0, 0}});
  asset->nodes = {body, wing};
  ModelClip clip; clip.id = "fly"; clip.name = "Flight"; clip.begin = 2; clip.end = 4;
  clip.nodes = {body.rest, wing.rest};
  clip.nodes[0].position = positions({2, 3, 4}, {{0,2,0,0}, {2,3,1,0}, {4,2,2,0}});
  clip.nodes[1].rotation = std::make_shared<const Track>(Property::Rotation, Interpolation::Smooth,
    std::vector<double>{2, 2.5, 4}, std::vector<Value>{{0,0,-.6,.8}, {0,0,.6,.8}, {0,0,-.6,.8}});
  asset->clips = {clip}; asset->Validate();
  return asset;
}
}

context("Shared model hierarchy") {
  test_that("instances share immutable channels but never share evaluated poses") {
    auto asset = bird_asset();
    ModelPlayer a; a.clip = "fly"; a.loop = Loop::Repeat;
    ModelPlayer b = a; b.offset = .5; b.speed = -2;
    ModelInstance first(asset, a), second(asset, b);
    expect_true(first.asset == second.asset);
    expect_true(first.nodes[1]->rotation == second.nodes[1]->rotation);
    expect_true(first.nodes[0]->position == asset->clips[0].nodes[0].position);
    expect_true(first.nodes[1] != second.nodes[1]);
    expect_true(first.clock != second.clock);
    auto original = first.At(.125);
    for (double time : {8., -1., 0., 2., .5, .125, 1e9, .125}) {
      auto independent = second.At(time);
      auto repeated = first.At(.125);
      expect_true(repeated.local == original.local);
      expect_true(repeated.global == original.global);
      expect_true(independent.local.data() != repeated.local.data());
    }
    auto pose = first.At(1);
    auto pivot = pose.global[1](point3f(0,0,0));
    for (int k = 0; k < 3; ++k) expect_true(std::abs(pivot[k] - (k == 2 ? 1 : 3)) < 1e-6);
    expect_true(pose.global[1] == pose.global[0]*pose.local[1]);
    auto placed = first.At(1, Translate(vec3f(5,0,0)));
    auto dirty = ComparePoses(pose, placed);
    expect_true((dirty.placement && !dirty.pose));
    expect_true(placed.global == pose.global);
    expect_true(!ComparePoses(first.At(0), first.At(2)).pose);
  }

  test_that("root removal retains bob and child articulation with conservative bounds") {
    auto asset = bird_asset();
    asset->clips[0].nodes[0].rotation = std::make_shared<const Track>(Property::Rotation,
      Interpolation::Linear, std::vector<double>{2,4}, std::vector<Value>{{0,0,0,1},{0,.6,0,.8}});
    ModelPlayer player; player.clip = "fly"; player.loop = Loop::Repeat;
    player.in_place = true; player.root = "body"; player.remove_translation = {true,false,true};
    player.remove_yaw = true;
    ModelInstance in_place(asset, player);
    auto pose = in_place.At(1);
    auto position = pose.local[0](point3f(0,0,0));
    for (int k = 0; k < 3; ++k) expect_true(std::abs(position[k] - (k == 1 ? 3 : 0)) < 1e-6);
    auto right = pose.local[0](vec3f(1,0,0));
    expect_true((std::abs(right[0]-1) < 1e-6 && std::abs(right[2]) < 1e-6));
    expect_true(!(pose.local[1] == asset->nodes[1].rest.At(0)));
    aabb box(point3f(-1,-.2,-.2), point3f(1,.2,.2));
    auto wing_bound = in_place.nodes[1]->Bounds(box,-5,7);
    auto bound = in_place.nodes[0]->Bounds(wing_bound,-5,7);
    for (int sample = 0; sample <= 300; ++sample) {
      auto current = in_place.At(-5+12.*sample/300);
      for (int corner = 0; corner < 8; ++corner) {
        point3f p((corner&1) ? 1 : -1, (corner&2) ? .2 : -.2, (corner&4) ? .2 : -.2);
        p = current.global[1](p);
        for (int k = 0; k < 3; ++k) expect_true((p[k] >= bound.min()[k] && p[k] <= bound.max()[k]));
      }
    }
    player.root = "wing";
    expect_error(ModelInstance(asset, player));
    player.root = "missing";
    expect_error(ModelInstance(asset, player));
  }

  test_that("yaw preserves swing and rejects singular trajectories before publication") {
    Track singular(Property::Rotation, Interpolation::Linear, {0,1}, {{1,0,0,0},{0,0,0,1}});
    expect_error(singular.ValidateYaw());
    Track crossing(Property::Rotation, Interpolation::Cubic, {0,1}, {{0,0,0,1},{0,0,0,1}},
      {{0,0,0,0},{-4,0,0,4}}, {{4,0,0,-4},{0,0,0,0}});
    // The full quaternion is valid, but its yaw projection vanishes at .5.
    expect_error(crossing.ValidateYaw());
    Track huge_swing(Property::Rotation, Interpolation::Cubic, {0,1}, {{0,0,0,1},{0,0,0,1}},
      {{0,0,0,0},{-1e10,0,0,0}}, {{1e10,0,0,0},{0,0,0,0}});
    expect_error(huge_swing.ValidateYaw());
    TransformTrack combined;
    double sx = std::sin(M_PI/12), cx = std::cos(M_PI/12);
    double sy = .5, cy = std::sqrt(.75);
    combined.rotation = std::make_shared<const Track>(Property::Rotation, Interpolation::Linear,
      std::vector<double>{0}, std::vector<Value>{{cy*sx,sy*cx,-sy*sx,cy*cx}});
    auto removal = std::make_shared<RootMotion>(); removal->yaw = true;
    combined.root_motion = removal;
    auto expected = RotateX(30).GetMatrix(), actual = combined.At(0).GetMatrix();
    for (int r = 0; r < 4; ++r) for (int c = 0; c < 4; ++c)
      expect_true(std::abs(expected.m[r][c]-actual.m[r][c]) < 1e-6);
  }

  test_that("invalid graphs and ambiguous labels never become target identities") {
    auto asset = bird_asset();
    asset->nodes[0].parent = 1;
    expect_error(asset->Validate());
    asset->nodes[0].parent = -1;
    asset->nodes[1].id = "body";
    expect_error(asset->Validate());
    asset = bird_asset();
    auto other = asset->clips[0]; other.id = "other";
    asset->clips.push_back(other); asset->Validate();
    expect_error(asset->ClipIndex("Flight"));
    expect_true(asset->ClipIndex("other") == 1);
    expect_error(asset->NodeIndex("Flight"));
  }
}
#endif
