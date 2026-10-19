#ifdef NOT_CRAN
#include "model_deformation.h"
#include "animatedtransform.h"
#include "../hitables/mesh3d.h"
#include "../materials/texturecache.h"
#include "../volumes/picking.h"
#include <testthat.h>

namespace {
using namespace animation;
DeformationMesh triangle_source() {
  DeformationMesh mesh;
  mesh.positions = {point3f(0,0,0), point3f(1,0,0), point3f(0,1,0)};
  mesh.normals = {normal3f(0,0,1), normal3f(0,0,1), normal3f(0,0,1)};
  mesh.tangents = {vec3f(1,0,0), vec3f(1,0,0), vec3f(1,0,0)};
  mesh.handedness = {1,1,1}; mesh.indices = {0,1,2};
  return mesh;
}
bool near(const point3f &a, const point3f &b) { return (a-b).length() < 1e-6; }
}
context("Exact model deformation reference") {
  test_that("frame publication validates every deformation before committing any") {
    SceneTimeContext clock;
    double first=0, second=0; int refits=0, notices=0;
    clock.has_model_poses = true;
    clock.frozen_pose_notice = [&] { ++notices; };
    clock.deformation_updates.push_back([&](double time) -> std::function<void()> {
      if (time==first) return {};
      return [&,time] { first=time; };
    });
    clock.deformation_updates.push_back([&](double time) -> std::function<void()> {
      if (time>2) throw std::runtime_error("invalid pose");
      if (2*time==second) return {};
      return [&,time] { second=2*time; };
    });
    clock.deformation_bounds.push_back([&](Float open, Float close) {
      expect_true(open<=close); ++refits;
    });
    clock.UpdateGeometry(1,1,1,1);
    expect_true((first==1 && second==2 && clock.epoch==1 && clock.geometry_revision==1 && refits==1));
    expect_error(clock.UpdateGeometry(3,3,3,3));
    expect_true((first==1 && second==2 && clock.epoch==1 && clock.geometry_revision==1 && refits==1));
    // An asymmetric exposure samples the frame pose, not either endpoint.
    clock.UpdateGeometry(2,2,2.25,2.5);
    expect_true((first==2 && second==4 && clock.epoch==2 && clock.pose_offset==0 && notices==1));
    clock.UpdateGeometry(2,2,1.5,2.5);
    expect_true((clock.geometry_revision==2 && notices==1));
    clock.UpdateGeometry(1,1,1,1);
    expect_true((clock.geometry_revision==3 && refits==4 && notices==1));
    clock.UpdateGeometry(1,1,1.25,1.25);
    expect_true((first==1.25 && second==2.5 && clock.pose_offset==.25 && notices==1));
    expect_error(clock.UpdateGeometry(2,2,2.5,2.5));
    expect_true((first==1.25 && clock.epoch==1 && clock.pose_offset==.25));
  }
  test_that("display picking and light PDFs follow committed deformation and refitted bounds") {
    auto source = std::make_shared<DeformationMesh>(triangle_source());
    MorphTarget target; target.positions = {vec3f(0,0,1),vec3f(0,0,1),vec3f(0,0,1)};
    source->targets = {target};
    auto asset = std::make_shared<ModelAsset>(); asset->id="pick-asset";
    ModelNode node; node.id="surface"; node.meshes={source}; node.weights={0};
    asset->nodes={node};
    ModelClip clip; clip.id="move"; clip.begin=0; clip.end=2; clip.nodes.resize(1);
    clip.weights={ {std::make_shared<Track>(Property::Scalar, Interpolation::Linear,
      std::vector<double>{0,2},std::vector<Value>{{0,0,0,0},{2,0,0,0}})} };
    asset->clips={clip}; asset->Validate();
    ModelPlayer player; player.clip="move";
    auto instance=std::make_shared<ModelInstance>(asset,player);
    Rcpp::NumericMatrix vertices(3,3), normals(3,3);
    Rcpp::IntegerMatrix indices(1,3);
    for (int i=0;i<3;++i) {
      indices(0,i)=i;
      for (int j=0;j<3;++j) {
        vertices(i,j)=source->positions[i][j]; normals(i,j)=source->normals[i][j];
      }
    }
    auto info=Rcpp::List::create(Rcpp::_["vertices"]=vertices, Rcpp::_["indices"]=indices,
      Rcpp::_["normals"]=normals, Rcpp::_["texcoords"]=Rcpp::NumericMatrix(0,2),
      Rcpp::_["color_vals"]=Rcpp::NumericMatrix(0,3), Rcpp::_["subdivision_levels"]=0,
      Rcpp::_["texture"]="", Rcpp::_["bump_texture"]="", Rcpp::_["bump_intensity"]=0);
    Transform identity; TextureCache textures; random_gen rng;
    auto mat=std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(1)));
    auto surface=std::make_shared<mesh3d>(info,mat,"",0,false,textures,false,false,
      0,0,BVHBuildOptions{},rng,&identity,&identity,false);
    surface->BindDeformation(instance,0,0);
    BVHAggregate world({surface},0,0,1,true);
    SceneTimeContext clock;
    clock.deformation_updates.push_back([&](double time) { return surface->PrepareDeformation(time); });
    clock.deformation_bounds.push_back([&](Float a, Float b) { surface->mesh_bvh->Refit(a,b); world.Refit(a,b); });
    for (double time : {0.,2.,.5,1.75,0.,2.}) {
      clock.UpdateGeometry(time,0,0,0);
      const Ray ray(point3f(.2,.2,5),vec3f(0,0,-1),Float(0));
      auto pick=PickRay(ray,&world,nullptr);
      expect_true(bool(pick));
      if (pick) expect_true(std::abs(pick->p[2]-time)<1e-6);
      hit_record hit;
      expect_true(world.hit(ray,.001,100,hit,rng));
      expect_true(std::abs(hit.t-(5-time))<1e-6);
      expect_true(std::abs(surface->pdf_value(ray.origin(),ray.direction(),rng,0)-
        2*(5-time)*(5-time))<1e-4);
    }
  }
  test_that("two joints and all influences reproduce bind and animated positions") {
    auto mesh = triangle_source();
    mesh.influences = {{{0,.5},{1,.5}}, {{0,1}}, {{1,.1},{1,.2},{1,.3},{1,.4},{0,0}}};
    mesh.Validate(2);
    auto rest = DeformMesh(mesh, {}, {Transform(),Transform()});
    for (size_t i=0;i<3;++i) expect_true(near(rest.positions[i],mesh.positions[i]));
    auto pose = DeformMesh(mesh, {}, {Translate(vec3f(2,0,0)),Translate(vec3f(0,2,0))});
    expect_true(near(pose.positions[0],point3f(1,1,0)));
    expect_true(near(pose.positions[1],point3f(3,0,0)));
    expect_true(near(pose.positions[2],point3f(0,3,0)));
    expect_true(near(mesh.positions[0],point3f(0,0,0)));
  }
  test_that("morphs precede skinning and arbitrary reevaluation never accumulates") {
    auto mesh = triangle_source();
    MorphTarget target; target.positions = {vec3f(1,0,0),vec3f(0,0,1),vec3f(0,0,2)};
    mesh.targets = {target}; mesh.Validate(0);
    auto only_morph = DeformMesh(mesh,{.5},{});
    expect_true(near(only_morph.positions[0],point3f(.5,0,0)));
    mesh.influences = {{{0,1}},{{0,1}},{{0,1}}}; mesh.Validate(1);
    auto combined = DeformMesh(mesh,{.5},{Scale(2,3,4)});
    expect_true(near(combined.positions[0],point3f(1,0,0)));
    expect_true(near(combined.positions[1],point3f(2,0,2)));
    for (double weight : {2.,0.,-1.,.2,.5}) DeformMesh(mesh,{weight},{RotateY(40)});
    auto again = DeformMesh(mesh,{.5},{Scale(2,3,4)});
    for (size_t i=0;i<3;++i) expect_true(near(again.positions[i],combined.positions[i]));
  }
  test_that("nonuniform and reflected skin maps retain orthonormal shading frames") {
    auto mesh = triangle_source();
    mesh.influences = {{{0,1}},{{0,1}},{{0,1}}};
    for (auto &normal : mesh.normals) normal = normal3f(0, std::sqrt(.5), std::sqrt(.5));
    auto pose = DeformMesh(mesh,{}, {Scale(-2,3,4)});
    for (size_t i=0;i<3;++i) {
      expect_true(std::abs(pose.normals[i][1]-.8)<1e-6);
      expect_true(std::abs(pose.normals[i][2]-.6)<1e-6);
      expect_true(std::abs(pose.tangents[i].length()-1)<1e-6);
      expect_true(std::abs(pose.tangents[i][1])<1e-6);
      expect_true(pose.handedness[i]==-1);
    }
  }
  test_that("posed surface normals include skin weight gradients and preserve seeking") {
    auto mesh = triangle_source();
    mesh.influences = {{{0,1}},{{1,1}},{{0,1}}};
    const std::vector<Transform> palette = {Transform(),Translate(vec3f(0,0,1))};
    auto authored = DeformMesh(mesh,{},palette);
    expect_true(authored.normals[0][2]==1);
    mesh.recalculate_normals = true;
    mesh.Validate(2);
    auto posed = DeformMesh(mesh,{},palette);
    for (int i=0;i<3;++i) {
      expect_true(std::abs(posed.normals[i][0]+std::sqrt(.5))<1e-6);
      expect_true(std::abs(posed.normals[i][2]-std::sqrt(.5))<1e-6);
      expect_true(std::abs(dot(posed.tangents[i],vec3f(posed.normals[i][0],
        posed.normals[i][1],posed.normals[i][2])))<1e-6);
    }
    DeformMesh(mesh,{}, {Transform(),Translate(vec3f(0,0,-2))});
    auto again = DeformMesh(mesh,{},palette);
    for (int i=0;i<3;++i) expect_true((again.normals[i]-posed.normals[i]).length()<1e-7);
    mesh.normals.clear(); mesh.tangents.clear(); mesh.handedness.clear();
    mesh.Validate(2);
    expect_true(DeformMesh(mesh,{},palette).normals.size()==3);
    mesh.positions[2]=mesh.positions[0];
    expect_error(DeformMesh(mesh,{},palette));
  }
  test_that("recomputed normals use corner angles and keep split vertices separate") {
    DeformationMesh mesh;
    mesh.recalculate_normals = true;
    mesh.positions={point3f(0,0,0),point3f(1,0,0),point3f(0,1,0),point3f(0,0,4),
      point3f(0,0,0),point3f(1,0,0),point3f(0,-1,0)};
    mesh.indices={0,1,2,0,3,1,4,5,6};
    mesh.Validate(0);
    const auto pose=DeformMesh(mesh,{},{});
    // Equal right angles at vertex 0, despite different incident face areas.
    expect_true(std::abs(pose.normals[0][1]-std::sqrt(.5))<1e-6);
    expect_true(std::abs(pose.normals[0][2]-std::sqrt(.5))<1e-6);
    expect_true(pose.normals[4][2]==-1);
  }
  test_that("malformed weights and singular blends fail explicitly") {
    auto mesh=triangle_source(); mesh.influences={{{0,0}},{{0,1}},{{0,1}}};
    expect_error(mesh.Validate(1));
    mesh.influences[0]={{2,1}}; expect_error(mesh.Validate(1));
    mesh.influences[0]={{0,-1}}; expect_error(mesh.Validate(1));
    mesh.influences={{{0,.5},{1,.5}},{{0,1}},{{0,1}}}; mesh.Validate(2);
    expect_error(DeformMesh(mesh,{}, {Scale(1,1,1),Scale(-1,1,1)}));
    expect_error(DeformMesh(mesh,{1},{}));
  }
}
#endif
