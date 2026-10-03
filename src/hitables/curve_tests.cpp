#ifdef NOT_CRAN
#include "curve.h"
#include <testthat.h>

context("Curve intersection and subdivision") {
  test_that("recursive intersection retains closest distance and the caller's interval") {
    Transform identity;
    auto mat = std::make_shared<hair>(point3f(.5f),1.55f,.3f,.3f,2);
    const point3f p[4] = {point3f(-1,0,4), point3f(2,0,4), point3f(-2,0,1), point3f(1,0,1)};
    auto common = std::make_shared<CurveCommon>(p,.025f,.025f,CurveType::Cylinder,nullptr);
    curve whole(0,1,common,mat,&identity,&identity,false);
    random_gen rng(99); RandomSampler sampler(rng);
    hit_record hit, other;
    Ray ray(point3f(0),vec3f(0,0,1),nullptr,.91f);
    expect_true(whole.hit(ray,0,10,hit,rng));
    expect_true((hit.t > 1 && hit.t < 2));
    expect_true(whole.hit(ray,2,10,hit,rng));
    expect_true((hit.t >= 2 && hit.t < 3));
    expect_true(whole.hit(ray,3,10,hit,rng));
    expect_true(hit.t >= 3);
    expect_false(whole.hit(ray,0,.9f,hit,rng));
    expect_false(whole.hit(ray,4.1f,10,hit,rng));
    expect_true(whole.hit(ray,0,10,hit,rng));
    expect_true(whole.hit(ray,0,10,other,&sampler));
    expect_true((hit.t == other.t && hit.u == other.u));
  }

  test_that("straight tapered split curves agree with the full primitive under transforms") {
    auto mat = std::make_shared<hair>(point3f(.5f),1.55f,.3f,.3f,2);
    const point3f p[4] = {point3f(0,-1,0),point3f(0,-1.f/3,0),point3f(0,1.f/3,0),point3f(0,1,0)};
    vec3f normals[2] = {vec3f(0,0,1),vec3f(0,0,1)};
    Transform object = Translate(vec3f(3,-2,5)) * RotateY(37) * Scale(2,1,.7f);
    Transform inverse = Inverse(object);
    random_gen rng(614);
    for (auto type : {CurveType::Flat,CurveType::Cylinder,CurveType::Ribbon}) {
      auto common = std::make_shared<CurveCommon>(p,.1f,.2f,type,normals);
      curve whole(0,1,common,mat,&object,&inverse,false);
      std::vector<std::shared_ptr<curve>> segments;
      for (int i=0;i<8;++i)
        segments.push_back(std::make_shared<curve>(i/8.f,(i+1)/8.f,common,mat,&object,&inverse,false));
      for (int j=0;j<512;++j) {
        Float y=-1.1f+2.2f*rng.unif_rand(), x=-.12f+.24f*rng.unif_rand();
        Ray local(point3f(x,y,-2),vec3f(0,0,2));
        Ray ray(object(local.origin()),object(local.direction()));
        hit_record a,b;
        bool full=whole.hit(ray,0,4,a,rng), split=false;
        Float closest=4;
        for (auto segment:segments) if(segment->hit(ray,0,closest,b,rng)) {
          split=true;closest=b.t;
        }
        expect_true(full==split);
        if(full && split) {
          expect_true(std::abs(a.t-closest)<2e-5);
          expect_true((std::isfinite(a.normal[0]) && std::isfinite(a.normal[1]) && std::isfinite(a.normal[2])));
        }
      }
    }
  }

  test_that("curved primitives agree with a finer independent subdivision away from silhouettes") {
    Transform identity;
    auto mat = std::make_shared<hair>(point3f(.5f),1.55f,.3f,.3f,2);
    const point3f p[4] = {point3f(-1,-1,0),point3f(1,-.5f,.4f),
                         point3f(-1,.5f,-.4f),point3f(1,1,0)};
    auto common = std::make_shared<CurveCommon>(p,.08f,.12f,CurveType::Cylinder,nullptr);
    curve whole(0,1,common,mat,&identity,&identity,false);
    std::vector<std::shared_ptr<curve>> fine;
    for(int i=0;i<64;++i)
      fine.push_back(std::make_shared<curve>(i/64.f,(i+1)/64.f,common,mat,&identity,&identity,false));
    random_gen rng(511);
    int common_hits=0, mismatches=0;
    Float largest_distance_error=0;
    for(int i=0;i<4096;++i) {
      Ray ray(point3f(-1.1f+2.2f*rng.unif_rand(),-1.1f+2.2f*rng.unif_rand(),-3),vec3f(0,0,1));
      hit_record a,b;
      bool coarse_hit=whole.hit(ray,0,6,a,rng), fine_hit=false;
      Float closest=6;
      for(auto part:fine) if(part->hit(ray,0,closest,b,rng)) { fine_hit=true;closest=b.t; }
      if(coarse_hit != fine_hit) ++mismatches;
      if(coarse_hit && fine_hit) {
        ++common_hits;
        largest_distance_error=std::max(largest_distance_error,std::abs(a.t-closest));
      }
    }
    // Ray-facing curve intersection is approximate to width/20. Subdivision
    // may move silhouette classifications by that tolerance; it is not a bitwise
    // mesh tessellation equivalence. This fixture has no overlapping segments.
    expect_true(common_hits > 100);
    expect_true(mismatches < 12);
    expect_true(largest_distance_error < .012f);
  }

  test_that("zero tangent and zero-width tips do not create invalid intersections") {
    Transform identity;
    auto mat = std::make_shared<hair>(point3f(0),1.55f,.3f,.3f,2);
    point3f p[4] = {point3f(0),point3f(0),point3f(0),point3f(0)};
    auto common=std::make_shared<CurveCommon>(p,.1f,.1f,CurveType::Cylinder,nullptr);
    curve zero(0,1,common,mat,&identity,&identity,false);
    hit_record hit; random_gen rng(1);
    expect_false(zero.hit(Ray(point3f(0,0,-1),vec3f(0,0,1)),0,2,hit,rng));
    expect_false(zero.hit(Ray(point3f(0),vec3f(0)),0,2,hit,rng));
  }
}
#endif
