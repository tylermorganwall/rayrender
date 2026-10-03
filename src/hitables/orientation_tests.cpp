#ifdef NOT_CRAN
#include "sphere.h"
#include "ellipsoid.h"
#include "disk.h"
#include "cylinder.h"
#include "instance.h"
#include <testthat.h>

context("Analytic primitive orientation") {
  test_that("mirroring and explicit reversal compose for every hit path") {
    auto mat = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(1)));
    random_gen rng(615); RandomSampler sampler(rng);
    for (bool mirror : {false, true}) for (bool reverse : {false, true}) {
      Transform transform = Translate(vec3f(3,-2,5)) * RotateZ(23) *
        Scale(mirror ? -2 : 2, 1.2f, .8f);
      Transform inverse = Inverse(transform);
      sphere ball(1,mat,nullptr,nullptr,&transform,&inverse,reverse);
      ellipsoid oval(point3f(0),1,vec3f(1,1.4f,.7f),mat,nullptr,nullptr,&transform,&inverse,reverse);
      disk plate(vec3f(0),1,0,mat,nullptr,nullptr,&transform,&inverse,reverse);
      cylinder tube(1,2,0,2*Float(M_PI),true,mat,nullptr,nullptr,&transform,&inverse,reverse);
      const Float sign = mirror != reverse ? -1 : 1;
      for (hitable* shape : std::vector<hitable*>{&ball,&oval,&plate,&tube}) {
        // Front/back roots, both cylinder caps, and rays starting inside.
        for (int side : {-1,1}) for (bool inside : {false,true}) {
          CATCH_INFO(shape->GetName() << " mirror=" << mirror << " reverse=" << reverse
                     << " side=" << side << " inside=" << inside);
          bool vertical = shape == &plate || shape == &tube;
          point3f origin = vertical ? point3f(.2f,side*(inside ? .1f : 3.f),.1f) :
            point3f(side*(inside ? .1f : 3.f),.1f,.1f);
          vec3f direction = vertical ? vec3f(0,inside ? side : -side,0) :
            vec3f(inside ? side : -side,0,0);
          if(shape == &plate && inside) direction = -direction;
          Ray ray(transform(origin), transform(direction));
          hit_record a,b;
          bool hit_a=shape->hit(ray, .0001f, 100, a, rng);
          bool hit_b=shape->hit(ray, .0001f, 100, b, &sampler);
          expect_true(hit_a);
          expect_true(hit_b);
          if(!hit_a || !hit_b) continue;
          expect_true(shape->HitP(ray,.0001f,100,rng));
          expect_true(shape->HitP(ray,.0001f,100,&sampler));
          point3f local = inverse(a.p);
          normal3f n;
          if (shape == &plate) n=normal3f(0,1,0);
          else if (shape == &tube) n=normal3f(0,local[1] > 0 ? 1 : -1,0);
          else if (shape == &oval) n=normal3f(local[0],local[1]/Sqr(1.4f),local[2]/Sqr(.7f));
          else n=convert_to_normal3(local);
          n.make_unit_vector();
          normal3f expected = transform(n * sign); expected.make_unit_vector();
          expect_true(dot(a.normal, expected) > .9999f);
          expect_true(dot(b.normal, expected) > .9999f);
          expect_true(dot(a.geometric_normal,expected) > .9999f);
          expect_true(dot(b.geometric_normal,expected) > .9999f);
          expect_true(dot(a.texture_object_normal,n*sign) > .9999f);
          expect_true(dot(b.texture_object_normal,n*sign) > .9999f);
        }
      }
      // Cylinder side hits and far roots use different branches than the caps.
      for(Float start : {3.f,0.f}) {
        Ray ray(transform(point3f(start,.2f,.1f)), transform(vec3f(-1,0,0)));
        hit_record a,b;
        expect_true(tube.hit(ray,.0001f,100,a,rng));
        expect_true(tube.hit(ray,.0001f,100,b,&sampler));
        point3f local = inverse(a.p);
        normal3f n(local[0],0,local[2]); n.make_unit_vector();
        normal3f expected = transform(n*sign); expected.make_unit_vector();
        expect_true(dot(a.normal,expected) > .9999f);
        expect_true(dot(b.normal,expected) > .9999f);
        expect_true(dot(a.texture_object_normal,n*sign) > .9999f);
        expect_true(dot(b.texture_object_normal,n*sign) > .9999f);
      }
      // Instance placement transforms the already oriented child interaction;
      // it must not apply primitive handedness a second time.
      Transform placement = Translate(vec3f(-4,2,1)) * RotateY(31) * Scale(-1,1,1);
      Transform undo = Inverse(placement);
      instance placed(&ball,&placement,&undo,nullptr);
      Ray ray(transform(point3f(3,.1f,.2f)),transform(vec3f(-1,0,0)));
      hit_record child, a, b;
      expect_true(ball.hit(ray,.0001f,100,child,rng));
      Ray placed_ray = placement(ray);
      expect_true(placed.hit(placed_ray,.0001f,100,a,rng));
      expect_true(placed.hit(placed_ray,.0001f,100,b,&sampler));
      auto expected = placement(child.normal); expected.make_unit_vector();
      expect_true(dot(a.normal,expected) > .9999f);
      expect_true(dot(b.normal,expected) > .9999f);
    }
  }
  test_that("flat cylinder bump normals follow orientation at either side root") {
    unsigned char pixels[16]; std::fill(pixels,pixels+16,128);
    auto bump=std::make_shared<bump_texture>(pixels,4,4,1,1);
    auto tex=std::make_shared<constant_texture>(point3f(1));
    auto mat=std::make_shared<MicrofacetTransmission>(tex,
      new TrowbridgeReitzDistribution(.3f,.3f,nullptr,false),point3f(1.5f),point3f(0));
    random_gen rng(97); RandomSampler sampler(rng);
    for(bool mirror : {false,true}) for(bool reverse : {false,true}) {
      Transform t=RotateY(17)*Scale(mirror ? -1 : 1,1,1), inv=Inverse(t);
      cylinder tube(1,2,0,2*Float(M_PI),true,mat,nullptr,bump,&t,&inv,reverse);
      for(Float origin : {0.f,3.f}) {
        Ray ray(t(point3f(origin,.2f,.1f)),t(vec3f(-1,0,0)));
        hit_record a,b;
        expect_true(tube.hit(ray,.0001f,100,a,rng));
        expect_true(tube.hit(ray,.0001f,100,b,&sampler));
        expect_true(dot(a.bump_normal,a.geometric_normal)>.9999f);
        expect_true(dot(b.bump_normal,b.geometric_normal)>.9999f);
      }
    }
  }
}
#endif
