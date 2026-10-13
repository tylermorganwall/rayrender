#include "preview_navigation.h"
#include "../volumes/picking.h"
#include <algorithm>
#include <vector>

void PreviewNavigation::SetDistance(Float value) {
  if (std::isfinite(value) && value > 0) distance = value;
}

void PreviewNavigation::Initialize(RayCamera& cam, hitable* world, const VolumeScene* volumes,
                                   const aabb& bounds, const std::function<bool()>& cancel) {
  // A central 3x3 stencil tolerates a center-ray miss or foreground outlier.
  // Use axial depth, not off-axis ray length, and ignore environment proxies.
  std::vector<Float> depths;
  for (Float u : {Float(.45), Float(.5), Float(.55)}) {
    for (Float v : {Float(.45), Float(.5), Float(.55)}) {
      if (cancel && cancel()) return;
      Ray ray;
      if (cam.get_fov() < 0) {
        if (!(cam.GenerateRay(CameraSample(point2f(1-u, 1-v), point2f(.5,.5), .5), &ray) > 0))
          continue;
      } else {
        ray = cam.get_ray(u, v, point3f(0), .5);
      }
      auto hit = PickRay(ray, world, volumes, .15, cancel);
      if (!hit || hit->background) continue;
      Float depth = dot(hit->p - cam.get_origin(), unit_vector(cam.get_w()));
      if (std::isfinite(depth) && depth > 0) depths.push_back(depth);
    }
  }
  SetDistance((cam.get_lookat() - cam.get_origin()).length());
  if (!depths.empty()) {
    std::sort(depths.begin(), depths.end());
    SetDistance(depths[depths.size()/2]);
  } else {
    // Bounds must exclude the environment sphere. They are only a fallback:
    // a huge floor should not control navigation when the subject is visible.
    bool valid = true;
    for (int axis = 0; axis < 3; ++axis)
      valid &= std::isfinite(bounds.min()[axis]) && std::isfinite(bounds.max()[axis]) &&
               bounds.max()[axis] >= bounds.min()[axis];
    if (valid) {
      Float depth = dot(bounds.Centroid() - cam.get_origin(), unit_vector(cam.get_w()));
      SetDistance(std::max(depth, bounds.Diag().length() * Float(.25)));
    }
  }
  initial_distance = distance;
  Reset(cam);
}

void PreviewNavigation::Reset(RayCamera& cam) {
  distance = initial_distance;
  Float focal = cam.get_focal_distance();
  cam.update_lookat(cam.get_origin() + unit_vector(cam.get_w()) * distance);
  cam.update_focal_absolute(focal);
}

Float PreviewNavigation::ViewHeight(RayCamera& cam) const {
  Float fov = cam.get_fov();
  if (fov == 0) return cam.get_ortho()[1];
  // A panoramic image has no single perspective plane. For it (and lens
  // cameras), use a stable 60-degree navigation cone without altering optics.
  if (fov < 0 || fov >= 180) fov = 60;
  return 2 * distance * std::tan(fov * Float(M_PI / 360));
}

bool PreviewNavigation::Move(RayCamera& cam, Float horizontal, Float vertical, Float forward,
                             bool shifted, bool orbit, Float multiplier, double seconds) {
  if (!(seconds > 0) || !std::isfinite(seconds) || !(multiplier > 0) ||
      (horizontal == 0 && vertical == 0 && forward == 0)) return false;
  const Float dt = seconds * multiplier;
  const Float focal = cam.get_focal_distance();
  if (shifted) {
    // Pitch and roll are angular velocities, independent of scene units/FOV.
    if (forward != 0) cam.rotate_forward(30 * dt * forward);
    if (horizontal != 0) cam.rotate_navigation_up(30 * dt * horizontal);
    cam.update_lookat(cam.get_origin() + unit_vector(cam.get_w()) * distance);
    if (vertical != 0) Move(cam, 0, vertical, 0, false, orbit, multiplier, seconds);
  } else if (orbit) {
    point3f pivot = cam.get_lookat();
    vec3f offset = cam.get_origin() - pivot;
    vec3f tangent = -horizontal * cam.get_navigation_u() + vertical * cam.get_v();
    Float strength = tangent.length();
    if (strength > 0 && offset.length() > 0) {
      vec3f axis = cross(offset, tangent);
      if (axis.length() > 0)
        offset = Rotate(30 * dt * std::min(strength, Float(1)), unit_vector(axis))(offset);
    }
    // Exponential dolly moves a fraction of the remaining radius and cannot
    // cross the pivot. Integrating several short intervals equals one long one.
    if (forward != 0) offset *= std::exp(-Float(.8) * dt * forward);
    point3f next_origin = pivot + offset;
    Float radius = (next_origin - pivot).length();
    // Stop before the radius rounds to zero in world coordinates: camera
    // frame construction cannot normalize an origin coincident with its target.
    if (radius > 0 && std::isfinite(radius)) {
      cam.update_position_absolute(next_origin);
      cam.update_lookat(pivot);
      SetDistance(radius);
    }
  } else {
    vec3f lateral = -horizontal * cam.get_navigation_u() + vertical * cam.get_v();
    if (lateral.length() > 1) lateral = unit_vector(lateral);
    vec3f delta = dt * (Float(.5) * ViewHeight(cam) * lateral +
                       Float(.5) * distance * forward * unit_vector(cam.get_w()));
    vec3f direction = unit_vector(cam.get_w());
    cam.update_position(delta, false, false);
    cam.update_lookat(cam.get_origin() + direction * distance);
  }
  cam.update_focal_absolute(focal);
  return true;
}

#ifdef NOT_CRAN
#include "../hitables/sphere.h"
#include "../hitables/hitablelist.h"
#include "../hitables/infinite_area_light.h"
#include "../materials/material.h"
#include <testthat.h>

context("Preview navigation") {
  test_that("visible geometry sets navigation depth without changing the image or lens") {
    Transform identity;
    auto material = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
    sphere subject(2, material, nullptr, nullptr, &identity, &identity, false);
    camera cam(point3f(0,0,-10), point3f(0,0,-9), vec3f(0,1,0), 30, 1, .2, 3, 0, 1, 1);
    cam.set_camera_flip_x(true);
    Ray before = cam.get_ray(.3,.6,point3f(.1,.2,0),.5);
    PreviewNavigation nav;
    nav.Initialize(cam, &subject, nullptr, aabb(point3f(-2),point3f(2)));
    Ray after = cam.get_ray(.3,.6,point3f(.1,.2,0),.5);
    expect_true((nav.Distance() > 7.9 && nav.Distance() < 8.2));
    expect_true((cam.get_lookat()-cam.get_origin()).length() == Approx(nav.Distance()));
    expect_true(cam.get_focal_distance() == 3);
    expect_true((before.d-after.d).length() < 1e-6);
    expect_true((before.o-after.o).length() < 1e-6);
    nav.SetDistance(20);
    cam.reset(); nav.Reset(cam);
    expect_true((nav.Distance() > 7.9 && nav.Distance() < 8.2));
  }
  test_that("bounds and empty scenes provide finite fallback scales") {
    camera cam(point3f(0,0,-10),point3f(0,0,-9),vec3f(0,1,0),60,1,0,1,0,1,1);
    PreviewNavigation nav;
    nav.Initialize(cam, nullptr, nullptr, aabb(point3f(-2),point3f(2)));
    expect_true(nav.Distance() == Approx(10));
    cam.reset();
    nav.Initialize(cam, nullptr, nullptr, aabb());
    expect_true(nav.Distance() == Approx(1));
  }
  test_that("environment proxies and isolated foreground hits do not set the scale") {
    Transform identity, foreground = Translate(vec3f(0,0,-8));
    Transform foreground_inverse = Inverse(foreground);
    auto material = std::make_shared<diffuse_material>(std::make_shared<constant_texture>(point3f(.5)));
    auto light = std::make_shared<ImageInfiniteLight>(
      std::make_shared<constant_texture>(point3f(1)), 4, 2, 0);
    hitable_list world;
    world.add(std::make_shared<InfiniteAreaLight>(light, 10000, point3f(0), &identity, &identity));
    camera cam(point3f(0,0,-10),point3f(0,0,-9),vec3f(0,1,0),30,1,0,1,0,1,1);
    PreviewNavigation nav;
    nav.Initialize(cam, &world, nullptr, aabb());
    expect_true(nav.Distance() == Approx(1));
    world.add(std::make_shared<sphere>(2, material, nullptr, nullptr, &identity, &identity, false));
    world.add(std::make_shared<sphere>(.01, material, nullptr, nullptr, &foreground, &foreground_inverse, false));
    nav.Initialize(cam, &world, nullptr, aabb());
    expect_true((nav.Distance() > 7.9 && nav.Distance() < 8.2));
  }
  test_that("movement scales with time, FOV and scene units while orbit uses angles") {
    for (bool flipped : {false,true}) {
      for (Float scale : {Float(.01),Float(1),Float(1000)}) {
        camera a(point3f(0,0,-10*scale),point3f(0),vec3f(0,1,0),60,1,0,3*scale,0,1,1);
        camera b = a;
        a.set_camera_flip_x(flipped); b.set_camera_flip_x(flipped);
        PreviewNavigation na, nb;
        na.SetDistance(10*scale); nb.SetDistance(10*scale);
        na.Move(a,1,0,0,false,false,1,1);
        for(int i=0;i<60;++i) nb.Move(b,1,0,0,false,false,1,1.0/60);
        expect_true((a.get_origin()-b.get_origin()).length()/scale < 1e-4);
        expect_true(dot(a.get_origin()-point3f(0,0,-10*scale), -a.get_navigation_u()) > 0);
        expect_true(a.get_focal_distance()/scale == Approx(3));
        a.reset(); b.reset();
        na.Move(a,1,0,0,false,true,1,1);
        for(int i=0;i<60;++i) nb.Move(b,1,0,0,false,true,1,1.0/60);
        expect_true((a.get_origin()-b.get_origin()).length()/scale < 1e-4);
        expect_true((a.get_origin()-a.get_lookat()).length()/scale == Approx(10));
        a.reset(); b.reset();
        na.Move(a,0,0,1,false,true,1,1);
        for(int i=0;i<60;++i) nb.Move(b,0,0,1,false,true,1,1.0/60);
        expect_true((a.get_origin()-b.get_origin()).length()/scale < 1e-4);
        expect_true(dot(a.get_w(),a.get_lookat()-a.get_origin()) > 0);
      }
    }
    camera wide(point3f(0,0,-10),point3f(0),vec3f(0,1,0),60,1,0,1,0,1,1);
    camera narrow(point3f(0,0,-10),point3f(0),vec3f(0,1,0),30,1,0,1,0,1,1);
    PreviewNavigation nav; nav.SetDistance(10);
    expect_true(nav.ViewHeight(wide)/nav.ViewHeight(narrow) == Approx(std::tan(M_PI/6)/std::tan(M_PI/12)));
    ortho_camera ortho(point3f(0,0,-10),point3f(0),vec3f(0,1,0),8,4,0,1,1);
    expect_true(nav.ViewHeight(ortho) == 4);
    nav.SetDistance(100);
    expect_true(nav.ViewHeight(ortho) == 4);
  }
}
#endif
