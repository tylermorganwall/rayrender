#ifndef RAYRENDER_POINT_LIGHT_H
#define RAYRENDER_POINT_LIGHT_H
#include "../core/ray.h"
#include <Rcpp.h>
#include <vector>

// Delta-position emitters have a discrete selection probability but no ordinary
// solid-angle PDF. They are sampled independently of the area/infinite mixture.
struct PointLightSample {
  point3f position{0}, radiance{0};
  vec3f wi{0};
  double pmf = 0;
};
struct PointLight {
  point3f position{0}, intensity{0};
  vec3f direction{0, -1, 0};
  bool spot = false;
  double cos_outer = 0, cos_inner = 1;
  bool projection = false;
  vec3f right{1, 0, 0}, up{0, 1, 0};
  double half_width = 1, half_height = 1, projection_power = 0;
  int width = 0, height = 0;
  std::vector<point3f> image; // Linear RGB, row-major, first row at image top.
  point3f Project(const vec3f &outgoing) const;
  double Falloff(const vec3f &outgoing) const;
  double Power() const;
  PointLightSample Sample(const point3f &p) const;
};
class PointLightSet {
  friend class WavefrontSceneCompiler;
public:
  PointLightSet() = default;
  explicit PointLightSet(const Rcpp::List &descriptions);
  bool Empty() const { return lights.empty(); }
  PointLightSample Sample(const point3f &p, double u) const;
private:
  std::vector<PointLight> lights;
  std::vector<double> cdf;
};
#endif
