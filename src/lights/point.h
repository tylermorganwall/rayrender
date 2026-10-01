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
  double Falloff(const vec3f &outgoing) const;
  double Power() const;
  PointLightSample Sample(const point3f &p) const;
};
class PointLightSet {
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
