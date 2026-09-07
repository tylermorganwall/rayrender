#ifndef RAYRENDER_VOLUME_INTERSECTIONS_H
#define RAYRENDER_VOLUME_INTERSECTIONS_H
#include "../core/ray.h"

// Zero-width tangencies do not cross a medium boundary. Evaluating the roots
// in double avoids discarding entries near a spawned Float ray origin.
inline bool VolumeSphereRoots(const Ray &ray, Float radius, Float &near, Float &far) {
  double a = 0, b = 0, c = -double(radius) * radius;
  for (int i = 0; i < 3; ++i) {
    a += double(ray.d[i]) * ray.d[i];
    b += 2 * double(ray.d[i]) * ray.o[i];
    c += double(ray.o[i]) * ray.o[i];
  }
  double discriminant = b * b - 4 * a * c;
  if (!(discriminant > 0))
    return false;
  double q = -.5 * (b + std::copysign(std::sqrt(discriminant), b));
  double t0 = q / a, t1 = c / q;
  if (t0 > t1)
    std::swap(t0, t1);
  near = Float(t0);
  far = Float(t1);
  return true;
}
#endif
