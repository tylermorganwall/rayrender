#ifndef RAYRENDER_VOLUME_INTERSECTIONS_H
#define RAYRENDER_VOLUME_INTERSECTIONS_H
#include "../core/ray.h"

// Inclusive intersections for ordered medium crossings. The plane equation
// returns exact zero on coordinate-aligned faces; a barycentric distance sum
// can round to either side and disagree between containment and transport rays.
inline bool VolumeTriangleIntersection(const Ray &r, const point3f &a, const point3f &b,
                                       const point3f &c, Float tmin, Float tmax,
                                       Float &t, Float &b0, Float &b1, Float &b2,
                                       double *precise_t = nullptr) {
  const point3f *vertices[3] = {&a, &b, &c};
  int z = MaxDimension(Abs(r.d)), x = (z + 1) % 3, y = (x + 1) % 3;
  double p[3][2], ab[3], ac[3];
  for (int i = 0; i < 3; ++i) {
    double dz = double((*vertices[i])[z]) - r.o[z];
    p[i][0] = double((*vertices[i])[x]) - r.o[x] - double(r.d[x]) / r.d[z] * dz;
    p[i][1] = double((*vertices[i])[y]) - r.o[y] - double(r.d[y]) / r.d[z] * dz;
    ab[i] = double(b[i]) - a[i];
    ac[i] = double(c[i]) - a[i];
  }
  double e[3] = {p[1][0] * p[2][1] - p[1][1] * p[2][0],
                 p[2][0] * p[0][1] - p[2][1] * p[0][0],
                 p[0][0] * p[1][1] - p[0][1] * p[1][0]};
  if ((e[0] < 0 || e[1] < 0 || e[2] < 0) && (e[0] > 0 || e[1] > 0 || e[2] > 0)) return false;
  double det = e[0] + e[1] + e[2];
  if (det == 0) return false;
  double n[3] = {ab[1] * ac[2] - ab[2] * ac[1],
                 ab[2] * ac[0] - ab[0] * ac[2],
                 ab[0] * ac[1] - ab[1] * ac[0]};
  double numerator = 0, denominator = 0;
  for (int i = 0; i < 3; ++i) {
    numerator += n[i] * (double(a[i]) - r.o[i]);
    denominator += n[i] * r.d[i];
  }
  if (denominator == 0) return false;
  double distance = numerator / denominator;
  if (distance < std::max({0.0, double(tmin), r.medium_t_min}) || distance > tmax) return false;
  t = Float(distance);
  if (precise_t) *precise_t = distance;
  b0 = Float(e[0] / det); b1 = Float(e[1] / det); b2 = Float(e[2] / det);
  return true;
}

// Zero-width tangencies do not cross a medium boundary. Evaluating the roots
// in double avoids discarding entries near a spawned Float ray origin.
inline bool VolumeSphereRoots(const Ray &ray, Float radius, Float &t_near, Float &t_far) {
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
  t_near = Float(t0);
  t_far = Float(t1);
  return true;
}
#endif
