#ifndef RAYRENDER_VOLUME_LIGHTS_H
#define RAYRENDER_VOLUME_LIGHTS_H
#include "../hitables/hitablelist.h"

// Direction proposals are evaluated as a mixture in solid angle. Visibility
// supplies the first emitted radiance along that direction, including textured
// emitters and the environment. Selection probability appears only in Pdf().
class VolumeLightAdapter {
public:
  virtual ~VolumeLightAdapter() = default;
  virtual vec3f Sample(const point3f &, Sampler *, Float time) const = 0;
  virtual Float Pdf(const point3f &, const vec3f &, random_gen &, Float time) const = 0;
};
class VolumeLightSampler {
public:
  explicit VolumeLightSampler(const hitable_list &);
  vec3f Sample(const point3f &, Sampler *, Float time) const;
  Float Pdf(const point3f &, const vec3f &, random_gen &, Float time) const;

private:
  std::unique_ptr<VolumeLightAdapter> lights;
};
#endif
