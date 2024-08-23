#ifndef AABBH
#define AABBH

#include "ray.h"
#include "rng.h"
#include "point3.h"
#include "mathinline.h"
#include "sampler.h"

#include <algorithm>
#include "simd.h" // Make sure this includes our SIMD abstraction layer


class aabb {
  public: 
    aabb() {
      Float minNum = std::numeric_limits<Float>::lowest();
      Float maxNum = std::numeric_limits<Float>::max();
      bounds[0] = point3f(maxNum, maxNum, maxNum);
      bounds[1] = point3f(minNum, minNum, minNum);
    }
    aabb(vec3f a) {
      bounds[0] = point3f(a.x(),a.y(),a.z());
      bounds[1] = point3f(a.x(),a.y(),a.z());
    }
    aabb(point3f a) {
      bounds[0] = a;
      bounds[1] = a;
    }
    aabb(const vec3f& a, const vec3f& b) { 
      bounds[0] = point3f(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = point3f(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
    }
    aabb(const point3f& a, const point3f& b) { 
      bounds[0] = point3f(fmin(a.x(), b.x()), fmin(a.y(), b.y()),fmin(a.z(), b.z()));
      bounds[1] = point3f(fmax(a.x(), b.x()), fmax(a.y(), b.y()),fmax(a.z(), b.z()));
    }
    aabb(const aabb &box) {
      bounds[0] = box.bounds[0]; 
      bounds[1] = box.bounds[1];
    } 
    
    point3f min() const {return(bounds[0]);}
    point3f max() const {return(bounds[1]);}
    
    
    const Float surface_area() const;
    const Float Volume() const;
    const point3f Centroid() const;
    const point3f Diag() const;
    
    const bool hit(const ray& r, Float tmin, Float tmax, random_gen& rng) const;
    const bool hit(const ray& r, Float tmin, Float tmax, Sampler* sampler) const;
    
    const point3f offset(const point3f o);
    const point3f offset(const vec3f o);
    
    const point3f Corner(int corner) const;
    const point3f Lerp(const point3f &t) const;
    const point2f Lerp(const point2f &t) const;
    
    void Expand(Float delta);
    point3f bounds[2];
};

inline aabb surrounding_box(aabb box0, aabb box1) {
  point3f small(fmin(box0.min().x(), box1.min().x()),
             fmin(box0.min().y(), box1.min().y()),
             fmin(box0.min().z(), box1.min().z()));
  point3f big(fmax(box0.max().x(), box1.max().x()),
           fmax(box0.max().y(), box1.max().y()),
           fmax(box0.max().z(), box1.max().z()));
  return(aabb(small,big));
}

inline aabb surrounding_box(aabb box0, point3f point1) {
  point3f small(fmin(box0.min().x(), point1.x()),
                fmin(box0.min().y(), point1.y()),
                fmin(box0.min().z(), point1.z()));
  point3f big(fmax(box0.max().x(), point1.x()),
              fmax(box0.max().y(), point1.y()),
              fmax(box0.max().z(), point1.z()));
  return(aabb(small,big));
}

inline aabb surrounding_box(aabb box0, vec3f point1) {
  point3f small(fmin(box0.min().x(), point1.x()),
                fmin(box0.min().y(), point1.y()),
                fmin(box0.min().z(), point1.z()));
  point3f big(fmax(box0.max().x(), point1.x()),
              fmax(box0.max().y(), point1.y()),
              fmax(box0.max().z(), point1.z()));
  return(aabb(small,big));
}

inline aabb Expand(aabb box, Float delta) {
  return(aabb(box.min() + -point3f(delta, delta, delta),
              box.max() + point3f(delta, delta, delta)));
}

inline aabb Expand(aabb box, vec3f delta) {
  return(aabb(box.min() - delta,
              box.max() + delta));
}


inline bool Inside(const point3f &p, const aabb &b) {
  return (p.x() >= b.min().x() && p.x() <= b.max().x() && 
          p.y() >= b.min().y() && p.y() <= b.max().y() && 
          p.z() >= b.min().z() && p.z() <= b.max().z());
}

inline bool InsideExclusive(const point3f &p, const aabb &b) {
  return (p.x() >= b.min().x() && p.x()  < b.max().x() && p.y() >= b.min().y() &&
          p.y() <  b.max().y() && p.z() >= b.min().z() && p.z()  < b.max().z());
}

inline std::ostream& operator<<(std::ostream &os, const aabb &t) {
  os << "Low: " << t.bounds[0] <<  " High: " << t.bounds[1];
  return os;
}


struct SimdAABB {
  SimdFloat min_x, min_y, min_z;
  SimdFloat max_x, max_y, max_z;
  
  bool intersect(const ray& r, Float tmin, Float tmax) const {
    // Load ray origin and direction
    SimdFloat origin_x = simd_set1(r.origin().x());
    SimdFloat origin_y = simd_set1(r.origin().y());
    SimdFloat origin_z = simd_set1(r.origin().z());
    
    SimdFloat inv_dir_x = simd_set1(r.inverse_dir().x());
    SimdFloat inv_dir_y = simd_set1(r.inverse_dir().y());
    SimdFloat inv_dir_z = simd_set1(r.inverse_dir().z());
    
    // Compute t0 and t1 for each axis
    SimdFloat tx1 = simd_mul(simd_sub(min_x, origin_x), inv_dir_x);
    SimdFloat tx2 = simd_mul(simd_sub(max_x, origin_x), inv_dir_x);
    
    SimdFloat ty1 = simd_mul(simd_sub(min_y, origin_y), inv_dir_y);
    SimdFloat ty2 = simd_mul(simd_sub(max_y, origin_y), inv_dir_y);
    
    SimdFloat tz1 = simd_mul(simd_sub(min_z, origin_z), inv_dir_z);
    SimdFloat tz2 = simd_mul(simd_sub(max_z, origin_z), inv_dir_z);
    
    // Compute tmin and tmax for each AABB
    SimdFloat tmin_x = simd_min(tx1, tx2);
    SimdFloat tmax_x = simd_max(tx1, tx2);
    SimdFloat tmin_y = simd_min(ty1, ty2);
    SimdFloat tmax_y = simd_max(ty1, ty2);
    SimdFloat tmin_z = simd_min(tz1, tz2);
    SimdFloat tmax_z = simd_max(tz1, tz2);
    
    SimdFloat tmin_aabb = simd_max(simd_max(tmin_x, tmin_y), simd_max(tmin_z, simd_set1(tmin)));
    SimdFloat tmax_aabb = simd_min(simd_min(tmax_x, tmax_y), simd_min(tmax_z, simd_set1(tmax)));
    
    // Check if there's an intersection
    SimdMask hit_mask = simd_less_equal(tmin_aabb, tmax_aabb);
    
    // Convert SIMD mask to boolean result
    return simd_any_true(hit_mask);
  }
  
  // Constructor taking an array of AABBs
  template<size_t N>
  SimdAABB(const std::array<aabb, N>& aabbs) {
    static_assert(N == SIMD_WIDTH, "Number of AABBs must match SIMD width");
    
    // Temporary arrays to hold the values before loading into SIMD registers
    float min_x_values[N], min_y_values[N], min_z_values[N];
    float max_x_values[N], max_y_values[N], max_z_values[N];
    
    // Fill the temporary arrays
    for (size_t i = 0; i < N; ++i) {
      const auto& aabb = aabbs[i];
      min_x_values[i] = aabb.min().x();
      min_y_values[i] = aabb.min().y();
      min_z_values[i] = aabb.min().z();
      max_x_values[i] = aabb.max().x();
      max_y_values[i] = aabb.max().y();
      max_z_values[i] = aabb.max().z();
    }
    
    // Load the values into SIMD registers
    min_x = simd_load(min_x_values);
    min_y = simd_load(min_y_values);
    min_z = simd_load(min_z_values);
    max_x = simd_load(max_x_values);
    max_y = simd_load(max_y_values);
    max_z = simd_load(max_z_values);
  }
};

template<typename... Args>
SimdAABB make_simd_aabb(Args&&... args);

#endif
