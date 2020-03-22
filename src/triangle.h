#ifndef TRIANGLEH
#define TRIANGLEH

#include "hitable.h"
#include "material.h"

class triangle : public hitable {
public:
  triangle() {}
  ~triangle() {
    // Rcpp::Rcout << "bvh delete " << typeid(*mp).name() << "\n";
    if(single) {
      delete mp;
      delete alpha_mask;
    }
  }
  triangle(vec3 _a, vec3 _b, vec3 _c, bool _single, material *mat, alpha_texture *alpha_mask) :
  a(_a), b(_b), c(_c), single(_single), mp(mat), alpha_mask(alpha_mask) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normal.make_unit_vector();
    normals_provided = false;
  };
  triangle(vec3 _a, vec3 _b, vec3 _c, vec3 _na, vec3 _nb, vec3 _nc, bool _single, 
           material *mat, alpha_texture *alpha_mask) :
    a(_a), b(_b), c(_c), na(_na), nb(_nb), nc(_nc), single(_single), mp(mat), alpha_mask(alpha_mask) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normals_provided = true;
  };
  virtual bool hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(Float t0, Float t1, aabb& box) const {
    vec3 min_v(fmin(fmin(a.x(), b.x()), c.x()), 
               fmin(fmin(a.y(), b.y()), c.y()), 
               fmin(fmin(a.z(), b.z()), c.z()));
    vec3 max_v(fmax(fmax(a.x(), b.x()), c.x()), 
               fmax(fmax(a.y(), b.y()), c.y()), 
               fmax(fmax(a.z(), b.z()), c.z()));
    
    vec3 difference = max_v - min_v;
    
    if (difference.x() < 1E-5) max_v.e[0] += 1E-5;
    if (difference.y() < 1E-5) max_v.e[1] += 1E-5;
    if (difference.z() < 1E-5) max_v.e[2] += 1E-5;
    
    box = aabb(min_v, max_v);
    return(true);
  }
  virtual Float pdf_value(const vec3& o, const vec3& v, random_gen& rng) { 
    hit_record rec;
    if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, rng)) {
      Float distance = rec.t * rec.t * v.squared_length();;
      Float cosine = dot(v, rec.normal);
      return(distance / (cosine * area));
    }
    return 0; 
  }
  
  virtual vec3 random(const vec3& origin, random_gen& rng) {
    Float r1 = rng.unif_rand();
    Float r2 = rng.unif_rand();
    Float sr1 = sqrt(r1);
    vec3 random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
    return(random_point - origin); 
  }
  vec3 normal;
  vec3 a, b, c, na, nb, nc;
  vec3 edge1, edge2;
  Float area;
  bool normals_provided;
  bool single;
  material *mp;
  alpha_texture *alpha_mask;
};

bool triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  vec3 pvec = cross(r.direction(), edge2);
  Float det = dot(pvec, edge1);

  // no culling
  if (std::fabs(det) < 1E-15) {
    return(false);
  }
  Float invdet = 1.0 / det;
  vec3 tvec = r.origin() - a;
  Float u = dot(pvec, tvec) * invdet;
  if (u < 0.0 || u > 1.0) {
    return(false);
  }

  vec3 qvec = cross(tvec, edge1);
  Float v = dot(qvec, r.direction()) * invdet;
  if (v < 0 || u + v > 1.0) {
    return(false);
  }
  Float t = dot(qvec, edge2) * invdet; 
  
  if (t < t_min || t > t_max) {
    return(false);
  }
  if(alpha_mask) {
    if(alpha_mask->channel_value(u, v, rec.p) < rng.unif_rand()) {
      return(false);
    }
  }
  Float w = 1 - u - v;
  rec.t = t;
  rec.p = r.point_at_parameter(t);
  rec.u = u;
  rec.v = v;
  if(normals_provided) {
    vec3 normal_temp = w * na + u * nb + v * nc;
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal_temp) < 0 ? normal_temp : -normal_temp;
    } else {
      rec.normal = normal_temp;
    }
  } else {
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal) < 0 ? normal : -normal;
    } else {
      rec.normal = normal;
    }
  }
  rec.mat_ptr = mp;
  return(true);
}

#endif
