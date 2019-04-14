#ifndef TRIANGLEH
#define TRIANGLEH

#include "hitable.h"

class triangle : public hitable {
public:
  triangle() {}
  triangle(vec3 _a, vec3 _b, vec3 _c, material *mat) :
  a(_a), b(_b), c(_c), mp(mat) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normal.make_unit_vector();
    normals_provided = false;
  };
  triangle(vec3 _a, vec3 _b, vec3 _c, vec3 _na, vec3 _nb, vec3 _nc, material *mat) :
    a(_a), b(_b), c(_c), na(_na), nb(_nb), nc(_nc), mp(mat) {
    edge1 = b-a;
    edge2 = c-a;
    normal = cross(edge1, edge2);
    area = normal.length()/2;
    normals_provided = true;
  };
  virtual bool hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng);
  virtual bool bounding_box(float t0, float t1, aabb& box) const {
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
  virtual float pdf_value(const vec3& o, rand_point& v, random_gen& rng) { 
    hit_record rec;
    if (this->hit(ray(o, v.dir), 0.001, FLT_MAX, rec, rng)) {
      float distance = rec.t * rec.t * v.dir.squared_length();;
      float cosine = dot(v.dir, normal);
      return(distance / (cosine * area));
    }
    return 0; 
  }
  
  virtual rand_point random(const vec3& origin, random_gen& rng) {
    rand_point temp;
    float r1 = rng.unif_rand();
    float r2 = rng.unif_rand();
    float sr1 = sqrt(r1);
    vec3 random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
    temp.dir = random_point - origin;
    temp.normal = normal;
    return(temp); 
  }
  vec3 normal;
  vec3 a, b, c, na, nb, nc;
  vec3 edge1, edge2;
  float area;
  bool normals_provided;
  material *mp;
};

bool triangle::hit(const ray& r, float t_min, float t_max, hit_record& rec, random_gen& rng) {
  vec3 pvec = cross(r.direction(), edge2);
  float det = dot(pvec, edge1);

  // no culling
  if (std::fabs(det) < 1E-15) {
    return(false);
  }
  float invdet = 1.0 / det;
  vec3 tvec = r.origin() - a;
  float u = dot(pvec, tvec) * invdet;
  if (u < 0.0 || u > 1.0) {
    return(false);
  }

  vec3 qvec = cross(tvec, edge1);
  float v = dot(qvec, r.direction()) * invdet;
  if (v < 0 || u + v > 1.0) {
    return(false);
  }
  float t = dot(qvec, edge2) * invdet; 
  
  if (t < t_min || t > t_max) {
    return(false);
  }
  float w = 1 - u - v;
  rec.t = t;
  rec.p = r.point_at_parameter(t);
  rec.u = u;
  rec.v = v;
  if(normals_provided) {
    rec.normal = w * na + u * nb + v * nc;
  } else {
    rec.normal = normal;
  }
  rec.mat_ptr = mp;
  return(true);
}

#endif