#include "triangle.h"
#include "RcppThread.h"

bool triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  point3f p0t = p0 - vec3f(r.origin());
  point3f p1t = p1 - vec3f(r.origin());
  point3f p2t = p2 - vec3f(r.origin());
  
  int kz = MaxDimension(Abs(r.direction()));
  int kx = kz + 1;
  if (kx == 3) kx = 0;
  int ky = kx + 1;
  if (ky == 3) ky = 0;
  vec3f d = Permute(r.direction(), kx, ky, kz);
  p0t = Permute(p0t, kx, ky, kz);
  p1t = Permute(p1t, kx, ky, kz);
  p2t = Permute(p2t, kx, ky, kz);

  // Apply shear transformation to translated vertex positions
  Float Sx = -d.x() / d.z();
  Float Sy = -d.y() / d.z();
  Float Sz = 1.f / d.z();
  p0t[0] += Sx * p0t.z();
  p0t[1] += Sy * p0t.z();
  p1t[0] += Sx * p1t.z();
  p1t[1] += Sy * p1t.z();
  p2t[0] += Sx * p2t.z();
  p2t[1] += Sy * p2t.z();
  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = DifferenceOfProducts(p1t.x(), p2t.y(), p1t.y(), p2t.x());
  Float e1 = DifferenceOfProducts(p2t.x(), p0t.y(), p2t.y(), p0t.x());
  Float e2 = DifferenceOfProducts(p0t.x(), p1t.y(), p0t.y(), p1t.x());

  // Fall back to double precision test at triangle edges
  if (sizeof(Float) == sizeof(float) &&
      (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
    double p2txp1ty = (double)p2t.x() * (double)p1t.y();
    double p2typ1tx = (double)p2t.y() * (double)p1t.x();
    e0 = (float)(p2typ1tx - p2txp1ty);
    double p0txp2ty = (double)p0t.x() * (double)p2t.y();
    double p0typ2tx = (double)p0t.y() * (double)p2t.x();
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.x() * (double)p0t.y();
    double p1typ0tx = (double)p1t.y() * (double)p0t.x();
    e2 = (float)(p1typ0tx - p1txp0ty);
  }

  if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
    return false;
  Float det = e0 + e1 + e2;
  if (det == 0) return false;

  p0t[2] *= Sz;
  p1t[2] *= Sz;
  p2t[2] *= Sz;
  Float tScaled = e0 * p0t.z() + e1 * p1t.z() + e2 * p2t.z();
  if (det < 0 && (tScaled >= 0 || tScaled < t_max * det)) {
    return false;
  } else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det)) {
    return false;
  }

  // Compute barycentric coordinates and $t$ value for triangle intersection
  Float invDet = 1 / det;
  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  Float t = tScaled * invDet;
  Float maxZt = MaxComponent(Abs(vec3f(p0t.z(), p1t.z(), p2t.z())));
  Float deltaZ = gamma(3) * maxZt;

  // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
  Float maxXt = MaxComponent(Abs(vec3f(p0t.x(), p1t.x(), p2t.x())));
  Float maxYt = MaxComponent(Abs(vec3f(p0t.y(), p1t.y(), p2t.y())));
  Float deltaX = gamma(5) * (maxXt + maxZt);
  Float deltaY = gamma(5) * (maxYt + maxZt);

  // Compute $\delta_e$ term for triangle $t$ error bounds
  Float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

  // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
  Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
  Float deltaT = 3 *
    (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
    std::abs(invDet);
  if (t <= deltaT) return false;

  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);

  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  Float determinant = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
  bool degenerateUV = std::fabs(determinant) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinant;
    rec.dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    rec.dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || cross(rec.dpdu, rec.dpdv).squared_length() == 0) {
    // Handle zero determinant for triangle partial derivative matrix
    vec3f ng = cross(p2 - p0, p1 - p0);
    if (ng.squared_length() == 0) {
      // The triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    }
    CoordinateSystem(unit_vector(ng), &rec.dpdu, &rec.dpdv);
  }
  //Add error calc
  Float xAbsSum = (std::abs(b0 * p0.x()) + std::abs(b1 * p1.x()) +
    std::abs(b2 * p2.x()));
  Float yAbsSum = (std::abs(b0 * p0.y()) + std::abs(b1 * p1.y()) +
    std::abs(b2 * p2.y()));
  Float zAbsSum = (std::abs(b0 * p0.z()) + std::abs(b1 * p1.z()) +
    std::abs(b2 * p2.z()));
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);

  point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

  Float uHit = uvHit[0];
  Float vHit = uvHit[1];

  bool alpha_miss = false;
  normal3f normal = normal3f(unit_vector(cross(dp02, dp12)));
  
  int mat_id = mesh->face_material_id[face_number];
  alpha_texture* alpha_mask = mesh->alpha_textures[mat_id].get();
  if(alpha_mask) {
    if(alpha_mask->value(uHit, vHit, rec.p) < rng.unif_rand()) {
      alpha_miss = true;
    }
  }
  rec.t = t;
  rec.p = pHit;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);

  rec.has_bump = false;

  // Use that to calculate normals
  if(mesh->has_normals) {
    normal3f n1 = mesh->n[n[0]];
    normal3f n2 = mesh->n[n[1]];
    normal3f n3 = mesh->n[n[2]];
    rec.normal = unit_vector(b0 * n1 + b1 * n2 + b2 * n3);
  } else {
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal) < 0 ? normal : -normal;
    } else {
      rec.normal = normal;
    }
  }
  bump_texture* bump_tex = mesh->bump_textures[mat_id].get();
  
  if(bump_tex) {
    vec3f norm_bump = dot(r.direction(), normal) < 0 ? rec.normal.convert_to_vec3() : -rec.normal.convert_to_vec3();

    point3f bvbu = bump_tex->value(uHit, vHit, rec.p);
    rec.bump_normal = cross(rec.dpdu + bvbu.x() * norm_bump ,
                            rec.dpdv - bvbu.y() * norm_bump);
    rec.bump_normal.make_unit_vector();
    rec.has_bump = true;
  }
  rec.u = uHit;
  rec.v = vHit;
  
  rec.mat_ptr = mesh->mtl_materials[mesh->face_material_id[face_number]].get();
  rec.alpha_miss = alpha_miss;
  rec.shape = this;
  
  return(true);
}


bool triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler *sampler) {
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  point3f p0t = p0 - vec3f(r.origin());
  point3f p1t = p1 - vec3f(r.origin());
  point3f p2t = p2 - vec3f(r.origin());
  
  int kz = MaxDimension(Abs(r.direction()));
  int kx = kz + 1;
  if (kx == 3) kx = 0;
  int ky = kx + 1;
  if (ky == 3) ky = 0;
  vec3f d = Permute(r.direction(), kx, ky, kz);
  p0t = Permute(p0t, kx, ky, kz);
  p1t = Permute(p1t, kx, ky, kz);
  p2t = Permute(p2t, kx, ky, kz);
  
  // Apply shear transformation to translated vertex positions
  Float Sx = -d.x() / d.z();
  Float Sy = -d.y() / d.z();
  Float Sz = 1.f / d.z();
  p0t[0] += Sx * p0t.z();
  p0t[1] += Sy * p0t.z();
  p1t[0] += Sx * p1t.z();
  p1t[1] += Sy * p1t.z();
  p2t[0] += Sx * p2t.z();
  p2t[1] += Sy * p2t.z();
  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = DifferenceOfProducts(p1t.x(), p2t.y(), p1t.y(), p2t.x());
  Float e1 = DifferenceOfProducts(p2t.x(), p0t.y(), p2t.y(), p0t.x());
  Float e2 = DifferenceOfProducts(p0t.x(), p1t.y(), p0t.y(), p1t.x());
  
  // Fall back to double precision test at triangle edges
  if (sizeof(Float) == sizeof(float) &&
      (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
    double p2txp1ty = (double)p2t.x() * (double)p1t.y();
    double p2typ1tx = (double)p2t.y() * (double)p1t.x();
    e0 = (float)(p2typ1tx - p2txp1ty);
    double p0txp2ty = (double)p0t.x() * (double)p2t.y();
    double p0typ2tx = (double)p0t.y() * (double)p2t.x();
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.x() * (double)p0t.y();
    double p1typ0tx = (double)p1t.y() * (double)p0t.x();
    e2 = (float)(p1typ0tx - p1txp0ty);
  }
  
  if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
    return false;
  Float det = e0 + e1 + e2;
  if (det == 0) return false;
  
  p0t[2] *= Sz;
  p1t[2] *= Sz;
  p2t[2] *= Sz;
  Float tScaled = e0 * p0t.z() + e1 * p1t.z() + e2 * p2t.z();
  if (det < 0 && (tScaled >= 0 || tScaled < t_max * det)) {
    return false;
  } else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det)) {
    return false;
  }
  
  // Compute barycentric coordinates and $t$ value for triangle intersection
  Float invDet = 1 / det;
  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  Float t = tScaled * invDet;
  Float maxZt = MaxComponent(Abs(vec3f(p0t.z(), p1t.z(), p2t.z())));
  Float deltaZ = gamma(3) * maxZt;
  
  // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
  Float maxXt = MaxComponent(Abs(vec3f(p0t.x(), p1t.x(), p2t.x())));
  Float maxYt = MaxComponent(Abs(vec3f(p0t.y(), p1t.y(), p2t.y())));
  Float deltaX = gamma(5) * (maxXt + maxZt);
  Float deltaY = gamma(5) * (maxYt + maxZt);
  
  // Compute $\delta_e$ term for triangle $t$ error bounds
  Float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);
  
  // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
  Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
  Float deltaT = 3 *
    (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
    std::abs(invDet);
  if (t <= deltaT) return false;
  
  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);
  
  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  Float determinant = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
  bool degenerateUV = std::fabs(determinant) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinant;
    rec.dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    rec.dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || cross(rec.dpdu, rec.dpdv).squared_length() == 0) {
    // Handle zero determinant for triangle partial derivative matrix
    vec3f ng = cross(p2 - p0, p1 - p0);
    if (ng.squared_length() == 0) {
      // The triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    }
    CoordinateSystem(unit_vector(ng), &rec.dpdu, &rec.dpdv);
  }
  //Add error calc
  Float xAbsSum = (std::abs(b0 * p0.x()) + std::abs(b1 * p1.x()) +
    std::abs(b2 * p2.x()));
  Float yAbsSum = (std::abs(b0 * p0.y()) + std::abs(b1 * p1.y()) +
    std::abs(b2 * p2.y()));
  Float zAbsSum = (std::abs(b0 * p0.z()) + std::abs(b1 * p1.z()) +
    std::abs(b2 * p2.z()));
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  
  point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
  
  Float uHit = uvHit[0];
  Float vHit = uvHit[1];
  
  bool alpha_miss = false;
  
  // if(alpha_mask) {
  //   if(alpha_mask->value(u, v, rec.p) < rng.unif_rand()) {
  //     alpha_miss = true;
  //   }
  // }
  rec.t = t;
  rec.p = pHit;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  normal3f normal = normal3f(unit_vector(cross(dp02, dp12)));
  
  rec.has_bump = false;
  
  // Use that to calculate normals
  if(mesh->n) {
    normal3f n1 = mesh->n[n[0]];
    normal3f n2 = mesh->n[n[1]];
    normal3f n3 = mesh->n[n[2]];
    rec.normal = unit_vector(b0 * n1 + b1 * n2 + b2 * n3);
  } else {
    // if(alpha_mask) {
    //   rec.normal = dot(r.direction(), normal) < 0 ? normal : -normal;
    // } else {
      rec.normal = normal;
    // }
  }
  
  // if(bump_tex) {
  //   vec3f norm_bump = dot(r.direction(), normal) < 0 ? rec.normal.convert_to_vec3() : -rec.normal.convert_to_vec3();
  // 
  //   point3f bvbu = bump_tex->value(u, v, rec.p);
  //   rec.bump_normal = cross(rec.dpdu + bvbu.x() * norm_bump ,
  //                           rec.dpdv - bvbu.y() * norm_bump);
  //   rec.bump_normal.make_unit_vector();
  //   rec.has_bump = true;
  // }
  rec.u = uHit;
  rec.v = vHit;
  
  rec.mat_ptr = mesh->mtl_materials[mesh->face_material_id[*v / 3]].get();
  rec.alpha_miss = alpha_miss;
  rec.shape = this;
  
  return(true);
}


bool triangle::bounding_box(Float t0, Float t1, aabb& box) const {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  point3f min_v(fmin(fmin(a.x(), b.x()), c.x()),
                fmin(fmin(a.y(), b.y()), c.y()),
                fmin(fmin(a.z(), b.z()), c.z()));
  point3f max_v(fmax(fmax(a.x(), b.x()), c.x()),
                fmax(fmax(a.y(), b.y()), c.y()),
                fmax(fmax(a.z(), b.z()), c.z()));

  point3f difference = max_v + -min_v;

  if (difference.x() < 1E-5) max_v.e[0] += 1E-5;
  if (difference.y() < 1E-5) max_v.e[1] += 1E-5;
  if (difference.z() < 1E-5) max_v.e[2] += 1E-5;

  box = aabb(min_v, max_v);
  return(true);
}

Float triangle::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) { 
  hit_record rec;
  if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, rng)) {
    Float distance = rec.t * rec.t * v.squared_length();;
    Float cosine = dot(v, rec.normal);
    return(distance / (cosine * Area()));
  }
  return 0; 
}

Float triangle::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) { 
  hit_record rec;
  if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, sampler)) {
    Float distance = rec.t * rec.t * v.squared_length();;
    Float cosine = dot(v, rec.normal);
    return(distance / (cosine * Area()));
  }
  return 0; 
}

vec3f triangle::random(const point3f& origin, random_gen& rng, Float time) {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  Float r1 = rng.unif_rand();
  Float r2 = rng.unif_rand();
  Float sr1 = sqrt(r1);
  point3f random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin);
}
vec3f triangle::random(const point3f& origin, Sampler* sampler, Float time) {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  vec2f u = sampler->Get2D();
  Float r1 = u.x();
  Float r2 = u.y();
  Float sr1 = sqrt(r1);
  point3f random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin);
}

void triangle::GetUVs(point2f uv[3]) const {
  if (mesh->uv) {
    uv[0] = mesh->uv[t[0]];
    uv[1] = mesh->uv[t[1]];
    uv[2] = mesh->uv[t[2]];
  } else {
    uv[0] = point2f(0, 0);
    uv[1] = point2f(1, 0);
    uv[2] = point2f(1, 1);
  }
}

Float triangle::Area() const {
  // Get triangle vertices in _p0_, _p1_, and _p2_
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  return 0.5 * cross(p1 - p0, p2 - p0).length();
}