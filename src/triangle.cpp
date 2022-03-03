#include "triangle.h"
#include "RcppThread.h"

bool triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) {
  point3f p0t = a - vec3f(r.origin());
  point3f p1t = b - vec3f(r.origin());
  point3f p2t = c - vec3f(r.origin());
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
  Float e0 = p1t.x() * p2t.y() - p1t.y() * p2t.x();
  Float e1 = p2t.x() * p0t.y() - p2t.y() * p0t.x();
  Float e2 = p0t.x() * p1t.y() - p0t.y() * p1t.x();
  
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
  uv[0] = point2f(0, 0);
  uv[1] = point2f(1, 0);
  uv[2] = point2f(1, 1);
  
  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = a - c, dp12 = b - c;
  Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
  bool degenerateUV = std::fabs(determinant) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinant;
    dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || cross(dpdu, dpdv).squared_length() == 0) {
    // Handle zero determinant for triangle partial derivative matrix
    vec3f ng = cross(c - a, b - a);
    if (ng.squared_length() == 0) {
      // The triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    }
    CoordinateSystem(unit_vector(ng), &dpdu, &dpdv);
  }
  //Add error calc
  Float xAbsSum = (std::abs(b0 * a.x()) + std::abs(b1 * b.x()) +
    std::abs(b2 * c.x()));
  Float yAbsSum = (std::abs(b0 * a.y()) + std::abs(b1 * b.y()) +
    std::abs(b2 * c.y()));
  Float zAbsSum = (std::abs(b0 * a.z()) + std::abs(b1 * b.z()) +
    std::abs(b2 * c.z()));
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  
  point3f pHit = b0 * a + b1 * b + b2 * c;
  // point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
  // Float u = uvHit[0];
  // Float v = uvHit[1];
  
  bool alpha_miss = false;
  
  
  //Calculate UV's using old method
  vec3f pvec = cross(r.direction(), edge2);
  det = dot(pvec, edge1);

  // no culling
  if (std::fabs(det) < 1E-15) {
    return(false);
  }
  Float invdet = 1.0 / det;
  vec3f tvec = vec3f(r.origin()) - a;
  Float u = dot(pvec, tvec) * invdet;
  if (u < 0.0 || u > 1.0) {
    return(false);
  }
  //
  vec3f qvec = cross(tvec, edge1);
  Float v = dot(qvec, r.direction()) * invdet;
  if (v < 0 || u + v > 1.0) {
    return(false);
  }
  // Float t = dot(qvec, edge2) * invdet;
  
  if (t < t_min || t > t_max) {
    return(false);
  }
  if(alpha_mask) {
    if(alpha_mask->channel_value(u, v, rec.p) < rng.unif_rand()) {
      alpha_miss = true;
    }
  }
  Float w = 1 - u - v;
  rec.t = t;
  rec.p = pHit;
  rec.u = u;
  rec.v = v;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);


  rec.has_bump = false;
  
  if(bump_tex) {
    //Get UV values + calculate dpdu/dpdv
    vec3f u_val, v_val;
    vec2f uv[3];
    u_val = bump_tex->u_vec;
    v_val = bump_tex->v_vec;
    uv[0] = vec2f(u_val.x(), v_val.x());
    uv[1] = vec2f(u_val.y(), v_val.y());
    uv[2] = vec2f(u_val.z(), v_val.z());
    vec2f duv02 = uv[0] - uv[2];
    vec2f duv12 = uv[1] - uv[2];
    vec3f dp02 = edge1;
    vec3f dp12 = edge2;
    
    Float determinant = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
    if (determinant == 0) {
      onb uvw;
      uvw.build_from_w(cross(edge2,edge1));
      rec.dpdu = uvw.u();
      rec.dpdv = uvw.v();
    } else {
      Float invdet = 1 / determinant;
      rec.dpdu = -( duv12[1] * dp02 - duv02[1] * dp12) * invdet;
      rec.dpdv = -(-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
  } else {
    rec.dpdu = vec3f(0,0,0);
    rec.dpdv = vec3f(0,0,0);
  }
  // Use that to calculate normals
  if(normals_provided) {
    normal3f normal_temp = w * na + u * nb + v * nc;
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
  if(bump_tex) {
    point3f bvbu = bump_tex->mesh_value(rec.u, rec.v, rec.p);
    rec.bump_normal = dot(r.direction(), normal) < 0 ? 
    rec.normal + normal3f( bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv) :
      rec.normal -  normal3f(bvbu.x() * rec.dpdu - bvbu.y() * rec.dpdv);
    rec.bump_normal.make_unit_vector();
    rec.has_bump = true;
  }
  rec.mat_ptr = mp.get();
  rec.alpha_miss = alpha_miss;
  rec.shape = this;
  
  return(true);
}


bool triangle::hit(const ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) {
  point3f p0t = a - vec3f(r.origin());
  point3f p1t = b - vec3f(r.origin());
  point3f p2t = c - vec3f(r.origin());
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
  Float e0 = p1t.x() * p2t.y() - p1t.y() * p2t.x();
  Float e1 = p2t.x() * p0t.y() - p2t.y() * p0t.x();
  Float e2 = p0t.x() * p1t.y() - p0t.y() * p1t.x();
  
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
  if (det < 0 && (tScaled >= 0 || tScaled < r.tMax * det)) {
    return false;
  } else if (det > 0 && (tScaled <= 0 || tScaled > r.tMax * det)) {
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
  uv[0] = point2f(0, 0);
  uv[1] = point2f(1, 0);
  uv[2] = point2f(1, 1);
  
  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = a - c, dp12 = b - c;
  Float determinant = duv02[0] * duv12[1] - duv02[1] * duv12[0];
  bool degenerateUV = std::fabs(determinant) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinant;
    dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || cross(dpdu, dpdv).squared_length() == 0) {
    // Handle zero determinant for triangle partial derivative matrix
    vec3f ng = cross(c - a, b - a);
    if (ng.squared_length() == 0) {
      // The triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    }
    CoordinateSystem(unit_vector(ng), &dpdu, &dpdv);
  }
  //Add error calc
  Float xAbsSum = (std::abs(b0 * a.x()) + std::abs(b1 * b.x()) +
    std::abs(b2 * c.x()));
  Float yAbsSum = (std::abs(b0 * a.y()) + std::abs(b1 * b.y()) +
    std::abs(b2 * c.y()));
  Float zAbsSum = (std::abs(b0 * a.z()) + std::abs(b1 * b.z()) +
    std::abs(b2 * c.z()));

  point3f pHit = b0 * a + b1 * b + b2 * c;
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
  Float u = uvHit[0];
  Float v = uvHit[1];
  
  bool alpha_miss = false;
  // vec3f pvec = cross(r.direction(), edge2);
  // Float det = dot(pvec, edge1);
  
  // no culling
  // if (std::fabs(det) < 1E-15) {
  //   return(false);
  // }
  // Float invdet = 1.0 / det;
  // vec3f tvec = vec3f(r.origin()) - a;
  // Float u = dot(pvec, tvec) * invdet;
  // if (u < 0.0 || u > 1.0) {
  //   return(false);
  // }
  // // 
  // vec3f qvec = cross(tvec, edge1);
  // Float v = dot(qvec, r.direction()) * invdet;
  // if (v < 0 || u + v > 1.0) {
  //   return(false);
  // }
  // Float t = dot(qvec, edge2) * invdet;
  // 
  if (t < t_min || t > t_max) {
    return(false);
  }
  if(alpha_mask) {
    if(alpha_mask->channel_value(u, v, rec.p) < sampler->Get1D()) {
      alpha_miss = true;
    }
  }
  Float w = 1 - u - v;
  rec.t = t;
  rec.p = pHit;
  rec.u = u;
  rec.v = v;
  rec.has_bump = false;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  
  if(bump_tex) {
    //Get UV values + calculate dpdu/dpdv
    vec3f u_val, v_val;
    vec2f uv[3];
    u_val = bump_tex->u_vec;
    v_val = bump_tex->v_vec;
    uv[0] = vec2f(u_val.x(), v_val.x());
    uv[1] = vec2f(u_val.y(), v_val.y());
    uv[2] = vec2f(u_val.z(), v_val.z());
    vec2f duv02 = uv[0] - uv[2];
    vec2f duv12 = uv[1] - uv[2];
    vec3f dp02 = edge1;
    vec3f dp12 = edge2;
    
    Float determinant = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
    if (determinant == 0) {
      onb uvw;
      uvw.build_from_w(cross(edge2,edge1));
      rec.dpdu = uvw.u();
      rec.dpdv = uvw.v();
    } else {
      Float invdet = 1 / determinant;
      rec.dpdu = -( duv12[1] * dp02 - duv02[1] * dp12) * invdet;
      rec.dpdv = -(-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
    }
  } else {
    rec.dpdu = vec3f(0,0,0);
    rec.dpdv = vec3f(0,0,0);
  }
  // Use that to calculate normals
  if(normals_provided) {
    normal3f normal_temp = w * na + u * nb + v * nc;
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
  if(bump_tex) {
    point3f bvbu = bump_tex->mesh_value(rec.u, rec.v, rec.p);
    rec.bump_normal = dot(r.direction(), normal) < 0 ? 
    rec.normal + normal3f( bvbu.x() * rec.dpdu + bvbu.y() * rec.dpdv) :
      rec.normal -  normal3f(bvbu.x() * rec.dpdu - bvbu.y() * rec.dpdv);
    rec.bump_normal.make_unit_vector();
    rec.has_bump = true;
  }
  rec.mat_ptr = mp.get();
  rec.alpha_miss = alpha_miss;
  rec.shape = this;
  
  return(true);
}

bool triangle::bounding_box(Float t0, Float t1, aabb& box) const {
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
    return(distance / (cosine * area));
  }
  return 0; 
}

Float triangle::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) { 
  hit_record rec;
  if (this->hit(ray(o, v), 0.001, FLT_MAX, rec, sampler)) {
    Float distance = rec.t * rec.t * v.squared_length();;
    Float cosine = dot(v, rec.normal);
    return(distance / (cosine * area));
  }
  return 0; 
}

vec3f triangle::random(const point3f& origin, random_gen& rng, Float time) {
  Float r1 = rng.unif_rand();
  Float r2 = rng.unif_rand();
  Float sr1 = sqrt(r1);
  point3f random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin); 
}
vec3f triangle::random(const point3f& origin, Sampler* sampler, Float time) {
  vec2f u = sampler->Get2D();
  Float r1 = u.x();
  Float r2 = u.y();
  Float sr1 = sqrt(r1);
  point3f random_point((1.0 - sr1) * a + sr1 * (1.0 - r2) * b + sr1 * r2 * c);
  return(random_point - origin); 
}
