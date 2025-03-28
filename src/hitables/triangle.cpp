#include "../hitables/triangle.h"
#include "RcppThread.h"
#include "../utils/raylog.h"
#include "../math/vectypes.h"

const bool triangle::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Triangle");

  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  vec3f p0t = p0 - r.origin();
  vec3f p1t = p1 - r.origin();
  vec3f p2t = p2 - r.origin();

  {
    int kx = r.kx;
    int ky = r.ky;
    int kz = r.kz;

    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);
  }
  vec3f Svec = r.Svec;
  const Float Sz = Svec.xyz.z;
  vec3f Svec_z = vec3f(1,1,Sz);
  const vec3f zero_z(1,1,0);

  Svec *= zero_z;
  p0t += Svec * p0t.xyz.z;
  p1t += Svec * p1t.xyz.z;
  p2t += Svec * p2t.xyz.z;

  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = DifferenceOfProducts(p1t.xyz.x, p2t.xyz.y, p1t.xyz.y, p2t.xyz.x);
  Float e1 = DifferenceOfProducts(p2t.xyz.x, p0t.xyz.y, p2t.xyz.y, p0t.xyz.x);
  Float e2 = DifferenceOfProducts(p0t.xyz.x, p1t.xyz.y, p0t.xyz.y, p1t.xyz.x);
  // __builtin_prefetch(&p0t[2],1);
  // __builtin_prefetch(&p1t[2],1);
  // __builtin_prefetch(&p2t[2],1);
  // __builtin_prefetch(&Sz,0);

  // Fall back to double precision test at triangle edges
  #ifndef RAY_FLOAT_AS_DOUBLE
  if (e0 == 0.f || e1 == 0.f || e2 == 0.f) [[unlikely]]	{
    double p2txp1ty = (double)p2t.xyz.x * (double)p1t.xyz.y;
    double p2typ1tx = (double)p2t.xyz.y * (double)p1t.xyz.x;
    e0 = (float)(p2typ1tx - p2txp1ty);
    double p0txp2ty = (double)p0t.xyz.x * (double)p2t.xyz.y;
    double p0typ2tx = (double)p0t.xyz.y * (double)p2t.xyz.x;
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.xyz.x * (double)p0t.xyz.y;
    double p1typ0tx = (double)p1t.xyz.y * (double)p0t.xyz.x;
    e2 = (float)(p1typ0tx - p1txp0ty);
  }
  #endif
  // __builtin_prefetch(&p0t[2],1);
  // __builtin_prefetch(&p1t[2]);
  // __builtin_prefetch(&p2t[2]);

  Float det = e0 + e1 + e2;

  // Check if all e0, e1, e2 have the same sign as det
  if ((e0 * det < 0) || (e1 * det < 0) || (e2 * det < 0)) {
      return false;
  }

  p0t *= Svec_z;
  p1t *= Svec_z;
  p2t *= Svec_z;
  // p0t[2] *= Sz;
  // p1t[2] *= Sz;
  // p2t[2] *= Sz;
  Float tScaled = e0 * p0t.xyz.z + e1 * p1t.xyz.z + e2 * p2t.xyz.z;
  if (det < 0 && (tScaled >= 0 || tScaled < t_max * det)) {
    return false;
  } else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det)) {
    return false;
  }

  // Compute barycentric coordinates and $t$ value for triangle intersection
  Float invDet = 1 / det;
  Float t = tScaled * invDet;
  {
    Float maxZt = MaxComponent(Abs(vec3f(p0t.xyz.z, p1t.xyz.z, p2t.xyz.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(vec3f(p0t.xyz.x, p1t.xyz.x, p2t.xyz.x)));
    Float maxYt = MaxComponent(Abs(vec3f(p0t.xyz.y, p1t.xyz.y, p2t.xyz.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
    Float deltaT = 3 *
      (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
      ffabs(invDet);
    if (t <= deltaT) {
      return false;
    }
  }

  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  point3f b0v(b0);
  point3f b1v(b1);
  point3f b2v(b2);

  vec3f bVec(b0, b1, b2);

  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);

  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  normal3f normal = convert_to_normal3(cross(dp02, dp12));

  Float determinantUV = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
  bool degenerateUV = ffabs(determinantUV) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinantUV;
    dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || parallelVectors(dpdu, dpdv)) {
    // Handle zero determinant for triangle partial derivative matrix
    vec3f ng = cross(p2 - p0, p1 - p0);
    if (ng.squared_length() == 0) {
      // The triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    }
    CoordinateSystem(unit_vector(ng), &dpdu, &dpdv);
  }
  //Everything after this could be delayed until later
  
  point3f pHit = b0v * p0 + b1v * p1 + b2v * p2;
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];

  point3f bSum0[3];
  bSum0[0] = Abs(b0v * p0);
  bSum0[1] = Abs(b1v * p1);
  bSum0[2] = Abs(b2v * p2);
  vec3f absSum = convert_to_vec3(bSum0[0] + bSum0[1] + bSum0[2]);

  
  Float uHit = uvHit[0];
  Float vHit = uvHit[1];
  if(mesh->has_vertex_colors) {
    uHit = b0;
    vHit = b1;
  } 
  
  bool alpha_miss = false;
  int mat_id = mesh->face_material_id[face_number];

  alpha_texture* alpha_mask = mesh->alpha_textures[mat_id].get();
  if(alpha_mask) {
    if(alpha_mask->value(uHit, vHit, rec.p) < rng.unif_rand()) {
      alpha_miss = true;
    }
  }
  rec.p = pHit;
  rec.t = t;
  
  // __builtin_prefetch(mesh->bump_textures[mat_id].get()); //WILL_READ_ONLY
  // __builtin_prefetch(mesh->mesh_materials[mat_id].get()); //WILL_READ_ONLY
  // Use that to calculate normals
  if(mesh->has_normals && n[0] != -1 && n[1] != -1 && n[2] != -1) {
    normal3f n1 = mesh->n[n[0]];
    normal3f n2 = mesh->n[n[1]];
    normal3f n3 = mesh->n[n[2]];
    
    vec3f np = convert_to_vec3(b0v * n1 + b1v * n2 + b2v * n3);
    if(np.squared_length() == 0) {
      rec.normal = normal;
    } else {
      np.make_unit_vector();
      if(mesh->has_consistent_normals) {
        bool flip = dot(np, r.direction()) > 0;
        Float af1 = mesh->alpha_v[n[0]];
        Float af2 = mesh->alpha_v[n[1]];
        Float af3 = mesh->alpha_v[n[2]];
        Float af = b0 * af1 + b1 * af2 + b2 * af3;
        vec3f i = (-r.direction()); //normalizing here really increases time--see if you can normalize once
        i *= flip ? -1 : 1;
        Float b = dot(i,np);
        Float q = (1 - 2 * M_1_PI * af) * (1 - (2 * M_1_PI) * af)/(1 + 2 * (1 - 2 * M_1_PI) * af);
        Float g = 1 + q * (b - 1);
        Float rho = sqrt(q * (1 + g) / (1 +b));
        vec3f r1 = (g + rho * b) * np - rho * i;
        rec.normal = convert_to_normal3(i + r1);
      } else {
        rec.normal = convert_to_normal3(np);
      }
    }
  } else {
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal) < 0 ? normal : -normal;
    } else {
      rec.normal = normal;
    }
  }
  rec.normal.make_unit_vector();

  rec.dpdu = dpdu;
  rec.dpdv = dpdv;
  rec.pError = vec3f(gamma(7)) * absSum;

  rec.u = uHit;
  rec.v = vHit;
  rec.has_bump = false;
  rec.alpha_miss = alpha_miss;

  bump_texture* bump_tex = mesh->bump_textures[mat_id].get();

  if(bump_tex) {
    vec3f norm_bump = dot(r.direction(), rec.normal) < 0 ? convert_to_vec3(rec.normal) : -convert_to_vec3(rec.normal);

    point3f bvbu = bump_tex->value(uHit, vHit, rec.p);
    rec.bump_normal = convert_to_normal3(cross(rec.dpdu + bvbu.xyz.x * norm_bump ,
                            rec.dpdv - bvbu.xyz.y * norm_bump));
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = Faceforward(rec.bump_normal,rec.normal);
    rec.has_bump = true;
  }

  rec.shape = this;
  rec.mat_ptr = mesh->mesh_materials[mat_id].get();
  return(true);
}


const bool triangle::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Triangle");
  
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  vec3f p0t = p0 - r.origin();
  vec3f p1t = p1 - r.origin();
  vec3f p2t = p2 - r.origin();
  
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
  Float Sx = -d.xyz.x / d.xyz.z;
  Float Sy = -d.xyz.y / d.xyz.z;
  Float Sz = 1.f / d.xyz.z;
  p0t[0] += Sx * p0t.xyz.z;
  p0t[1] += Sy * p0t.xyz.z;
  p1t[0] += Sx * p1t.xyz.z;
  p1t[1] += Sy * p1t.xyz.z;
  p2t[0] += Sx * p2t.xyz.z;
  p2t[1] += Sy * p2t.xyz.z;
  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = DifferenceOfProducts(p1t.xyz.x, p2t.xyz.y, p1t.xyz.y, p2t.xyz.x);
  Float e1 = DifferenceOfProducts(p2t.xyz.x, p0t.xyz.y, p2t.xyz.y, p0t.xyz.x);
  Float e2 = DifferenceOfProducts(p0t.xyz.x, p1t.xyz.y, p0t.xyz.y, p1t.xyz.x);
  
  // Fall back to double precision test at triangle edges
  if (sizeof(Float) == sizeof(float) &&
      (e0 == 0.0f || e1 == 0.0f || e2 == 0.0f)) {
    double p2txp1ty = (double)p2t.xyz.x * (double)p1t.xyz.y;
    double p2typ1tx = (double)p2t.xyz.y * (double)p1t.xyz.x;
    e0 = (float)(p2typ1tx - p2txp1ty); 
    double p0txp2ty = (double)p0t.xyz.x * (double)p2t.xyz.y;
    double p0typ2tx = (double)p0t.xyz.y * (double)p2t.xyz.x;
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.xyz.x * (double)p0t.xyz.y;
    double p1typ0tx = (double)p1t.xyz.y * (double)p0t.xyz.x;
    e2 = (float)(p1typ0tx - p1txp0ty);
  }
  __builtin_prefetch(&p0t[2], 0, 1);
  __builtin_prefetch(&p1t[2], 0, 1);
  __builtin_prefetch(&p2t[2], 0, 1);

  if ((e0 < 0 || e1 < 0 || e2 < 0) && (e0 > 0 || e1 > 0 || e2 > 0))
    return false;
  Float det = e0 + e1 + e2;
  if (det == 0) return false;
  
  p0t[2] *= Sz;
  p1t[2] *= Sz;
  p2t[2] *= Sz;
  Float tScaled = e0 * p0t.xyz.z + e1 * p1t.xyz.z + e2 * p2t.xyz.z;
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
  Float maxZt = MaxComponent(Abs(vec3f(p0t.xyz.z, p1t.xyz.z, p2t.xyz.z)));
  Float deltaZ = gamma(3) * maxZt;
  
  // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
  Float maxXt = MaxComponent(Abs(vec3f(p0t.xyz.x, p1t.xyz.x, p2t.xyz.x)));
  Float maxYt = MaxComponent(Abs(vec3f(p0t.xyz.y, p1t.xyz.y, p2t.xyz.y)));
  Float deltaX = gamma(5) * (maxXt + maxZt);
  Float deltaY = gamma(5) * (maxYt + maxZt);
  
  // Compute $\delta_e$ term for triangle $t$ error bounds
  Float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);
  
  // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
  Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
  Float deltaT = 3 *
    (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
    ffabs(invDet);
  if (t <= deltaT) return false;
  
  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);
  
  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  Float determinant = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
  bool degenerateUV = std::abs(determinant) < 1e-8;
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
  Float xAbsSum = (ffabs(b0 * p0.xyz.x) + ffabs(b1 * p1.xyz.x) +
    ffabs(b2 * p2.xyz.x));
  Float yAbsSum = (ffabs(b0 * p0.xyz.y) + ffabs(b1 * p1.xyz.y) +
    ffabs(b2 * p2.xyz.y));
  Float zAbsSum = (ffabs(b0 * p0.xyz.z) + ffabs(b1 * p1.xyz.z) +
    ffabs(b2 * p2.xyz.z));
  // rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  
  point3f pHit = b0 * p0 + b1 * p1 + b2 * p2;
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
  
  if(mesh->has_vertex_colors) {
    uvHit = point2f(b0,b1);
  } 
  
  bool alpha_miss = false;
  normal3f normal = convert_to_normal3(unit_vector(cross(dp02, dp12)));
  int mat_id = mesh->face_material_id[face_number];
  
  alpha_texture* alpha_mask = mesh->alpha_textures[mat_id].get();
  if(alpha_mask) {
    if(alpha_mask->value(uvHit[0], uvHit[1], rec.p) < sampler->Get1D()) {
      alpha_miss = true;
    }
  }
  rec.t = t;
  rec.p = pHit;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  rec.has_bump = false;
  
  // Use that to calculate normals
  if(mesh->has_normals && n[0] != -1 && n[1] != -1 && n[2] != -1) {
    normal3f n1 = mesh->n[n[0]];
    normal3f n2 = mesh->n[n[1]];
    normal3f n3 = mesh->n[n[2]];
    normal3f np = (b0 * n1 + b1 * n2 + b2 * n3);
    if(np.squared_length() == 0) {
      rec.normal = normal;
    } else {
      np.make_unit_vector();
      if(mesh->has_consistent_normals) {
        bool flip = dot(np, r.direction()) > 0;
        Float af1 = mesh->alpha_v[n[0]];
        Float af2 = mesh->alpha_v[n[1]];
        Float af3 = mesh->alpha_v[n[2]];
        Float af = b0 * af1 + b1 * af2 + b2 * af3;
        vec3f i = -unit_vector(r.direction());
        i *= flip ? -1 : 1;
        Float b = dot(i,np);
        Float q = (1 - 2 * M_1_PI * af) * (1 - (2 * M_1_PI) * af)/(1 + 2 * (1 - 2 * M_1_PI) * af);
        Float g = 1 + q * (b - 1);
        Float rho = sqrt(q * (1 + g) / (1 +b));
        normal3f r1 = (g + rho * b) * np - rho * normal3f(i.xyz.x,i.xyz.y,i.xyz.z);
        rec.normal = unit_vector(normal3f(i.xyz.x,i.xyz.y,i.xyz.z) + r1);
      } else {
        rec.normal = np;
      }
    }
  } else {
    if(alpha_mask) {
      rec.normal = dot(r.direction(), normal) < 0 ? normal : -normal;
    } else {
      rec.normal = normal;
    }
  }
  bump_texture* bump_tex = mesh->bump_textures[mat_id].get();

  if(bump_tex) {
    vec3f norm_bump = dot(r.direction(), rec.normal) < 0 ? convert_to_vec3(rec.normal) : -convert_to_vec3(rec.normal);
    
    point3f bvbu = bump_tex->value(uvHit[0], uvHit[1], rec.p);
    rec.bump_normal = convert_to_normal3(cross(rec.dpdu + bvbu.xyz.x * norm_bump ,
                            rec.dpdv - bvbu.xyz.y * norm_bump));
    rec.bump_normal.make_unit_vector();
    rec.bump_normal = Faceforward(rec.bump_normal,rec.normal);
    rec.has_bump = true;
  }
  rec.u = uvHit[0];
  rec.v = uvHit[1];
  
  rec.mat_ptr = mesh->mesh_materials[mat_id].get();
  rec.alpha_miss = alpha_miss;
  rec.shape = this;
  
  return(true);
}

bool triangle::HitP(const Ray& r, Float t_min, Float t_max, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Triangle");

  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  vec3f p0t = p0 - r.origin();
  vec3f p1t = p1 - r.origin();
  vec3f p2t = p2 - r.origin();

  {
    int kx = r.kx;
    int ky = r.ky;
    int kz = r.kz;

    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);
  }
  vec3f Svec = r.Svec;
  const Float Sz = Svec.xyz.z;
  vec3f Svec_z = vec3f(1,1,Sz);
  const vec3f zero_z(1,1,0);

  Svec *= zero_z;
  p0t += Svec * p0t.xyz.z;
  p1t += Svec * p1t.xyz.z;
  p2t += Svec * p2t.xyz.z;

  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = DifferenceOfProducts(p1t.xyz.x, p2t.xyz.y, p1t.xyz.y, p2t.xyz.x);
  Float e1 = DifferenceOfProducts(p2t.xyz.x, p0t.xyz.y, p2t.xyz.y, p0t.xyz.x);
  Float e2 = DifferenceOfProducts(p0t.xyz.x, p1t.xyz.y, p0t.xyz.y, p1t.xyz.x);

  // Fall back to double precision test at triangle edges
  #ifndef RAY_FLOAT_AS_DOUBLE
  if (e0 == 0.f || e1 == 0.f || e2 == 0.f) [[unlikely]]	{
    double p2txp1ty = (double)p2t.xyz.x * (double)p1t.xyz.y;
    double p2typ1tx = (double)p2t.xyz.y * (double)p1t.xyz.x;
    e0 = (float)(p2typ1tx - p2txp1ty);
    double p0txp2ty = (double)p0t.xyz.x * (double)p2t.xyz.y;
    double p0typ2tx = (double)p0t.xyz.y * (double)p2t.xyz.x;
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.xyz.x * (double)p0t.xyz.y;
    double p1typ0tx = (double)p1t.xyz.y * (double)p0t.xyz.x;
    e2 = (float)(p1typ0tx - p1txp0ty);
  }
  #endif

  Float det = e0 + e1 + e2;

  // Check if all e0, e1, e2 have the same sign as det
  if ((e0 * det < 0) || (e1 * det < 0) || (e2 * det < 0)) {
      return false;
  }

  p0t *= Svec_z;
  p1t *= Svec_z;
  p2t *= Svec_z;

  Float tScaled = e0 * p0t.xyz.z + e1 * p1t.xyz.z + e2 * p2t.xyz.z;
  if (det < 0 && (tScaled >= 0 || tScaled < t_max * det)) {
    return false;
  } else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det)) {
    return false;
  }

  // Compute barycentric coordinates and $t$ value for triangle intersection
  Float invDet = 1 / det;
  Float t = tScaled * invDet;
  {
    Float maxZt = MaxComponent(Abs(vec3f(p0t.xyz.z, p1t.xyz.z, p2t.xyz.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(vec3f(p0t.xyz.x, p1t.xyz.x, p2t.xyz.x)));
    Float maxYt = MaxComponent(Abs(vec3f(p0t.xyz.y, p1t.xyz.y, p2t.xyz.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
    Float deltaT = 3 *
      (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
      ffabs(invDet);
    if (t <= deltaT) {
      return false;
    }
  }

  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  point3f b0v(b0);
  point3f b1v(b1);
  point3f b2v(b2);

  vec3f bVec(b0, b1, b2);

  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);

  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;

  Float determinantUV = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
  bool degenerateUV = ffabs(determinantUV) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinantUV;
    dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || parallelVectors(dpdu, dpdv)) {
    // Handle zero determinant for triangle partial derivative matrix
    vec3f ng = cross(p2 - p0, p1 - p0);
    if (ng.squared_length() == 0) {
      // The triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    }
  }
  return(true);
}

bool triangle::HitP(const Ray& r, Float t_min, Float t_max, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Triangle");

  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  vec3f p0t = p0 - r.origin();
  vec3f p1t = p1 - r.origin();
  vec3f p2t = p2 - r.origin();

  {
    int kx = r.kx;
    int ky = r.ky;
    int kz = r.kz;

    p0t = Permute(p0t, kx, ky, kz);
    p1t = Permute(p1t, kx, ky, kz);
    p2t = Permute(p2t, kx, ky, kz);
  }
  vec3f Svec = r.Svec;
  const Float Sz = Svec.xyz.z;
  vec3f Svec_z = vec3f(1,1,Sz);
  const vec3f zero_z(1,1,0);

  Svec *= zero_z;
  p0t += Svec * p0t.xyz.z;
  p1t += Svec * p1t.xyz.z;
  p2t += Svec * p2t.xyz.z;

  // Compute edge function coefficients _e0_, _e1_, and _e2_
  Float e0 = DifferenceOfProducts(p1t.xyz.x, p2t.xyz.y, p1t.xyz.y, p2t.xyz.x);
  Float e1 = DifferenceOfProducts(p2t.xyz.x, p0t.xyz.y, p2t.xyz.y, p0t.xyz.x);
  Float e2 = DifferenceOfProducts(p0t.xyz.x, p1t.xyz.y, p0t.xyz.y, p1t.xyz.x);

  // Fall back to double precision test at triangle edges
  #ifndef RAY_FLOAT_AS_DOUBLE
  if (e0 == 0.f || e1 == 0.f || e2 == 0.f) [[unlikely]]	{
    double p2txp1ty = (double)p2t.xyz.x * (double)p1t.xyz.y;
    double p2typ1tx = (double)p2t.xyz.y * (double)p1t.xyz.x;
    e0 = (float)(p2typ1tx - p2txp1ty);
    double p0txp2ty = (double)p0t.xyz.x * (double)p2t.xyz.y;
    double p0typ2tx = (double)p0t.xyz.y * (double)p2t.xyz.x;
    e1 = (float)(p0typ2tx - p0txp2ty);
    double p1txp0ty = (double)p1t.xyz.x * (double)p0t.xyz.y;
    double p1typ0tx = (double)p1t.xyz.y * (double)p0t.xyz.x;
    e2 = (float)(p1typ0tx - p1txp0ty);
  }
  #endif

  Float det = e0 + e1 + e2;

  // Check if all e0, e1, e2 have the same sign as det
  if ((e0 * det < 0) || (e1 * det < 0) || (e2 * det < 0)) {
      return false;
  }

  p0t *= Svec_z;
  p1t *= Svec_z;
  p2t *= Svec_z;

  Float tScaled = e0 * p0t.xyz.z + e1 * p1t.xyz.z + e2 * p2t.xyz.z;
  if (det < 0 && (tScaled >= 0 || tScaled < t_max * det)) {
    return false;
  } else if (det > 0 && (tScaled <= 0 || tScaled > t_max * det)) {
    return false;
  }

  // Compute barycentric coordinates and $t$ value for triangle intersection
  Float invDet = 1 / det;
  Float t = tScaled * invDet;
  {
    Float maxZt = MaxComponent(Abs(vec3f(p0t.xyz.z, p1t.xyz.z, p2t.xyz.z)));
    Float deltaZ = gamma(3) * maxZt;

    // Compute $\delta_x$ and $\delta_y$ terms for triangle $t$ error bounds
    Float maxXt = MaxComponent(Abs(vec3f(p0t.xyz.x, p1t.xyz.x, p2t.xyz.x)));
    Float maxYt = MaxComponent(Abs(vec3f(p0t.xyz.y, p1t.xyz.y, p2t.xyz.y)));
    Float deltaX = gamma(5) * (maxXt + maxZt);
    Float deltaY = gamma(5) * (maxYt + maxZt);

    // Compute $\delta_e$ term for triangle $t$ error bounds
    Float deltaE = 2 * (gamma(2) * maxXt * maxYt + deltaY * maxXt + deltaX * maxYt);

    // Compute $\delta_t$ term for triangle $t$ error bounds and check _t_
    Float maxE = MaxComponent(Abs(vec3f(e0, e1, e2)));
    Float deltaT = 3 *
      (gamma(3) * maxE * maxZt + deltaE * maxZt + deltaZ * maxE) *
      ffabs(invDet);
    if (t <= deltaT) {
      return false;
    }
  }

  Float b0 = e0 * invDet;
  Float b1 = e1 * invDet;
  Float b2 = e2 * invDet;
  point3f b0v(b0);
  point3f b1v(b1);
  point3f b2v(b2);

  vec3f bVec(b0, b1, b2);

  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);

  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;

  Float determinantUV = DifferenceOfProducts(duv02[0],duv12[1],duv02[1],duv12[0]);
  bool degenerateUV = ffabs(determinantUV) < 1e-8;
  if (!degenerateUV) {
    Float invdet = 1 / determinantUV;
    dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * invdet;
    dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * invdet;
  }
  if (degenerateUV || parallelVectors(dpdu, dpdv)) {
    // Handle zero determinant for triangle partial derivative matrix
    vec3f ng = cross(p2 - p0, p1 - p0);
    if (ng.squared_length() == 0) {
      // The triangle is actually degenerate; the intersection is
      // bogus.
      return false;
    }
  }
  return(true);
}


bool triangle::bounding_box(Float t0, Float t1, aabb& box) const {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  point3f min_v(ffmin(ffmin(a.xyz.x, b.xyz.x), c.xyz.x),
                ffmin(ffmin(a.xyz.y, b.xyz.y), c.xyz.y),
                ffmin(ffmin(a.xyz.z, b.xyz.z), c.xyz.z));
  point3f max_v(ffmax(ffmax(a.xyz.x, b.xyz.x), c.xyz.x),
                ffmax(ffmax(a.xyz.y, b.xyz.y), c.xyz.y),
                ffmax(ffmax(a.xyz.z, b.xyz.z), c.xyz.z));

  point3f difference = max_v + -min_v;

  if (difference.xyz.x < 1E-5) max_v.e[0] += 1E-5;
  if (difference.xyz.y < 1E-5) max_v.e[1] += 1E-5;
  if (difference.xyz.z < 1E-5) max_v.e[2] += 1E-5;

  box = aabb(min_v, max_v);
  return(true);
}

Float triangle::pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) { 
  if (this->HitP(Ray(o, v), 0.001, FLT_MAX, rng)) {
    return(1 / SolidAngle(o));
  }
  return 0; 
}

Float triangle::pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) { 
  if (this->HitP(Ray(o, v), 0.001, FLT_MAX, sampler)) {
    return(1 / SolidAngle(o));
  }
  return 0; 
}

Float triangle::SolidAngle(point3f p) const {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  return SphericalTriangleArea(unit_vector(a - p), 
                               unit_vector(b - p),
                               unit_vector(c - p));
}

vec3f triangle::random(const point3f& origin, random_gen& rng, Float time) {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  Float r1 = sqrt(rng.unif_rand());
  Float r2 = rng.unif_rand();
  Float u1 = 1.0f-r1;
  Float u2 = r2*r1;
  point3f random_point = (u1 * a + u2 * b + (1 - u1 - u2) * c);
  return(random_point - origin);
}
vec3f triangle::random(const point3f& origin, Sampler* sampler, Float time) {
  const point3f &a = mesh->p[v[0]];
  const point3f &b = mesh->p[v[1]];
  const point3f &c = mesh->p[v[2]];
  vec2f u = sampler->Get2D();
  Float r1 = sqrt(u.xy.x);
  Float r2 = u.xy.y;
  Float u1 = 1.0f-r1;
  Float u2 = r2*r1;
  point3f random_point = (u1 * a + u2 * b + (1 - u1 - u2) * c);
  return(random_point - origin);
}

void triangle::GetUVs(point2f uv[3]) const {
  if (mesh->has_tex && t[0] != -1 && t[1] != -1 && t[2] != -1) {
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
