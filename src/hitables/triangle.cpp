#include "../materials/normalmap.h"
#include "../hitables/triangle.h"
#include "RcppThread.h"
#include "../utils/raylog.h"
#include "../math/vectypes.h"
#include "../volumes/intersections.h"


OpaqueShadowType triangle::ShadowType() const {
  const int id = mesh->face_material_id[face_number];
  // Alpha masks require the ordered stochastic walk; opaque lamps block light too.
  return mesh->alpha_textures[id] ? OpaqueShadowType::Unsupported
      : opaque_shadow_material(mesh->mesh_materials[id].get());
}

bool triangle::OpaqueHit(const Ray& r, Float t_min, Float t_max, random_gen& rng) const {
  if (!r.segment_absorption) {
    hit_record rec;
    return hit(r, t_min, t_max, rec, rng);
  }
  const point3f &p0 = mesh->p[v[0]], &p1 = mesh->p[v[1]], &p2 = mesh->p[v[2]];
  Float distance, b0, b1, b2;
  if (!VolumeTriangleIntersection(r, p0, p1, p2, t_min, t_max, distance, b0, b1, b2))
    return false;

  // Match hit()'s degenerate-geometry rejection before omitting the shading
  // record: no interpolated normals, textures, bump map, or position error.
  point2f uv[3];
  GetUVs(uv);
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  Float determinant = DifferenceOfProducts(duv02[0], duv12[1], duv02[1], duv12[0]);
  bool degenerate = ffabs(determinant) < 1e-8;
  if (!degenerate) {
    Float inv = 1 / determinant;
    vec3f dpdu = (duv12[1] * dp02 - duv02[1] * dp12) * inv;
    vec3f dpdv = (-duv12[0] * dp02 + duv02[0] * dp12) * inv;
    degenerate = parallelVectors(dpdu, dpdv);
  }
  return !degenerate || cross(p2 - p0, p1 - p0).squared_length() != 0;
}

const bool triangle::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, random_gen& rng) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Triangle");

  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  Float t,b0,b1,b2;
  double precise_t = INFINITY;
  if(r.segment_absorption) {
    if(!VolumeTriangleIntersection(r,p0,p1,p2,t_min,t_max,t,b0,b1,b2,&precise_t)) return false;
  } else {
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
  t = tScaled * invDet;
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

  b0 = e0 * invDet;
  b1 = e1 * invDet;
  b2 = e2 * invDet;


  }
  point3f b0v(b0), b1v(b1), b2v(b2);
  vec3f dpdu, dpdv;
  point2f uv[3];
  GetUVs(uv);

  // Compute deltas for triangle partial derivatives
  vec2f duv02 = uv[0] - uv[2], duv12 = uv[1] - uv[2];
  vec3f dp02 = p0 - p2, dp12 = p1 - p2;
  normal3f normal = convert_to_normal3(cross(dp02, dp12));
  // Vertices are already transformed, so their cross product contains the
  // transform's handedness. Remove that sign, then apply explicit reversal.
  if (reverseOrientation ^ transformSwapsHandedness) normal = -normal;

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
  if (r.segment_absorption)
    for (int a = 0; a < 3; ++a)
      if (p0[a] == p1[a] && p0[a] == p2[a]) pHit[a] = p0[a];
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
  rec.geometric_normal = unit_vector(normal);
  rec.physical_shading_normal = rec.geometric_normal;
  rec.p = pHit;
  rec.t = t;
  rec.precise_t = precise_t;
  
  // __builtin_prefetch(mesh->bump_textures[mat_id].get()); //WILL_READ_ONLY
  // __builtin_prefetch(mesh->mesh_materials[mat_id].get()); //WILL_READ_ONLY
  // Use that to calculate normals
  if(mesh->has_normals && n[0] != -1 && n[1] != -1 && n[2] != -1) {
    normal3f n1 = mesh->n[n[0]];
    normal3f n2 = mesh->n[n[1]];
    normal3f n3 = mesh->n[n[2]];
    
    vec3f np = convert_to_vec3(b0v * n1 + b1v * n2 + b2v * n3);
    // Vertex normals were inverse-transpose transformed, without the cross
    // product's handedness sign. Only explicit reversal remains to apply.
    if (reverseOrientation) np = -np;
    if(np.squared_length() == 0) {
      rec.normal = normal;
    } else {
      np.make_unit_vector();
      rec.physical_shading_normal = normal3f(np[0],np[1],np[2]);
      if(mesh->has_consistent_normals) {
        bool flip = dot(np, r.direction()) > 0;
        Float af1 = mesh->alpha_v[n[0]];
        Float af2 = mesh->alpha_v[n[1]];
        Float af3 = mesh->alpha_v[n[2]];
        Float af = b0 * af1 + b1 * af2 + b2 * af3;
        vec3f i = -unit_vector(r.direction()); // Match the Sampler overload without changing ray t.
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
  rec.base_shading_normal = rec.physical_shading_normal;
  rec.normal.make_unit_vector();

  rec.dpdu = dpdu;
  rec.dpdv = dpdv;
  rec.pError = vec3f(gamma(7)) * absSum;

  rec.u = uHit;
  rec.v = vHit;
  rec.has_bump = false;
  rec.alpha_miss = alpha_miss;

  bump_texture* bump_tex = mesh->bump_textures[mat_id].get();

  rec.dndu = rec.dndv = normal3f(0);
  const bool smooth_normals = mesh->has_normals && n[0] != -1 && n[1] != -1 && n[2] != -1;
  if (smooth_normals && !degenerateUV) {
    const Float normal_det = DifferenceOfProducts(duv02[0], duv12[1], duv02[1], duv12[0]);
    const auto dn02=mesh->n[n[0]]-mesh->n[n[2]], dn12=mesh->n[n[1]]-mesh->n[n[2]];
    rec.dndu = (duv12[1] * dn02 - duv02[1] * dn12) / normal_det;
    rec.dndv = (-duv12[0] * dn02 + duv02[0] * dn12) / normal_det;
    if (reverseOrientation) {
      rec.dndu = -rec.dndu;
      rec.dndv = -rec.dndv;
    }
  }
  rec.shape = this;
  rec.texture_object_p = (*WorldToObject)(rec.p);
  rec.texture_object_normal = unit_vector((*WorldToObject)(rec.geometric_normal));
  rec.ComputeDifferentials(r);
  if (bump_tex) {
    const auto ns=rec.physical_shading_normal;
    // PBRT's triangle shading frame uses the interpolated normal and projected
    // u tangent. UV footprints above use geometric tangents, before this frame.
    if (smooth_normals) {
      auto bitangent = cross(convert_to_vec3(ns), rec.dpdu);
      if (bitangent.squared_length() > 0) {
        rec.dpdu = cross(bitangent, convert_to_vec3(ns));
        rec.dpdv = bitangent;
      } else
        CoordinateSystem(convert_to_vec3(ns), &rec.dpdu, &rec.dpdv);
    }
    const TextureFootprint footprint{rec.dudx, rec.dvdx, rec.dudy, rec.dvdy, rec.has_differentials};
    TextureEvalContext context;
    if (bump_tex->height_texture) {
      context = TextureEvalContext::FromHit(rec);
      context.u = uvHit[0]; context.v = uvHit[1];
    }
    rec.bump_normal = bump_tex->perturb(uvHit[0], uvHit[1], rec.p, ns, rec.dpdu, rec.dpdv, rec.dndu,
                                        rec.dndv, footprint, bump_tex->height_texture ? &context : nullptr);
    rec.physical_shading_normal=rec.bump_normal;
    rec.has_bump = true;
  }

  rec.mat_ptr = mesh->mesh_materials[mat_id].get();
  return(true);
}


const bool triangle::hit(const Ray& r, Float t_min, Float t_max, hit_record& rec, Sampler* sampler) const {
  SCOPED_CONTEXT("Hit");
  SCOPED_TIMER_COUNTER("Triangle");
  
  const point3f &p0 = mesh->p[v[0]];
  const point3f &p1 = mesh->p[v[1]];
  const point3f &p2 = mesh->p[v[2]];
  
  Float t,b0,b1,b2;
  double precise_t = INFINITY;
  if(r.segment_absorption) {
    if(!VolumeTriangleIntersection(r,p0,p1,p2,t_min,t_max,t,b0,b1,b2,&precise_t)) return false;
  } else {
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
  b0 = e0 * invDet;
  b1 = e1 * invDet;
  b2 = e2 * invDet;
  t = tScaled * invDet;
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
  
  }
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
  if (r.segment_absorption)
    for (int a = 0; a < 3; ++a)
      if (p0[a] == p1[a] && p0[a] == p2[a]) pHit[a] = p0[a];
  point2f uvHit = b0 * uv[0] + b1 * uv[1] + b2 * uv[2];
  
  if(mesh->has_vertex_colors) {
    uvHit = point2f(b0,b1);
  } 
  
  bool alpha_miss = false;
  normal3f normal = convert_to_normal3(unit_vector(cross(dp02, dp12)));
  if (reverseOrientation ^ transformSwapsHandedness) normal = -normal;
  int mat_id = mesh->face_material_id[face_number];
  
  alpha_texture* alpha_mask = mesh->alpha_textures[mat_id].get();
  if(alpha_mask) {
    if(alpha_mask->value(uvHit[0], uvHit[1], rec.p) < sampler->Get1D()) {
      alpha_miss = true;
    }
  }
  rec.t = t;
  rec.geometric_normal = unit_vector(normal);
  rec.physical_shading_normal = rec.geometric_normal;
  rec.precise_t = precise_t;
  rec.p = pHit;
  rec.pError = gamma(7) * vec3f(xAbsSum, yAbsSum, zAbsSum);
  rec.has_bump = false;
  
  // Use that to calculate normals
  if(mesh->has_normals && n[0] != -1 && n[1] != -1 && n[2] != -1) {
    normal3f n1 = mesh->n[n[0]];
    normal3f n2 = mesh->n[n[1]];
    normal3f n3 = mesh->n[n[2]];
    normal3f np = (b0 * n1 + b1 * n2 + b2 * n3);
    if (reverseOrientation) np = -np;
    if(np.squared_length() == 0) {
      rec.normal = normal;
    } else {
      np.make_unit_vector();
      rec.physical_shading_normal = normal3f(np[0],np[1],np[2]);
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
  rec.base_shading_normal = rec.physical_shading_normal;
  bump_texture* bump_tex = mesh->bump_textures[mat_id].get();

  rec.dndu = rec.dndv = normal3f(0);
  const bool smooth_normals = mesh->has_normals && n[0] != -1 && n[1] != -1 && n[2] != -1;
  if (smooth_normals && !degenerateUV) {
    const Float normal_det = DifferenceOfProducts(duv02[0], duv12[1], duv02[1], duv12[0]);
    const auto dn02=mesh->n[n[0]]-mesh->n[n[2]], dn12=mesh->n[n[1]]-mesh->n[n[2]];
    rec.dndu = (duv12[1] * dn02 - duv02[1] * dn12) / normal_det;
    rec.dndv = (-duv12[0] * dn02 + duv02[0] * dn12) / normal_det;
    if (reverseOrientation) {
      rec.dndu = -rec.dndu;
      rec.dndv = -rec.dndv;
    }
  }
  rec.shape = this;
  rec.texture_object_p = (*WorldToObject)(rec.p);
  rec.texture_object_normal = unit_vector((*WorldToObject)(rec.geometric_normal));
  rec.ComputeDifferentials(r);
  if (bump_tex) {
    const auto ns=rec.physical_shading_normal;
    // PBRT's triangle shading frame uses the interpolated normal and projected
    // u tangent. UV footprints above use geometric tangents, before this frame.
    if (smooth_normals) {
      auto bitangent = cross(convert_to_vec3(ns), rec.dpdu);
      if (bitangent.squared_length() > 0) {
        rec.dpdu = cross(bitangent, convert_to_vec3(ns));
        rec.dpdv = bitangent;
      } else
        CoordinateSystem(convert_to_vec3(ns), &rec.dpdu, &rec.dpdv);
    }
    const TextureFootprint footprint{rec.dudx, rec.dvdx, rec.dudy, rec.dvdy, rec.has_differentials};
    TextureEvalContext context;
    if (bump_tex->height_texture) {
      context = TextureEvalContext::FromHit(rec);
      context.u = uvHit[0]; context.v = uvHit[1];
    }
    rec.bump_normal = bump_tex->perturb(uvHit[0], uvHit[1], rec.p, ns, rec.dpdu, rec.dpdv, rec.dndu,
                                        rec.dndv, footprint, bump_tex->height_texture ? &context : nullptr);
    rec.physical_shading_normal=rec.bump_normal;
    rec.has_bump = true;
  }
  rec.u = mesh->has_vertex_colors ? b0 : uvHit[0];
  rec.v = mesh->has_vertex_colors ? b1 : uvHit[1];
  
  rec.mat_ptr = mesh->mesh_materials[mat_id].get();
  rec.alpha_miss = alpha_miss;
  
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

#ifdef NOT_CRAN
#include <testthat.h>
#include "../volumes/lights.h"

context("Triangle emission orientation") {
  test_that("reversal and mirrored transforms agree across hit paths and shading normals") {
    float vertices[] = {-1,-1,0, 1,-1,0, 0,1,0};
    float normals[] = {.2,0,.98, .2,0,.98, .2,0,.98};
    int indices[] = {0,1,2};
    auto lamp = std::make_shared<diffuse_light>(
        std::make_shared<constant_texture>(point3f(1)), 2, false);
    for (Float sx : {Float(-1), Float(1)}) {
      for (Float sz : {Float(-1), Float(1)}) {
        Transform transform = Scale(sx, 1, sz), inverse = Inverse(transform);
        for (bool reversed : {false, true}) {
          for (bool smooth : {false, true}) {
            TriangleMesh mesh(vertices, indices, smooth ? normals : nullptr, nullptr,
                              3, 3, nullptr, nullptr, lamp, &transform, &inverse, reversed);
            triangle tri(&mesh, mesh.vertexIndices.data(), mesh.normalIndices.data(),
                         mesh.texIndices.data(), 0, &transform, &inverse, reversed);
            const Float sign = reversed ? -1 : 1;
            const normal3f expected = unit_vector(transform(normal3f(0,0,1))) * sign;
            const normal3f shading = smooth ? unit_vector(transform(normal3f(.2,0,.98))) * sign : expected;
            for (bool segment : {false, true}) {
              for (Float side : {Float(-1), Float(1)}) {
                Ray ray(point3f(0,0,side*2), vec3f(0,0,-side));
                ray.segment_absorption = segment;
                random_gen rng(13); RandomSampler sampler(rng);
                hit_record a, b;
                expect_true(tri.hit(ray, 0, 10, a, rng));
                expect_true(tri.hit(ray, 0, 10, b, &sampler));
                for (const auto &hit : {a, b}) {
                  for (int c = 0; c < 3; ++c) {
                    expect_true(std::abs(hit.geometric_normal[c] - expected[c]) < 1e-6);
                    expect_true(std::abs(hit.physical_shading_normal[c] - shading[c]) < 1e-6);
                    expect_true(std::abs(hit.normal[c] - shading[c]) < 1e-6);
                  }
                  bool invisible = false;
                  const Float emitted = lamp->emitted(ray, hit, hit.u, hit.v, hit.p, invisible)[0];
                  expect_true(emitted == (dot(expected, ray.d) < 0 ? 2 : 0));
                }
              }
            }
          }
        }
      }
    }
  }

  test_that("light selection and sampled endpoints use the reversed emission side") {
    float vertices[] = {-1,-1,0, 1,-1,0, 0,1,0};
    int indices[] = {0,1,2};
    auto lamp = std::make_shared<diffuse_light>(
        std::make_shared<constant_texture>(point3f(1)), 2, false);
    for (Float sx : {Float(-1), Float(1)}) {
      Transform transform = Scale(sx,1,1), inverse = Inverse(transform);
      TriangleMesh mesh(vertices, indices, nullptr, nullptr, 3, 3, nullptr, nullptr,
                        lamp, &transform, &inverse, false);
      hitable_list lights;
      for (bool reversed : {false, true})
        lights.add(std::make_shared<triangle>(&mesh, mesh.vertexIndices.data(),
            mesh.normalIndices.data(), mesh.texIndices.data(), 0, &transform, &inverse, reversed));
      VolumeLightSampler selection(lights);
      for (Float side : {Float(-1), Float(1)}) {
        VolumeLightSampler::Context context{point3f(0,0,side*3), normal3f(0), 0};
        const size_t front = side > 0 ? 0 : 1;
        expect_true(selection.SelectionPmf(context, front) > 5 * selection.SelectionPmf(context, 1-front));
        random_gen rng(9); RandomSampler sampler(rng);
        for (int i = 0; i < 16; ++i) {
          auto sample = selection.SampleEmitter(context, &sampler, rng);
          Ray ray(context.p, sample.wi); ray.segment_absorption = true;
          hit_record hit;
          expect_true(selection.Endpoint(sample, ray, hit, rng));
          bool invisible = false;
          const Float emitted = hit.mat_ptr->emitted(ray, hit, hit.u, hit.v, hit.p, invisible)[0];
          expect_true(emitted == (sample.index == front ? 2 : 0));
        }
      }
    }
  }
}
#endif
