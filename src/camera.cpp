#include "camera.h"
#include "low_discrepancy.h"

camera::camera(point3f lookfrom, point3f _lookat, vec3f _vup, Float vfov, 
               Float _aspect, Float aperture, Float _focus_dist,
               Float t0, Float t1) {
  time0 = t0;
  time1 = t1;
  aspect = _aspect;
  fov = vfov;
  start_fov = fov;
  lens_radius = aperture / 2;
  start_lens_radius = lens_radius;
  Float theta = vfov * static_cast<Float>(M_PI)/180;
  half_height = tan(theta/2);
  half_width = aspect * half_height;
  origin = lookfrom;
  start_origin = origin;
  lookat = _lookat;
  start_lookat = lookat;
  vup = _vup;
  focus_dist = _focus_dist;
  start_focus_dist = focus_dist;
  w = unit_vector(lookfrom - lookat);
  u = unit_vector(cross(vup, w)); //unit_same
  v = cross(w, u);
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
}

ray camera::get_ray(Float s, Float t, point3f u3, Float u1) {
  point3f rd = lens_radius * u3;
  vec3f offset = u * rd.x() + v * rd.y();
  Float time = time0 + u1 * (time1 - time0);
  return(ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset, time)); 
}

void camera::update_position(vec3f delta, bool update_uvw, bool update_focal) {
  origin += delta;
  if(update_uvw) {
    w = unit_vector(origin - lookat);
    if(update_focal) {
      origin += w * (focus_dist - (origin - lookat).length());
    } else {
      focus_dist = (origin - lookat).length();
    }
    u = unit_vector(cross(vup, w));
    v = cross(w, u);
  }
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void camera::update_fov(Float delta_fov) {
  fov += delta_fov;
  fov = std::fmax(std::fmin(fov,179.9f),0.1f);
  Float theta = fov * static_cast<Float>(M_PI)/180;
  half_height = tan(theta/2);
  half_width = aspect * half_height;
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
}

void camera::update_aperture(Float delta_aperture) {
  lens_radius += delta_aperture / 2;
  lens_radius = std::fmax(lens_radius,0);
}

void camera::update_focal_distance(Float delta_focus) {
  focus_dist += delta_focus;
  focus_dist = std::fmax(0.001,focus_dist);
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
}

void camera::update_look_direction(vec3f dir) {
  w = unit_vector(dir);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void camera::update_lookat(point3f point) {
  lookat = point;
  w = unit_vector(origin-lookat);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void camera::update_position_absolute(point3f point) {
  origin = point;
  w = unit_vector(origin - lookat);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void camera::update_ortho_absolute(vec2f o_size) {
}

void camera::update_aperture_absolute(Float aperture) {
  lens_radius = aperture;
}

void camera::update_focal_absolute(Float focal_length) {
  focus_dist = focal_length;
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void camera::update_fov_absolute(Float fov_new) {
  fov = fov_new;
  fov = std::fmax(std::fmin(fov,179.9f),0.1f);
  Float theta = fov * static_cast<Float>(M_PI)/180;
  half_height = tan(theta/2);
  half_width = aspect * half_height;
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
}



void camera::reset() {
  origin = start_origin;
  lookat = start_lookat;
  focus_dist = start_focus_dist;
  fov = start_fov;
  Float theta = start_fov * static_cast<Float>(M_PI)/180;
  half_height = tan(theta/2);
  half_width = aspect * half_height;
  lens_radius = start_lens_radius;
  w = unit_vector(origin - lookat);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - half_width * focus_dist *  u - half_height * focus_dist * v - focus_dist * w;
  horizontal = 2.0f * half_width * focus_dist * u;
  vertical = 2.0f * half_height * focus_dist * v;
}


ortho_camera::ortho_camera(point3f lookfrom, point3f _lookat, vec3f _vup, 
             Float _cam_width, Float _cam_height, 
             Float t0, Float t1) {
  time0 = t0;
  time1 = t1;
  origin = lookfrom;
  start_origin = origin;
  cam_width = _cam_width;
  cam_height = _cam_height;
  start_cam_width= cam_width;
  start_cam_height= cam_height;
  lookat = _lookat;
  start_lookat = lookat;
  vup = _vup;
  focus_dist = (lookfrom-lookat).length();
  w = unit_vector(lookfrom - lookat);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
  initial_ratio = cam_width / cam_height;
}

ray ortho_camera::get_ray(Float s, Float t, point3f u3, Float u) {
  Float time = time0 + u * (time1 - time0);
  return(ray(lower_left_corner + s * horizontal + t * vertical, -w, time)); 
}

void ortho_camera::update_position(vec3f delta, bool update_uvw, bool update_focal) {
  origin += delta;
  if(update_uvw) {
    w = unit_vector(origin - lookat);
    if(update_focal) {
      origin += w * (focus_dist - (origin - lookat).length());
    } 
    u = unit_vector(cross(vup, w));
    v = cross(w, u);
  }
  focus_dist = (origin - lookat).length();
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void ortho_camera::update_fov(Float delta_fov) {
  cam_width += delta_fov * initial_ratio;
  cam_height += delta_fov;
  cam_width = std::fmax(cam_width,0.001);
  cam_height = std::fmax(cam_height,0.001);
  
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
}

void ortho_camera::update_aperture(Float delta_aperture) {
  Rprintf("Cannot change aperture while using an orthographic camera.\n");
}

void ortho_camera::update_focal_distance(Float delta_focus) {
  Rprintf("Cannot change focal distance while using an orthographic camera.\n");
}

void ortho_camera::update_look_direction(vec3f dir) {
  w = unit_vector(dir);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void ortho_camera::update_lookat(point3f point) {
  lookat = point;
  focus_dist = (origin-lookat).length();
}

void ortho_camera::update_position_absolute(point3f point) {
  origin = point;
  w = unit_vector(origin - lookat);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void ortho_camera::update_ortho_absolute(vec2f o_size) {
  cam_width = std::fmax(o_size.x(),0.001);
  cam_height = std::fmax(o_size.y(),0.001);
  
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
}

void ortho_camera::update_aperture_absolute(Float aperture) {
}

void ortho_camera::update_focal_absolute(Float focal_length) {
  focus_dist = focal_length;
}

void ortho_camera::update_fov_absolute(Float fov_new) {
}



void ortho_camera::reset() {
  origin = start_origin;
  lookat = start_lookat;
  cam_width = start_cam_width;
  cam_height = start_cam_height;
  w = unit_vector(origin - lookat);
  u = unit_vector(cross(vup, w));
  v = cross(w, u);
  lower_left_corner = origin - cam_width/2 *  u - cam_height/2 * v;
  horizontal = cam_width * u;
  vertical = cam_height * v;
}


environment_camera::environment_camera(point3f lookfrom, point3f lookat, vec3f _vup, 
                   Float t0, Float t1) {
  time0 = t0;
  time1 = t1;
  origin = lookfrom;
  start_origin = lookfrom;
  start_lookat = lookat;
  vup = _vup;
  w = unit_vector(lookfrom - lookat);
  v = unit_vector(-cross(vup, w));
  u = cross(w, v);
  uvw = onb(w,v,u);
}


ray environment_camera::get_ray(Float s, Float t, point3f u3, Float u1) {
  Float time = time0 + u1 * (time1 - time0);
  Float theta = static_cast<Float>(M_PI) * t;
  Float phi = 2 * static_cast<Float>(M_PI) * s;
  vec3f dir(std::sin(theta) * std::cos(phi), 
            std::sin(theta) * std::sin(phi),
            std::cos(theta));
  dir = uvw.local_to_world(dir);
  return(ray(origin, dir, time)); 
}

void environment_camera::update_position(vec3f delta, bool update_uvw, bool update_focal) {
  Float focus_dist = (origin - lookat).length();
  origin += delta;
  if(update_uvw) {
    w = unit_vector(origin - lookat);
    if(update_focal) {
      origin += w * (focus_dist - (origin - lookat).length());
    } else {
      focus_dist = (origin - lookat).length();
    }
    v = unit_vector(-cross(vup, w));
    u = cross(w, v);
    uvw = onb(w,v,u);
  }
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void environment_camera::update_fov(Float delta_fov) {
  Rprintf("Cannot change FOV while using an environment camera.\n");
}

void environment_camera::update_aperture(Float delta_aperture) {
  Rprintf("Cannot change aperture while using an environment camera.\n");
}

void environment_camera::update_focal_distance(Float delta_focus) {
  Rprintf("Cannot change focal distance while using an environment camera.\n");
}

void environment_camera::update_look_direction(vec3f dir) {
  w = unit_vector(dir);
  v = unit_vector(-cross(vup, w));
  u = cross(w, v);
  uvw = onb(w,v,u);
}

void environment_camera::update_lookat(point3f point) {
  lookat = point;
}

void environment_camera::update_position_absolute(point3f point) {
  origin = point;
  w = unit_vector(origin - lookat);
  v = unit_vector(-cross(vup, w));
  u = cross(w, v);
  uvw = onb(w,v,u);
  if(w.length() == 0 && u.length() == 0) {
    reset();
  } 
}

void environment_camera::update_ortho_absolute(vec2f o_size) {
}

void environment_camera::update_aperture_absolute(Float aperture) {
}

void environment_camera::update_focal_absolute(Float focal_length) {
}

void environment_camera::update_fov_absolute(Float fov_new) {
}

vec3f environment_camera::get_w() {return(w);}
vec3f environment_camera::get_u() {return(-v);}
vec3f environment_camera::get_v() {return(-u);}

void environment_camera::reset() {
  origin = start_origin;
  lookat = start_lookat;
  w = unit_vector(origin - lookat);
  v = unit_vector(-cross(vup, w));
  u = cross(w, v);
  uvw = onb(w,v,u);
}

RealisticCamera::RealisticCamera(const AnimatedTransform &CameraToWorld,
                                 Float shutterOpen, Float shutterClose, 
                                 Float apertureDiameter, 
                                 Float cam_width, Float cam_height,
                                 Float _focusDistance,
                                 bool simpleWeighting,
                                 std::vector<Float> &lensData,
                                 Float film_size, Float camera_scale,
                                 Float _iso,
                                 vec3f _camera_up, Transform _CamTransform, 
                                 point3f _lookat)
  : CameraToWorld(CameraToWorld), 
    shutterOpen(shutterOpen), shutterClose(shutterClose), 
    simpleWeighting(simpleWeighting), cam_width(cam_width), cam_height(cam_height),
    diag(film_size * camera_scale), iso(_iso), camera_up(_camera_up), CamTransform(_CamTransform),
    focusDistance(_focusDistance), start_focusDistance(_focusDistance), start_lookat(_lookat),
    lookat(_lookat) {
  CameraMovement = Transform(Matrix4x4(1,0,0,0,
                                       0,1,0,0,
                                       0,0,1,0,
                                       0,0,0,1));
  
  if(lensData.size() == 0) {
    init = false;
  } else {
    init = true;
  }

  if(init) {
    for (int i = 0; i < (int)lensData.size(); i += 4) {
      if (lensData[i] == 0) {
        if (apertureDiameter > lensData[i + 3]*camera_scale) {
          Rcpp::Rcout << "Specified aperture diameter is greater than maximum possible.  Clamping it.\n";
        } else {
          lensData[i + 3] = apertureDiameter;
        }
      }
      elementInterfaces.push_back(LensElementInterface(
      {lensData[i] * (Float).001 * camera_scale, lensData[i + 1] * (Float).001 * camera_scale,
       lensData[i + 2], lensData[i + 3] * Float(.001) / Float(2.) * camera_scale}));
    }
    min_aperture = elementInterfaces[0].apertureRadius * camera_scale;
    for(unsigned int i = 1; i < elementInterfaces.size(); i++) {
      min_aperture = elementInterfaces[i].apertureRadius < min_aperture ? elementInterfaces[i].apertureRadius : min_aperture;
    }
    elementInterfaces.back().thickness = FocusThickLens(focusDistance);
    
// #ifdef REALDEBUG
//     std::ofstream myfile;
//     myfile.open("cam_info_debug.txt", std::ios::app | std::ios::out);
//     
//     for (int i = 0; i < elementInterfaces.size(); i++) {
//       const LensElementInterface &element = elementInterfaces[i];
//       myfile      << element.curvatureRadius << "\t"
//                   << element.thickness <<  "\t" 
//                   << element.eta << "\t"
//                   << element.apertureRadius*2 <<"\n";
//     }
//     myfile.close();
//     
// #endif
    // Compute exit pupil bounds at sampled points on the film
    int nSamples = 64;
    exitPupilBounds.resize(nSamples);
    // ParallelFor([&](int i) {
    for(int i = 0; i < nSamples; i++) {
      Float r0 = (Float)i / nSamples * diag / 2;
      Float r1 = (Float)(i + 1) / nSamples * diag / 2;
      exitPupilBounds[i] = BoundExitPupil(r0, r1);
    }
  }
  // }, nSamples);
  
  // if (simpleWeighting) {
  //   Warning("\"simpleweighting\" option with RealisticCamera no longer "
  //             "necessarily matches regular camera images. Further, pixel "
  //             "values will vary a bit depending on the aperture size. See "
  //             "this discussion for details: "
  //             "https://github.com/mmp/pbrt-v3/issues/162#issuecomment-348625837");
  // }
}

Float RealisticCamera::FocusDistance(Float filmDistance) {
  // Find offset ray from film center through lens
  Bounds2f bounds = BoundExitPupil(0, min_aperture/10);
  
  const std::array<Float, 3> scaleFactors = {0.1f, 0.01f, 0.001f};
  Float lu = 0.0f;
  
  ray ray2;
  
  // Try some different and decreasing scaling factor to find focus ray
  // more quickly when `aperturediameter` is too small.
  // (e.g. 2 [mm] for `aperturediameter` with wide.22mm.dat),
  bool foundFocusRay = false;
  for (Float scale : scaleFactors) {
    lu = scale * bounds.pMax[0];
    if (TraceLensesFromFilm(ray(point3f(0, 0, LensRearZ() - filmDistance),
                                vec3f(lu, 0, filmDistance)),
                                &ray2)) {
      foundFocusRay = true;
      break;
    }
  }
  
  if (!foundFocusRay) {
    Rcpp::Rcout << "Focus ray at lens pos(" << lu << ",0) didn't make it through the lenses at distance " <<  filmDistance << "\n";
    return Infinity;
  }
  
  // Compute distance _zFocus_ where ray intersects the principal axis
  Float tFocus = -ray2.origin().x() / ray2.direction().x();
  Float zFocus = ray2(tFocus).z();
  if (zFocus < 0) zFocus = Infinity;
  return zFocus;
}

Float RealisticCamera::FocusBinarySearch(Float focusDistance) {
  Float filmDistanceLower, filmDistanceUpper;
  // Find _filmDistanceLower_, _filmDistanceUpper_ that bound focus distance
  filmDistanceLower = filmDistanceUpper = FocusThickLens(focusDistance);
  while (FocusDistance(filmDistanceLower) > focusDistance)
    filmDistanceLower *= 1.005f;
  while (FocusDistance(filmDistanceUpper) < focusDistance)
    filmDistanceUpper /= 1.005f;
  
  // Do binary search on film distances to focus
  for (int i = 0; i < 20; ++i) {
    Float fmid = 0.5f * (filmDistanceLower + filmDistanceUpper);
    Float midFocus = FocusDistance(fmid);
    if (midFocus < focusDistance)
      filmDistanceLower = fmid;
    else
      filmDistanceUpper = fmid;
  }
  return 0.5f * (filmDistanceLower + filmDistanceUpper);
}

inline bool Refract(const vec3f &wi, normal3f &n, Float eta, vec3f *wt) {
  Float cosTheta_i = dot(n, wi);
  // Potentially flip interface orientation for Snell's law
  if (cosTheta_i < 0) {
    eta = 1 / eta;
    cosTheta_i = -cosTheta_i;
    n = -n;
  }
  
  // Compute $\cos\,\theta_\roman{t}$ using Snell's law
  Float sin2Theta_i = std::max<Float>(0, 1 - Sqr(cosTheta_i));
  Float sin2Theta_t = sin2Theta_i / Sqr(eta);
  // Handle total internal reflection case
  if (sin2Theta_t >= 1)
    return false;
  
  Float cosTheta_t = SafeSqrt(1 - sin2Theta_t);
  
  *wt = -wi / eta + (cosTheta_i / eta - cosTheta_t) * vec3f(n.x(),n.y(),n.z());
  return(true);
}


inline bool Quadratic(Float A, Float B, Float C, Float *t0, Float *t1) {
  if (A == 0) {
    if (B == 0) {
      return false;
    }
    *t0 = *t1 = -C / B;
    return true;
  }
  // Find quadratic discriminant
  Float discrim = DifferenceOfProducts(B, B, 4*A, C);
  if (discrim < 0.) return false;
  Float rootDiscrim = std::sqrt(discrim);
  
  // Compute quadratic _t_ values
  Float q;
  if ((float)B < 0) {
    q = -.5 * (B - rootDiscrim);
  } else {
    q = -.5 * (B + rootDiscrim);
  }
  *t0 = q / A;
  *t1 = C / q;
  if ((float)*t0 > (float)*t1) {
    std::swap(*t0, *t1);
  }
  return true;
}


bool RealisticCamera::IntersectSphericalElement(Float radius,Float zCenter, 
                                                const ray &ray2, Float *t, normal3f *n) {
  const vec3f& rd = ray2.direction();
  const point3f& ro = ray2.origin();
  
  point3f o = ro - vec3f(0, 0, zCenter);
  Float A = rd.x()*rd.x() + rd.y()*rd.y() + rd.z()*rd.z();
  Float B = 2 * (rd.x()*o.x() + rd.y()*o.y() + rd.z()*o.z());
  Float C = o.x()*o.x() + o.y()*o.y() + o.z()*o.z() - radius*radius;
  Float t0, t1;
  if (!Quadratic(A, B, C, &t0, &t1)) {
    return false;
  }
  bool useCloserT = (rd.z() > 0) ^ (radius < 0);
  *t = useCloserT ? std::fmin(t0, t1) : std::fmax(t0, t1);
  if (*t < 0) {
    return false;
  }
  *n = Faceforward(unit_vector(convert_to_normal3(o + *t * rd)), -rd);
  return true;
}

bool RealisticCamera::TraceLensesFromScene(const ray &rCamera,
                                           ray* rOut) const {
  Float elementZ = -LensFrontZ();
  // Transform _rCamera_ from camera to lens system space
  
  static const Transform CameraToLens = Scale(1, 1, -1);
  ray rLens = CameraToLens(rCamera);
  for (size_t i = 0; i < elementInterfaces.size(); ++i) {
    const LensElementInterface &element = elementInterfaces[i];
    // Compute intersection of ray with lens element
    Float t;
    normal3f n;
    bool isStop = (element.curvatureRadius == 0);
    if (isStop) {
      t = (elementZ - rLens.origin().z()) / rLens.direction().z();
      
      //Debug
    } else {
      Float radius = element.curvatureRadius;
      Float zCenter = elementZ + element.curvatureRadius;

      if (!IntersectSphericalElement(radius, zCenter, rLens, &t, &n)) {
        return false;
      }
    }
    // Test intersection point against element aperture
    point3f pHit = rLens(t);
    Float r2 = pHit.x() * pHit.x() + pHit.y() * pHit.y();
    if (r2 > element.apertureRadius * element.apertureRadius) {
      return false;
    }
    rLens.A = pHit;
    
    // Update ray path for from-scene element interface interaction
    if (!isStop) {
      vec3f wt;
      Float etaI = (i == 0 || elementInterfaces[i - 1].eta == 0)
        ? 1
      : elementInterfaces[i - 1].eta;
      Float etaT =
        (elementInterfaces[i].eta != 0) ? elementInterfaces[i].eta : 1;
      if (!Refract(unit_vector(-rLens.direction()), n, etaT / etaI, &wt)) {
        return false;
      }
      rLens.B = wt;
    }
    elementZ += element.thickness;
  }
  // Transform _rLens_ from lens system space back to camera space
  if (rOut) {
    static const Transform LensToCamera = Scale(1, 1, -1);
    *rOut = LensToCamera(rLens);
  }
  
  return true;
}

bool RealisticCamera::TraceLensesFromFilm(const ray &rCamera, ray *rOut) const {
  Float elementZ = 0;
  // Transform _rCamera_ from camera to lens system space
  static const Transform CameraToLens = Scale(1, 1, -1);
  ray rLens = CameraToLens(rCamera);
  
  for (int i = elementInterfaces.size() - 1; i >= 0; --i) {
    const LensElementInterface &element = elementInterfaces[i];
    // Update ray from film accounting for interaction with _element_
    elementZ -= element.thickness;
    
    // Compute intersection of ray with lens element
    Float t ;
    normal3f n;
    bool isStop = (element.curvatureRadius == 0);
    if (isStop) {
      // The refracted ray computed in the previous lens element
      // interface may be pointed towards film plane(+z) in some
      // extreme situations; in such cases, 't' becomes negative.
      if (rLens.direction().z() >= 0.0) {
        return false;
      }
      t = (elementZ - rLens.origin().z()) / rLens.direction().z();
    } else {
      Float radius = element.curvatureRadius;
      Float zCenter = elementZ + element.curvatureRadius;
      if (!IntersectSphericalElement(radius, zCenter, rLens, &t, &n)) {
        return false;
      }
    }
    
    // Test intersection point against element aperture
    point3f pHit = rLens(t);
    
    Float r2 = pHit.x() * pHit.x() + pHit.y() * pHit.y();
    if (r2 > element.apertureRadius * element.apertureRadius) {
      return false;
    }
    rLens.A = pHit;
    
    // Update ray path for element interface interaction
    if (!isStop) {
      vec3f w;
      Float etaI = element.eta;
      Float etaT = (i > 0 && elementInterfaces[i - 1].eta != 0)
        ? elementInterfaces[i - 1].eta
        : 1;
      if (!Refract(unit_vector(-rLens.direction()), n, etaT / etaI, &w)) {
        return false;
      }
      rLens.B = w;
    }
  }
  rLens = ray(rLens.origin(), rLens.direction());
  // Transform _rLens_ from lens system space back to camera space
  if (rOut) {
    static const Transform LensToCamera = Scale(1, 1, -1);
    *rOut = LensToCamera(rLens);
  }
  return true;
}


void RealisticCamera::ComputeCardinalPoints(const ray &rIn,
                                            const ray &rOut, Float *pz, Float *fz) {
  Float tf = -rOut.origin().x() / rOut.direction().x();
  *fz = -rOut(tf).z();
  Float tp = (rIn.origin().x() - rOut.origin().x()) / rOut.direction().x();
  *pz = -rOut(tp).z();
}

//This is all working correctly (checked with PBRT)
void RealisticCamera::ComputeThickLensApproximation(Float pz[2],
                                                    Float fz[2]) const {
  Float x =  min_aperture/10;// * film->diagonal;
  ray rScene(point3f(x, 0, LensFrontZ() + 1), vec3f(0, 0, -1));
  ray rFilm;
  if(!TraceLensesFromScene(rScene, &rFilm)) {
    throw std::runtime_error("Unable to trace ray from scene to film for thick lens approximation. Is aperture stop extremely small?");
  }
  ComputeCardinalPoints(rScene, rFilm, &pz[0], &fz[0]);
  rFilm = ray(point3f(x, 0, LensRearZ() - 1), vec3f(0, 0, 1));
  if(!TraceLensesFromFilm(rFilm, &rScene)) {
    throw std::runtime_error("Unable to trace ray from film to scene for thick lens approximation. Is aperture stop extremely small?");
  }
  ComputeCardinalPoints(rFilm, rScene, &pz[1], &fz[1]);
}


//I believe this is working correctly 
Float RealisticCamera::FocusThickLens(Float focusDistance, bool throw_error) {
  Float pz[2], fz[2];
  ComputeThickLensApproximation(pz, fz);
  // LOG(INFO) << StringPrintf("Cardinal points: p' = %f f' = %f, p = %f f = %f.\n",
      // pz[0], fz[0], pz[1], fz[1]);
  // LOG(INFO) << StringPrintf("Effective focal length %f\n", fz[0] - pz[0]);
  // Compute translation of lens, _delta_, to focus at _focusDistance_
  Float f = fz[0] - pz[0];
  // Rcpp::Rcout << " Effective focal length " << fz[0] - pz[0] << "\n";
  Float z = -focusDistance;
  Float c = (pz[1] - z - pz[0]) * (pz[1] - z - 4 * f - pz[0]);
  if(c < 0) {
    if(throw_error) {
      throw std::runtime_error("Coefficient must be positive. It looks focusDistance is too short for a given lenses configuration");
    } else {
      return(0.f);
    }
  }
  Float delta = 0.5f * (pz[1] - z + pz[0] - std::sqrt(c));
  return elementInterfaces.back().thickness + delta;
}

Bounds2f RealisticCamera::BoundExitPupil(Float pFilmX0, Float pFilmX1) const {
  Bounds2f pupilBounds;
  // Sample a collection of points on the rear lens to find exit pupil
  const int nSamples = 1024 * 1024;

  int nExitingRays = 0;
  
  // Compute bounding box of projection of rear element on sampling plane
  Float rearRadius = RearElementRadius();
  Bounds2f projRearBounds(point2f(-1.5f * rearRadius, -1.5f * rearRadius),
                          point2f(1.5f * rearRadius, 1.5f * rearRadius));
  // Rcpp::Rcout << projRearBounds << "\n";
  
  for (int i = 0; i < nSamples; ++i) {
    // Find location of sample points on $x$ segment and rear lens element
    point3f pFilm(lerp(Float(i + 0.5f) / (Float)nSamples, pFilmX0, pFilmX1), 0, 0); 
    Float u[2] = {spacefillr::RadicalInverse(0, i), spacefillr::RadicalInverse(1, i)};
    point3f pRear(lerp(u[0], projRearBounds.pMin.x(), projRearBounds.pMax.x()),
                  lerp(u[1], projRearBounds.pMin.y(), projRearBounds.pMax.y()),
                  LensRearZ());
    
    // Expand pupil bounds if ray makes it through the lens system
    if (Inside(point2f(pRear.x(), pRear.y()), pupilBounds) || 
        TraceLensesFromFilm(ray(pFilm, pRear - pFilm), nullptr)) {
      pupilBounds = UnionB(pupilBounds, point2f(pRear.x(), pRear.y()));
      ++nExitingRays;
    }
  }
  
  // Return entire element bounds if no rays made it through the lens system
  if (nExitingRays == 0) {
    // Rcpp::Rcout << "Unable to find exit pupil in x = [" << pFilmX0 << " , "  <<pFilmX1 <<  "]  on film.\n";
    return projRearBounds;
  }
  
  // Expand bounds to account for sample spacing
  pupilBounds = Expand(pupilBounds, 2 * projRearBounds.Diagonal().length() / std::sqrt(nSamples));
  return pupilBounds;
}

point3f RealisticCamera::SampleExitPupil(const point2f &pFilm,
                                         const point2f &lensSample,
                                         Float *sampleBoundsArea) const {
  // Find exit pupil bound for sample distance from film center
  Float rFilm = std::sqrt(pFilm.x() * pFilm.x() + pFilm.y() * pFilm.y());
  int rIndex = rFilm / (diag / 2) * exitPupilBounds.size();

  rIndex = std::min((int)exitPupilBounds.size() - 1, rIndex);
  Bounds2f pupilBounds = exitPupilBounds[rIndex];

  if (sampleBoundsArea) {
    *sampleBoundsArea = pupilBounds.Area();
  }
  
  // Generate sample point inside exit pupil bound
  point2f pLens = pupilBounds.Lerp(lensSample);
  
  // Return sample point rotated by angle of _pFilm_ with $+x$ axis
  Float sinTheta = (rFilm != 0) ? pFilm.y() / rFilm : 0;
  Float cosTheta = (rFilm != 0) ? pFilm.x() / rFilm : 1;
  return point3f(cosTheta * pLens.x() - sinTheta * pLens.y(),
                 sinTheta * pLens.x() + cosTheta * pLens.y(), LensRearZ());
}


Bounds2f  RealisticCamera::GetPhysicalExtent() const {
  Float aspect = (Float)cam_height / (Float)cam_width;
  Float x = std::sqrt(diag * diag / (1 + aspect * aspect));
  Float y = aspect * x;
  return Bounds2f(point2f(-x / 2, -y / 2), point2f(x / 2, y / 2));
}


Float RealisticCamera::GenerateRay(const CameraSample &sample, ray *ray2) const {
  // Find point on film, _pFilm_, corresponding to _sample.pFilm_
  point2f pFilm2 = GetPhysicalExtent().Lerp(sample.pFilm);
  point3f pFilm(-pFilm2.x(), pFilm2.y(), 0);

  // Trace ray from _pFilm_ through lens system
  Float exitPupilBoundsArea;
  point3f pRear = SampleExitPupil(point2f(pFilm.x(), pFilm.y()), sample.pLens,
                                  &exitPupilBoundsArea);

  ray rFilm(pFilm, unit_vector(pRear - pFilm), 
            lerp(sample.time, shutterOpen, shutterClose));
  if (!TraceLensesFromFilm(rFilm, ray2)) {
    return 0;
  }
  // Finish initialization of _RealisticCamera_ ray
  *ray2 = (CamTransform(*ray2));
  ray2->B = unit_vector(ray2->direction());
  
  // Return weighting for _RealisticCamera_ ray
  Float cosTheta = unit_vector(rFilm.direction()).z();
  Float cos4Theta = (cosTheta * cosTheta) * (cosTheta * cosTheta);
  if (simpleWeighting) {
    return cos4Theta * exitPupilBoundsArea / exitPupilBounds[0].Area();
  } else {
    return (shutterClose - shutterOpen) *
      (cos4Theta * exitPupilBoundsArea) / (LensRearZ() * LensRearZ());
  }
  return(0);
}

point3f RealisticCamera::get_origin() {
  return(CamTransform(point3f(0.f)));
}

void RealisticCamera::update_position(vec3f delta, bool update_uvw, bool update_focal) {
  if(update_uvw) {
    point3f old_lookfrom =  CamTransform(point3f(0.f));
    try{
      CamTransform = Inverse(LookAt(old_lookfrom + delta, lookat, camera_up));
    } catch (const std::runtime_error& error) {
      Rprintf("Can't move there, singular inverse matrix--moving back to previous position.\n");
      CamTransform = Inverse(LookAt(old_lookfrom, lookat, camera_up));
    }
  } else {
    CamTransform = Translate(delta) * CamTransform;
  }
}

void RealisticCamera::update_fov(Float delta_fov) {
  Rprintf("Cannot adjust FOV using a realistic camera.\n");
}

void RealisticCamera::update_aperture(Float delta_aperture) {
  Rprintf("Cannot interactively adjust aperture using a realistic camera.\n");
}

void RealisticCamera::update_focal_distance(Float delta_focus) {
  Float new_fd = FocusThickLens(focusDistance + delta_focus, false);
  if(new_fd > 0) {
    focusDistance = focusDistance + delta_focus;
    elementInterfaces.back().thickness = new_fd;
  } else {
    Rprintf("Cannot focus to distance %.2f; Maintaining focal distance of %.2f\n",focusDistance + delta_focus,focusDistance);
  }
}

//This actually takes the position of the intersection
void RealisticCamera::update_look_direction(vec3f dir) {
  
  CamTransform = Transform(LookAt(get_origin(), convert_to_point3(dir), camera_up).GetInverseMatrix());
}

void RealisticCamera::update_lookat(point3f point) {
  lookat = point;
  CamTransform = Transform(LookAt(get_origin(), point, camera_up).GetInverseMatrix());
}

void RealisticCamera::update_position_absolute(point3f point) {
}

void RealisticCamera::update_ortho_absolute(vec2f o_size) {
}

void RealisticCamera::update_aperture_absolute(Float aperture) {
}

void RealisticCamera::update_focal_absolute(Float focal_length) {
}
void RealisticCamera::update_fov_absolute(Float fov_new) {
}

void RealisticCamera::reset() {
  CamTransform = (*CameraToWorld.GetStartTransform());
  focusDistance = start_focusDistance;
  lookat = start_lookat;
  elementInterfaces.back().thickness =  FocusThickLens(focusDistance);
}


vec3f RealisticCamera::get_w() {return(-CamTransform.w());}
vec3f RealisticCamera::get_u() {return(-CamTransform.u());}
vec3f RealisticCamera::get_v() {return(CamTransform.v());}


