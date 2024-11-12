#ifndef CAMERAH
#define CAMERAH

#include "vec3.h"
#include "point3.h"
#include "onbh.h"
#include "animatedtransform.h"
#include "bounds.h"


struct CameraSample {
  point2f pFilm;
  point2f pLens;
  Float time;
  CameraSample(point2f pFilm_, point2f pLens_, Float time_) : 
    pFilm(pFilm_), pLens(pLens_), time(time_) {};
  CameraSample(vec2f pFilm_, vec2f pLens_, Float time_) : 
    pFilm(point2f(pFilm_.x(),pFilm_.y())), pLens(point2f(pLens_.x(),pLens_.y())), time(time_) {};
};

class RayCamera {
  public:
    RayCamera() {};
    virtual ~RayCamera() {};
    virtual ray get_ray(Float s, Float t, point3f u3, Float u) {
      return(ray());
    };
    virtual void update_position(vec3f delta, bool update_uvw, bool update_focal = true) = 0;
    virtual void update_fov(Float delta_fov)  = 0;
    virtual void update_aperture(Float delta_aperture)  = 0;
    virtual void update_focal_distance(Float delta_focus)  = 0;
    virtual void update_look_direction(vec3f dir) = 0;
    virtual void update_lookat(point3f point) = 0;
    virtual void update_position_absolute(point3f point) = 0;
    virtual void update_ortho_absolute(vec2f o_size) = 0;
    virtual void update_aperture_absolute(Float aperture) = 0;
    virtual void update_focal_absolute(Float focal_length) = 0;
    virtual void update_fov_absolute(Float fov) = 0;
    
    virtual void reset()  = 0;
    virtual Float GenerateRay(const CameraSample &sample, ray* ray2) const {
      return(0.0);
    };
    virtual vec3f get_w() = 0;
    virtual vec3f get_u() = 0;
    virtual vec3f get_v() = 0;
    virtual Float get_fov() = 0;
    virtual Float get_aperture() = 0;
    virtual Float get_focal_distance() = 0;
    virtual point3f get_origin() = 0;
    virtual Float get_iso() {return(1.f);}
    virtual vec3f get_up() {return(vec3f(0,1,0));}
    virtual point3f get_lookat() {return(point3f(0,0,0));}
    virtual point2f get_ortho() {return(point2f(1.f,1.f));}
    
    
};

class camera : public RayCamera {
  public:
    camera(point3f lookfrom, point3f _lookat, vec3f _vup, Float vfov, Float aspect, Float aperture, Float _focus_dist,
           Float t0, Float t1);
    ray get_ray(Float s, Float t, point3f u3, Float u);

    void update_position(vec3f delta, bool update_uvw, bool update_focal = true);
    void update_fov(Float delta_fov);
    void update_aperture(Float delta_aperture);
    void update_focal_distance(Float delta_focus);
    void update_look_direction(vec3f dir);
    void update_lookat(point3f point);
    void update_position_absolute(point3f point);
    void update_ortho_absolute(vec2f o_size);
    void update_aperture_absolute(Float aperture);
    void update_focal_absolute(Float focal_length);
    void update_fov_absolute(Float fov_new);
    
    void reset();
    vec3f get_w() {return(w);}
    vec3f get_u() {return(u);}
    vec3f get_v() {return(v);}
    Float get_fov() {return(fov);}
    Float get_aperture() {return(lens_radius * 2);}
    Float get_focal_distance() {return(focus_dist);}

    point3f get_origin() {return(origin);}
    vec3f get_up() {return(vup);}
    point3f get_lookat() {return(lookat);}
    point2f get_ortho() {return(point2f(1.f,1.f));}
    
    Float half_height;
    Float half_width;
    point3f origin;
    point3f lookat;
    Float focus_dist;
    vec3f vup;
    point3f lower_left_corner;
    vec3f horizontal;
    vec3f vertical;
    vec3f u, v, w;
    Float time0, time1;
    Float lens_radius;
    Float start_lens_radius;
    point3f start_origin;
    Float start_focus_dist;
    Float aspect;
    Float fov;
    Float start_fov;
    point3f start_lookat;
    
};

class ortho_camera : public RayCamera {
public:
  ortho_camera(point3f lookfrom, point3f _lookat, vec3f _vup, 
               Float _cam_width, Float _cam_height, 
               Float t0, Float t1);
  ray get_ray(Float s, Float t, point3f u3, Float u);
  void update_position(vec3f delta, bool update_uvw, bool update_focal = true);
  void update_fov(Float delta_fov);
  void update_aperture(Float delta_aperture);
  void update_focal_distance(Float delta_focus);
  void update_look_direction(vec3f dir);
  void update_lookat(point3f point);
  void update_position_absolute(point3f point);
  void update_ortho_absolute(vec2f o_size);
  void update_aperture_absolute(Float aperture);
  void update_focal_absolute(Float focal_length);
  void update_fov_absolute(Float fov_new);
  
  void reset();
  vec3f get_w() {return(w);}
  vec3f get_u() {return(u);}
  vec3f get_v() {return(v);}
  Float get_fov() {return(0);} 
  Float get_aperture() {return(0);} 
  Float get_focal_distance() {return(0);}
  point3f get_origin() {return(origin);}
  vec3f get_up() {return(vup);}
  point3f get_lookat() {return(lookat);}
  point2f get_ortho() {return(point2f(cam_width,cam_height));}
  
  point3f origin;
  point3f lower_left_corner;
  point3f start_origin;
  point3f lookat;
  vec3f vup;
  vec3f horizontal;
  vec3f vertical;
  vec3f u, v, w;
  Float time0, time1;
  Float cam_width, cam_height;
  Float start_cam_width, start_cam_height;
  point3f start_lookat;
  Float focus_dist;
  Float initial_ratio;
};


class environment_camera : public RayCamera {
  public:
    environment_camera(point3f lookfrom, point3f lookat, vec3f _vup, 
                       Float t0, Float t1);
    ray get_ray(Float s, Float t, point3f u3, Float u);
    void update_position(vec3f delta, bool update_uvw, bool update_focal = true);
    void update_fov(Float delta_fov);
    void update_aperture(Float delta_aperture);
    void update_focal_distance(Float delta_focus);
    void update_look_direction(vec3f dir);
    void update_lookat(point3f point);
    void update_position_absolute(point3f point);
    void update_ortho_absolute(vec2f o_size);
    void update_aperture_absolute(Float aperture);
    void update_focal_absolute(Float focal_length);
    void update_fov_absolute(Float fov_new);
    
    void reset();
    vec3f get_w();
    vec3f get_u();
    vec3f get_v();
    Float get_fov() {return(360);}
    Float get_aperture() {return(0);}
    Float get_focal_distance() {return(0);}
    point3f get_origin() {return(origin);}
    vec3f get_up() {return(vup);}
    point3f get_lookat() {return(lookat);}
    point2f get_ortho() {return(point2f(1.f,1.f));}
    
    point3f origin;
    point3f start_origin;
    vec3f u, v, w;
    Float nx, ny;
    Float time0, time1;
    onb uvw;
    vec3f vup;
    point3f lookat;
    point3f start_lookat;
    
};

class RealisticCamera  : public RayCamera {
public:
  // RealisticCamera Public Methods
  RealisticCamera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
                  Float shutterClose, Float apertureDiameter, Float cam_width, Float cam_height,
                  Float focusDistance, bool simpleWeighting,
                  std::vector<Float> &lensData,
                  Float film_size, Float camera_scale, Float _iso,
                  vec3f _camera_up, Transform _CamTransform, 
                  point3f _lookat);
  Float GenerateRay(const CameraSample &sample, ray* ray2) const;
  void update_position(vec3f delta, bool update_uvw, bool update_focal = true);
  void update_fov(Float delta_fov);
  void update_aperture(Float delta_aperture);
  void update_focal_distance(Float delta_focus);
  void update_look_direction(vec3f dir);
  void update_lookat(point3f point);
  void update_position_absolute(point3f point);
  void update_ortho_absolute(vec2f o_size);
  void update_aperture_absolute(Float aperture);
  void update_focal_absolute(Float focal_length);
  void update_fov_absolute(Float fov_new);
  
  void reset();
  vec3f get_w();
  vec3f get_u();
  vec3f get_v();
  Float get_fov() {return(-1);}
  Float get_aperture() {return(0);}
  Float get_focal_distance() {return(focusDistance);}
  point3f get_origin();
  Float get_iso() {return(iso);}
  vec3f get_up() {return(camera_up);}
  point3f get_lookat() {return(lookat);}
  point2f get_ortho() {return(point2f(1.f,1.f));}
  
private:
  // RealisticCamera Private Declarations
  struct LensElementInterface {
    Float curvatureRadius;
    Float thickness;
    Float eta;
    Float apertureRadius;
  };
  
  // RealisticCamera Private Data
  std::vector<LensElementInterface> elementInterfaces;
  std::vector<Bounds2f> exitPupilBounds;
  
  // RealisticCamera Private Methods
  Float LensRearZ() const { 
    return elementInterfaces.back().thickness; 
  }
  Float LensFrontZ() const {
    Float zSum = 0;
    for (const LensElementInterface &element : elementInterfaces)
      zSum += element.thickness;
    return zSum;
  }
  Float RearElementRadius() const {
    return elementInterfaces.back().apertureRadius;
  }
  
  bool TraceLensesFromFilm(const ray &r, ray *rOut) const;
  static bool IntersectSphericalElement(Float radius, Float zCenter,
                                        const ray &ray, Float *t,
                                        normal3f *n);
  bool TraceLensesFromScene(const ray &rCamera, ray* rOut) const;
  static void ComputeCardinalPoints(const ray &rIn, const ray &rOut, Float *p,
                                    Float *f);
  void ComputeThickLensApproximation(Float pz[2], Float f[2]) const;
  Float FocusThickLens(Float focusDistance, bool throw_error = true);
  Bounds2f BoundExitPupil(Float pFilmX0, Float pFilmX1) const;
  point3f SampleExitPupil(const point2f &pFilm, const point2f &lensSample,
                          Float *sampleBoundsArea) const;
  Float FocusBinarySearch(Float focusDistance);
  Float FocusDistance(Float filmDistance);
  Bounds2f GetPhysicalExtent() const;
  
  AnimatedTransform CameraToWorld;
  Transform CameraMovement;
  
  Float shutterOpen;
  Float shutterClose;
  const bool simpleWeighting;
  Float cam_width;
  Float cam_height;
  Float diag;
  Float min_aperture;
  bool init;
  Float iso;
  vec3f camera_up;
  Transform CamTransform;
  point3f origin;
  Float focusDistance;
  Float start_focusDistance;
  point3f start_lookat;
  point3f lookat;
};


  
#endif
