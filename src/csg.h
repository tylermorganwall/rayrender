#ifndef CSGH
#define CSGH

#include "hitable.h"
#include "material.h"
#include "mathinline.h"
#include <algorithm>
#include <cmath>
#include "point3.h"

class ImplicitShape { 
  public: 
    virtual Float getDistance(const point3f& from) const = 0; 
    virtual bool bbox(Float t0, Float t1, aabb& box) const = 0; 
    virtual ~ImplicitShape() {} 
}; 

class csg_sphere : public ImplicitShape {
  public: 
    csg_sphere(const vec3f& c, const float& r) : 
      center(c), radius(r) {} 
    Float getDistance(const point3f& from) const { 
      return((from - center).length() - radius); 
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = aabb(center - vec3f(radius,radius,radius), center + vec3f(radius,radius,radius));
      return(true);
    }
    vec3f center;
    float radius; 
}; 

class csg_plane : public ImplicitShape { 
  public: 
    csg_plane(const vec3f nn, const point3f pp, Float width_x, Float width_z) : 
      n(nn), pointOnPlane(pp), width_x(width_x), width_z(width_z) {
      axis.build_from_w(n);
      axis.swap_yz();
    } 
    Float getDistance(const point3f& from) const {
      Float dist   = dot(axis.v(),from - pointOnPlane);
      Float dist_x = dot(axis.u(),from - pointOnPlane);
      Float dist_z = dot(axis.w(),from - pointOnPlane);
      Float diff_x = std::fabs(dist_x) - width_x/2;
      Float diff_z = std::fabs(dist_z) - width_z/2;
      dist = diff_x > 0 ? std::sqrt(diff_x * diff_x + dist * dist) : dist;
      dist = diff_z > 0 ? std::sqrt(diff_z * diff_z + dist * dist) : dist;
      return(dist); 
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = aabb(pointOnPlane - vec3f(width_x,0.1,width_z),
                 pointOnPlane + vec3f(width_x,0.1,width_z));
      return(true);
    }
    vec3f n;
    point3f pointOnPlane;
    Float width_x, width_z;
    onb axis;
}; 

class csg_box : public ImplicitShape { 
  public: 
    csg_box(const point3f& c, point3f width_) :  center(c), width(width_) {} 
    Float getDistance(const point3f& from_old) const {
      point3f from = from_old + -center;
      point3f q = Abs(from) + -width/2;
      const static point3f zeros(0,0,0);
      return(Max(q, zeros).length() + 
             fmin(MaxComponent(q),0.0)); 
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = aabb(center + -Abs(width)/2, center + Abs(width)/2);
      return(true);
    }
    point3f center;
    point3f width; 
}; 

class csg_rounded_box : public ImplicitShape { 
  public: 
    csg_rounded_box(const vec3f& c, vec3f width_, Float radius) :  center(c), width(width_), radius(radius){} 
    Float getDistance(const point3f& from_old) const {
      vec3f from = from_old  - center;
      vec3f q = Abs(from)  - width/2;
      const static vec3f zeros(0,0,0);
      return(Max(q, zeros).length() + 
             fmin(MaxComponent(q),0.0)-radius); 
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = aabb(center-Abs(width)/2, center+Abs(width)/2);
      return(true);
    }
    point3f center;
    vec3f width; 
    Float radius;
}; 

class csg_list : public ImplicitShape { 
  public: 
    csg_list(std::vector<std::shared_ptr<ImplicitShape> > shapes) : shapes(shapes) {} 
    Float getDistance(const point3f& from) const {
      float min_dist = INFINITY;
      float temp;
      for(const auto& shape: shapes) {
        temp = shape->getDistance(from);
        min_dist = temp < min_dist ? temp : min_dist;
      }
      return(min_dist);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      const static vec3f zeros(0,0,0);
      aabb temp1(zeros,zeros), temp2(zeros,zeros);
      for(const auto& shape: shapes) {
        shape->bbox(t0,t1,temp2);
        temp1 = surrounding_box(temp1, temp2);
      }
      box = temp1;
      return(true);
    }
    std::vector<std::shared_ptr<ImplicitShape> > shapes;
}; 

class csg_torus : public ImplicitShape { 
  public: 
    csg_torus(const point3f& c, float ring_radius, float cross_radius) : center(c), 
      ring_radius(ring_radius), cross_radius(cross_radius) {} 
    Float getDistance(const point3f& from_old) const {
      vec3f from = from_old - center;
      vec2f q = vec2f(std::sqrt(from.x()*from.x() + from.z()*from.z()) - ring_radius, from.y());
      return(q.length()-cross_radius);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = aabb(center-vec3f(ring_radius+cross_radius,cross_radius,ring_radius+cross_radius), 
                 center+vec3f(ring_radius+cross_radius,cross_radius,ring_radius+cross_radius));
      return(true);
    }
    point3f center;
    float ring_radius;
    float cross_radius;
    
}; 

class csg_capsule : public ImplicitShape { 
  public: 
    csg_capsule(point3f start, point3f end, Float radius) :  start(start), end(end), radius(radius){} 
    Float getDistance(const point3f& from) const {
      vec3f pa = from - start; 
      vec3f ba = end - start;
      float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
      return((pa - ba*h).length() - radius);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      vec3f min = vec3f(ffmin(start.x(),end.x()),ffmin(start.y(),end.y()),ffmin(start.z(),end.z()));
      vec3f max = vec3f(ffmax(start.x(),end.x()),ffmax(start.y(),end.y()),ffmax(start.z(),end.z()));
      box = aabb(min-radius, max+radius);
      return(true);
    }
    point3f start, end; 
    Float radius;
}; 

class csg_cylinder : public ImplicitShape { 
  public: 
    csg_cylinder(point3f start, point3f end, Float radius, Float corner_radius) :  start(start), end(end), 
      radius(radius), corner_radius(corner_radius) {
      ba = end - start;
      baba = dot(ba,ba);
      inv_baba = 1.0/baba;
    } 
    Float getDistance(const point3f& from) const {
      vec3f pa = from - start; 
      float paba = dot(pa,ba);
      float x = (pa*baba-ba*paba).length() - radius*baba;
      float y = std::fabs(paba-baba*0.5)-baba*0.5;
      float x2 = x*x;
      float y2 = y*y*baba;
      float d = (std::fmax(x,y)<0.0)?-std::fmin(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
      return(d > 0 ? std::sqrt(std::fabs(d)) * inv_baba - corner_radius : -std::sqrt(std::fabs(d)) * inv_baba - corner_radius);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      vec3f min = vec3f(ffmin(start.x(),end.x()),ffmin(start.y(),end.y()),ffmin(start.z(),end.z()));
      vec3f max = vec3f(ffmax(start.x(),end.x()),ffmax(start.y(),end.y()),ffmax(start.z(),end.z()));
      box = aabb(min-radius, max+radius);
      return(true);
    }
    point3f start, end;
    vec3f ba; 
    Float radius, corner_radius, baba, inv_baba;
}; 

class csg_ellipsoid : public ImplicitShape { 
  public: 
    csg_ellipsoid(const point3f& c, point3f axes) : center(c), axes(axes)  {
      inv_axes = point3f(1.0/axes.x(),1.0/axes.y(),1.0/axes.z());
    } 
    Float getDistance(const point3f& from_old) const {
      point3f from = point3f(from_old - center); 
      float k0 = (from * inv_axes).length();
      float k1 = (from * (inv_axes*inv_axes)).length();
      return(k0 < 1.0 ? (k0 - 1.0) * MinComponent(axes) : k0*(k0-1.0)/k1);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = aabb(center+-axes, center+axes);
      return(true);
    }
    point3f center, axes, inv_axes; 
}; 

class csg_rounded_cone : public ImplicitShape { 
  public: 
    csg_rounded_cone(point3f start, point3f end, Float r1, Float r2) :  
      start(start), end(end), r1(r1), r2(r2) {} 
    Float getDistance(const point3f& from) const {
      vec3f pa = from - start; 
      vec3f ba = end - start;
      float l2 = dot(ba,ba);
      float rr = r1 - r2;
      float a2 = l2 - rr*rr;
      float il2 = 1.0/l2;
      
      float y = dot(pa,ba);
      float z = y - l2;
      float x2 = dot(pa*l2 - ba*y,pa*l2 - ba*y);
      float y2 = y*y*l2;
      float z2 = z*z*l2;
      
      float k = sgn(rr) * rr*rr*x2;
      if(sgn(z)*a2*z2 > k ) {
        return(std::sqrt(x2 + z2)*il2 - r2);
      }
      if(sgn(y)*a2*y2 < k ) {
        return(std::sqrt(x2 + y2)*il2 - r1);
      }
      return((std::sqrt(x2*a2*il2)+y*rr)*il2 - r1);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      point3f min = point3f(ffmin(start.x(),end.x()),ffmin(start.y(),end.y()),ffmin(start.z(),end.z()));
      point3f max = point3f(ffmax(start.x(),end.x()),ffmax(start.y(),end.y()),ffmax(start.z(),end.z()));
      box = aabb(min-r1-r2, max+r1+r2);
      return(true);
    }
    point3f start, end; 
    Float r1, r2;
}; 

class csg_cone : public ImplicitShape { 
  public: 
    csg_cone(point3f start, point3f end, Float radius) :  start(start), end(end), 
      radius(radius){
      height = (end - start).length();
      axis.build_from_w(end - start);
      axis.swap_yz();
    } 
    Float getDistance(const point3f& from_old) const {
      vec3f from_trans = axis.world_to_local(from_old - start);
      vec3f from = from_trans - vec3f(0,height,0);
      vec2f q = vec2f(radius,-height);
      
      vec2f w = vec2f( std::sqrt(from.x() * from.x() +  from.z() *  from.z()), from.y() );
      vec2f a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
      vec2f b = w - q*vec2f( clamp( w.x()/q.x(), 0.0, 1.0 ), 1.0 );
      float k = sgn( q.y() );
      float d = std::min(dot( a, a ),dot(b, b));
      float s = std::fmax( k*(w.x()*q.y()-w.y()*q.x()),k*(w.y()-q.y()));
      return(std::sqrt(d)*sgn(s));
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      point3f min = point3f(ffmin(start.x(),end.x()),ffmin(start.y(),end.y()),ffmin(start.z(),end.z()));
      point3f max = point3f(ffmax(start.x(),end.x()),ffmax(start.y(),end.y()),ffmax(start.z(),end.z()));
      box = aabb(min-radius, max+radius);
      return(true);
    }
    point3f start, end; 
    Float radius, height;
    onb axis;
}; 

class csg_pyramid : public ImplicitShape { 
  public: 
    csg_pyramid(const point3f& c, float h, float base) : center_bottom(c), h(h), base(base) {
      base_inv = vec3f(1.0/base,1,1.0/base);
      m2 = h*h + 0.25;
      m2_inv = 1/m2;
      m2_inv_buff = 1 / (m2+0.25);
    } 
    Float getDistance(const point3f& from_old) const {
      vec3f from = from_old - center_bottom;
      from = from * base_inv;
      from.e[0] = std::fabs(from.e[0]); 
      from.e[2] = std::fabs(from.e[2]);
      from = (from.z() > from.x()) ? vec3f(from.z(),from.y(),from.x()) : from;
      from -= vec3f(0.5,0,0.5);
      
      vec3f q = vec3f(from.z(), h*from.y() - 0.5*from.x(), h*from.x() + 0.5*from.y());
      
      float s = std::fmax(-q.x(),0.0);
      float t = clamp((q.y()-0.5*from.z())*m2_inv_buff, 0.0, 1.0 );
      
      float a = m2*(q.x()+s)*(q.x()+s) + q.y()*q.y();
      float b = m2*(q.x()+0.5*t)*(q.x()+0.5*t) + (q.y()-m2*t)*(q.y()-m2*t);
      
      float d2 = std::fmin(q.y(),-q.x()*m2-q.y()*0.5) > 0.0 ? 0.0 : std::fmin(a,b);
      
      return(std::sqrt((d2+q.z()*q.z()) * m2_inv ) * sgn(std::fmax(q.z(),-from.y())));
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = aabb(center_bottom-vec3f(base,0,base), center_bottom+vec3f(base,h,base));
      return(true);
    }
    point3f center_bottom;
    vec3f base_inv;
    float h, base, m2, m2_inv, m2_inv_buff;
}; 

class csg_triangle : public ImplicitShape { 
  public: 
    csg_triangle(point3f v1,point3f v2,point3f v3) : a(v1), b(v2), c(v3) {
      ba = b - a; 
      cb = c - b; 
      ac = a - c; 
      nor = cross( ba, ac );
    } 
    Float getDistance(const point3f& from) const {
      vec3f pa = from - a;
      vec3f pb = from - b;
      vec3f pc = from - c;
      return(std::sqrt((sgn(dot(cross(ba,nor),pa)) +
                   sgn(dot(cross(cb,nor),pb)) +
                   sgn(dot(cross(ac,nor),pc)) < 2.0) ?
                   std::min(std::min(
            (ba*clamp(dot(ba,pa)/ba.squared_length(),0.0,1.0)-pa).squared_length(),
            (cb*clamp(dot(cb,pb)/cb.squared_length(),0.0,1.0)-pb).squared_length() ),
            (ac*clamp(dot(ac,pc)/ac.squared_length(),0.0,1.0)-pc).squared_length() ) :
              dot(nor,pa)*dot(nor,pa)/nor.squared_length() ));
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      point3f min = point3f(ffmin(ffmin(a.x(),b.x()),c.x()),
                      ffmin(ffmin(a.y(),b.y()),c.y()),
                      ffmin(ffmin(a.z(),b.z()),c.z()));
      
      point3f max = point3f(ffmax(ffmax(a.x(),b.x()),c.x()),
                      ffmax(ffmax(a.y(),b.y()),c.y()),
                      ffmax(ffmax(a.z(),b.z()),c.z()));
      box = Expand(aabb(min, max),0.01);
      return(true);
    }
    point3f a,b,c;
    vec3f ba, cb, ac, nor;
}; 

class csg_elongate : public ImplicitShape { 
  public: 
    csg_elongate(std::shared_ptr<ImplicitShape> shape, point3f center,vec3f elongate) : 
      shape(shape),center(center), elongate(elongate) {} 
    Float getDistance(const point3f& from_old) const {
      vec3f from = from_old - center;
      point3f q = from - clamp(from, -elongate, elongate);
      return(shape->getDistance(q + center));
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      shape->bbox(t0,t1,box);
      box = Expand(box, elongate);
      return(true);
    }
    std::shared_ptr<ImplicitShape> shape;
    point3f center;
    vec3f elongate;
}; 

class csg_elongate_robust : public ImplicitShape { 
  public: 
    csg_elongate_robust(std::shared_ptr<ImplicitShape> shape, point3f center, vec3f elongate) : 
      shape(shape),center(center),elongate(elongate) {} 
    Float getDistance(const point3f& from_old) const {
      vec3f from = from_old - center;
      const static vec3f zeros(0,0,0);
      const static vec3f inf(INFINITY,INFINITY,INFINITY);
      vec3f q = Abs(from)-elongate;
      return(shape->getDistance(sgn(from) * clamp(q,zeros,inf) + vec3f(center)) + 
             std::fmin(MaxComponent(q),0.0));
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      shape->bbox(t0,t1,box);
      box = Expand(box, elongate);
      return(true);
    }
    std::shared_ptr<ImplicitShape> shape;
    point3f center;
    vec3f elongate;
}; 

class csg_round : public ImplicitShape { 
  public: 
    csg_round(std::shared_ptr<ImplicitShape> shape, Float r) : 
      shape(shape),r(r) {} 
    Float getDistance(const point3f& from) const {
      return(shape->getDistance(from) - r);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      shape->bbox(t0,t1,box);
      box = Expand(box, r);
      return(true);
    }
    std::shared_ptr<ImplicitShape> shape;
    Float r;
}; 

class csg_onion : public ImplicitShape { 
  public: 
    csg_onion(std::shared_ptr<ImplicitShape> shape, Float thickness) : 
    shape(shape),thickness(thickness) {} 
    Float getDistance(const point3f& from) const {
      return(fabs(shape->getDistance(from)) - thickness);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      shape->bbox(t0,t1,box);
      box = Expand(box,thickness);
      return(true);
    }
    std::shared_ptr<ImplicitShape> shape;
    Float thickness;
}; 

class csg_scale : public ImplicitShape { 
  public: 
    csg_scale(std::shared_ptr<ImplicitShape> shape, Float scale) : 
      shape(shape),scale(scale) {} 
    Float getDistance(const point3f& from) const {
      return(shape->getDistance(from/scale)*scale);
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      shape->bbox(t0,t1,box);
      box = aabb(box.min() * scale, box.max() * scale);
      return(true);
    }
    std::shared_ptr<ImplicitShape> shape;
    Float scale;
}; 

class csg_rotate : public ImplicitShape { 
  public: 
    csg_rotate(std::shared_ptr<ImplicitShape> shape, vec3f pivot_point, vec3f up) : 
      shape(shape),pivot_point(pivot_point), up(up) {
      axis.build_from_w(up);
      axis.swap_yz();
      
      aabb box2;
      shape->bbox(0,0,box2); //Pre-compute rotated bbox
      vec3f corners[8];
      corners[0] = axis.local_to_world(vec3f(box2.min()) - pivot_point) + pivot_point;
      corners[1] = axis.local_to_world(vec3f(box2.min().x(), box2.min().y(), box2.max().z())- pivot_point) + pivot_point;
      corners[2] = axis.local_to_world(vec3f(box2.min().x(), box2.max().y(), box2.min().z())- pivot_point) + pivot_point;
      corners[3] = axis.local_to_world(vec3f(box2.max().x(), box2.min().y(), box2.min().z())- pivot_point) + pivot_point;
      corners[4] = axis.local_to_world(vec3f(box2.min().x(), box2.max().y(), box2.max().z())- pivot_point) + pivot_point;
      corners[5] = axis.local_to_world(vec3f(box2.max().x(), box2.min().y(), box2.max().z())- pivot_point) + pivot_point;
      corners[6] = axis.local_to_world(vec3f(box2.max().x(), box2.max().y(), box2.min().z())- pivot_point) + pivot_point;
      corners[7] = axis.local_to_world(vec3f(box2.max()) - pivot_point) + pivot_point;
      vec3f temp_min = corners[0];
      vec3f temp_max = corners[7];
      
      for(int i = 1; i < 7; i++) {
        temp_min = vec3f(ffmin(temp_min.x(), corners[i].x()),
                        ffmin(temp_min.y(), corners[i].y()),
                        ffmin(temp_min.z(), corners[i].z()));
        temp_max = vec3f(ffmax(temp_max.x(), corners[i].x()),
                        ffmax(temp_max.y(), corners[i].y()),
                        ffmax(temp_max.z(), corners[i].z()));
      }
      aabb new_box(temp_min,temp_max);
      box_cache = new_box;
    } 
    csg_rotate(std::shared_ptr<ImplicitShape> shape, vec3f pivot_point, vec3f up, vec3f axis_x, vec3f axis_z) : 
      shape(shape),pivot_point(pivot_point), up(up) {
      axis.axis[0] = axis_x;
      axis.axis[1] = up;
      axis.axis[2] = axis_z;
      
      aabb box2;
      shape->bbox(0,0,box2); //Pre-compute rotated bbox
      vec3f corners[8];
      corners[0] = axis.local_to_world(vec3f(box2.min()) - pivot_point) + pivot_point;
      corners[1] = axis.local_to_world(vec3f(box2.min().x(), box2.min().y(), box2.max().z())- pivot_point) + pivot_point;
      corners[2] = axis.local_to_world(vec3f(box2.min().x(), box2.max().y(), box2.min().z())- pivot_point) + pivot_point;
      corners[3] = axis.local_to_world(vec3f(box2.max().x(), box2.min().y(), box2.min().z())- pivot_point) + pivot_point;
      corners[4] = axis.local_to_world(vec3f(box2.min().x(), box2.max().y(), box2.max().z())- pivot_point) + pivot_point;
      corners[5] = axis.local_to_world(vec3f(box2.max().x(), box2.min().y(), box2.max().z())- pivot_point) + pivot_point;
      corners[6] = axis.local_to_world(vec3f(box2.max().x(), box2.max().y(), box2.min().z())- pivot_point) + pivot_point;
      corners[7] = axis.local_to_world(vec3f(box2.max()) - pivot_point) + pivot_point;
      vec3f temp_min = corners[0];
      vec3f temp_max = corners[7];
      
      for(int i = 1; i < 7; i++) {
        temp_min = vec3f(ffmin(temp_min.x(), corners[i].x()),
                        ffmin(temp_min.y(), corners[i].y()),
                        ffmin(temp_min.z(), corners[i].z()));
        temp_max = vec3f(ffmax(temp_max.x(), corners[i].x()),
                        ffmax(temp_max.y(), corners[i].y()),
                        ffmax(temp_max.z(), corners[i].z()));
      }
      aabb new_box(temp_min,temp_max);
      box_cache = new_box;
    } 
    Float getDistance(const point3f& from) const {
      return(shape->getDistance(axis.world_to_local(from - pivot_point) + vec3f(pivot_point)));
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      box = box_cache;
      return(true);
    }
    std::shared_ptr<ImplicitShape> shape;
    point3f pivot_point;
    vec3f up;
    onb axis; 
    aabb box_cache;
}; 

class csg_translate : public ImplicitShape { 
  public: 
    csg_translate(std::shared_ptr<ImplicitShape> shape, vec3f translate) : 
      shape(shape),translate(translate) {} 
    Float getDistance(const point3f& from) const {
      return(shape->getDistance(from - translate));
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      shape->bbox(t0,t1,box);
      box.bounds[0] += translate; 
      box.bounds[1] += translate; 
      return(true);
    }
    std::shared_ptr<ImplicitShape> shape;
    vec3f translate;
}; 

struct unionFunc { 
  float operator() (float a, float b) const { return std::min(a, b); } 
}; 

struct subtractFunc { 
  float operator() (float a, float b) const { return std::fmax(a, -b); } 
}; 

struct intersectionFunc { 
  float operator() (float a, float b) const { return std::fmax(a, b); } 
}; 

struct blendFunc { 
  blendFunc(const float &_k) : k(_k) {} 
  float operator() (float a, float b) const { 
    float h = std::fmax(k-std::fabs(a-b), 0.0 )/k;
    return(std::fmin( a, b ) - h*h*k*(1.0/4.0));
  } 
  float k; 
}; 

struct blendFuncMinus { 
  blendFuncMinus(const float &_k) : k(_k) {} 
  float operator() (float a, float b) const { 
    float h = std::fmax(k-std::fabs(-b-a), 0.0)/k;
    return(std::fmax(-b, a) + h*h*k*(1.0/4.0));
  }
  float k; 
}; 

struct mixFunc { 
  mixFunc(const float &_t) : t(_t) {} 
  float operator() (float a, float b) const { 
    return a * (1 - t) + b * t; 
  } 
  float t; 
}; 

template<typename Op, typename ... Args> 
class CSG : public ImplicitShape { 
  public: 
    CSG( 
      const std::shared_ptr<ImplicitShape> s1, 
      const std::shared_ptr<ImplicitShape> s2, 
      Args&& ... args) : op(std::forward<Args>(args) ...), shape1(s1), shape2(s2) {} 
    Float getDistance(const point3f& from) const { 
      return(op(shape1->getDistance(from), shape2->getDistance(from))); 
    } 
    virtual bool bbox(Float t0, Float t1, aabb& box) const {
      aabb temp1, temp2;
      shape1->bbox(t0,t1,temp1);
      shape2->bbox(t0,t1,temp2);
      box = surrounding_box(temp1,temp2);
      return(true); 
    } 
    Op op; 
    const std::shared_ptr<ImplicitShape> shape1, shape2; 
}; 

using Union = CSG<unionFunc>;
using Subtract = CSG<subtractFunc>;
using Intersect = CSG<intersectionFunc>;
using Blend = CSG<blendFunc, float>;
using BlendMinus = CSG<blendFuncMinus, float>;
using Mix = CSG<mixFunc, float>;

class csg: public hitable {
  public:
    csg() {}
    csg(std::shared_ptr<material> mat, std::shared_ptr<ImplicitShape> shapes,
        std::shared_ptr<Transform> ObjectToWorld, std::shared_ptr<Transform> WorldToObject, bool reverseOrientation) : 
          hitable(ObjectToWorld, WorldToObject, reverseOrientation), 
          mat_ptr(mat), shapes(shapes) {
      aabb box;
      bool temp = shapes->bbox(0,1,box);
      if(temp) {
        max_dist = fmax(100,(box.max()-box.min()).length());
      }
      if(std::isinf(max_dist)) {
        Rcpp::Rcout << "min: " << box.min() << "\n";
        Rcpp::Rcout << "max: " << box.max() << "\n";
        throw std::runtime_error("error");
      }
    };
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, random_gen& rng);
    virtual bool hit(const ray& r, Float tmin, Float tmax, hit_record& rec, Sampler* sampler);
    
    virtual bool bounding_box(Float t0, Float t1, aabb& box) const;
    virtual Float pdf_value(const point3f& o, const vec3f& v, random_gen& rng, Float time) {
      return(1);
    }
    virtual Float pdf_value(const point3f& o, const vec3f& v, Sampler* sampler, Float time) {
      return(1);
    }
    virtual vec3f random(const point3f& o, random_gen& rng, Float time) {
      return(vec3f(0,1,0));
    }
    virtual vec3f random(const point3f& o, Sampler* sampler, Float time) {
      return(vec3f(0,1,0));
    }
    virtual std::string GetName() const {
      return(std::string("CSG"));
    }
    std::shared_ptr<material> mat_ptr;
    std::shared_ptr<ImplicitShape> shapes;
    Float max_dist;
    vec3f last_intersection;
};

inline std::shared_ptr<ImplicitShape> parse_csg(Rcpp::List csg_info) {
  if(Rcpp::as<int>(csg_info["csg_type"]) == 1) { //Operation
    std::shared_ptr<ImplicitShape> obj1 = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object1"]));
    std::shared_ptr<ImplicitShape> obj2 = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object2"]));
    int op_type = Rcpp::as<int>(csg_info["operation"]);
    if(op_type == 1) {
      return(std::make_shared<Union>(obj1, obj2));
    } else if(op_type == 2) {
      return(std::make_shared<Subtract>(obj1, obj2));
    } else if(op_type == 3) {
      return(std::make_shared<Intersect>(obj1, obj2));
    } else if(op_type == 4) {
      double blend_val = Rcpp::as<double>(csg_info["blend"]);
      return(std::make_shared<Blend>(obj1, obj2, blend_val));
    } else if(op_type == 5) {
      double blend_val = Rcpp::as<double>(csg_info["blend"]);
      return(std::make_shared<Mix>(obj1, obj2, blend_val));
    } else if(op_type == 6) {
      double blend_val = Rcpp::as<double>(csg_info["blend"]);
      return(std::make_shared<BlendMinus>(obj1, obj2, blend_val));
    }
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 2) { //Sphere
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    double r = Rcpp::as<double>(csg_info["radius"]);
    return(std::make_shared<csg_sphere>(vec3f(x,y,z), r));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 3) { //Plane
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    double width_x = Rcpp::as<double>(csg_info["width_x"]);
    double width_z = Rcpp::as<double>(csg_info["width_z"]);
    
    Rcpp::NumericVector n = Rcpp::as<Rcpp::NumericVector>(csg_info["normal"]);
    return(std::make_shared<csg_plane>(vec3f(n(0),n(1),n(2)),
                                       vec3f(x,y,z), width_x, width_z));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 4) { //Box
    Float x = Rcpp::as<Float>(csg_info["x"]);
    Float y = Rcpp::as<Float>(csg_info["y"]);
    Float z = Rcpp::as<Float>(csg_info["z"]);
    Rcpp::NumericVector w = Rcpp::as<Rcpp::NumericVector>(csg_info["width"]);
    return(std::make_shared<csg_box>(vec3f(x,y,z), vec3f(w(0),w(1),w(2))));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 5) { //Curved Box
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    Rcpp::NumericVector w = Rcpp::as<Rcpp::NumericVector>(csg_info["width"]);
    double r = Rcpp::as<double>(csg_info["radius"]);
    return(std::make_shared<csg_rounded_box>(vec3f(x,y,z), vec3f(w(0),w(1),w(2)), r));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 6) { //List of shapes
    std::vector<std::shared_ptr<ImplicitShape> > shapes;
    Rcpp::List entries = Rcpp::as<Rcpp::List>(csg_info["shapes"]);
    for(int i = 0; i < entries.length(); i++) {
      shapes.push_back(parse_csg(entries[i]));
    }
    return(std::make_shared<csg_list>(shapes));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 7) { //Torus
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    double ring_radius = Rcpp::as<double>(csg_info["ring_radius"]);
    double cross_radius = Rcpp::as<double>(csg_info["cross_radius"]);
    return(std::make_shared<csg_torus>(vec3f(x,y,z), ring_radius, cross_radius));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 8) { //Capsule
    Rcpp::NumericVector start = Rcpp::as<Rcpp::NumericVector>(csg_info["start"]);
    Rcpp::NumericVector end = Rcpp::as<Rcpp::NumericVector>(csg_info["end"]);
    double r = Rcpp::as<double>(csg_info["radius"]);
    return(std::make_shared<csg_capsule>(vec3f(start(0),start(1),start(2)), 
                                       vec3f(end(0),end(1),end(2)), r));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 9) { //Cylinder
    Rcpp::NumericVector start = Rcpp::as<Rcpp::NumericVector>(csg_info["start"]);
    Rcpp::NumericVector end = Rcpp::as<Rcpp::NumericVector>(csg_info["end"]);
    double r = Rcpp::as<double>(csg_info["radius"]);
    double cr = Rcpp::as<double>(csg_info["corner_radius"]);
    return(std::make_shared<csg_cylinder>(vec3f(start(0),start(1),start(2)), 
                                         vec3f(end(0),end(1),end(2)), r, cr));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 10) { //Ellipsoid
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    Rcpp::NumericVector a = Rcpp::as<Rcpp::NumericVector>(csg_info["axes"]);
    return(std::make_shared<csg_ellipsoid>(vec3f(x,y,z), vec3f(a(0),a(1),a(2))));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 11) { //Rounded Cone
    Rcpp::NumericVector start = Rcpp::as<Rcpp::NumericVector>(csg_info["start"]);
    Rcpp::NumericVector end = Rcpp::as<Rcpp::NumericVector>(csg_info["end"]);
    double r = Rcpp::as<double>(csg_info["radius"]);
    double cr = Rcpp::as<double>(csg_info["upper_radius"]);
    return(std::make_shared<csg_rounded_cone>(vec3f(start(0),start(1),start(2)), 
                                      vec3f(end(0),end(1),end(2)), r, cr));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 12) { //Cone
    Rcpp::NumericVector start = Rcpp::as<Rcpp::NumericVector>(csg_info["start"]);
    Rcpp::NumericVector end = Rcpp::as<Rcpp::NumericVector>(csg_info["end"]);
    double r = Rcpp::as<double>(csg_info["radius"]);
    return(std::make_shared<csg_cone>(vec3f(start(0),start(1),start(2)), 
                                      vec3f(end(0),end(1),end(2)), r));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 13) { //Pyramid
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    double h = Rcpp::as<double>(csg_info["height"]);
    double base = Rcpp::as<double>(csg_info["base"]);
    return(std::make_shared<csg_pyramid>(vec3f(x,y,z), h, base));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 14) { //Triangle
    Rcpp::NumericVector v1 = Rcpp::as<Rcpp::NumericVector>(csg_info["v1"]);
    Rcpp::NumericVector v2 = Rcpp::as<Rcpp::NumericVector>(csg_info["v2"]);
    Rcpp::NumericVector v3 = Rcpp::as<Rcpp::NumericVector>(csg_info["v3"]);
    return(std::make_shared<csg_triangle>(vec3f(v1(0),v1(1),v1(2)),
                                          vec3f(v2(0),v2(1),v2(2)),
                                          vec3f(v3(0),v3(1),v3(2))));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 15) { //Elongate
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    Rcpp::NumericVector e = Rcpp::as<Rcpp::NumericVector>(csg_info["elongate"]);
    
    std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
    return(std::make_shared<csg_elongate>(obj, vec3f(x,y,z),vec3f(e(0),e(1),e(2))));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 16) { //Elongate (2)
    double x = Rcpp::as<double>(csg_info["x"]);
    double y = Rcpp::as<double>(csg_info["y"]);
    double z = Rcpp::as<double>(csg_info["z"]);
    Rcpp::NumericVector e = Rcpp::as<Rcpp::NumericVector>(csg_info["elongate"]);
    
    std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
    return(std::make_shared<csg_elongate_robust>(obj, vec3f(x,y,z),vec3f(e(0),e(1),e(2))));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 17) { //Round 
    double r = Rcpp::as<double>(csg_info["radius"]);
    std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
    return(std::make_shared<csg_round>(obj,r));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 18) { //Onion 
    double thick = Rcpp::as<double>(csg_info["thickness"]);
    std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
    return(std::make_shared<csg_onion>(obj,thick));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 19) { //Scale 
    double s = Rcpp::as<double>(csg_info["scale"]);
    std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
    return(std::make_shared<csg_scale>(obj,s));
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 20) { //Onion 
    Rcpp::NumericVector up = Rcpp::as<Rcpp::NumericVector>(csg_info["up"]);    
    Rcpp::NumericVector axis_x = Rcpp::as<Rcpp::NumericVector>(csg_info["axis_x"]);    
    Rcpp::NumericVector axis_z = Rcpp::as<Rcpp::NumericVector>(csg_info["axis_z"]);    
    Rcpp::NumericVector pp = Rcpp::as<Rcpp::NumericVector>(csg_info["pivot_point"]);    
    if(axis_x(0) == 0 && axis_x(1) == 0 && axis_x(2) == 0 && 
       axis_z(0) == 0 && axis_z(1) == 0 && axis_z(2) == 0) {
      std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
      return(std::make_shared<csg_rotate>(obj,vec3f(pp(0),pp(1),pp(2)),
                                          vec3f(up(0),up(1),up(2))));
    } else {
      std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
      return(std::make_shared<csg_rotate>(obj,vec3f(pp(0),pp(1),pp(2)),
                                          vec3f(up(0),up(1),up(2)),
                                          vec3f(axis_x(0),axis_x(1),axis_x(2)),
                                          vec3f(axis_z(0),axis_z(1),axis_z(2))));
    }
  } else if (Rcpp::as<int>(csg_info["csg_type"]) == 21) { //Onion 
    Rcpp::NumericVector t = Rcpp::as<Rcpp::NumericVector>(csg_info["translate"]);   
    std::shared_ptr<ImplicitShape> obj = parse_csg(Rcpp::as<Rcpp::List>(csg_info["object"]));
    return(std::make_shared<csg_translate>(obj,vec3f(t(0),t(1),t(2))));
  }
  return(nullptr);
}

#endif
