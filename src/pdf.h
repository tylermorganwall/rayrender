#ifndef PDFH
#define PDFH

#include "onbh.h"
#include "hitable.h"

inline vec3 random_cosine_direction() {
  float r1 = drand48();
  float r2 = drand48();
  float z = sqrt(1-r2);
  float phi = 2*M_PI*r1;
  float x = cos(phi)*2*sqrt(r2);
  float y = sin(phi)*2*sqrt(r2);
  return vec3(x, y, z);
}

class pdf {
public: 
  virtual float value(const vec3& direction) const = 0;
  virtual vec3 generate() const = 0;
};

class cosine_pdf : public pdf {
public:
  cosine_pdf(const vec3& w) {
    uvw.build_from_w(w);
  }
  virtual float value(const vec3& direction) const {
    float cosine = dot(unit_vector(direction), uvw.w());
    if(cosine > 0) {
      return(cosine/M_PI);
    } else {
      return(0);
    }
  } 
  virtual vec3 generate() const {
    return(uvw.local(random_cosine_direction()));
  }
  onb uvw;
};

class hitable_pdf : public pdf {
public:
  hitable_pdf(hitable *p, const vec3& origin) : ptr(p), o(origin) {}
  virtual float value(const vec3& direction) const {
    return(ptr->pdf_value(o,direction));
  }
  virtual vec3 generate() const {
    return(ptr->random(o)); 
  }
  hitable *ptr;
  vec3 o;
};

class mixture_pdf : public pdf {
public:
  mixture_pdf(pdf *p0, pdf *p1) {
    p[0] = p0;
    p[1] = p1;
  }
  virtual float value(const vec3& direction) const {
    return(0.5 * p[0]->value(direction) + 0.5 * p[1]->value(direction));
  }
  virtual vec3 generate() const {
    if(drand48() < 0.5) {
      return(p[0]->generate());
    } else {
      return(p[1]->generate());
    } 
  }
  pdf *p[2];
};

#endif