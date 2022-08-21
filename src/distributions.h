#ifndef DISTRIBUTIONSH
#define DISTRIBUTIONSH

#include "float.h"
#include "vec2.h"
#include "mathinline.h"
#include <vector>
#include <memory>


struct Distribution1D {
  Distribution1D(const Float *f, int n);
  int Count() const;
  Float SampleContinuous(Float u, Float *pdf, int *off = nullptr) const;
  int SampleDiscrete(Float u, Float *pdf = nullptr,
                     Float *uRemapped = nullptr) const;
  Float DiscretePDF(int index) const;
  size_t GetSize();
  // Distribution1D Public Data
  std::vector<Float> func, cdf;
  Float funcInt;
};

class Distribution2D {
public:
  // Distribution2D Public Methods
  Distribution2D(const Float *data, int nu, int nv);
  vec2f SampleContinuous(const vec2f &u, Float *pdf) const;
  Float Pdf(const vec2f &p) const;
  size_t GetSize();
private:
  // Distribution2D Private Data
  std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
  std::unique_ptr<Distribution1D> pMarginal;
};


#endif
