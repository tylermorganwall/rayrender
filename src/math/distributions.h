#ifndef DISTRIBUTIONSH
#define DISTRIBUTIONSH

#include "../math/float.h"
#include "../math/vec2.h"
#include "../math/mathinline.h"
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

// Piecewise-constant UV density with Walker/Vose lookup on both axes. Keeping
// row and conditional-column tables avoids spending most of a float variate's
// mantissa on a single large 2D cell index. Each remapped coin jitters its axis.
class AliasDistribution2D {
public:
  AliasDistribution2D(const Float *data, int width, int height);
  vec2f SampleContinuous(const vec2f &u, Float *pdf) const;
  Float Pdf(const vec2f &uv) const;
  size_t GetSize() const;

private:
  struct Entry {
    Float cutoff = 1, density = 0;
    uint32_t alias = 0;
  };
  static void Build(const std::vector<double> &weights, Entry *entries);
  static int Sample(const Entry *entries, int count, Float u, Float &residual);
  int width, height;
  std::vector<Entry> columns, rows;
};


#endif
