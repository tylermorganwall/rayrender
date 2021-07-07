#include "distributions.h"

Distribution1D::Distribution1D(const Float *f, int n) : func(f, f + n), cdf(n + 1) {
  // Compute integral of step function at $x_i$
  cdf[0] = 0;
  for (int i = 1; i < n + 1; ++i) {
    cdf[i] = cdf[i - 1] + func[i - 1] / n;
  }
  
  // Transform step function integral into CDF
  funcInt = cdf[n];
  if (funcInt == 0) {
    for (int i = 1; i < n + 1; ++i) {
      cdf[i] = Float(i) / Float(n);
    }
  } else {
    for (int i = 1; i < n + 1; ++i) {
      cdf[i] /= funcInt;
    }
  }
}

int Distribution1D::Count() const { 
  return((int)func.size()); 
}

Float Distribution1D::SampleContinuous(Float u, Float *pdf, int *off) const {
  // Find surrounding CDF segments and _offset_
  int offset = FindInterval((int)cdf.size(),
                            [&](int index) { return cdf[index] <= u;});
  if (off) {
    *off = offset;
  }
  // Compute offset along CDF segment
  Float du = u - cdf[offset];
  if ((cdf[offset + 1] - cdf[offset]) > 0) {
    du /= (cdf[offset + 1] - cdf[offset]);
  }
  // Compute PDF for sampled offset
  if (pdf) {
    *pdf = (funcInt > 0) ? func[offset] / funcInt : 0;
  }
  // Return $x\in{}[0,1)$ corresponding to sample
  return (offset + du) / Count();
}

int Distribution1D::SampleDiscrete(Float u, Float *pdf,
                   Float *uRemapped) const {
  // Find surrounding CDF segments and _offset_
  int offset = FindInterval((int)cdf.size(),
                            [&](int index) { return cdf[index] <= u; });
  if (pdf) {
    *pdf = (funcInt > 0) ? func[offset] / (funcInt * Count()) : 0;
  }
  if (uRemapped) {
    *uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
  }
  return offset;
}

Float Distribution1D::DiscretePDF(int index) const {
  return(func[index] / (funcInt * Count()));
}



vec2f Distribution2D::SampleContinuous(const vec2f &u, Float *pdf) const {
  Float pdfs[2];
  int v;
  Float d1 = pMarginal->SampleContinuous(u[1], &pdfs[1], &v);
  Float d0 = pConditionalV[v]->SampleContinuous(u[0], &pdfs[0]);
  *pdf = pdfs[0] * pdfs[1];
  return(vec2f(d0, d1));
}

Float Distribution2D::Pdf(const vec2f &p) const {
  int iu = clamp(int(p[0] * pConditionalV[0]->Count()), 0,
                 pConditionalV[0]->Count() - 1);
  int iv = clamp(int(p[1] * pMarginal->Count()), 0, pMarginal->Count() - 1);
  return(pConditionalV[iv]->func[iu] / pMarginal->funcInt);
}

Distribution2D::Distribution2D(const Float *func, int nu, int nv) {
  pConditionalV.reserve(nv);
  for (int v = 0; v < nv; ++v) {
    // Compute conditional sampling distribution for $\tilde{v}$
    pConditionalV.emplace_back(new Distribution1D(&func[v * nu], nu));
  }
  // Compute marginal sampling distribution $p[\tilde{v}]$
  std::vector<Float> marginalFunc;
  marginalFunc.reserve(nv);
  for (int v = 0; v < nv; ++v)
    marginalFunc.push_back(pConditionalV[v]->funcInt);
  pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
}
