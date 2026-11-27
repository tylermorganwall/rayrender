#ifndef RAYRENDER_INTEL_FUNCTIONAL_H
#define RAYRENDER_INTEL_FUNCTIONAL_H

// icpx 2024 diagnoses libstdc++ 13's own use of std::result_of in std::bind
// when RcppThread instantiates parallelFor. Scope the workaround to that header.
#if defined(__INTEL_LLVM_COMPILER) && __INTEL_LLVM_COMPILER < 20250000
#include <cstddef>
#if defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE == 13
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <functional>
#pragma clang diagnostic pop
#endif
#endif

#endif
