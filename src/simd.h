#ifndef SIMDH
#define SIMDH

#include <cstdint>
#include <algorithm>
#include <stdio.h>
#include "point3.h"
#include <cstddef> // For alignas

// #define HAS_NEON

// SIMD vector size (4 for SSE, 8 for AVX)
#ifdef HAS_AVX
#define SIMD_WIDTH 8
#else
#define SIMD_WIDTH 4
#endif


#ifdef HAS_NEON
    #include <arm_neon.h>
#elif defined(HAS_SSE)
    #include <xmmintrin.h>
#else
    #error "No SIMD support available"
#endif

#ifdef HAS_SSE
typedef __m128 fvec4;
typedef __m128i ivec4;
#elif defined(HAS_NEON)
typedef uint32x4_t uivec4;
typedef int32x4_t ivec4;
typedef float32x4_t fvec4;
#endif

typedef struct alignas(16) FVec4 {
public:
    union {
        fvec4 v;
        float xyzw[4];
    };
    FVec4(fvec4 f4) : v(f4) {}
    FVec4()  {
      xyzw[0] = 0.f;
      xyzw[1] = 0.f;
      xyzw[2] = 0.f;
      xyzw[3] = 0.f;
    }
    FVec4(float x_, float y_, float z_, float w_) {
      xyzw[0] = x_;
      xyzw[1] = y_;
      xyzw[2] = z_;
      xyzw[3] = w_;
    }
    FVec4(float x_, float y_, float z_) {
      xyzw[0] = x_;
      xyzw[1] = y_;
      xyzw[2] = z_;
      xyzw[3] = 0.f;
    }
    FVec4(point3f v) {
      xyzw[0] = v.x();
      xyzw[1] = v.y();
      xyzw[2] = v.z();
      xyzw[3] = 0;
    }

    float operator[](int i) const { return xyzw[i]; }
    float& operator[](int i) { return xyzw[i]; }
} FVec4;


typedef struct alignas(16) IVec4 {
    union {
#ifdef HAS_SSE
        __m128i v;
#elif defined(HAS_NEON)
        int32x4_t v;
        uint32x4_t uv;
#endif
        alignas(16) int xyzw[4];
    };

    IVec4()  {
      xyzw[0] = 0;
      xyzw[1] = 0;
      xyzw[2] = 0;
      xyzw[3] = 0;
    }
    IVec4(int x_, int y_, int z_, int w_)  {
      xyzw[0] = x_;
      xyzw[1] = y_;
      xyzw[2] = z_;
      xyzw[3] = w_;
    }
    IVec4(int x_, int y_, int z_) {
      xyzw[0] = x_;
      xyzw[1] = y_;
      xyzw[2] = z_;
      xyzw[3] = 0;
    }
#ifdef HAS_SSE
    IVec4(__m128i m128i_) : v(m128i_) {};
#elif defined(HAS_NEON)
    IVec4(int32x4_t i32x4_) : v(i32x4_) {};
#endif

  int operator[](int i) const { return xyzw[i]; }
  int& operator[](int i) { return xyzw[i]; }

  int i0() const {
  #ifdef HAS_NEON
      return vgetq_lane_s32(v, 0);
  #elif defined(HAS_SSE)
      return _mm_extract_epi32(v, 0);
  #else
      return xyzw[0];
  #endif
  }

  int i1() const {
  #ifdef HAS_NEON
      return vgetq_lane_s32(v, 1);
  #elif defined(HAS_SSE)
      return _mm_extract_epi32(v, 1);
  #else
      return xyzw[1];
  #endif
  }

  int i2() const {
  #ifdef HAS_NEON
      return vgetq_lane_s32(v, 2);
  #elif defined(HAS_SSE)
      return _mm_extract_epi32(v, 2);
  #else
      return xyzw[2];
  #endif
  }

  int i3() const {
  #ifdef HAS_NEON
      return vgetq_lane_s32(v, 3);
  #elif defined(HAS_SSE)
      return _mm_extract_epi32(v, 3);
  #else
      return xyzw[3];
  #endif
  }
} IVec4;


// SIMD operations
inline FVec4 simd_load(const float* ptr) {
  FVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_load_ps(ptr);
#elif defined(HAS_SSE)
  result.v = _mm_load_ps(ptr);
#elif defined(HAS_NEON)
  result.v = vld1q_f32(ptr);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = ptr[i];
  }
#endif
  return result;
}


// SIMD operations
inline IVec4 simd_load_int(const int* ptr) {
  IVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_load_ps(ptr);
#elif defined(HAS_SSE)
  result.v = _mm_load_ps(ptr);
#elif defined(HAS_NEON)
  result.v[0] = ptr[0];
  result.v[1] = ptr[1];
  result.v[2] = ptr[2];
  result.v[3] = ptr[3];
  // result.v = vld1q_f32(ptr);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = ptr[i];
  }
#endif
  return result;
}

typedef FVec4 SimdMask;

inline FVec4 simd_set1(float value) {
  FVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_set1_ps(value);
#elif defined(HAS_SSE)
  result.v = _mm_set1_ps(value);
#elif defined(HAS_NEON)
  result.v = vdupq_n_f32(value);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = value;
  }
#endif
  return result;
}


inline FVec4 simd_sub(FVec4 a, FVec4 b) {
  FVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_sub_ps(a.v, b.v);
#elif defined(HAS_SSE)
  result.v = _mm_sub_ps(a.v, b.v);
#elif defined(HAS_NEON)
  result.v = vsubq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = a.v[i] - b.v[i];
  }
#endif
  return result;
}


inline FVec4 simd_mul(FVec4 a, FVec4 b) {
  FVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_mul_ps(a.v, b.v);
#elif defined(HAS_SSE)
  result.v = _mm_mul_ps(a.v, b.v);
#elif defined(HAS_NEON)
  result.v = vmulq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = a.v[i] * b.v[i];
  }
#endif
  return result;
}


inline FVec4 simd_min(FVec4 a, FVec4 b) {
  FVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_min_ps(a.v, b.v);
#elif defined(HAS_SSE)
  result.v = _mm_min_ps(a.v, b.v);
#elif defined(HAS_NEON)
  result.v = vminq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = std::min(a.v[i], b.v[i]);
  }
#endif
  return result;
}


inline FVec4 simd_max(FVec4 a, FVec4 b) {
  FVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_max_ps(a.v, b.v);
#elif defined(HAS_SSE)
  result.v = _mm_max_ps(a.v, b.v);
#elif defined(HAS_NEON)
  result.v = vmaxq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = std::max(a.v[i], b.v[i]);
  }
#endif
  return result;
}


inline FVec4 simd_less_equal(FVec4 a, FVec4 b) {
  FVec4 result;
#ifdef HAS_AVX
  result.v = _mm256_cmp_ps(a.v, b.v, _CMP_LE_OQ);
#elif defined(HAS_SSE)
  result.v = _mm_cmple_ps(a.v, b.v);
#elif defined(HAS_NEON)
  result.v = vcleq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = (a.v[i] <= b.v[i]) ? 0xFFFFFFFF : 0;
  }
#endif
  return result;
}


inline bool simd_any_true(SimdMask mask) {
#ifdef HAS_AVX
  return _mm256_movemask_ps(mask.v) != 0;
#elif defined(HAS_SSE)
  return _mm_movemask_ps(mask.v) != 0;
#elif defined(HAS_NEON)
  return vmaxvq_u32(vreinterpretq_u32_f32(mask.v)) != 0;
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    if (mask.v[i] != 0) return true;
  }
  return false;
#endif
}

inline IVec4 simd_cast_float_to_int(FVec4 mask) {
#ifdef HAS_SSE
    IVec4 result;
    result.v = _mm_castps_si128(mask.v);
    return result;
#elif defined(HAS_NEON)
    IVec4 result;
    result.v = vreinterpretq_s32_f32(mask.v);
    return result;
#else
    // Fallback for non-SIMD
    IVec4 result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = (mask.xyzw[i] != 0.0f) ? -1 : 0;
    }
    return result;
#endif
}


inline IVec4 simd_cast_to_int(SimdMask mask) {
#ifdef HAS_SSE
    IVec4 result;
    result.v = _mm_castps_si128(mask.v);
    return result;
#elif defined(HAS_NEON)
    IVec4 result;
    result.v = vreinterpretq_s32_f32(mask.v);
    return result;
#else
    // Fallback for non-SIMD
    IVec4 result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = (mask.xyzw[i] != 0.0f) ? -1 : 0;
    }
    return result;
#endif
}

inline SimdMask simd_setmask(bool m0, bool m1, bool m2, bool m3) {
    SimdMask result;
#ifdef HAS_SSE
    __m128i mask = _mm_set_epi32(
        m3 ? -1 : 0,
        m2 ? -1 : 0,
        m1 ? -1 : 0,
        m0 ? -1 : 0
    );
    result.v = _mm_castsi128_ps(mask);
#elif defined(HAS_NEON)
    int32x4_t mask_int = {
        m0 ? -1 : 0,
        m1 ? -1 : 0,
        m2 ? -1 : 0,
        m3 ? -1 : 0
    };
    result.v = vreinterpretq_f32_s32(mask_int);
#else
    for (int i = 0; i < 4; ++i) {
        uint32_t bits = ((i == 0 && m0) || (i == 1 && m1) || (i == 2 && m2) || (i == 3 && m3)) ? 0xFFFFFFFF : 0x00000000;
        result.xyzw[i] = reinterpret_cast<float&>(bits);
    }
#endif
    return result;
}


inline FVec4 simd_blend(SimdMask mask, FVec4 a, FVec4 b) {
#ifdef HAS_SSE
    // For SSE, use bitwise operations
    FVec4 result;
    result.v = _mm_or_ps(_mm_and_ps(mask.v, a.v), _mm_andnot_ps(mask.v, b.v));
    return result;
#elif defined(HAS_NEON)
    // For NEON, use vbslq_f32
    FVec4 result;
    result.v = vbslq_f32(vreinterpretq_u32_f32(mask.v), a.v, b.v);
    return result;
#else
    // Fallback for non-SIMD
    FVec4 result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = mask.xyzw[i] ? a.xyzw[i] : b.xyzw[i];
    }
    return result;
#endif
}

inline SimdMask simd_cmpge(FVec4 a, FVec4 b) {
#ifdef HAS_SSE
    SimdMask result;
    result.v = _mm_cmpge_ps(a.v, b.v);
    return result;
#elif defined(HAS_NEON)
    SimdMask result;
    result.v = vcgeq_f32(a.v, b.v);
    return result;
#else
    // Fallback for non-SIMD
    SimdMask result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = (a.xyzw[i] >= b.xyzw[i]) ? -1.0f : 0.0f;
    }
    return result;
#endif
}

// SSE Implementation
#ifdef HAS_SSE

inline FVec4 simd_swap_pairs(FVec4 a) {
    // Swap (0,1) and (2,3)
    FVec4 result;
    // Shuffle mask: _MM_SHUFFLE(z, y, x, w)
    // For swapping (0,1) and (2,3), we use _MM_SHUFFLE(2, 3, 0, 1)
    result.v = _mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(2, 3, 0, 1));
    return result;
}

inline FVec4 simd_reverse(FVec4 a) {
    // Reverse the vector to (3,2,1,0)
    FVec4 result;
    result.v = _mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(0, 1, 2, 3));
    return result;
}

// Similarly for IVec4
inline IVec4 simd_swap_pairs(IVec4 a) {
    // Swap (0,1) and (2,3)
    IVec4 result;
    result.v = _mm_shuffle_epi32(a.v, _MM_SHUFFLE(2, 3, 0, 1));
    return result;
}

inline IVec4 simd_reverse(IVec4 a) {
    // Reverse the vector to (3,2,1,0)
    IVec4 result;
    result.v = _mm_shuffle_epi32(a.v, _MM_SHUFFLE(0, 1, 2, 3));
    return result;
}

#endif // HAS_SSE

// NEON Implementation
#ifdef HAS_NEON

inline FVec4 simd_swap_pairs(FVec4 a) {
    // Swap (0,1) and (2,3)
    FVec4 result;
    // Swap the first two elements
    float32x2_t low = vget_low_f32(a.v);
    low = vrev64_f32(low); // Swap elements in low 64 bits
    // Swap the last two elements
    float32x2_t high = vget_high_f32(a.v);
    high = vrev64_f32(high); // Swap elements in high 64 bits
    // Recombine
    result.v = vcombine_f32(low, high);
    return result;
}

inline FVec4 simd_reverse(FVec4 a) {
    // Reverse the vector to (3,2,1,0)
    FVec4 result;
    // Use vrev64q_f32 to reverse pairs, then swap the pairs
    result.v = vrev64q_f32(a.v); // Swap elements within 64-bit lanes: (1,0,3,2)
    result.v = vcombine_f32(vget_high_f32(result.v), vget_low_f32(result.v)); // Swap the 64-bit lanes
    return result;
}

// Similarly for IVec4
inline IVec4 simd_swap_pairs(IVec4 a) {
    // Swap (0,1) and (2,3)
    IVec4 result;
    // Swap the first two elements
    int32x2_t low = vget_low_s32(a.v);
    low = vrev64_s32(low); // Swap elements in low 64 bits
    // Swap the last two elements
    int32x2_t high = vget_high_s32(a.v);
    high = vrev64_s32(high); // Swap elements in high 64 bits
    // Recombine
    result.v = vcombine_s32(low, high);
    return result;
}

inline IVec4 simd_reverse(IVec4 a) {
    // Reverse the vector to (3,2,1,0)
    IVec4 result;
    // Use vrev64q_s32 to reverse pairs, then swap the pairs
    result.v = vrev64q_s32(a.v); // Swap elements within 64-bit lanes: (1,0,3,2)
    result.v = vcombine_s32(vget_high_s32(result.v), vget_low_s32(result.v)); // Swap the 64-bit lanes
    return result;
}

#endif // HAS_NEON

// Fallback Implementation
#ifndef HAS_SSE
#ifndef HAS_NEON

inline FVec4 simd_swap_pairs(FVec4 a) {
    FVec4 result;
    result.v[0] = a.v[1];
    result.v[1] = a.v[0];
    result.v[2] = a.v[3];
    result.v[3] = a.v[2];
    return result;
}

inline FVec4 simd_reverse(FVec4 a) {
    FVec4 result;
    result.v[0] = a.v[3];
    result.v[1] = a.v[2];
    result.v[2] = a.v[1];
    result.v[3] = a.v[0];
    return result;
}

// Similarly for IVec4
inline IVec4 simd_swap_pairs(IVec4 a) {
    IVec4 result;
    result.v[0] = a.v[1];
    result.v[1] = a.v[0];
    result.v[2] = a.v[3];
    result.v[3] = a.v[2];
    return result;
}

inline IVec4 simd_reverse(IVec4 a) {
    IVec4 result;
    result.v[0] = a.v[3];
    result.v[1] = a.v[2];
    result.v[2] = a.v[1];
    result.v[3] = a.v[0];
    return result;
}

#endif // HAS_NEON
#endif // HAS_SSE


inline FVec4 simd_shuffle(FVec4 a, int idx0, int idx1, int idx2, int idx3) {
    FVec4 result;
    result.v[0] = a.v[idx0];
    result.v[1] = a.v[idx1];
    result.v[2] = a.v[idx2];
    result.v[3] = a.v[idx3];
    return result;
}

inline IVec4 simd_shuffle(IVec4 a, int idx0, int idx1, int idx2, int idx3) {
    IVec4 result;
    result.v[0] = a.v[idx0];
    result.v[1] = a.v[idx1];
    result.v[2] = a.v[idx2];
    result.v[3] = a.v[idx3];
    return result;
}

// inline SimdMask simd_cmpneq(IVec4 a, IVec4 b) {
// #ifdef HAS_AVX
//     SimdMask result;
//     result.v = _mm256_castsi256_ps(_mm256_cmp_epi32(a.v, b.v, _MM_CMPINT_NE));
//     return result;
// #elif defined(HAS_SSE)
//     SimdMask result;
//     result.v = _mm_castsi128_ps(_mm_cmpneq_epi32(a.v, b.v));
//     return result;
// #elif defined(HAS_NEON)
//     SimdMask result;
//     result.v = vreinterpretq_f32_u32(vmvnq_u32(vceqq_s32(a.v, b.v)));
//     return result;
// #else
//     SimdMask result;
//     for (int i = 0; i < SIMD_WIDTH; ++i) {
//         result.v[i] = (a.v[i] != b.v[i]) ? 0xFFFFFFFF : 0;
//     }
//     return result;
// #endif
// }

inline IVec4 simd_cmpneq(IVec4 a, IVec4 b) {
#ifdef HAS_AVX
    SimdMask result;
    result.v = _mm256_castsi256_ps(_mm256_cmp_epi32(a.v, b.v, _MM_CMPINT_NE));
    return result;
#elif defined(HAS_SSE)
    SimdMask result;
    result.v = _mm_castsi128_ps(_mm_cmpneq_epi32(a.v, b.v));
    return result;
#elif defined(HAS_NEON)
    IVec4 result;
    result.v = vmvnq_u32(vceqq_s32(a.v, b.v));
    return result;
#else
    SimdMask result;
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.v[i] = (a.v[i] != b.v[i]) ? 0xFFFFFFFF : 0;
    }
    return result;
#endif
}


// inline SimdMask simd_and(SimdMask a, SimdMask b) {
// #ifdef HAS_AVX
//     SimdMask result;
//     result.v = _mm256_and_ps(a.v, b.v);
//     return result;
// #elif defined(HAS_SSE)
//     SimdMask result;
//     result.v = _mm_and_ps(a.v, b.v);
//     return result;
// #elif defined(HAS_NEON)
//     SimdMask result;
//     result.v = vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(a.v), vreinterpretq_u32_f32(b.v)));
//     return result;
// #else
//     SimdMask result;
//     for (int i = 0; i < SIMD_WIDTH; ++i) {
//         result.v[i] = a.v[i] & b.v[i];
//     }
//     return result;
// #endif
// }

inline IVec4 simd_and(IVec4 a, IVec4 b) {
#ifdef HAS_AVX
    SimdMask result;
    result.v = _mm256_and_ps(a.v, b.v);
    return result;
#elif defined(HAS_SSE)
    SimdMask result;
    result.v = _mm_and_ps(a.v, b.v);
    return result;
#elif defined(HAS_NEON)
    IVec4 result;
    result.v = (vandq_u32((a.v), (b.v)));
    return result;
#else
    SimdMask result;
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.v[i] = a.v[i] & b.v[i];
    }
    return result;
#endif
}

inline void simd_extract_fvec4(FVec4 src, float* dest) {
#ifdef HAS_NEON
    vst1q_f32(dest, src.v); // Store NEON float32x4 to array
#elif defined(HAS_SSE)
    _mm_store_ps(dest, src.v); // Store SSE __m128 to array
#else
    // Fallback for non-SIMD
    for (int i = 0; i < 4; ++i) {
        dest[i] = src.xyzw[i];
    }
#endif
}

inline IVec4 simd_not_equals_minus_one(IVec4 a) {
#ifdef HAS_AVX2
    // AVX2 supports integer operations on 256-bit registers
    IVec4 result;
    __m256i minus_one = _mm256_set1_epi32(-1);
    __m256i cmp_eq = _mm256_cmpeq_epi32(a.v, minus_one);
    __m256i cmp_neq = _mm256_xor_si256(cmp_eq, _mm256_set1_epi32(-1));
    result.v = _mm256_srli_epi32(cmp_neq, 31);
    return result;
#elif defined(HAS_SSE2)
    // SSE2 supports integer operations on 128-bit registers
    IVec4 result;
    __m128i minus_one = _mm_set1_epi32(-1);
    __m128i cmp_eq = _mm_cmpeq_epi32(a.v, minus_one);
    __m128i cmp_neq = _mm_xor_si128(cmp_eq, _mm_set1_epi32(-1));
    result.v = _mm_srli_epi32(cmp_neq, 31);
    return result;
#elif defined(HAS_NEON)
    // NEON supports integer comparisons
    IVec4 result;
    int32x4_t minus_one = vdupq_n_s32(-1);
    uint32x4_t cmp_eq = vceqq_s32(a.v, minus_one);
    uint32x4_t cmp_neq = vmvnq_u32(cmp_eq);
    result.v = vshrq_n_u32(cmp_neq, 31);
    return result;
#else
    // Scalar fallback
    IVec4 result;
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.v[i] = (a.v[i] != -1) ? 1 : 0;
    }
    return result;
#endif
}

inline SimdMask simd_cmpgt(FVec4 a, FVec4 b) {
#ifdef HAS_AVX
    SimdMask result;
    result.v = _mm256_cmp_ps(a.v, b.v, _CMP_GT_OQ);
    return result;
#elif defined(HAS_SSE)
    SimdMask result;
    result.v = _mm_cmpgt_ps(a.v, b.v);
    return result;
#elif defined(HAS_NEON)
    SimdMask result;
    result.v = vcgtq_f32(a.v, b.v);
    return result;
#else
    SimdMask result;
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.v[i] = (a.v[i] > b.v[i]) ? 0xFFFFFFFF : 0;
    }
    return result;
#endif
}

inline SimdMask simd_cmplt(FVec4 a, FVec4 b) {
    #ifdef HAS_SSE
        SimdMask result;
        result.v = _mm_cmplt_ps(a.v, b.v);
        return result;
    #elif defined(HAS_NEON)
        SimdMask result;
        result.v = vcltq_f32(a.v, b.v);
        return result;
    #else
        // Fallback for non-SIMD
        SimdMask result;
        for (int i = 0; i < 4; ++i) {
            result.xyzw[i] = (a.xyzw[i] < b.xyzw[i]) ? -1.0f : 0.0f;
        }
        return result;
    #endif
}

inline IVec4 simd_blend_int(SimdMask mask, IVec4 a, IVec4 b) {
    IVec4 mask_int = simd_cast_to_int(mask);

    #ifdef HAS_SSE
        IVec4 result;
        result.v = _mm_or_si128(_mm_and_si128(mask_int.v, a.v), _mm_andnot_si128(mask_int.v, b.v));
        return result;
    #elif defined(HAS_NEON)
        IVec4 result;
        result.v = vbslq_s32(vreinterpretq_u32_s32(mask_int.v), a.v, b.v);
        return result;
    #else
        // Fallback for non-SIMD
        IVec4 result;
        for (int i = 0; i < 4; ++i) {
            result.xyzw[i] = mask_int.xyzw[i] ? a.xyzw[i] : b.xyzw[i];
        }
        return result;
    #endif
}

inline SimdMask simd_not(SimdMask mask) {
#ifdef HAS_SSE
    SimdMask result;
    result.v = _mm_xor_ps(mask.v, _mm_castsi128_ps(_mm_set1_epi32(-1)));
    return result;
#elif defined(HAS_NEON)
    SimdMask result;
    result.v = vreinterpretq_f32_u32(vmvnq_u32(vreinterpretq_u32_f32(mask.v)));
    return result;
#else
    SimdMask result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = ~reinterpret_cast<unsigned int&>(mask.xyzw[i]);
    }
    return result;
#endif
}

inline SimdMask simd_set1(unsigned int value) {
    SimdMask result;
#ifdef HAS_SSE
    result.v = _mm_castsi128_ps(_mm_set1_epi32(value));
#elif defined(HAS_NEON)
    result.v = vreinterpretq_f32_u32(vdupq_n_u32(value));
#else
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = *reinterpret_cast<float*>(&value);
    }
#endif
    return result;
}

#include <cstring>

// Helper functions to convert between float and uint32_t
inline uint32_t float_to_uint32(float f) {
    uint32_t i;
    memcpy(&i, &f, sizeof(float));
    return i;
}

inline float uint32_to_float(uint32_t i) {
    float f;
    memcpy(&f, &i, sizeof(float));
    return f;
}

// Pack index into the lowest 2 bits of the float
inline float pack_index(float value, uint32_t index) {
    uint32_t int_value = float_to_uint32(value);
    int_value = (int_value & ~0x3) | (index & 0x3); // Set the lowest 2 bits
    return uint32_to_float(int_value);
}

// Extract index from the lowest 2 bits of the float
inline uint32_t extract_index(float value) {
    uint32_t int_value = float_to_uint32(value);
    return int_value & 0x3;
}

inline int simd_extract_hitmask(const IVec4& vec) {
#ifdef HAS_SSE
    // SSE implementation

    // Mask the least significant bits (elements are 0 or 1)
    __m128i masked = _mm_and_si128(vec.v, _mm_set1_epi32(1));

    // Multiply each element by its corresponding power of 2
    // Powers of 2: [3]=8, [2]=4, [1]=2, [0]=1
    __m128i powers = _mm_set_epi32(8, 4, 2, 1);
    __m128i weighted = _mm_mullo_epi32(masked, powers);

    // Sum the weighted bits
    // Since SSE2 lacks horizontal add for integers, we sum manually
    __m128i temp1 = _mm_add_epi32(weighted, _mm_srli_si128(weighted, 8)); // Sum pairs: [0]+[2], [1]+[3]
    __m128i temp2 = _mm_add_epi32(temp1, _mm_srli_si128(temp1, 4));       // Sum all elements

    int mask = _mm_cvtsi128_si32(temp2); // Extract the lowest 32 bits

    return mask & 0xF; // Ensure only the lower 4 bits are used

#elif defined(HAS_NEON)
    // NEON implementation

    // Mask the least significant bits (elements are 0 or 1)
    uint32x4_t masked = vandq_u32(vec.v, vdupq_n_u32(1));

    // Multiply each element by its corresponding power of 2
    uint32x4_t powers = {1, 2, 4, 8}; // Powers of 2: 2^0, 2^1, 2^2, 2^3
    uint32x4_t weighted = vmulq_u32(masked, powers);

    // Sum the weighted bits
    uint32_t mask = vaddvq_u32(weighted); // Requires ARMv8.1-A

    return static_cast<int>(mask & 0xF);

#else
    // Scalar fallback
    int mask = ((vec[0] & 1) << 0) |
               ((vec[1] & 1) << 1) |
               ((vec[2] & 1) << 2) |
               ((vec[3] & 1) << 3);
    return mask;
#endif
}


inline IVec4 sort_simd_4_floats(FVec4 values);

#endif