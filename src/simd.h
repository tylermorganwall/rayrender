#ifndef SIMDH
#define SIMDH

#include <cstdint>
#include <algorithm>
#include <stdio.h>
#include <cstddef> // For alignas

// SIMD vector size (4 for SSE, 8 for AVX)
#ifdef HAS_AVX
#define SIMD_WIDTH 8
#else
#define SIMD_WIDTH 4
#endif

#ifndef HAS_SSE
#ifndef HAS_NEON
#define NO_SSE
#endif
#endif


#ifdef HAS_NEON
    #include <arm_neon.h>
#elif defined(HAS_SSE)
    #include <xmmintrin.h>
#endif

#ifdef HAS_SSE
    typedef __m128 fvec4;
    typedef __m128i ivec4;
#elif defined(HAS_NEON)
    typedef float32x4_t fvec4;
    typedef int32x4_t ivec4;
    typedef uint32x4_t uivec4;
#else
    typedef float fvec4[SIMD_WIDTH]; 
    typedef int ivec4[SIMD_WIDTH];
    typedef unsigned int uivec4[SIMD_WIDTH];
#endif

typedef struct alignas(16) FVec4 {
public:
    union {
        fvec4 v;
        float xyzw[4];
    };
#ifndef NO_SSE
    FVec4(fvec4 f4) : v(f4) {}
#else
    FVec4(fvec4 f4) {
      for (int i = 0; i < SIMD_WIDTH; ++i) {
          v[i] = f4[i];
      }
    }
#endif
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
    FVec4(float t_) {
#if defined(HAS_SSE)
        v = _mm_set1_ps(t_);
#elif defined(HAS_NEON)
        v = vdupq_n_f32(t_);
#else
        for (int i = 0; i < SIMD_WIDTH; ++i) {
            xyzw[i] = t_;
        }
#endif
    }

    float operator[](int i) const { return xyzw[i]; }
    float& operator[](int i) { return xyzw[i]; }
} FVec4;


typedef struct alignas(16) IVec4 {
    union {
        ivec4 v;
#ifdef HAS_NEON
        uint32x4_t uv;
#endif
        alignas(16) int xyzw[4];
    };
#ifndef NO_SSE
    IVec4(ivec4 i4) : v(i4) {}
#else
    IVec4(ivec4 i4) {
      for (int i = 0; i < SIMD_WIDTH; ++i) {
          v[i] = i4[i];
      }
    }
#endif
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
    
  int operator[](int i) const { return xyzw[i]; }
  int& operator[](int i) { return xyzw[i]; }

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


inline FVec4 simd_add(FVec4 a, FVec4 b) {
    FVec4 result;
#ifdef HAS_AVX
    result.v = _mm256_add_ps(a.v, b.v);
#elif defined(HAS_SSE)
    result.v = _mm_add_ps(a.v, b.v);
#elif defined(HAS_NEON)
    result.v = vaddq_f32(a.v, b.v);
#else
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.v[i] = a.v[i] + b.v[i];
    }
#endif
    return result;
}

inline FVec4 simd_div(FVec4 a, FVec4 b) {
    FVec4 result;
#ifdef HAS_AVX
    result.v = _mm256_div_ps(a.v, b.v);
#elif defined(HAS_SSE)
    result.v = _mm_div_ps(a.v, b.v);
#elif defined(HAS_NEON)
    // NEON doesn't have a divide instruction, so we approximate
    float32x4_t reciprocal = vrecpeq_f32(b.v);
    reciprocal = vmulq_f32(vrecpsq_f32(b.v, reciprocal), reciprocal);
    reciprocal = vmulq_f32(vrecpsq_f32(b.v, reciprocal), reciprocal);
    result.v = vmulq_f32(a.v, reciprocal);
#else
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.v[i] = a.v[i] / b.v[i];
    }
#endif
    return result;
}


inline FVec4 simd_set(float e0, float e1, float e2, float e3) {
    FVec4 result;
#ifdef HAS_SSE
    result.v = _mm_set_ps(e3, e2, e1, e0);  // Note the reverse order
#elif defined(HAS_NEON)
    float32_t data[4] = { e0, e1, e2, e3 };
    result.v = vld1q_f32(data);
#else
    result.xyzw[0] = e0;
    result.xyzw[1] = e1;
    result.xyzw[2] = e2;
    result.xyzw[3] = e3;
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
    result.v[i] = (a.v[i] <= b.v[i]) ? -1.0f : 0.0f;
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
    IVec4 result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = reinterpret_cast<int&>(mask.xyzw[i]);
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
    result.xyzw[0] = a.xyzw[idx0];
    result.xyzw[1] = a.xyzw[idx1];
    result.xyzw[2] = a.xyzw[idx2];
    result.xyzw[3] = a.xyzw[idx3];
    return result;
}

inline IVec4 simd_shuffle(IVec4 a, int idx0, int idx1, int idx2, int idx3) {
    IVec4 result;
#ifdef HAS_SSE
    // For SSE2, use _mm_set_epi32 to set elements
    result.v = _mm_set_epi32(
        a.xyzw[idx3], // Note: _mm_set_epi32 sets elements in order [3,2,1,0]
        a.xyzw[idx2],
        a.xyzw[idx1],
        a.xyzw[idx0]
    );
#elif defined(HAS_NEON)
    // For NEON, create an array and load it into the vector
    int32_t values[4] = {
        a.xyzw[idx0],
        a.xyzw[idx1],
        a.xyzw[idx2],
        a.xyzw[idx3]
    };
    result.v = vld1q_s32(values);
#else
    // Scalar fallback
    result.xyzw[0] = a.xyzw[idx0];
    result.xyzw[1] = a.xyzw[idx1];
    result.xyzw[2] = a.xyzw[idx2];
    result.xyzw[3] = a.xyzw[idx3];
#endif
    return result;
}

inline IVec4 simd_cmpneq(IVec4 a, IVec4 b) {
#ifdef HAS_AVX
    IVec4 result;
    __m256i cmp = _mm256_cmpeq_epi32(a.v, b.v);
    result.v = _mm256_xor_si256(cmp, _mm256_set1_epi32(-1)); // Invert bits
    return result;
#elif defined(HAS_SSE)
    IVec4 result;
    __m128i cmp = _mm_cmpeq_epi32(a.v, b.v);
    result.v = _mm_xor_si128(cmp, _mm_set1_epi32(-1)); // Invert bits
    return result;
#elif defined(HAS_NEON)
    IVec4 result;
    result.v = vmvnq_s32(vceqq_s32(a.v, b.v));
    return result;
#else
    IVec4 result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = (a.xyzw[i] != b.xyzw[i]) ? -1 : 0;
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

static inline float sgn_local(float val) {
  return (float(0) < val) - (val < float(0));
}

inline int sgn_local(int val) {
  return (0 < val) - (val < 0);
}

inline FVec4 simd_sgn(const FVec4& a) {
  FVec4 result;
#ifdef HAS_SSE2
    __m128i zero = _mm_setzero_si128();
    __m128i one = _mm_set1_epi32(1);
    __m128i neg_one = _mm_set1_epi32(-1);

    __m128i gt_mask = _mm_cmpgt_epi32(a.v, zero);      // x > 0
    __m128i lt_mask = _mm_cmpgt_epi32(zero, a.v);      // x < 0

    __m128i pos_result = _mm_and_si128(gt_mask, one);     // 1 where x > 0
    __m128i neg_result = _mm_and_si128(lt_mask, neg_one); // -1 where x < 0

    result.v = _mm_add_epi32(pos_result, neg_result);

#elif defined(HAS_NEON)
    int32x4_t zero = vdupq_n_s32(0);
    int32x4_t one = vdupq_n_s32(1);
    int32x4_t neg_one = vdupq_n_s32(-1);

    uint32x4_t gt_mask = vcgtq_s32(a.v, zero); // x > 0
    uint32x4_t lt_mask = vcltq_s32(a.v, zero); // x < 0

    int32x4_t pos_result = vandq_s32(vreinterpretq_s32_u32(gt_mask), one);
    int32x4_t neg_result = vandq_s32(vreinterpretq_s32_u32(lt_mask), neg_one);

    result.v = vaddq_s32(pos_result, neg_result);

#else
    // Fallback to scalar implementation
    result.xyzw[0] = sgn_local(a.xyzw[0]);
    result.xyzw[1] = sgn_local(a.xyzw[1]);
    result.xyzw[2] = sgn_local(a.xyzw[2]);
    result.xyzw[3] = sgn_local(a.xyzw[3]);
#endif
    return result;
}


inline IVec4 simd_sgn(const IVec4& a) {
  IVec4 result;
#ifdef HAS_SSE2
    __m128i zero = _mm_setzero_si128();
    __m128i one = _mm_set1_epi32(1);
    __m128i neg_one = _mm_set1_epi32(-1);

    __m128i gt_mask = _mm_cmpgt_epi32(p.e.v, zero);      // x > 0
    __m128i lt_mask = _mm_cmpgt_epi32(zero, p.e.v);      // x < 0

    __m128i pos_result = _mm_and_si128(gt_mask, one);     // 1 where x > 0
    __m128i neg_result = _mm_and_si128(lt_mask, neg_one); // -1 where x < 0

    result.v = _mm_add_epi32(pos_result, neg_result);

#elif defined(HAS_NEON)
    int32x4_t zero = vdupq_n_s32(0);
    int32x4_t one = vdupq_n_s32(1);
    int32x4_t neg_one = vdupq_n_s32(-1);

    uint32x4_t gt_mask = vcgtq_s32(a.v, zero); // x > 0
    uint32x4_t lt_mask = vcltq_s32(a.v, zero); // x < 0

    int32x4_t pos_result = vandq_s32(vreinterpretq_s32_u32(gt_mask), one);
    int32x4_t neg_result = vandq_s32(vreinterpretq_s32_u32(lt_mask), neg_one);

    result.v = vaddq_s32(pos_result, neg_result);

#else
    // Fallback to scalar implementation
    result.xyzw[0] = sgn_local(a.xyzw[0]);
    result.xyzw[1] = sgn_local(a.xyzw[1]);
    result.xyzw[2] = sgn_local(a.xyzw[2]);
    result.xyzw[3] = sgn_local(a.xyzw[3]);
#endif
    return result;
}

inline IVec4 simd_and(IVec4 a, IVec4 b) {
#ifdef HAS_AVX2
    IVec4 result;
    result.v = _mm256_and_si256(a.v, b.v);
    return result;
#elif defined(HAS_SSE)
    IVec4 result;
    result.v = _mm_and_si128(a.v, b.v);
    return result;
#elif defined(HAS_NEON)
    IVec4 result;
    result.v = vandq_s32(a.v, b.v);
    return result;
#else
    IVec4 result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = a.xyzw[i] & b.xyzw[i];
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
#elif defined(HAS_SSE)
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

inline FVec4 simd_abs(const FVec4& b) {
    FVec4 result;
#ifdef HAS_SSE
    result.v = _mm_andnot_ps(_mm_set1_ps(-0.0f), b.v);
#elif defined(HAS_NEON)
    result.v = vabsq_f32(b.v);
#else
#ifdef RAY_FLOAT_AS_DOUBLE
    result.e[0] = std::fabs(b.e[0]);
    result.e[1] = std::fabs(b.e[1]);
    result.e[2] = std::fabs(b.e[2]);
    result.e[3] = 0.0f;
#else
    result.e[0] = std::fabsf(b.e[0]);
    result.e[1] = std::fabsf(b.e[1]);
    result.e[2] = std::fabsf(b.e[2]);
    result.e[3] = 0.0f;
#endif
#endif
    return result;
}

inline float simd_dot(FVec4 a, FVec4 b) {
#ifdef HAS_AVX
    // For AVX, use _mm256_dp_ps or manually compute
    __m256 mul = _mm256_mul_ps(a.v, b.v);
    // Sum the elements manually
    __m128 hi = _mm256_extractf128_ps(mul, 1); // Upper 128 bits
    __m128 lo = _mm256_castps256_ps128(mul);   // Lower 128 bits
    __m128 sum = _mm_add_ps(lo, hi);           // Sum lower and upper parts
    sum = _mm_hadd_ps(sum, sum);
    sum = _mm_hadd_ps(sum, sum);
    return _mm_cvtss_f32(sum);
#elif defined(HAS_SSE)
    // For SSE, use dot product intrinsics or manual horizontal addition
    __m128 mul = _mm_mul_ps(a.v, b.v);
#if defined(__SSE4_1__)
    // If SSE4.1 is available, use _mm_dp_ps
    __m128 sum = _mm_dp_ps(a.v, b.v, 0x71); // 0x71 masks elements and specifies which elements to sum
    return _mm_cvtss_f32(sum);
#else
    // Without SSE4.1, perform horizontal adds
    __m128 shuf = _mm_shuffle_ps(mul, mul, _MM_SHUFFLE(2, 3, 0, 1)); // Shuffle
    __m128 sums = _mm_add_ps(mul, shuf); // Add pairs
    shuf = _mm_movehl_ps(shuf, sums);    // Move high part
    sums = _mm_add_ss(sums, shuf);       // Final sum
    return _mm_cvtss_f32(sums);
#endif
#elif defined(HAS_NEON)
    // For NEON, use vmlaq_f32 and vaddvq_f32 (vaddvq_f32 requires ARMv8)
#if defined(__aarch64__)
    // ARMv8 and above
    float32x4_t mul = vmulq_f32(a.v, b.v);
    float result = vaddvq_f32(mul); // Sum across vector
    return result;
#else
    // For ARMv7, sum manually
    float32x4_t mul = vmulq_f32(a.v, b.v);
    float32x2_t sum_pair = vadd_f32(vget_low_f32(mul), vget_high_f32(mul));
    float32x2_t sum = vpadd_f32(sum_pair, sum_pair);
    return vget_lane_f32(sum, 0);
#endif
#else
    // Scalar fallback
    float result = 0.0f;
    for (int i = 0; i < 4; ++i) { // Only sum the first 3 elements
        result += a.xyzw[i] * b.xyzw[i];
    }
    return result;
#endif
}


inline FVec4 simd_cross(FVec4 a, FVec4 b) {
    FVec4 result;
#ifdef HAS_SSE
    // For SSE, use shuffle and multiply operations
    __m128 a_yzx = _mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(3, 0, 2, 1)); // (a.z, a.x, a.y)
    __m128 b_yzx = _mm_shuffle_ps(b.v, b.v, _MM_SHUFFLE(3, 0, 2, 1)); // (b.z, b.x, b.y)

    __m128 c = _mm_sub_ps(
        _mm_mul_ps(a.v, b_yzx),
        _mm_mul_ps(a_yzx, b.v)
    );

    result.v = _mm_shuffle_ps(c, c, _MM_SHUFFLE(3, 0, 2, 1)); // Shuffle back to (c.x, c.y, c.z)
#elif defined(HAS_NEON)
    // Allocate aligned memory for arrays
    alignas(16) float a_values[4];
    alignas(16) float b_values[4];

    // Store vectors into arrays
    vst1q_f32(a_values, a.v);
    vst1q_f32(b_values, b.v);

     // ARM code adaptation
    // Load elements into NEON registers
    float32x2_t vec_a_1 = vld1_f32(a_values + 1); // [a_y, a_z]
    float32x2_t vec_a_2 = vld1_f32(a_values);     // [a_x, a_y]

    float32x2_t vec_b_1 = vld1_f32(b_values + 1); // [b_y, b_z]
    float32x2_t vec_b_2 = vld1_f32(b_values);     // [b_x, b_y]

    // Combine to create 128-bit vectors
    float32x4_t vec_a = vcombine_f32(vec_a_2, vec_a_1); // [a_x, a_y, a_y, a_z]
    float32x4_t vec_b = vcombine_f32(vec_b_2, vec_b_1); // [b_x, b_y, b_y, b_z]

    // Rotate vectors
    float32x4_t vec_a_rot = vextq_f32(vec_a, vec_a, 1); // [a_y, a_y, a_z, a_x]
    float32x4_t vec_b_rot = vextq_f32(vec_b, vec_b, 1); // [b_y, b_y, b_z, b_x]

    // Compute cross product
    float32x4_t prod = vmulq_f32(vec_a, vec_b_rot); // [a_x * b_y, a_y * b_y, a_y * b_z, a_z * b_x]
    // [a_x * b_y - a_y * b_x, a_y * b_y - a_y * b_y, a_y * b_z - a_z * b_y, a_z * b_x - a_x * b_z]
    prod = vmlsq_f32(prod, vec_a_rot, vec_b); 
    // So this is [c_z, 0, c_x, c_y]
    prod = vextq_f32(prod, prod, 2);
    
    result.v = prod;
#else
    // Scalar fallback
    result.xyzw[0] = DifferenceOfProducts(a.xyzw[1], b.xyzw[2], a.xyzw[2], b.xyzw[1]);
    result.xyzw[1] = DifferenceOfProducts(a.xyzw[2], b.xyzw[0], a.xyzw[0], b.xyzw[2]);
    result.xyzw[2] = DifferenceOfProducts(a.xyzw[0], b.xyzw[1], a.xyzw[1], b.xyzw[0]);
    result.xyzw[3] = 0.0f;
#endif
    // Ensure the fourth element is zero
    result.xyzw[3] = 0.0f;
    return result;
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


inline IVec4 simd_set1(int value) {
    IVec4 result;
#ifdef HAS_AVX2
    result.v = _mm256_set1_epi32(value);
#elif defined(HAS_SSE)
    result.v = _mm_set1_epi32(value);
#elif defined(HAS_NEON)
    result.v = vdupq_n_s32(value);
#else
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.xyzw[i] = value;
    }
#endif
    return result;
}

inline float simd_squared_length(FVec4 a) {
    float result;
#if defined(HAS_SSE)
    // SSE implementation using _mm intrinsics
    __m128 v = a.v;
    __m128 v_squared = _mm_mul_ps(v, v);          // Element-wise square
#if defined(__SSE3__)
    // Use _mm_hadd_ps for horizontal addition (requires SSE3)
    __m128 sum = _mm_hadd_ps(v_squared, v_squared);
    sum = _mm_hadd_ps(sum, sum);
    result = _mm_cvtss_f32(sum);                 // Extract the result
#else
    // Manual horizontal addition for SSE2
    __m128 shuf = _mm_shuffle_ps(v_squared, v_squared, _MM_SHUFFLE(2, 3, 0, 1)); // Shuffle
    __m128 sums = _mm_add_ps(v_squared, shuf);  // Add pairs
    shuf = _mm_movehl_ps(shuf, sums);           // Move high part
    sums = _mm_add_ss(sums, shuf);              // Add final pair
    result = _mm_cvtss_f32(sums);               // Extract the result
#endif                          // Extract the result
#elif defined(HAS_NEON)
    // NEON implementation using ARM intrinsics
    float32x4_t v = a.v;
    float32x4_t v_squared = vmulq_f32(v, v);
    // Sum all elements of v_squared
    #if defined(__aarch64__)
        // AArch64 provides vaddvq_f32 to sum across vector
        result = vaddvq_f32(v_squared);
    #else
        // For ARMv7, sum manually
        float32x2_t sum_pair = vadd_f32(vget_low_f32(v_squared), vget_high_f32(v_squared));
        result = vget_lane_f32(vpadd_f32(sum_pair, sum_pair), 0);
    #endif
#else
    // Scalar fallback
    result = a.xyzw[0] * a.xyzw[0] +
             a.xyzw[1] * a.xyzw[1] +
             a.xyzw[2] * a.xyzw[2] +
             a.xyzw[3] * a.xyzw[3];
#endif
    return result;
}

inline IVec4 simd_add(IVec4 a, IVec4 b) {
    IVec4 result;
#ifdef HAS_AVX2
    result.v = _mm256_add_epi32(a.v, b.v);
#elif defined(HAS_SSE)
    result.v = _mm_add_epi32(a.v, b.v);
#elif defined(HAS_NEON)
    result.v = vaddq_s32(a.v, b.v);
#else
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.xyzw[i] = a.xyzw[i] + b.xyzw[i];
    }
#endif
    return result;
}

inline IVec4 simd_sub(IVec4 a, IVec4 b) {
    IVec4 result;
#ifdef HAS_AVX2
    result.v = _mm256_sub_epi32(a.v, b.v);
#elif defined(HAS_SSE)
    result.v = _mm_sub_epi32(a.v, b.v);
#elif defined(HAS_NEON)
    result.v = vsubq_s32(a.v, b.v);
#else
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.xyzw[i] = a.xyzw[i] - b.xyzw[i];
    }
#endif
    return result;
}

inline IVec4 simd_mul(IVec4 a, IVec4 b) {
    IVec4 result;
#ifdef HAS_AVX2
    result.v = _mm256_mullo_epi32(a.v, b.v);
#elif defined(HAS_SSE4_1)
    result.v = _mm_mullo_epi32(a.v, b.v);
#elif defined(HAS_NEON)
    result.v = vmulq_s32(a.v, b.v);
#else
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.xyzw[i] = a.xyzw[i] * b.xyzw[i];
    }
#endif
    return result;
}

// Integer division is not directly supported in SIMD; we can approximate or use scalar code
inline IVec4 simd_div(IVec4 a, IVec4 b) {
    IVec4 result;
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.xyzw[i] = a.xyzw[i] / b.xyzw[i];
    }
    return result;
}

// inline IVec4 simd_blend_int(SimdMask mask, IVec4 a, IVec4 b) {
//     IVec4 mask_int = simd_cast_to_int(mask);

//     #ifdef HAS_SSE
//         IVec4 result;
//         result.v = _mm_or_si128(_mm_and_si128(mask_int.v, a.v), _mm_andnot_si128(mask_int.v, b.v));
//         return result;
//     #elif defined(HAS_NEON)
//         IVec4 result;
//         result.v = vbslq_s32(vreinterpretq_u32_s32(mask_int.v), a.v, b.v);
//         return result;
//     #else
//         // Fallback for non-SIMD
//         IVec4 result;
//         for (int i = 0; i < 4; ++i) {
//             result.xyzw[i] = mask_int.xyzw[i] ? a.xyzw[i] : b.xyzw[i];
//         }
//         return result;
//     #endif
// }

inline IVec4 simd_blend_int(SimdMask mask, IVec4 a, IVec4 b) {
#ifdef HAS_AVX
    IVec4 result;
    result.v = _mm256_blendv_epi8(b.v, a.v, _mm256_castps_si256(mask.v));
    return result;
#elif defined(HAS_SSE)
    IVec4 result;
    result.v = _mm_or_si128(_mm_and_si128(_mm_castps_si128(mask.v), a.v),
                            _mm_andnot_si128(_mm_castps_si128(mask.v), b.v));
    return result;
#elif defined(HAS_NEON)
    IVec4 result;
    result.v = vbslq_s32(vreinterpretq_u32_f32(mask.v), a.v, b.v);
    return result;
#else
    IVec4 result;
    for (int i = 0; i < 4; ++i) {
        result.xyzw[i] = (reinterpret_cast<uint32_t&>(mask.xyzw[i])) ? a.xyzw[i] : b.xyzw[i];
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
        uint32_t bits = ~reinterpret_cast<uint32_t&>(mask.xyzw[i]);
        result.xyzw[i] = reinterpret_cast<float&>(bits);
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
    __m128i masked = _mm_and_si128(vec.v, _mm_set1_epi32(1)); // Now masked contains 0 or 1 in each 32-bit element

    // Shift left by 31 bits to place the bit in the sign position
    __m128i shifted = _mm_slli_epi32(masked, 31);

    // Cast to float to use _mm_movemask_ps
    __m128 float_vals = _mm_castsi128_ps(shifted);

    // Use _mm_movemask_ps to extract the sign bits
    int mask = _mm_movemask_ps(float_vals); // Extracts the sign bits of each float

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