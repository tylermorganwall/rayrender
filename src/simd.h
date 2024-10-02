#ifndef SIMDH
#define SIMDH

#include <cstdint>
#include <algorithm>
#include <stdio.h>

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

typedef struct FVec4 {
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


typedef struct IVec4 {
    union {
#ifdef HAS_SSE
        __m128i v;
#elif defined(HAS_NEON)
        int32x4_t v;
        uint32x4_t uv;
#endif
        int xyzw[4];
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

#endif