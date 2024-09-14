#ifndef SIMDH
#define SIMDH

#include <cstdint>
#include <algorithm>
#include <stdio.h>

// SIMD vector size (4 for SSE, 8 for AVX)
#ifdef __AVX__
#define SIMD_WIDTH 8
#else
#define SIMD_WIDTH 4
#endif



#ifdef __ARM_NEON
    #include <arm_neon.h>
#elif defined(__SSE__)
    #include <xmmintrin.h>
#else
    #error "No SIMD support available"
#endif


typedef struct FVec4 {
public:
    union {
#if defined(__x86_64__)
        __m128 m128;
#elif defined(__aarch64__)
        float32x4_t f32x4;
#endif
        float xyzw[4];
        // struct {
        //   float x;
        //   float y;
        //   float z;
        //   float w;
        // };
    };
#if defined(__x86_64__)
    FVec4(__m128 f4) : m128(f4) {}
#elif defined(__aarch64__)
    FVec4(float32x4_t f4) : f32x4(f4) {}
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
#if defined(__x86_64__)
        __m128i m128i;
#elif defined(__aarch64__)
        int32x4_t i32x4;
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
#if defined(__x86_64__)
    IVec4(__m128i m128i_) : m128i(m128i_) {};
#elif defined(__aarch64__)
    IVec4(int32x4_t i32x4_) : i32x4(i32x4_) {};
#endif

    int operator[](int i) const { return xyzw[i]; }
    int& operator[](int i) { return xyzw[i]; }
} IVec4;


// SIMD vector type
struct SimdFloat {
#ifdef __AVX__
  __m256 v;
#elif defined(__SSE__)
  __m128 v;
#elif defined(__ARM_NEON)
  float32x4_t v;
#else
  float v[SIMD_WIDTH];
#endif
};


// SIMD operations
inline SimdFloat simd_load(const float* ptr) {
  SimdFloat result;
#ifdef __AVX__
  result.v = _mm256_load_ps(ptr);
#elif defined(__SSE__)
  result.v = _mm_load_ps(ptr);
#elif defined(__ARM_NEON)
  result.v = vld1q_f32(ptr);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = ptr[i];
  }
#endif
  return result;
}

typedef SimdFloat SimdMask;

inline SimdFloat simd_set1(float value) {
  SimdFloat result;
#ifdef __AVX__
  result.v = _mm256_set1_ps(value);
#elif defined(__SSE__)
  result.v = _mm_set1_ps(value);
#elif defined(__ARM_NEON)
  result.v = vdupq_n_f32(value);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = value;
  }
#endif
  return result;
}

inline SimdFloat simd_sub(SimdFloat a, SimdFloat b) {
  SimdFloat result;
#ifdef __AVX__
  result.v = _mm256_sub_ps(a.v, b.v);
#elif defined(__SSE__)
  result.v = _mm_sub_ps(a.v, b.v);
#elif defined(__ARM_NEON)
  result.v = vsubq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = a.v[i] - b.v[i];
  }
#endif
  return result;
}


inline SimdFloat simd_mul(SimdFloat a, SimdFloat b) {
  SimdFloat result;
#ifdef __AVX__
  result.v = _mm256_mul_ps(a.v, b.v);
#elif defined(__SSE__)
  result.v = _mm_mul_ps(a.v, b.v);
#elif defined(__ARM_NEON)
  result.v = vmulq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = a.v[i] * b.v[i];
  }
#endif
  return result;
}


inline SimdFloat simd_min(SimdFloat a, SimdFloat b) {
  SimdFloat result;
#ifdef __AVX__
  result.v = _mm256_min_ps(a.v, b.v);
#elif defined(__SSE__)
  result.v = _mm_min_ps(a.v, b.v);
#elif defined(__ARM_NEON)
  result.v = vminq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = std::min(a.v[i], b.v[i]);
  }
#endif
  return result;
}


inline SimdFloat simd_max(SimdFloat a, SimdFloat b) {
  SimdFloat result;
#ifdef __AVX__
  result.v = _mm256_max_ps(a.v, b.v);
#elif defined(__SSE__)
  result.v = _mm_max_ps(a.v, b.v);
#elif defined(__ARM_NEON)
  result.v = vmaxq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = std::max(a.v[i], b.v[i]);
  }
#endif
  return result;
}


inline SimdMask simd_less_equal(SimdFloat a, SimdFloat b) {
  SimdMask result;
#ifdef __AVX__
  result.v = _mm256_cmp_ps(a.v, b.v, _CMP_LE_OQ);
#elif defined(__SSE__)
  result.v = _mm_cmple_ps(a.v, b.v);
#elif defined(__ARM_NEON)
  result.v = vcleq_f32(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = (a.v[i] <= b.v[i]) ? 0xFFFFFFFF : 0;
  }
#endif
  return result;
}


inline bool simd_any_true(SimdMask mask) {
#ifdef __AVX__
  return _mm256_movemask_ps(mask.v) != 0;
#elif defined(__SSE__)
  return _mm_movemask_ps(mask.v) != 0;
#elif defined(__ARM_NEON)
  return vmaxvq_u32(vreinterpretq_u32_f32(mask.v)) != 0;
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    if (mask.v[i] != 0) return true;
  }
  return false;
#endif
}

inline IVec4 simd_cast_to_int(FVec4 mask) {
#ifdef __ARM_NEON
    // NEON: reinterpret the result of comparison from float to int
    return vreinterpretq_s32_u32(vreinterpretq_u32_f32(mask.f32x4));
#elif defined(__SSE__)
    // SSE: cast floating-point mask to integer mask
    return _mm_castps_si128(mask.v);
#else
    // Fallback for non-SIMD
    IVec4 result;
    for (int i = 0; i < SIMD_WIDTH; ++i) {
        result.v[i] = (mask.v[i] != 0.0f) ? 0xFFFFFFFF : 0;
    }
    return result;
#endif
}

#endif