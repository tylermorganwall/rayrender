#ifndef SIMDH
#define SIMDH

#include <cstdint>
#include <algorithm>

// SIMD vector size (4 for SSE, 8 for AVX)
#ifdef __AVX__
#define SIMD_WIDTH 8
#else
#define SIMD_WIDTH 4
#endif

// SIMD vector type
struct SimdFloat {
#ifdef __AVX__
  __m256 v;
#elif defined(__SSE__)
  __m128 v;
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
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = ptr[i];
  }
#endif
  return result;
}


inline SimdFloat simd_add(SimdFloat a, SimdFloat b) {
  SimdFloat result;
#ifdef __AVX__
  result.v = _mm256_add_ps(a.v, b.v);
#elif defined(__SSE__)
  result.v = _mm_add_ps(a.v, b.v);
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    result.v[i] = a.v[i] + b.v[i];
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
#else
  for (int i = 0; i < SIMD_WIDTH; ++i) {
    if (mask.v[i] != 0) return true;
  }
  return false;
#endif
}

#endif