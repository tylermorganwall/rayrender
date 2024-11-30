#include "simd.h"
#include "Rcpp.h"

// const static SimdMask first_mask = simd_setmask(true, false, true, false);
// const static SimdMask second_mask = simd_setmask(true, true, false, false);
const static SimdMask first_mask = simd_setmask(false, true, false, true);
const static SimdMask second_mask = simd_setmask(false, false, true, true);

inline IVec4 sort_simd_4_floats(FVec4 values) {
    values.v[0] = pack_index(values.v[0], 0);
    values.v[1] = pack_index(values.v[1], 1);
    values.v[2] = pack_index(values.v[2], 2);
    values.v[3] = pack_index(values.v[3], 3);

    // Level 1: Swap pairs and compare
    {
        FVec4 tmp_values = simd_swap_pairs(values);

        SimdMask lower_mask = simd_cmplt(tmp_values, values);
        SimdMask higher_mask = simd_cmplt(values, tmp_values);

        // Blend values and indices based on the mask
        FVec4 low_values = simd_blend(lower_mask, tmp_values, values);
        FVec4 high_values = simd_blend(higher_mask, tmp_values, values);
        values = simd_blend(first_mask, low_values, high_values);
    }

    // Level 2: Reverse and compare
    {
        FVec4 tmp_values = simd_reverse(values);

        SimdMask lower_mask = simd_cmplt(tmp_values, values);
        SimdMask higher_mask = simd_cmplt(values, tmp_values);

        // Blend values and indices based on the mask
        FVec4 low_values = simd_blend(lower_mask, tmp_values, values);
        FVec4 high_values = simd_blend(higher_mask, tmp_values, values);
        values = simd_blend(second_mask, low_values, high_values);
    }

    // Level 1: Swap pairs and compare
    {
        FVec4 tmp_values = simd_swap_pairs(values);

        SimdMask lower_mask = simd_cmplt(tmp_values, values);
        SimdMask higher_mask = simd_cmplt(values, tmp_values);

        FVec4 low_values = simd_blend(lower_mask, tmp_values, values);
        FVec4 high_values = simd_blend(higher_mask, tmp_values, values);
        values = simd_blend(first_mask, low_values, high_values);
    }
    return (IVec4(extract_index(values[0]),
                  extract_index(values[1]), 
                  extract_index(values[2]),
                  extract_index(values[3])));
}

#include <testthat.h>

context("simd_load loads values correctly") {
  test_that("[simd_load]") {
    float values[SIMD_WIDTH] = {1.0f, 2.0f, 3.0f, 4.0f};
    FVec4 vec = simd_load(values);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        expect_true(vec.xyzw[i] == Approx(values[i]));
    }
  }
}


context("simd_sub subtracts vectors correctly") {
  test_that("[simd_sub]") {
    float a_values[SIMD_WIDTH] = {5.0f, 6.0f, 7.0f, 8.0f};
    float b_values[SIMD_WIDTH] = {1.0f, 2.0f, 3.0f, 4.0f};

    FVec4 a = simd_load(a_values);
    FVec4 b = simd_load(b_values);

    FVec4 result = simd_sub(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        expect_true(result.xyzw[i] == Approx(a_values[i] - b_values[i]));
    }
  }
}

context("simd_div performs element-wise division correctly") {
    test_that("[simd_div]") {
        float a_values[SIMD_WIDTH] = {10.0f, 20.0f, 30.0f, 40.0f};
        float b_values[SIMD_WIDTH] = {2.0f, 4.0f, 5.0f, 8.0f};

        FVec4 a = simd_load(a_values);
        FVec4 b = simd_load(b_values);

        FVec4 result = simd_div(a, b);

        for (int i = 0; i < SIMD_WIDTH; ++i) {
            float expected = a_values[i] / b_values[i];
            expect_true(fabs(result.xyzw[i] - expected) < 1e-6f);
        }
    }
}


context("simd_mul multiplies vectors correctly") {
  test_that("[simd_mul]") {
    float a_values[SIMD_WIDTH] = {2.0f, 3.0f, 4.0f, 5.0f};
    float b_values[SIMD_WIDTH] = {6.0f, 7.0f, 8.0f, 9.0f};

    FVec4 a = simd_load(a_values);
    FVec4 b = simd_load(b_values);

    FVec4 result = simd_mul(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        expect_true(result.xyzw[i] == Approx(a_values[i] * b_values[i]));
    }
  }
}


context("simd_min computes minimum correctly") {
  test_that("[simd_min]") {
    float a_values[SIMD_WIDTH] = {5.0f, 2.0f, 7.0f, 1.0f};
    float b_values[SIMD_WIDTH] = {3.0f, 4.0f, 6.0f, 2.0f};

    FVec4 a = simd_load(a_values);
    FVec4 b = simd_load(b_values);

    FVec4 result = simd_min(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        expect_true(result.xyzw[i] == Approx(std::min(a_values[i], b_values[i])));
    }
  }
}

context("simd_max computes maximum correctly") {
  test_that("[simd_max]") {
    float a_values[SIMD_WIDTH] = {5.0f, 2.0f, 7.0f, 1.0f};
    float b_values[SIMD_WIDTH] = {3.0f, 4.0f, 6.0f, 2.0f};

    FVec4 a = simd_load(a_values);
    FVec4 b = simd_load(b_values);

    FVec4 result = simd_max(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        expect_true(result.xyzw[i] == Approx(std::max(a_values[i], b_values[i])));
    }
  }
}


context("simd_less_equal compares vectors correctly") {
  test_that("[simd_less_equal]") {
    float a_values[SIMD_WIDTH] = {1.0f, 4.0f, 6.0f, 8.0f};
    float b_values[SIMD_WIDTH] = {2.0f, 3.0f, 6.0f, 7.0f};

    FVec4 a = simd_load(a_values);
    FVec4 b = simd_load(b_values);

    FVec4 mask = simd_less_equal(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        bool expected = a_values[i] <= b_values[i];
        expect_true((mask.xyzw[i] != 0.0f) == expected);
    }
  }
}


context("simd_any_true detects any true values") {
  test_that("[simd_any_true]") {
    float mask_values[SIMD_WIDTH] = {-1.0f, 0.0f, 0.0f, 0.0f};
    SimdMask mask = simd_load(mask_values);

    expect_true(simd_any_true(mask) == true);

    mask_values[0] = 0.0f;
    mask = simd_load(mask_values);

    expect_true(simd_any_true(mask) == false);
  }
}


context("simd_cast_float_to_int casts floats to ints") {
  test_that("[simd_cast_float_to_int]") {
    float values[SIMD_WIDTH] = {0.0f, -1.0f, 3.14f, -2.71f};
    FVec4 float_vec = simd_load(values);

    IVec4 int_vec = simd_cast_float_to_int(float_vec);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        int expected = *reinterpret_cast<int*>(&values[i]);
        expect_true(int_vec.xyzw[i] == expected);
    }
  }
}


context("simd_swap_pairs swaps pairs correctly") {
  test_that("[simd_swap_pairs]") {
    float values[SIMD_WIDTH] = {1.0f, 2.0f, 3.0f, 4.0f};
    FVec4 vec = simd_load(values);

    FVec4 result = simd_swap_pairs(vec);

    expect_true(result.xyzw[0] == Approx(2.0f));
    expect_true(result.xyzw[1] == Approx(1.0f));
    expect_true(result.xyzw[2] == Approx(4.0f));
    expect_true(result.xyzw[3] == Approx(3.0f));
  }
}

context("simd_reverse reverses elements correctly") {
  test_that("[simd_reverse]") {
    float values[SIMD_WIDTH] = {1.0f, 2.0f, 3.0f, 4.0f};
    FVec4 vec = simd_load(values);

    FVec4 result = simd_reverse(vec);

    expect_true(result.xyzw[0] == Approx(4.0f));
    expect_true(result.xyzw[1] == Approx(3.0f));
    expect_true(result.xyzw[2] == Approx(2.0f));
    expect_true(result.xyzw[3] == Approx(1.0f));
  }
}

context("simd_blend blends vectors correctly") {
  test_that("[simd_blend]") {
    float a_values[SIMD_WIDTH] = {1.0f, 2.0f, 3.0f, 4.0f};
    float b_values[SIMD_WIDTH] = {5.0f, 6.0f, 7.0f, 8.0f};
    bool mask_bools[SIMD_WIDTH] = {true, false, true, false};

    FVec4 a = simd_load(a_values);
    FVec4 b = simd_load(b_values);
    SimdMask mask = simd_setmask(mask_bools[0], mask_bools[1], mask_bools[2], mask_bools[3]);

    FVec4 result = simd_blend(mask, a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        float expected = mask_bools[i] ? a_values[i] : b_values[i];
        expect_true(result.xyzw[i] == Approx(expected));
    }
  }
}

context("simd_cmpge compares vectors correctly") {
  test_that("[simd_cmpge]") {
    float a_values[SIMD_WIDTH] = {5.0f, 2.0f, 7.0f, 1.0f};
    float b_values[SIMD_WIDTH] = {3.0f, 4.0f, 7.0f, 2.0f};

    FVec4 a = simd_load(a_values);
    FVec4 b = simd_load(b_values);

    SimdMask mask = simd_cmpge(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        bool expected = a_values[i] >= b_values[i];
        expect_true((mask.xyzw[i] != 0.0f) == expected);
    }
  }
}

context("simd_cmpneq compares integers correctly") {
  test_that("[simd_cmpneq]") {
    int a_values[SIMD_WIDTH] = {1, 2, 3, 4};
    int b_values[SIMD_WIDTH] = {1, 3, 3, 5};

    IVec4 a = IVec4(a_values[0], a_values[1], a_values[2], a_values[3]);
    IVec4 b = IVec4(b_values[0], b_values[1], b_values[2], b_values[3]);

    IVec4 result = simd_cmpneq(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        int expected = (a_values[i] != b_values[i]) ? -1 : 0;
        expect_true(result.xyzw[i] == expected);
    }
  }
}

context("simd_and performs bitwise AND correctly") {
  test_that("[simd_and]") {
    unsigned a_values[SIMD_WIDTH] = {0xFFFFFFFF, 0x00000000, 0xAAAAAAAA, 0x55555555};
    unsigned b_values[SIMD_WIDTH] = {0xFFFFFFFF, 0xFFFFFFFF, 0x55555555, 0xAAAAAAAA};

    IVec4 a = IVec4(a_values[0], a_values[1], a_values[2], a_values[3]);
    IVec4 b = IVec4(b_values[0], b_values[1], b_values[2], b_values[3]);

    IVec4 result = simd_and(a, b);

    for (int i = 0; i < SIMD_WIDTH; ++i) {
        int expected = a_values[i] & b_values[i];
        expect_true(result.xyzw[i] == expected);
    }
  }
}

context("simd_extract_hitmask extracts mask correctly") {
  test_that("[simd_extract_hitmask]") {
    int values[SIMD_WIDTH] = {1, 0, 1, 0}; // LSBs: 1, 0, 1, 0
    IVec4 vec = IVec4(values[0], values[1], values[2], values[3]);

    int mask = simd_extract_hitmask(vec);

    int expected_mask = ((values[0] & 1) << 0) |
                        ((values[1] & 1) << 1) |
                        ((values[2] & 1) << 2) |
                        ((values[3] & 1) << 3);

    expect_true(mask == expected_mask);
  }
}

context("simd_dot computes dot product correctly") {
    test_that("[simd_dot]") {
        float a_values[SIMD_WIDTH] = {1.0f, 3.0f, -5.0f, 0.0f};
        float b_values[SIMD_WIDTH] = {4.0f, -2.0f, -1.0f, 0.0f};

        FVec4 a = simd_load(a_values);
        FVec4 b = simd_load(b_values);

        float result = simd_dot(a, b);

        float expected = a_values[0] * b_values[0] + a_values[1] * b_values[1] + a_values[2] * b_values[2];

        expect_true(fabs(result - expected) < 1e-6f);
    }
}

context("simd_cross computes cross product correctly") {
  test_that("[simd_cross]") {
      float a_values[SIMD_WIDTH] = {1.0f, 2.0f, 3.0f, 0.0f};
      float b_values[SIMD_WIDTH] = {4.0f, 5.0f, 6.0f, 0.0f};

      FVec4 a = simd_load(a_values);
      FVec4 b = simd_load(b_values);

      FVec4 result = simd_cross(a, b);

      float expected_values[3];
      expected_values[0] = a_values[1] * b_values[2] - a_values[2] * b_values[1];
      expected_values[1] = a_values[2] * b_values[0] - a_values[0] * b_values[2];
      expected_values[2] = a_values[0] * b_values[1] - a_values[1] * b_values[0];

      for (int i = 0; i < 3; ++i) {
          expect_true(fabs(result.xyzw[i] - expected_values[i]) < 1e-6f);
      }

      // Ensure the fourth element is zero
      expect_true(result.xyzw[3] == 0.0f);
  }
}

context("simd_add performs vector addition correctly") {
    test_that("[simd_add]") {
        float a_values[SIMD_WIDTH] = {1.0f, 2.0f, 3.0f, 4.0f};
        float b_values[SIMD_WIDTH] = {5.0f, 6.0f, 7.0f, 8.0f};

        FVec4 a = simd_load(a_values);
        FVec4 b = simd_load(b_values);

        FVec4 result = simd_add(a, b);

        for (int i = 0; i < SIMD_WIDTH; ++i) {
            float expected = a_values[i] + b_values[i];
            expect_true(fabs(result.xyzw[i] - expected) < 1e-6f);
        }
    }
}