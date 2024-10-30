#include "simd.h"

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