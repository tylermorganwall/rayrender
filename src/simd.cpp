#include "simd.h"

// const static SimdMask first_mask = simd_setmask(true, false, true, false);
// const static SimdMask second_mask = simd_setmask(true, true, false, false);
const static SimdMask first_mask = simd_setmask(false, true, false, true);
const static SimdMask second_mask = simd_setmask(false, false, true, true);

IVec4 sort_simd_4_floats(FVec4 values) {
    IVec4 indices(0,1,2,3);

    // Level 1: Swap pairs and compare
    {
        FVec4 tmp_values = simd_swap_pairs(values);
        IVec4 tmp_indices = simd_swap_pairs(indices);

        SimdMask lower_mask = simd_cmplt(tmp_values, values);
        SimdMask higher_mask = simd_cmplt(values, tmp_values);

        // Blend values and indices based on the mask
        FVec4 low_values = simd_blend(lower_mask, tmp_values, values);
        FVec4 high_values = simd_blend(higher_mask, tmp_values, values);
        values = simd_blend(first_mask, low_values, high_values);

        IVec4 low_indices = simd_blend_int(lower_mask, tmp_indices, indices);
        IVec4 higher_indices = simd_blend_int(higher_mask, tmp_indices, indices);
        indices = simd_blend_int(first_mask, low_indices, higher_indices);
    }

    // Level 2: Reverse and compare
    {
        FVec4 tmp_values = simd_reverse(values);
        IVec4 tmp_indices = simd_reverse(indices);

        SimdMask lower_mask = simd_cmplt(tmp_values, values);
        SimdMask higher_mask = simd_cmplt(values, tmp_values);

        // Blend values and indices based on the mask
        FVec4 low_values = simd_blend(lower_mask, tmp_values, values);
        FVec4 high_values = simd_blend(higher_mask, tmp_values, values);
        values = simd_blend(second_mask, low_values, high_values);

        IVec4 low_indices = simd_blend_int(lower_mask, tmp_indices, indices);
        IVec4 higher_indices = simd_blend_int(higher_mask, tmp_indices, indices);
        indices = simd_blend_int(second_mask, low_indices, higher_indices);
    }

    // Level 1: Swap pairs and compare
    {
        FVec4 tmp_values = simd_swap_pairs(values);
        IVec4 tmp_indices = simd_swap_pairs(indices);

        SimdMask lower_mask = simd_cmplt(tmp_values, values);
        SimdMask higher_mask = simd_cmplt(values, tmp_values);

        FVec4 low_values = simd_blend(lower_mask, tmp_values, values);
        FVec4 high_values = simd_blend(higher_mask, tmp_values, values);
        values = simd_blend(first_mask, low_values, high_values);

        IVec4 low_indices = simd_blend_int(lower_mask, tmp_indices, indices);
        IVec4 higher_indices = simd_blend_int(higher_mask, tmp_indices, indices);
        indices = simd_blend_int(first_mask, low_indices, higher_indices);
    }

    return indices;
}