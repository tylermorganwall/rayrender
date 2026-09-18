#include "src/math/simd.h"
#include <cassert>
#include <climits>
int main() {
  const int values[] = {0, 1, -1, INT_MIN, INT_MAX, INT_MIN + 1, INT_MAX - 1};
  for (int x : values) {
    IVec4 signs = simd_sgn(IVec4(x, 0, -1, 1));
    assert(signs[0] == ((x > 0) - (x < 0)));
    assert(signs[1] == 0 && signs[2] == -1 && signs[3] == 1);
  }
  IVec4 product = simd_mul(IVec4(-2, 0, 1, 7), IVec4(3, 5, -4, 6));
  assert(product[0] == -6 && product[1] == 0 && product[2] == -4 && product[3] == 42);
}
