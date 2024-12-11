#include "Rcpp.h"
#include <testthat.h>
#include "vec3.h"
#include "point3.h"
#include "normal.h"
#include "simd.h"
#include "vectypes.h"  // Your provided header with inline functions

context("Conversion functions") {
    test_that("[convert_to_point3(vec3f)]") {
        vec3f v; v.e[0] = 1.0f; v.e[1] = 2.0f; v.e[2] = 3.0f;
        point3f p = convert_to_point3(v);

        expect_true(p.x() == Approx(1.0f));
        expect_true(p.y() == Approx(2.0f));
        expect_true(p.z() == Approx(3.0f));
    }

    test_that("[convert_to_point3(normal3f)]") {
        normal3f n; n.e[0] = -1.0f; n.e[1] = -2.0f; n.e[2] = -3.0f;
        point3f p = convert_to_point3(n);

        expect_true(p.x() == Approx(-1.0f));
        expect_true(p.y() == Approx(-2.0f));
        expect_true(p.z() == Approx(-3.0f));
    }

    test_that("[convert_to_vec3(point3f)]") {
        point3f p; p.e[0] = 4.0f; p.e[1] = 5.0f; p.e[2] = 6.0f;
        vec3f v = convert_to_vec3(p);

        expect_true(v.x() == Approx(4.0f));
        expect_true(v.y() == Approx(5.0f));
        expect_true(v.z() == Approx(6.0f));
    }

    test_that("[convert_to_vec3(normal3f)]") {
        normal3f n; n.e[0] = 7.0f; n.e[1] = 8.0f; n.e[2] = 9.0f;
        vec3f v = convert_to_vec3(n);

        expect_true(v.x() == Approx(7.0f));
        expect_true(v.y() == Approx(8.0f));
        expect_true(v.z() == Approx(9.0f));
    }

    test_that("[convert_to_normal3(vec3f)]") {
        vec3f v; v.e[0] = 10.0f; v.e[1] = -10.0f; v.e[2] = 20.0f;
        normal3f n = convert_to_normal3(v);

        expect_true(n.x() == Approx(10.0f));
        expect_true(n.y() == Approx(-10.0f));
        expect_true(n.z() == Approx(20.0f));
    }

    test_that("[convert_to_normal3(point3f)]") {
        point3f p; p.e[0] = -4.5f; p.e[1] = 0.5f; p.e[2] = 2.5f;
        normal3f n = convert_to_normal3(p);

        expect_true(n.x() == Approx(-4.5f));
        expect_true(n.y() == Approx(0.5f));
        expect_true(n.z() == Approx(2.5f));
    }
}

context("Arithmetic operators with point3, vec3, normal3") {
    test_that("[point3 * vec3]") {
        point3f p; p.e[0] = 1.0f; p.e[1] = 2.0f; p.e[2] = 3.0f;
        vec3f v; v.e[0] = 2.0f; v.e[1] = -1.0f; v.e[2] = 0.5f;

        point3f result = p * v;
        expect_true(result.x() == Approx(2.0f));   // 1*2
        expect_true(result.y() == Approx(-2.0f));  // 2*(-1)
        expect_true(result.z() == Approx(1.5f));   // 3*0.5
    }

    test_that("[point3 * normal3f]") {
        point3f p; p.e[0] = 4.0f; p.e[1] = -2.0f; p.e[2] = 10.0f;
        normal3f n; n.e[0] = 0.5f; n.e[1] = 2.0f; n.e[2] = -1.0f;

        point3f result = p * n;
        expect_true(result.x() == Approx(2.0f));   // 4*0.5
        expect_true(result.y() == Approx(-4.0f));  // -2*2
        expect_true(result.z() == Approx(-10.0f)); // 10*(-1)
    }

    test_that("[normal3f * vec3]") {
        normal3f n; n.e[0] = 1.0f; n.e[1] = 2.0f; n.e[2] = 3.0f;
        vec3f v; v.e[0] = 2.0f; v.e[1] = 2.0f; v.e[2] = -1.0f;

        normal3f result = n * v;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(4.0f));
        expect_true(result.z() == Approx(-3.0f));
    }


    test_that("[point3f *= point3f]") {
        point3f p1; p1.e[0] = 2.0f; p1.e[1] = 2.0f; p1.e[2] = 2.0f;
        point3f p2; p2.e[0] = 3.0f; p2.e[1] = -1.0f; p2.e[2] = 4.0f;

        point3f result = p1 *= p2;
        expect_true(result.x() == Approx(6.0f));   // 2*3
        expect_true(result.y() == Approx(-2.0f));  // 2*(-1)
        expect_true(result.z() == Approx(8.0f));   // 2*4
    }

    test_that("[point3f *= normal3f]") {
        point3f p; p.e[0] = 1.0f; p.e[1] = -2.0f; p.e[2] = 5.0f;
        normal3f n; n.e[0] = 2.0f; n.e[1] = 2.0f; n.e[2] = 0.5f;

        point3f result = p *= n;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(-4.0f));
        expect_true(result.z() == Approx(2.5f));
    }
}

context("vec3<Float> arithmetic operators") {
    test_that("[vec3 + vec3]") {
        vec3f v1, v2;
        v1.e[0] = 1.0f; v1.e[1] = 2.0f; v1.e[2] = 3.0f;
        v2.e[0] = 4.0f; v2.e[1] = -1.0f; v2.e[2] = 0.0f;

        vec3f result = v1 + v2;
        expect_true(result.x() == Approx(5.0f));
        expect_true(result.y() == Approx(1.0f));
        expect_true(result.z() == Approx(3.0f));
    }

    test_that("[vec3 - vec3]") {
        vec3f v1, v2;
        v1.e[0] = 5.0f; v1.e[1] = 5.0f; v1.e[2] = 5.0f;
        v2.e[0] = 3.0f; v2.e[1] = -2.0f; v2.e[2] = 2.0f;

        vec3f result = v1 - v2;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(7.0f));
        expect_true(result.z() == Approx(3.0f));
    }

    test_that("[vec3 * vec3]") {
        vec3f v1, v2;
        v1.e[0] = 2.0f; v1.e[1] = 3.0f; v1.e[2] = 4.0f;
        v2.e[0] = 1.0f; v2.e[1] = 2.0f; v2.e[2] = -1.0f;

        vec3f result = v1 * v2;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(6.0f));
        expect_true(result.z() == Approx(-4.0f));
    }

    test_that("[vec3 / vec3]") {
        vec3f v1, v2;
        v1.e[0] = 10.0f; v1.e[1] = 20.0f; v1.e[2] = -4.0f;
        v2.e[0] = 2.0f;  v2.e[1] = 5.0f;  v2.e[2] = -2.0f;

        vec3f result = v1 / v2;
        expect_true(result.x() == Approx(5.0f));
        expect_true(result.y() == Approx(4.0f));
        expect_true(result.z() == Approx(2.0f));
    }

    test_that("[vec3 + Float]") {
        vec3f v;
        v.e[0] = 1.0f; v.e[1] = 2.0f; v.e[2] = 3.0f;

        vec3f result = v + 1.0f;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(3.0f));
        expect_true(result.z() == Approx(4.0f));
    }

    test_that("[vec3 - Float]") {
        vec3f v;
        v.e[0] = 4.0f; v.e[1] = 6.0f; v.e[2] = 8.0f;

        vec3f result = v - 2.0f;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(4.0f));
        expect_true(result.z() == Approx(6.0f));
    }

    test_that("[Float * vec3]") {
        vec3f v;
        v.e[0] = -1.0f; v.e[1] = 0.5f; v.e[2] = 2.0f;

        vec3f result = 2.0f * v;
        expect_true(result.x() == Approx(-2.0f));
        expect_true(result.y() == Approx(1.0f));
        expect_true(result.z() == Approx(4.0f));
    }

    test_that("[vec3 * Float]") {
        vec3f v;
        v.e[0] = 1.0f; v.e[1] = -2.0f; v.e[2] = 3.0f;

        vec3f result = v * 3.0f;
        expect_true(result.x() == Approx(3.0f));
        expect_true(result.y() == Approx(-6.0f));
        expect_true(result.z() == Approx(9.0f));
    }

    test_that("[vec3 / Float]") {
        vec3f v;
        v.e[0] = 6.0f; v.e[1] = -6.0f; v.e[2] = 12.0f;

        vec3f result = v / 2.0f;
        expect_true(result.x() == Approx(3.0f));
        expect_true(result.y() == Approx(-3.0f));
        expect_true(result.z() == Approx(6.0f));
    }

    test_that("[Float / vec3]") {
        vec3f v;
        v.e[0] = 2.0f; v.e[1] = -2.0f; v.e[2] = 4.0f;

        vec3f result = 8.0f / v;
        expect_true(result.x() == Approx(4.0f));
        expect_true(result.y() == Approx(-4.0f));
        expect_true(result.z() == Approx(2.0f));
    }

    test_that("[cross(vec3, vec3)]") {
        vec3f v1, v2;
        v1.e[0] = 1.0f; v1.e[1] = 2.0f; v1.e[2] = 3.0f;
        v2.e[0] = 4.0f; v2.e[1] = 5.0f; v2.e[2] = 6.0f;

        vec3f result = cross(v1, v2);
        // Cross product: (2*6 - 3*5, 3*4 - 1*6, 1*5 - 2*4) = (12-15,12-6,5-8) = (-3,6,-3)
        expect_true(result.x() == Approx(-3.0f));
        expect_true(result.y() == Approx(6.0f));
        expect_true(result.z() == Approx(-3.0f));
    }
}

context("point3<Float> arithmetic operators") {
    test_that("[point3 + point3]") {
        point3f p1, p2;
        p1.e[0] = 1.0f; p1.e[1] = 2.0f; p1.e[2] = 3.0f;
        p2.e[0] = -1.0f; p2.e[1] = 0.0f; p2.e[2] = 5.0f;

        point3f result = p1 + p2;
        expect_true(result.x() == Approx(0.0f));
        expect_true(result.y() == Approx(2.0f));
        expect_true(result.z() == Approx(8.0f));
    }

    test_that("[point3 + vec3]") {
        point3f p; p.e[0] = 2.0f; p.e[1] = 2.0f; p.e[2] = 2.0f;
        vec3f v; v.e[0] = 1.0f; v.e[1] = -3.0f; v.e[2] = 4.0f;

        point3f result = p + v;
        expect_true(result.x() == Approx(3.0f));
        expect_true(result.y() == Approx(-1.0f));
        expect_true(result.z() == Approx(6.0f));
    }

    test_that("[vec3 + point3]") {
        vec3f v; v.e[0] = 3.0f; v.e[1] = 1.0f; v.e[2] = 0.0f;
        point3f p; p.e[0] = -1.0f; p.e[1] = 2.0f; p.e[2] = 5.0f;

        point3f result = v + p;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(3.0f));
        expect_true(result.z() == Approx(5.0f));
    }

    test_that("[point3 += vec3]") {
        point3f p; p.e[0] = 0.0f; p.e[1] = 1.0f; p.e[2] = 2.0f;
        vec3f v; v.e[0] = 1.0f; v.e[1] = 1.0f; v.e[2] = 1.0f;

        p += v;
        expect_true(p.x() == Approx(1.0f));
        expect_true(p.y() == Approx(2.0f));
        expect_true(p.z() == Approx(3.0f));
    }

    test_that("[point3 -= vec3]") {
        point3f p; p.e[0] = 4.0f; p.e[1] = 4.0f; p.e[2] = 4.0f;
        vec3f v; v.e[0] = 1.0f; v.e[1] = 2.0f; v.e[2] = 3.0f;

        p -= v;
        expect_true(p.x() == Approx(3.0f));
        expect_true(p.y() == Approx(2.0f));
        expect_true(p.z() == Approx(1.0f));
    }

    test_that("[point3 - point3 = vec3]") {
        point3f p1, p2;
        p1.e[0] = 5.0f; p1.e[1] = 3.0f; p1.e[2] = 2.0f;
        p2.e[0] = 2.0f; p2.e[1] = 1.0f; p2.e[2] = -1.0f;

        vec3f v = p1 - p2;
        expect_true(v.x() == Approx(3.0f));
        expect_true(v.y() == Approx(2.0f));
        expect_true(v.z() == Approx(3.0f));
    }

    test_that("[point3 - vec3]") {
        point3f p; p.e[0] = 1.0f; p.e[1] = 1.0f; p.e[2] = 1.0f;
        vec3f v; v.e[0] = -1.0f; v.e[1] = 2.0f; v.e[2] = 0.5f;

        point3f result = p - v;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(-1.0f));
        expect_true(result.z() == Approx(0.5f));
    }

    test_that("[point3 * point3]") {
        point3f p1, p2;
        p1.e[0] = 2.0f; p1.e[1] = 2.0f; p1.e[2] = 2.0f;
        p2.e[0] = 3.0f; p2.e[1] = 1.0f; p2.e[2] = 2.0f;

        point3f result = p1 * p2;
        expect_true(result.x() == Approx(6.0f));
        expect_true(result.y() == Approx(2.0f));
        expect_true(result.z() == Approx(4.0f));
    }

    test_that("[point3 / point3]") {
        point3f p1, p2;
        p1.e[0] = 6.0f; p1.e[1] = 9.0f; p1.e[2] = -4.0f;
        p2.e[0] = 3.0f; p2.e[1] = 3.0f; p2.e[2] = -2.0f;

        point3f result = p1 / p2;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(3.0f));
        expect_true(result.z() == Approx(2.0f));
    }

    test_that("[point3 + Float]") {
        point3f p; p.e[0] = 1.0f; p.e[1] = 0.0f; p.e[2] = -1.0f;

        point3f result = p + 2.0f;
        expect_true(result.x() == Approx(3.0f));
        expect_true(result.y() == Approx(2.0f));
        expect_true(result.z() == Approx(1.0f));
    }

    test_that("[point3 - Float]") {
        point3f p; p.e[0] = 4.0f; p.e[1] = 2.0f; p.e[2] = 0.0f;

        point3f result = p - 1.0f;
        expect_true(result.x() == Approx(3.0f));
        expect_true(result.y() == Approx(1.0f));
        expect_true(result.z() == Approx(-1.0f));
    }

    test_that("[Float * point3]") {
        point3f p; p.e[0] = 3.0f; p.e[1] = -3.0f; p.e[2] = 0.5f;

        point3f result = 2.0f * p;
        expect_true(result.x() == Approx(6.0f));
        expect_true(result.y() == Approx(-6.0f));
        expect_true(result.z() == Approx(1.0f));
    }

    test_that("[point3 * Float]") {
        point3f p; p.e[0] = -1.0f; p.e[1] = 4.0f; p.e[2] = 2.0f;

        point3f result = p * 3.0f;
        expect_true(result.x() == Approx(-3.0f));
        expect_true(result.y() == Approx(12.0f));
        expect_true(result.z() == Approx(6.0f));
    }

    test_that("[point3 / Float]") {
        point3f p; p.e[0] = 9.0f; p.e[1] = -6.0f; p.e[2] = 3.0f;

        point3f result = p / 3.0f;
        expect_true(result.x() == Approx(3.0f));
        expect_true(result.y() == Approx(-2.0f));
        expect_true(result.z() == Approx(1.0f));
    }

    test_that("[Float / point3]") {
        point3f p; p.e[0] = 3.0f; p.e[1] = -1.0f; p.e[2] = 1.0f;

        point3f result = 6.0f / p;
        expect_true(result.x() == Approx(2.0f));
        expect_true(result.y() == Approx(-6.0f));
        expect_true(result.z() == Approx(6.0f));
    }

    test_that("[cross(point3, point3)]") {
        point3f p1, p2;
        p1.e[0] = 1.0f; p1.e[1] = 2.0f; p1.e[2] = 3.0f;
        p2.e[0] = 4.0f; p2.e[1] = 5.0f; p2.e[2] = 6.0f;

        point3f result = cross(p1, p2);
        // Same as vec cross: (-3,6,-3)
        expect_true(result.x() == Approx(-3.0f));
        expect_true(result.y() == Approx(6.0f));
        expect_true(result.z() == Approx(-3.0f));
    }

    test_that("[Distance and DistanceSquared]") {
        point3f p1; p1.e[0] = 0.0f; p1.e[1] = 0.0f; p1.e[2] = 0.0f;
        point3f p2; p2.e[0] = 3.0f; p2.e[1] = 4.0f; p2.e[2] = 0.0f;

        Float dist = Distance(p1, p2);
        expect_true(dist == Approx(5.0f));

        Float dist_sq = DistanceSquared(p1, p2);
        expect_true(dist_sq == Approx(25.0f));
    }

    test_that("[Lerp between points]") {
        point3f p0, p1;
        p0.e[0] = 0.0f; p0.e[1] = 0.0f; p0.e[2] = 0.0f;
        p1.e[0] = 10.0f; p1.e[1] = 10.0f; p1.e[2] = 10.0f;

        Float t = 0.5f;
        point3f result = Lerp(t, p0, p1);
        expect_true(result.x() == Approx(5.0f));
        expect_true(result.y() == Approx(5.0f));
        expect_true(result.z() == Approx(5.0f));
    }

    test_that("[MinComponent, MaxComponent, MaxDimension]") {
        vec3f v;
        v.e[0] = -2.0f; v.e[1] = 3.0f; v.e[2] = 1.0f;

        Float min_c = MinComponent(v);
        Float max_c = MaxComponent(v);
        int max_d = MaxDimension(v);

        expect_true(min_c == Approx(-2.0f));
        expect_true(max_c == Approx(3.0f));
        expect_true(max_d == 1); // y is the largest dimension
    }

    test_that("[Min and Max for vec3]") {
        vec3f p1, p2;
        p1.e[0] = 1.0f; p1.e[1] = 4.0f; p1.e[2] = -1.0f;
        p2.e[0] = 0.0f; p2.e[1] = 5.0f; p2.e[2] = -2.0f;

        vec3f mn = Min(p1, p2);
        vec3f mx = Max(p1, p2);

        expect_true(mn.x() == Approx(0.0f));
        expect_true(mn.y() == Approx(4.0f));
        expect_true(mn.z() == Approx(-2.0f));

        expect_true(mx.x() == Approx(1.0f));
        expect_true(mx.y() == Approx(5.0f));
        expect_true(mx.z() == Approx(-1.0f));
    }

    test_that("[Permute and PermuteInPlace]") {
        vec3f v;
        v.e[0] = 1.0f; v.e[1] = 2.0f; v.e[2] = 3.0f;

        vec3f perm = Permute(v, 2, 1, 0);
        expect_true(perm.x() == Approx(3.0f));
        expect_true(perm.y() == Approx(2.0f));
        expect_true(perm.z() == Approx(1.0f));

        PermuteInPlace(v, 1, 2, 0);
        expect_true(v.x() == Approx(2.0f));
        expect_true(v.y() == Approx(3.0f));
        expect_true(v.z() == Approx(1.0f));
    }
}
