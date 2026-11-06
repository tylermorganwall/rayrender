#include <Rcpp.h>
#include <libopenusd.h>
#include <cmath>
#include <algorithm>

// [[Rcpp::export(rng = false)]]
SEXP usd_inspect_native(SEXP filename, SEXP selectors, SEXP properties) {
    return libopenusd_inspect(filename, selectors, properties);
}

// [[Rcpp::export(rng = false)]]
void usd_author_native(SEXP document, SEXP filename) {
    libopenusd_write_stage(document, filename);
}

// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix usd_instance_quaternions_native(Rcpp::NumericMatrix angles, Rcpp::IntegerVector order) {
    if (angles.ncol() != 3 || order.size() != 3) Rcpp::stop("Invalid instance rotation dimensions.");
    bool used[3] = {false, false, false};
    for (int axis : order) {
        if (axis < 1 || axis > 3 || used[axis - 1]) Rcpp::stop("Rotation order must be a permutation of 1:3.");
        used[axis - 1] = true;
    }
    Rcpp::NumericMatrix result(angles.nrow(), 4);
    const double radians = std::acos(-1.) / 360.;
    for (int row = 0; row < angles.nrow(); ++row) {
        double q[4] = {1, 0, 0, 0};
        // Row-vector rotations compose in reverse Hamilton-product order.
        for (int axis : order) {
            double theta = angles(row, axis - 1) * radians;
            if (!std::isfinite(theta)) Rcpp::stop("Nonfinite instance rotation.");
            double r[4] = {std::cos(theta), 0, 0, 0};
            r[axis] = std::sin(theta);
            double next[4] = {
                r[0]*q[0] - r[1]*q[1] - r[2]*q[2] - r[3]*q[3],
                r[0]*q[1] + q[0]*r[1] + r[2]*q[3] - r[3]*q[2],
                r[0]*q[2] + q[0]*r[2] + r[3]*q[1] - r[1]*q[3],
                r[0]*q[3] + q[0]*r[3] + r[1]*q[2] - r[2]*q[1]
            };
            std::copy(next, next + 4, q);
        }
        for (int column = 0; column < 4; ++column) result(row, column) = q[column];
    }
    return result;
}
