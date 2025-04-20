// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include "../util/print.h"

#include "../util/check.h"

#include <charconv>   // for std::to_chars
#include <string>
#include <limits>     // for std::numeric_limits
#include <cmath>      // for std::isinf, std::isnan

// #include <double-conversion/double-conversion.h>

namespace pbrt {

namespace detail {

// std::string FloatToString(float v) {
//     double_conversion::DoubleToStringConverter doubleConverter(
//         double_conversion::DoubleToStringConverter::NO_FLAGS, "Inf", "NaN", 'e',
//         -6 /* decimal_in_shortest_low */, 9 /* decimal_in_shortest_high */,
//         5 /* max_leading_padding_zeroes_in_precision_mode */,
//         5 /*  max_trailing_padding_zeroes_in_precision_mode */);
//     char buf[64];
//     double_conversion::StringBuilder result(buf, PBRT_ARRAYSIZE(buf));
//     doubleConverter.ToShortestSingle(v, &result);
//     int length = result.position();
//     return std::string(buf, length);
// }

// std::string DoubleToString(double v) {
//     double_conversion::DoubleToStringConverter doubleConverter(
//         double_conversion::DoubleToStringConverter::NO_FLAGS, "Inf", "NaN", 'e',
//         -6 /* decimal_in_shortest_low */, 9 /* decimal_in_shortest_high */,
//         5 /* max_leading_padding_zeroes_in_precision_mode */,
//         5 /*  max_trailing_padding_zeroes_in_precision_mode */);
//     char buf[64];
//     double_conversion::StringBuilder result(buf, PBRT_ARRAYSIZE(buf));
//     doubleConverter.ToShortest(v, &result);
//     int length = result.position();
//     return std::string(buf, length);
// }

std::string FloatToString(float v) {
    // Check special cases first
    if (std::isnan(v)) {
        return "NaN";
    }
    if (std::isinf(v)) {
        // signbit(v) can check for negative vs positive infinity
        // but a direct comparison is fine too:
        return (v > 0) ? "Inf" : "-Inf";
    }

    // Choose a buffer size. For float round-trip safety, 9 or 10 digits is typical.
    // We'll aim for the round-trip guarantee: "max_digits10".
    // That’s enough digits so that the string -> float -> string yields the same float.
    char buf[64];
    auto result = std::to_chars(
        buf,
        buf + sizeof(buf),
        v,
        std::chars_format::general,
        std::numeric_limits<float>::max_digits10
    );

    if (result.ec != std::errc{}) {
        // Handle error if needed; e.g. return an empty string or some fallback.
        return {};
    }

    // `result.ptr` is where the conversion wrote up to; use that to size the std::string.
    return std::string(buf, result.ptr);
}

std::string DoubleToString(double v) {
    if (std::isnan(v)) {
        return "NaN";
    }
    if (std::isinf(v)) {
        return (v > 0) ? "Inf" : "-Inf";
    }

    char buf[64];
    auto result = std::to_chars(
        buf,
        buf + sizeof(buf),
        v,
        std::chars_format::general,
        std::numeric_limits<double>::max_digits10
    );

    if (result.ec != std::errc{}) {
        // Handle error
        return {};
    }
    return std::string(buf, result.ptr);
}

void stringPrintfRecursive(std::string *s, const char *fmt) {
    const char *c = fmt;
    // No args left; make sure there aren't any extra formatting
    // specifiers.
    while (*c) {
        if (*c == '%') {
            if (c[1] != '%')
                LOG_FATAL("Not enough optional values passed to Printf.");
            ++c;
        }
        *s += *c++;
    }
}

// 1. Copy from fmt to *s, up to the next formatting directive.
// 2. Advance fmt past the next formatting directive and return the
//    formatting directive as a string.
std::string copyToFormatString(const char **fmt_ptr, std::string *s) {
    const char *&fmt = *fmt_ptr;
    while (*fmt) {
        if (*fmt != '%') {
            *s += *fmt;
            ++fmt;
        } else if (fmt[1] == '%') {
            // "%%"; let it pass through
            *s += '%';
            *s += '%';
            fmt += 2;
        } else
            // fmt is at the start of a formatting directive.
            break;
    }

    std::string nextFmt;
    while (*fmt) {
        char c = *fmt;
        nextFmt += c;
        ++fmt;
        // Is it a conversion specifier?
        if (c == 'd' || c == 'i' || c == 'o' || c == 'u' || c == 'x' || c == 'e' ||
            c == 'E' || c == 'f' || c == 'F' || c == 'g' || c == 'G' || c == 'a' ||
            c == 'A' || c == 'c' || c == 'C' || c == 's' || c == 'S' || c == 'p')
            break;
    }

    return nextFmt;
}

}  // namespace detail

}  // namespace pbrt
