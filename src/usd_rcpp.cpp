#include <Rcpp.h>
#include <libopenusd.h>
#include <cmath>
#include <algorithm>
#include <climits>
#include <iomanip>
#include <locale>
#include <memory>
#include <sstream>
#include <fstream>
#include <filesystem>

namespace {
// USDA uses UTF-8 quoted strings. Escape delimiters/control characters before
// writing any scene-supplied names or paths into the layer grammar.
void quote(std::ostream &out, const char *text) {
    out << '"';
    for (const unsigned char *p = reinterpret_cast<const unsigned char *>(text); *p; ++p) {
        switch (*p) {
            case '"': out << "\\\""; break;
            case '\\': out << "\\\\"; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default:
                if (*p < 32 || *p == 127) {
                    const char *hex = "0123456789abcdef";
                    out << "\\x" << hex[*p >> 4] << hex[*p & 15];
                } else out << *p;
        }
    }
    out << '"';
}

void pack_value(std::ostream &out, SEXP value, int depth, bool template_only = false);

void pack_fields(std::ostream &out, SEXP fields, int depth, bool templates) {
    Rcpp::List values(fields);
    Rcpp::CharacterVector names = values.names();
    if (names.size() != values.size()) Rcpp::stop("USD dictionary fields must be named.");
    out << "{\n";
    for (R_xlen_t i = 0; i < values.size(); ++i) {
        out << "dictionary ";
        quote(out, Rf_translateCharUTF8(names[i]));
        out << " = ";
        pack_value(out, values[i], depth + 1, templates);
        out << '\n';
    }
    out << '}';
}

void pack_value(std::ostream &out, SEXP value, int depth, bool template_only) {
    // Column templates omit data while retaining attributes for reconstruction.
    if (depth > 512) Rcpp::stop("USD metadata nesting exceeds 512 levels.");
    const char *kind = nullptr, *usd_type = nullptr;
    switch (TYPEOF(value)) {
        case NILSXP: kind = "NULL"; break;
        case VECSXP: kind = "list"; break;
        case LGLSXP: kind = "logical"; usd_type = "bool[]"; break;
        case INTSXP: kind = "integer"; usd_type = "int[]"; break;
        case REALSXP: kind = "double"; usd_type = "double[]"; break;
        case STRSXP: kind = "character"; usd_type = "string[]"; break;
        case RAWSXP: kind = "raw"; usd_type = "uchar[]"; break;
        default: Rcpp::stop("USD export cannot preserve an R %s; use a file-backed or numeric descriptor.", Rf_type2char(TYPEOF(value)));
    }
    const R_xlen_t size = (template_only || value == R_NilValue) ? 0 : XLENGTH(value);
    out << "{\nstring type = \"" << kind << "\"\n";
    if (TYPEOF(value) == VECSXP) {
        out << "dictionary items = {\n";
        for (R_xlen_t i = 0; i < size; ++i) {
            if (i % 65536 == 0) Rcpp::checkUserInterrupt();
            out << "dictionary n" << std::setfill('0') << std::setw(8) << i << " = ";
            pack_value(out, VECTOR_ELT(value, i), depth + 1);
            out << '\n';
        }
        out << "}\n";
    } else if (value != R_NilValue) {
        // Native arrays preserve precision and avoid one R allocation per
        // element. NA/NaN/Inf live in explicit masks, matching schema version 1.
        std::vector<R_xlen_t> missing;
        std::vector<const char *> codes;
        out << usd_type << " data = [";
        for (R_xlen_t i = 0; i < size; ++i) {
            if (i % 65536 == 0) Rcpp::checkUserInterrupt();
            if (i) out << ", ";
            const char *code = nullptr;
            switch (TYPEOF(value)) {
                case LGLSXP:
                    if (LOGICAL(value)[i] == NA_LOGICAL) code = "NA";
                    out << (LOGICAL(value)[i] == TRUE ? "true" : "false"); break;
                case INTSXP:
                    if (INTEGER(value)[i] == NA_INTEGER) code = "NA";
                    out << (code ? 0 : INTEGER(value)[i]); break;
                case REALSXP: {
                    double number = REAL(value)[i];
                    if (R_IsNA(number)) code = "NA";
                    else if (std::isnan(number)) code = "NaN";
                    else if (std::isinf(number)) code = number > 0 ? "Inf" : "-Inf";
                    out << (code ? 0 : number); break;
                }
                case STRSXP:
                    if (STRING_ELT(value, i) == NA_STRING) code = "NA";
                    quote(out, code ? "" : Rf_translateCharUTF8(STRING_ELT(value, i))); break;
                case RAWSXP: out << static_cast<unsigned int>(RAW(value)[i]); break;
            }
            if (code) { missing.push_back(i + 1); codes.push_back(code); }
        }
        out << "]\n";
        if (!missing.empty()) {
            out << "int[] missingIndices = [";
            for (size_t i = 0; i < missing.size(); ++i) {
                if (missing[i] > INT_MAX) Rcpp::stop("USD missing-value index exceeds schema version 1 limits.");
                if (i) out << ", ";
                out << missing[i];
            }
            out << "]\nstring[] missingCodes = [";
            for (size_t i = 0; i < codes.size(); ++i) {
                if (i) out << ", ";
                quote(out, codes[i]);
            }
            out << "]\n";
        }
    }
    Rcpp::RObject object(value);
    const auto attribute_names = object.attributeNames();
    if (!attribute_names.empty()) {
        out << "dictionary attributes = {\n";
        for (const auto &name : attribute_names) {
            out << "dictionary ";
            quote(out, name.c_str());
            out << " = ";
            Rcpp::RObject attribute = object.attr(name);
            pack_value(out, attribute, depth + 1);
            out << '\n';
        }
        out << "}\n";
    }
    out << '}';
}
} // namespace

// [[Rcpp::export(rng = false)]]
SEXP usd_inspect_native(SEXP filename, SEXP selectors, SEXP properties) {
    return libopenusd_inspect(filename, selectors, properties);
}

namespace {
void write_literal(std::ostream &out, SEXP value, const std::string &type) {
    bool array = type.size() >= 2 && type.substr(type.size() - 2) == "[]";
    std::string base = array ? type.substr(0, type.size() - 2) : type;
    bool strings = base == "string" || base == "token" || base == "asset";
    int width = 1;
    if (base == "matrix4d") width = 16;
    else if (base == "quatf" || base == "quath" || base == "quatd") width = 4;
    else if (base == "float2" || base == "texCoord2f" || base == "int2") width = 2;
    else if (base == "float3" || base == "double3" || base == "point3f" || base == "normal3f" || base == "color3f" || base == "vector3f") width = 3;
    R_xlen_t count = XLENGTH(value) / width;
    if (XLENGTH(value) % width || (!array && count != 1)) Rcpp::stop("Invalid USD %s value length.", type);
    bool matrix = Rf_isMatrix(value);
    const int storage = TYPEOF(value);
    if (!strings && storage != REALSXP && storage != INTSXP && storage != LGLSXP) {
        Rcpp::stop("USD numeric properties require numeric data.");
    }
    if (array) out << '[';
    for (R_xlen_t row = 0; row < count; ++row) {
        if (row % 65536 == 0) Rcpp::checkUserInterrupt();
        if (row) out << ", ";
        if (width > 1) out << '(';
        for (int column = 0; column < width; ++column) {
            if (column) out << ", ";
            if (base == "matrix4d" && column % 4 == 0) out << '(';
            R_xlen_t index = matrix && base != "matrix4d" ? row + count * column : row * width + column;
            if (base == "matrix4d" && matrix) index = (column / 4) + 4 * (column % 4);
            if (strings) {
                if (TYPEOF(value) != STRSXP || STRING_ELT(value, index) == NA_STRING) Rcpp::stop("Invalid USD string value.");
                const char *text = Rf_translateCharUTF8(STRING_ELT(value, index));
                if (base == "asset") {
                    std::string asset(text);
                    if (asset.find("@@@") != std::string::npos || asset.find('\n') != std::string::npos) Rcpp::stop("Asset path contains an unsupported USD delimiter.");
                    out << "@@@" << asset << "@@@";
                } else quote(out, text);
            } else {
                double number = storage == REALSXP ? REAL(value)[index] :
                    storage == INTSXP ? INTEGER(value)[index] : LOGICAL(value)[index];
                if (storage != REALSXP && number == NA_INTEGER) {
                    Rcpp::stop("Missing value in USD %s interchange property.", type);
                }
                if (!std::isfinite(number)) Rcpp::stop("Nonfinite value in USD %s interchange property.", type);
                if (base == "bool") out << (number ? "true" : "false");
                else out << number;
            }
            if (base == "matrix4d" && column % 4 == 3) out << ')';
        }
        if (width > 1) out << ')';
    }
    if (array) out << ']';

}

// A document is a tree of small syntax strings and references to R values.
// Typed nodes are recognized before generic lists so their data is never
// flattened or formatted into an intermediate R character string.
void write_document(std::ostream &out, SEXP fragment, unsigned depth = 0) {
    if (depth > 512) Rcpp::stop("USD document nesting exceeds 512 levels.");
    if (Rf_isNull(fragment)) return;

    if (Rf_inherits(fragment, "rayrender_usd_metadata")) {
        Rcpp::List node(fragment);
        if (Rcpp::as<bool>(node["fields"])) {
            pack_fields(out, node["value"], 0, Rcpp::as<bool>(node["templates"]));
        } else {
            pack_value(out, node["value"], 0);
        }
        return;
    }
    if (Rf_inherits(fragment, "rayrender_usd_literal")) {
        Rcpp::List node(fragment);
        write_literal(out, node["value"], Rcpp::as<std::string>(node["type"]));
        return;
    }
    if (TYPEOF(fragment) == VECSXP) {
        for (R_xlen_t i = 0; i < XLENGTH(fragment); ++i) {
            Rcpp::checkUserInterrupt();
            write_document(out, VECTOR_ELT(fragment, i), depth + 1);
        }
        return;
    }
    if (TYPEOF(fragment) == STRSXP) {
        for (R_xlen_t i = 0; i < XLENGTH(fragment); ++i) {
            if (STRING_ELT(fragment, i) == NA_STRING) {
                Rcpp::stop("Missing text in USD document.");
            }
            out << Rf_translateCharUTF8(STRING_ELT(fragment, i));
        }
        return;
    }
    Rcpp::stop("Unexpected USD document fragment type.");
}
} // namespace

// Small scalar formatting remains useful for names and camera transforms.
// Arrays and metadata use write_document() directly instead.
// [[Rcpp::export(rng = false)]]
std::string usd_literal_native(SEXP value, std::string type) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    write_literal(out, value, type);
    return out.str();
}

// [[Rcpp::export(rng = false)]]
void usd_stream_native(SEXP document, std::string filename) {
    // Filesystem paths are UTF-8, including on Windows. Stream exceptions report
    // write/flush failures, such as a full disk, before the stage is published.
    std::ofstream out;
    out.exceptions(std::ios::badbit | std::ios::failbit);
#if defined(__cpp_char8_t)
    const std::filesystem::path path(std::u8string(filename.begin(), filename.end()));
#else
    const auto path = std::filesystem::u8path(filename);
#endif
    out.open(path, std::ios::binary | std::ios::trunc);
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    write_document(out, document);
    out.close();
}

// [[Rcpp::export(rng = false)]]
void usd_export_file_native(std::string source, std::string destination) {
    const auto *api = libopenusd_get_api();
    using Stage = std::unique_ptr<libopenusd_stage, void (*)(libopenusd_stage *)>;
    Stage stage(api->stage_open(source.c_str()), api->stage_release);
    if (!stage) Rcpp::stop("Cannot read streamed USD stage: %s", api->last_error());

    // USDA is already written: opening it validates the file. Binary targets
    // are converted by OpenUSD without passing their source text back through R.
    if (source == destination) return;
    if (!api->stage_export(stage.get(), destination.c_str(), 0)) {
        Rcpp::stop("%s", api->last_error());
    }
    stage.reset();
    Stage restored(api->stage_open(destination.c_str()), api->stage_release);
    if (!restored) Rcpp::stop("Cannot reopen exported stage: %s", api->last_error());
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
