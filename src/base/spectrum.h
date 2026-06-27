#ifndef RAYRENDER_BASE_SPECTRUM_H
#define RAYRENDER_BASE_SPECTRUM_H

#include "color_types.h"
#include "sampled_spectrum.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace rayrender {
namespace base {

enum class SpectrumType {
  Albedo,
  Illuminant,
  Unbounded
};

enum class SpectrumSemantic {
  ColorMatchingFunction,
  Illuminant,
  OpticalEta,
  OpticalK,
  GlassEta,
  CameraSensorResponse,
  Unbounded
};

enum class SpectrumOrderPolicy {
  RequireSorted,
  Sort
};

enum class SpectrumInterpolationPolicy {
  Linear
};

enum class SpectrumExtrapolationPolicy {
  Zero,
  Constant,
  Error
};

enum class SpectrumValueValidation {
  Finite,
  NonNegative,
  Bounded01
};

struct SpectrumDataPolicy {
  SpectrumOrderPolicy order = SpectrumOrderPolicy::RequireSorted;
  SpectrumInterpolationPolicy interpolation = SpectrumInterpolationPolicy::Linear;
  SpectrumExtrapolationPolicy extrapolation = SpectrumExtrapolationPolicy::Zero;
  SpectrumValueValidation validation = SpectrumValueValidation::NonNegative;
};

struct SpectrumMetadata {
  std::string name;
  SpectrumSemantic semantic = SpectrumSemantic::Unbounded;
  std::string wavelengthUnit = "nm";
  std::string valueUnit;
  std::string normalization;
  std::string sourceCitation;
  std::string sourceLicense;
  std::string assetGenerationScriptVersion;
  std::string sha256;
  std::string storage = "piecewise_linear";
};

inline std::string SpectrumSemanticName(SpectrumSemantic semantic) {
  switch (semantic) {
  case SpectrumSemantic::ColorMatchingFunction:
    return "color_matching_function";
  case SpectrumSemantic::Illuminant:
    return "illuminant";
  case SpectrumSemantic::OpticalEta:
    return "optical_eta";
  case SpectrumSemantic::OpticalK:
    return "optical_k";
  case SpectrumSemantic::GlassEta:
    return "glass_eta";
  case SpectrumSemantic::CameraSensorResponse:
    return "camera_sensor_response";
  case SpectrumSemantic::Unbounded:
    return "unbounded";
  }
  return "unbounded";
}

inline SpectrumSemantic ParseSpectrumSemantic(const std::string& value) {
  if (value == "color_matching_function") {
    return SpectrumSemantic::ColorMatchingFunction;
  }
  if (value == "illuminant") {
    return SpectrumSemantic::Illuminant;
  }
  if (value == "optical_eta") {
    return SpectrumSemantic::OpticalEta;
  }
  if (value == "optical_k") {
    return SpectrumSemantic::OpticalK;
  }
  if (value == "glass_eta") {
    return SpectrumSemantic::GlassEta;
  }
  if (value == "camera_sensor_response") {
    return SpectrumSemantic::CameraSensorResponse;
  }
  if (value == "unbounded") {
    return SpectrumSemantic::Unbounded;
  }
  throw std::invalid_argument("Unknown spectrum semantic: " + value);
}

inline void ValidateSpectrumValue(Float value, SpectrumValueValidation validation) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument("Spectrum values must be finite");
  }
  if (validation == SpectrumValueValidation::NonNegative && value < 0) {
    throw std::invalid_argument("Spectrum values must be non-negative");
  }
  if (validation == SpectrumValueValidation::Bounded01 && (value < 0 || value > 1)) {
    throw std::invalid_argument("Spectrum values must be in [0, 1]");
  }
}

class ConstantSpectrum {
public:
  explicit ConstantSpectrum(Float c = 0) : c_(c) {
    if (!std::isfinite(c_)) {
      throw std::invalid_argument("ConstantSpectrum value must be finite");
    }
  }

  Float operator()(Float) const { return c_; }

  SampledSpectrum Sample(const SampledWavelengths&) const { return SampledSpectrum(c_); }

  Float MaxValue() const { return c_; }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ ConstantSpectrum c: " << c_ << " ]";
    return out.str();
  }

private:
  Float c_ = 0;
};

class PiecewiseLinearSpectrum {
public:
  PiecewiseLinearSpectrum() = default;

  PiecewiseLinearSpectrum(
    std::vector<Float> lambdas,
    std::vector<Float> values,
    SpectrumDataPolicy policy = {},
    SpectrumMetadata metadata = {}
  )
    : policy_(policy), metadata_(std::move(metadata)) {
    Assign(std::move(lambdas), std::move(values));
  }

  static PiecewiseLinearSpectrum FromInterleaved(
    const std::vector<Float>& samples,
    bool extendVisibleInterval,
    SpectrumDataPolicy policy = {},
    SpectrumMetadata metadata = {}
  ) {
    if (samples.size() % 2 != 0) {
      throw std::invalid_argument("Interleaved spectrum samples must have wavelength/value pairs");
    }

    std::vector<Float> lambdas;
    std::vector<Float> values;
    lambdas.reserve(samples.size() / 2 + 2);
    values.reserve(samples.size() / 2 + 2);

    if (extendVisibleInterval && !samples.empty() && samples[0] > LambdaMin) {
      lambdas.push_back(LambdaMin - 1);
      values.push_back(samples[1]);
    }

    for (std::size_t i = 0; i < samples.size() / 2; ++i) {
      lambdas.push_back(samples[2 * i]);
      values.push_back(samples[2 * i + 1]);
    }

    if (extendVisibleInterval && !lambdas.empty() && lambdas.back() < LambdaMax) {
      lambdas.push_back(LambdaMax + 1);
      values.push_back(values.back());
    }

    return PiecewiseLinearSpectrum(std::move(lambdas), std::move(values), policy, metadata);
  }

  Float operator()(Float lambda) const {
    if (!std::isfinite(lambda)) {
      throw std::invalid_argument("Spectrum wavelength must be finite");
    }
    if (lambdas_.empty()) {
      return 0;
    }

    if (lambda < lambdas_.front()) {
      return Extrapolate(lambda, values_.front());
    }
    if (lambda > lambdas_.back()) {
      return Extrapolate(lambda, values_.back());
    }
    if (lambda == lambdas_.back()) {
      return values_.back();
    }

    auto upper = std::upper_bound(lambdas_.begin(), lambdas_.end(), lambda);
    std::size_t upperIndex = static_cast<std::size_t>(upper - lambdas_.begin());
    std::size_t lowerIndex = upperIndex - 1;
    Float t = (lambda - lambdas_[lowerIndex]) / (lambdas_[upperIndex] - lambdas_[lowerIndex]);
    return Lerp(t, values_[lowerIndex], values_[upperIndex]);
  }

  SampledSpectrum Sample(const SampledWavelengths& lambda) const {
    SampledSpectrum result;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      result[i] = (*this)(lambda[i]);
    }
    return result;
  }

  Float MaxValue() const {
    if (values_.empty()) {
      return 0;
    }
    return *std::max_element(values_.begin(), values_.end());
  }

  void Scale(Float scale) {
    if (!std::isfinite(scale)) {
      throw std::invalid_argument("Spectrum scale must be finite");
    }
    for (Float& value : values_) {
      value *= scale;
    }
  }

  const std::vector<Float>& Lambdas() const { return lambdas_; }
  const std::vector<Float>& Values() const { return values_; }
  const SpectrumMetadata& Metadata() const { return metadata_; }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ PiecewiseLinearSpectrum n: " << values_.size() << " ]";
    return out.str();
  }

private:
  void Assign(std::vector<Float> lambdas, std::vector<Float> values) {
    if (lambdas.size() != values.size()) {
      throw std::invalid_argument("Spectrum wavelengths and values must have equal length");
    }
    if (lambdas.size() < 2) {
      throw std::invalid_argument("PiecewiseLinearSpectrum requires at least two samples");
    }
    for (Float lambda : lambdas) {
      if (!std::isfinite(lambda)) {
        throw std::invalid_argument("Spectrum wavelengths must be finite");
      }
    }
    for (Float value : values) {
      ValidateSpectrumValue(value, policy_.validation);
    }

    if (policy_.order == SpectrumOrderPolicy::Sort) {
      std::vector<std::pair<Float, Float>> pairs;
      pairs.reserve(lambdas.size());
      for (std::size_t i = 0; i < lambdas.size(); ++i) {
        pairs.emplace_back(lambdas[i], values[i]);
      }
      std::sort(pairs.begin(), pairs.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
      });
      for (std::size_t i = 0; i < pairs.size(); ++i) {
        lambdas[i] = pairs[i].first;
        values[i] = pairs[i].second;
      }
    }

    for (std::size_t i = 1; i < lambdas.size(); ++i) {
      if (!(lambdas[i - 1] < lambdas[i])) {
        throw std::invalid_argument("Spectrum wavelengths must be strictly increasing");
      }
    }

    lambdas_ = std::move(lambdas);
    values_ = std::move(values);
  }

  Float Extrapolate(Float, Float endpointValue) const {
    switch (policy_.extrapolation) {
    case SpectrumExtrapolationPolicy::Zero:
      return 0;
    case SpectrumExtrapolationPolicy::Constant:
      return endpointValue;
    case SpectrumExtrapolationPolicy::Error:
      throw std::out_of_range("Spectrum wavelength is outside the sampled domain");
    }
    return 0;
  }

  std::vector<Float> lambdas_;
  std::vector<Float> values_;
  SpectrumDataPolicy policy_;
  SpectrumMetadata metadata_;
};

class DenselySampledSpectrum {
public:
  DenselySampledSpectrum(
    int lambdaMin = static_cast<int>(LambdaMin),
    int lambdaMax = static_cast<int>(LambdaMax),
    SpectrumMetadata metadata = {}
  )
    : lambdaMin_(lambdaMin),
      lambdaMax_(lambdaMax),
      values_(static_cast<std::size_t>(lambdaMax - lambdaMin + 1), 0),
      metadata_(std::move(metadata)) {
    if (lambdaMax_ < lambdaMin_) {
      throw std::invalid_argument("Dense spectrum lambdaMax must be >= lambdaMin");
    }
  }

  DenselySampledSpectrum(
    int lambdaMin,
    int lambdaMax,
    std::vector<Float> values,
    SpectrumMetadata metadata = {}
  )
    : lambdaMin_(lambdaMin),
      lambdaMax_(lambdaMax),
      values_(std::move(values)),
      metadata_(std::move(metadata)) {
    if (lambdaMax_ < lambdaMin_) {
      throw std::invalid_argument("Dense spectrum lambdaMax must be >= lambdaMin");
    }
    if (values_.size() != static_cast<std::size_t>(lambdaMax_ - lambdaMin_ + 1)) {
      throw std::invalid_argument("Dense spectrum value count must match wavelength interval");
    }
    for (Float value : values_) {
      ValidateSpectrumValue(value, SpectrumValueValidation::NonNegative);
    }
  }

  template <typename F>
  static DenselySampledSpectrum SampleFunction(
    F&& function,
    int lambdaMin = static_cast<int>(LambdaMin),
    int lambdaMax = static_cast<int>(LambdaMax),
    SpectrumMetadata metadata = {}
  ) {
    DenselySampledSpectrum spectrum(lambdaMin, lambdaMax, std::move(metadata));
    for (int lambda = lambdaMin; lambda <= lambdaMax; ++lambda) {
      Float value = static_cast<Float>(function(static_cast<Float>(lambda)));
      ValidateSpectrumValue(value, SpectrumValueValidation::NonNegative);
      spectrum.values_[static_cast<std::size_t>(lambda - lambdaMin)] = value;
    }
    return spectrum;
  }

  Float operator()(Float lambda) const {
    if (!std::isfinite(lambda)) {
      throw std::invalid_argument("Spectrum wavelength must be finite");
    }
    int offset = static_cast<int>(std::lround(lambda)) - lambdaMin_;
    if (offset < 0 || offset >= static_cast<int>(values_.size())) {
      return 0;
    }
    return values_[static_cast<std::size_t>(offset)];
  }

  SampledSpectrum Sample(const SampledWavelengths& lambda) const {
    SampledSpectrum result;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      result[i] = (*this)(lambda[i]);
    }
    return result;
  }

  Float MaxValue() const {
    if (values_.empty()) {
      return 0;
    }
    return *std::max_element(values_.begin(), values_.end());
  }

  void Scale(Float scale) {
    if (!std::isfinite(scale)) {
      throw std::invalid_argument("Spectrum scale must be finite");
    }
    for (Float& value : values_) {
      value *= scale;
    }
  }

  int LambdaMinValue() const { return lambdaMin_; }
  int LambdaMaxValue() const { return lambdaMax_; }
  const std::vector<Float>& Values() const { return values_; }
  const SpectrumMetadata& Metadata() const { return metadata_; }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ DenselySampledSpectrum lambda_min: " << lambdaMin_
        << " lambda_max: " << lambdaMax_ << " ]";
    return out.str();
  }

private:
  int lambdaMin_ = static_cast<int>(LambdaMin);
  int lambdaMax_ = static_cast<int>(LambdaMax);
  std::vector<Float> values_;
  SpectrumMetadata metadata_;
};

inline Float Blackbody(Float lambda, Float temperature) {
  if (temperature <= 0) {
    return 0;
  }
  constexpr Float c = static_cast<Float>(299792458);
  constexpr Float h = static_cast<Float>(6.62606957e-34);
  constexpr Float kb = static_cast<Float>(1.3806488e-23);
  Float l = lambda * static_cast<Float>(1e-9);
  Float denominator =
    std::pow(l, static_cast<Float>(5)) *
    (std::exp((h * c) / (l * kb * temperature)) - static_cast<Float>(1));
  Float le = (static_cast<Float>(2) * h * c * c) / denominator;
  return std::isfinite(le) ? le : 0;
}

class BlackbodySpectrum {
public:
  explicit BlackbodySpectrum(Float temperature) : temperature_(temperature) {
    if (!(temperature_ > 0) || !std::isfinite(temperature_)) {
      throw std::invalid_argument("BlackbodySpectrum temperature must be finite and positive");
    }
    Float lambdaMax = static_cast<Float>(2.8977721e-3) / temperature_;
    normalizationFactor_ = static_cast<Float>(1) / Blackbody(lambdaMax * static_cast<Float>(1e9), temperature_);
  }

  Float operator()(Float lambda) const { return Blackbody(lambda, temperature_) * normalizationFactor_; }

  SampledSpectrum Sample(const SampledWavelengths& lambda) const {
    SampledSpectrum result;
    for (int i = 0; i < NSpectrumSamples; ++i) {
      result[i] = (*this)(lambda[i]);
    }
    return result;
  }

  Float MaxValue() const { return 1; }

  Float Temperature() const { return temperature_; }

  std::string ToString() const {
    std::ostringstream out;
    out << "[ BlackbodySpectrum T: " << temperature_ << " ]";
    return out.str();
  }

private:
  Float temperature_ = 0;
  Float normalizationFactor_ = 1;
};

class Spectrum {
public:
  Spectrum() = default;
  Spectrum(ConstantSpectrum spectrum) : spectrum_(std::move(spectrum)) {}
  Spectrum(PiecewiseLinearSpectrum spectrum) : spectrum_(std::move(spectrum)) {}
  Spectrum(DenselySampledSpectrum spectrum) : spectrum_(std::move(spectrum)) {}
  Spectrum(BlackbodySpectrum spectrum) : spectrum_(std::move(spectrum)) {}

  bool IsValid() const { return !std::holds_alternative<std::monostate>(spectrum_); }
  explicit operator bool() const { return IsValid(); }

  Float operator()(Float lambda) const {
    if (!IsValid()) {
      return 0;
    }
    return std::visit(
      [lambda](const auto& spectrum) -> Float {
        using T = std::decay_t<decltype(spectrum)>;
        if constexpr (std::is_same<T, std::monostate>::value) {
          return 0;
        } else {
          return spectrum(lambda);
        }
      },
      spectrum_
    );
  }

  SampledSpectrum Sample(const SampledWavelengths& lambda) const {
    if (!IsValid()) {
      return SampledSpectrum(0);
    }
    return std::visit(
      [&lambda](const auto& spectrum) -> SampledSpectrum {
        using T = std::decay_t<decltype(spectrum)>;
        if constexpr (std::is_same<T, std::monostate>::value) {
          return SampledSpectrum(0);
        } else {
          return spectrum.Sample(lambda);
        }
      },
      spectrum_
    );
  }

  Float MaxValue() const {
    if (!IsValid()) {
      return 0;
    }
    return std::visit(
      [](const auto& spectrum) -> Float {
        using T = std::decay_t<decltype(spectrum)>;
        if constexpr (std::is_same<T, std::monostate>::value) {
          return 0;
        } else {
          return spectrum.MaxValue();
        }
      },
      spectrum_
    );
  }

  std::string ToString() const {
    if (!IsValid()) {
      return "(nullptr)";
    }
    return std::visit(
      [](const auto& spectrum) -> std::string {
        using T = std::decay_t<decltype(spectrum)>;
        if constexpr (std::is_same<T, std::monostate>::value) {
          return "(nullptr)";
        } else {
          return spectrum.ToString();
        }
      },
      spectrum_
    );
  }

private:
  std::variant<std::monostate, ConstantSpectrum, PiecewiseLinearSpectrum, DenselySampledSpectrum, BlackbodySpectrum>
    spectrum_;
};

inline Float InnerProduct(const Spectrum& f, const Spectrum& g) {
  Float integral = 0;
  for (int lambda = static_cast<int>(LambdaMin); lambda <= static_cast<int>(LambdaMax); ++lambda) {
    integral += f(static_cast<Float>(lambda)) * g(static_cast<Float>(lambda));
  }
  return integral;
}

inline Float IntegrateSpectrum(const Spectrum& spectrum) {
  return InnerProduct(spectrum, Spectrum(ConstantSpectrum(1)));
}

inline void NormalizeToY(PiecewiseLinearSpectrum& spectrum, const Spectrum& y) {
  Float yIntegral = InnerProduct(Spectrum(spectrum), y);
  if (!(yIntegral > 0) || !std::isfinite(yIntegral)) {
    throw std::invalid_argument("Cannot normalize a spectrum with non-positive Y integral");
  }
  spectrum.Scale(CIEYIntegral / yIntegral);
}

inline XYZ SpectrumToXYZ(const Spectrum& spectrum, const Spectrum& x, const Spectrum& y, const Spectrum& z) {
  return XYZ(
    InnerProduct(spectrum, x) / CIEYIntegral,
    InnerProduct(spectrum, y) / CIEYIntegral,
    InnerProduct(spectrum, z) / CIEYIntegral
  );
}

namespace detail {

inline std::string TrimCarriageReturn(std::string value) {
  if (!value.empty() && value.back() == '\r') {
    value.pop_back();
  }
  return value;
}

inline std::vector<std::string> SplitTabLine(const std::string& line) {
  std::vector<std::string> fields;
  std::size_t start = 0;
  while (true) {
    std::size_t tab = line.find('\t', start);
    if (tab == std::string::npos) {
      fields.push_back(TrimCarriageReturn(line.substr(start)));
      break;
    }
    fields.push_back(line.substr(start, tab - start));
    start = tab + 1;
  }
  return fields;
}

inline Float ParseFloatField(const std::string& value, const std::string& fieldName) {
  char* end = nullptr;
  double parsedDouble = std::strtod(value.c_str(), &end);
  Float parsed = static_cast<Float>(parsedDouble);
  if (end == value.c_str() || *end != '\0' || !std::isfinite(parsed)) {
    throw std::invalid_argument("Invalid floating-point field " + fieldName + ": " + value);
  }
  return parsed;
}

inline std::string JoinPath(const std::string& lhs, const std::string& rhs) {
  if (lhs.empty()) {
    return rhs;
  }
  char last = lhs[lhs.size() - 1];
  if (last == '/' || last == '\\') {
    return lhs + rhs;
  }
  return lhs + "/" + rhs;
}

inline bool FileExists(const std::string& path) {
  std::ifstream input(path.c_str(), std::ios::in);
  return input.good();
}

inline bool IsSpectralAssetDirectory(const std::string& path) {
  return FileExists(JoinPath(path, "named-spectra-v1.tsv")) &&
         FileExists(JoinPath(path, "named-spectra-v1-metadata.tsv"));
}

} // namespace detail

inline std::string FindSpectralAssetDirectory(const std::string& preferred = "") {
  std::vector<std::string> candidates;
  if (!preferred.empty()) {
    candidates.push_back(preferred);
  }
  if (const char* env = std::getenv("RAYRENDER_SPECTRAL_ASSET_DIR")) {
    candidates.emplace_back(env);
  }
  if (const char* env = std::getenv("RAYRENDER_PACKAGE_ROOT")) {
    candidates.push_back(detail::JoinPath(env, "extdata/spectral"));
  }
  candidates.emplace_back("inst/extdata/spectral");
  candidates.emplace_back("../inst/extdata/spectral");
  candidates.emplace_back("../../inst/extdata/spectral");

  for (const std::string& candidate : candidates) {
    if (detail::IsSpectralAssetDirectory(candidate)) {
      return candidate;
    }
  }
  throw std::runtime_error("Unable to find rayrender spectral asset directory");
}

class NamedSpectrumRegistry {
public:
  static NamedSpectrumRegistry LoadFromDirectory(const std::string& assetDirectory) {
    NamedSpectrumRegistry registry;
    registry.LoadMetadata(detail::JoinPath(assetDirectory, "named-spectra-v1-metadata.tsv"));
    registry.LoadData(detail::JoinPath(assetDirectory, "named-spectra-v1.tsv"));
    return registry;
  }

  static NamedSpectrumRegistry LoadDefault() {
    return LoadFromDirectory(FindSpectralAssetDirectory());
  }

  const Spectrum* Get(const std::string& name) const {
    auto iter = spectra_.find(name);
    if (iter == spectra_.end()) {
      return nullptr;
    }
    return &iter->second.spectrum;
  }

  const SpectrumMetadata* Metadata(const std::string& name) const {
    auto iter = spectra_.find(name);
    if (iter == spectra_.end()) {
      return nullptr;
    }
    return &iter->second.metadata;
  }

  const Spectrum& GetOrThrow(const std::string& name) const {
    const Spectrum* spectrum = Get(name);
    if (spectrum == nullptr) {
      throw std::out_of_range("Unknown named spectrum: " + name);
    }
    return *spectrum;
  }

  std::vector<std::string> Names() const {
    std::vector<std::string> names;
    names.reserve(spectra_.size());
    for (const auto& item : spectra_) {
      names.push_back(item.first);
    }
    return names;
  }

  std::size_t Size() const { return spectra_.size(); }

private:
  struct Entry {
    Spectrum spectrum;
    SpectrumMetadata metadata;
  };

  void LoadMetadata(const std::string& path) {
    std::ifstream input(path.c_str());
    if (!input) {
      throw std::runtime_error("Unable to open spectral metadata file: " + path);
    }

    std::string line;
    if (!std::getline(input, line)) {
      throw std::runtime_error("Spectral metadata file is empty: " + path);
    }
    std::vector<std::string> header = detail::SplitTabLine(line);
    const std::vector<std::string> expected = {
      "name",
      "semantic_type",
      "wavelength_unit",
      "value_unit",
      "normalization",
      "source_citation",
      "source_license",
      "asset_generation_script_version",
      "sha256",
      "storage"
    };
    if (header != expected) {
      throw std::runtime_error("Unexpected spectral metadata header in " + path);
    }

    while (std::getline(input, line)) {
      if (line.empty()) {
        continue;
      }
      std::vector<std::string> fields = detail::SplitTabLine(line);
      if (fields.size() != expected.size()) {
        throw std::runtime_error("Malformed spectral metadata row in " + path);
      }
      SpectrumMetadata metadata;
      metadata.name = fields[0];
      metadata.semantic = ParseSpectrumSemantic(fields[1]);
      metadata.wavelengthUnit = fields[2];
      metadata.valueUnit = fields[3];
      metadata.normalization = fields[4];
      metadata.sourceCitation = fields[5];
      metadata.sourceLicense = fields[6];
      metadata.assetGenerationScriptVersion = fields[7];
      metadata.sha256 = fields[8];
      metadata.storage = fields[9];
      if (metadata.wavelengthUnit != "nm") {
        throw std::runtime_error("Spectral asset wavelengths must be stored in nanometers: " + metadata.name);
      }
      metadataByName_[metadata.name] = metadata;
    }
  }

  void LoadData(const std::string& path) {
    std::ifstream input(path.c_str());
    if (!input) {
      throw std::runtime_error("Unable to open spectral data file: " + path);
    }

    std::string line;
    if (!std::getline(input, line)) {
      throw std::runtime_error("Spectral data file is empty: " + path);
    }
    std::vector<std::string> header = detail::SplitTabLine(line);
    const std::vector<std::string> expected = {"name", "lambda_nm", "value"};
    if (header != expected) {
      throw std::runtime_error("Unexpected spectral data header in " + path);
    }

    std::map<std::string, std::pair<std::vector<Float>, std::vector<Float>>> rows;
    while (std::getline(input, line)) {
      if (line.empty()) {
        continue;
      }
      std::vector<std::string> fields = detail::SplitTabLine(line);
      if (fields.size() != expected.size()) {
        throw std::runtime_error("Malformed spectral data row in " + path);
      }
      const std::string& name = fields[0];
      rows[name].first.push_back(detail::ParseFloatField(fields[1], "lambda_nm"));
      rows[name].second.push_back(detail::ParseFloatField(fields[2], "value"));
    }

    for (const auto& row : rows) {
      auto metadataIter = metadataByName_.find(row.first);
      if (metadataIter == metadataByName_.end()) {
        throw std::runtime_error("Spectral data row has no metadata: " + row.first);
      }

      SpectrumMetadata metadata = metadataIter->second;
      SpectrumDataPolicy policy;
      policy.order = SpectrumOrderPolicy::RequireSorted;
      policy.extrapolation = SpectrumExtrapolationPolicy::Zero;
      policy.validation = SpectrumValueValidation::NonNegative;

      Spectrum spectrum;
      if (metadata.storage == "dense_1nm") {
        int lambdaMin = static_cast<int>(std::lround(row.second.first.front()));
        int lambdaMax = static_cast<int>(std::lround(row.second.first.back()));
        std::vector<Float> values = row.second.second;
        for (std::size_t i = 0; i < row.second.first.size(); ++i) {
          Float expectedLambda = static_cast<Float>(lambdaMin + static_cast<int>(i));
          if (std::fabs(row.second.first[i] - expectedLambda) > static_cast<Float>(1e-4)) {
            throw std::runtime_error("Dense spectral asset is not sampled at 1 nm intervals: " + row.first);
          }
        }
        spectrum = Spectrum(DenselySampledSpectrum(lambdaMin, lambdaMax, std::move(values), metadata));
      } else if (metadata.storage == "piecewise_linear") {
        spectrum = Spectrum(PiecewiseLinearSpectrum(row.second.first, row.second.second, policy, metadata));
      } else {
        throw std::runtime_error("Unknown spectral asset storage class: " + metadata.storage);
      }

      spectra_[row.first] = Entry{spectrum, metadata};
    }

    if (spectra_.size() != metadataByName_.size()) {
      throw std::runtime_error("Spectral data and metadata row counts do not match");
    }
  }

  std::map<std::string, SpectrumMetadata> metadataByName_;
  std::map<std::string, Entry> spectra_;
};

inline const Spectrum* GetNamedSpectrum(const NamedSpectrumRegistry& registry, const std::string& name) {
  return registry.Get(name);
}

} // namespace base
} // namespace rayrender

#endif
