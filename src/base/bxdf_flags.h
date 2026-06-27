#ifndef RAYRENDER_BASE_BXDF_FLAGS_H
#define RAYRENDER_BASE_BXDF_FLAGS_H

#include <cstdint>
#include <type_traits>

namespace rayrender {
namespace base {

enum class BxDFFlags : std::uint32_t {
  Unset = 0,
  Reflection = 1u << 0,
  Transmission = 1u << 1,
  Diffuse = 1u << 2,
  Glossy = 1u << 3,
  Specular = 1u << 4,
  DiffuseReflection = (1u << 0) | (1u << 2),
  DiffuseTransmission = (1u << 1) | (1u << 2),
  GlossyReflection = (1u << 0) | (1u << 3),
  GlossyTransmission = (1u << 1) | (1u << 3),
  SpecularReflection = (1u << 0) | (1u << 4),
  SpecularTransmission = (1u << 1) | (1u << 4),
  All = (1u << 0) | (1u << 1) | (1u << 2) | (1u << 3) | (1u << 4)
};

enum class TransportMode {
  Radiance,
  Importance
};

enum class BxDFReflTransFlags : std::uint32_t {
  Unset = 0,
  Reflection = 1u << 0,
  Transmission = 1u << 1,
  All = (1u << 0) | (1u << 1)
};

template <typename Enum>
constexpr typename std::underlying_type<Enum>::type FlagBits(Enum value) {
  return static_cast<typename std::underlying_type<Enum>::type>(value);
}

constexpr BxDFFlags operator|(BxDFFlags lhs, BxDFFlags rhs) {
  return static_cast<BxDFFlags>(FlagBits(lhs) | FlagBits(rhs));
}

constexpr BxDFFlags operator&(BxDFFlags lhs, BxDFFlags rhs) {
  return static_cast<BxDFFlags>(FlagBits(lhs) & FlagBits(rhs));
}

inline BxDFFlags& operator|=(BxDFFlags& lhs, BxDFFlags rhs) {
  lhs = lhs | rhs;
  return lhs;
}

constexpr bool HasFlag(BxDFFlags flags, BxDFFlags flag) {
  return (FlagBits(flags) & FlagBits(flag)) == FlagBits(flag);
}

constexpr bool HasAny(BxDFFlags flags, BxDFFlags mask) {
  return (FlagBits(flags) & FlagBits(mask)) != 0;
}

constexpr bool HasAny(BxDFFlags flags, BxDFReflTransFlags mask) {
  return (FlagBits(flags) & FlagBits(mask)) != 0;
}

constexpr bool IsReflective(BxDFFlags flags) {
  return HasFlag(flags, BxDFFlags::Reflection);
}

constexpr bool IsTransmissive(BxDFFlags flags) {
  return HasFlag(flags, BxDFFlags::Transmission);
}

constexpr bool IsDiffuse(BxDFFlags flags) {
  return HasFlag(flags, BxDFFlags::Diffuse);
}

constexpr bool IsGlossy(BxDFFlags flags) {
  return HasFlag(flags, BxDFFlags::Glossy);
}

constexpr bool IsSpecular(BxDFFlags flags) {
  return HasFlag(flags, BxDFFlags::Specular);
}

constexpr bool IsNonSpecular(BxDFFlags flags) {
  return HasAny(flags, BxDFFlags::Diffuse | BxDFFlags::Glossy);
}

constexpr BxDFReflTransFlags operator|(BxDFReflTransFlags lhs, BxDFReflTransFlags rhs) {
  return static_cast<BxDFReflTransFlags>(FlagBits(lhs) | FlagBits(rhs));
}

constexpr BxDFReflTransFlags operator&(BxDFReflTransFlags lhs, BxDFReflTransFlags rhs) {
  return static_cast<BxDFReflTransFlags>(FlagBits(lhs) & FlagBits(rhs));
}

inline BxDFReflTransFlags& operator|=(BxDFReflTransFlags& lhs, BxDFReflTransFlags rhs) {
  lhs = lhs | rhs;
  return lhs;
}

constexpr bool HasFlag(BxDFReflTransFlags flags, BxDFReflTransFlags flag) {
  return (FlagBits(flags) & FlagBits(flag)) == FlagBits(flag);
}

constexpr bool HasAny(BxDFReflTransFlags flags, BxDFReflTransFlags mask) {
  return (FlagBits(flags) & FlagBits(mask)) != 0;
}

constexpr TransportMode operator!(TransportMode mode) {
  return mode == TransportMode::Radiance ? TransportMode::Importance : TransportMode::Radiance;
}

} // namespace base
} // namespace rayrender

#endif
