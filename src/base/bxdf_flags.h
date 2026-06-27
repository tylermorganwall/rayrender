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
  Specular = 1u << 4
};

enum class TransportMode {
  Radiance,
  Importance
};

enum class BxDFReflTransFlags : std::uint32_t {
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

} // namespace base
} // namespace rayrender

#endif
