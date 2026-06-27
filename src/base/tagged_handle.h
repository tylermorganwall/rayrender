#ifndef RAYRENDER_BASE_TAGGED_HANDLE_H
#define RAYRENDER_BASE_TAGGED_HANDLE_H

#include <cstdint>
#include <limits>

namespace rayrender {
namespace base {

enum class DispatchTag : std::uint8_t {
  Spectrum,
  Texture,
  Material,
  BxDF,
  Light,
  Medium,
  PhaseFunction
};

inline const char* DispatchTagName(DispatchTag tag) {
  switch (tag) {
  case DispatchTag::Spectrum:
    return "Spectrum";
  case DispatchTag::Texture:
    return "Texture";
  case DispatchTag::Material:
    return "Material";
  case DispatchTag::BxDF:
    return "BxDF";
  case DispatchTag::Light:
    return "Light";
  case DispatchTag::Medium:
    return "Medium";
  case DispatchTag::PhaseFunction:
    return "PhaseFunction";
  }
  return "Unknown";
}

template <DispatchTag TagValue>
class TaggedHandle {
public:
  using IndexType = std::uint32_t;
  using GenerationType = std::uint32_t;

  static constexpr DispatchTag Tag = TagValue;

  TaggedHandle() = default;

  static constexpr TaggedHandle Invalid() {
    return TaggedHandle();
  }

  static constexpr TaggedHandle FromIndex(IndexType index, GenerationType generation = 1) {
    return TaggedHandle(index, generation);
  }

  static constexpr IndexType InvalidIndex() {
    return std::numeric_limits<IndexType>::max();
  }

  constexpr bool IsValid() const {
    return index_ != InvalidIndex();
  }

  constexpr IndexType Index() const {
    return index_;
  }

  constexpr GenerationType Generation() const {
    return generation_;
  }

  constexpr bool operator==(const TaggedHandle& other) const {
    return index_ == other.index_ && generation_ == other.generation_;
  }

  constexpr bool operator!=(const TaggedHandle& other) const {
    return !(*this == other);
  }

private:
  constexpr TaggedHandle(IndexType index, GenerationType generation) : index_(index), generation_(generation) {}

  IndexType index_ = InvalidIndex();
  GenerationType generation_ = 0;
};

using SpectrumHandle = TaggedHandle<DispatchTag::Spectrum>;
using TextureHandle = TaggedHandle<DispatchTag::Texture>;
using MaterialHandle = TaggedHandle<DispatchTag::Material>;
using BxDFHandle = TaggedHandle<DispatchTag::BxDF>;
using LightHandle = TaggedHandle<DispatchTag::Light>;
using MediumHandle = TaggedHandle<DispatchTag::Medium>;
using PhaseFunctionHandle = TaggedHandle<DispatchTag::PhaseFunction>;

} // namespace base
} // namespace rayrender

#endif
