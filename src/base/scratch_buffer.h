#ifndef RAYRENDER_BASE_SCRATCH_BUFFER_H
#define RAYRENDER_BASE_SCRATCH_BUFFER_H

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <new>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace rayrender {
namespace base {

class ScratchBuffer {
public:
  static constexpr std::size_t MaxSupportedAlignment = 64;

  explicit ScratchBuffer(std::size_t capacityBytes)
      : storage_(capacityBytes + MaxSupportedAlignment - 1), capacity_(capacityBytes),
        base_(AlignPointer(storage_.data(), MaxSupportedAlignment)) {}

  ~ScratchBuffer() {
    Reset();
  }

  ScratchBuffer(const ScratchBuffer&) = delete;
  ScratchBuffer& operator=(const ScratchBuffer&) = delete;

  ScratchBuffer(ScratchBuffer&&) = delete;
  ScratchBuffer& operator=(ScratchBuffer&&) = delete;

  void* AllocateRaw(std::size_t size, std::size_t alignment) {
    if (alignment > MaxSupportedAlignment) {
      throw std::invalid_argument("ScratchBuffer requested alignment exceeds MaxSupportedAlignment");
    }
    std::size_t alignedOffset = AlignUp(offset_, alignment);
    if (alignedOffset > capacity_ || size > capacity_ - alignedOffset) {
      throw std::bad_alloc();
    }
    offset_ = alignedOffset + size;
    highWater_ = std::max(highWater_, offset_);
    return base_ + alignedOffset;
  }

  template <typename T, typename... Args>
  T* Create(Args&&... args) {
    static_assert(!std::is_array<T>::value, "ScratchBuffer cannot allocate arrays");
    static_assert(std::is_nothrow_destructible<T>::value, "ScratchBuffer objects must be nothrow destructible");

    std::size_t previousOffset = offset_;
    void* memory = AllocateRaw(sizeof(T), alignof(T));
    try {
      T* object = new (memory) T(std::forward<Args>(args)...);
      if (!std::is_trivially_destructible<T>::value) {
        try {
          destructors_.push_back({object, &Destroy<T>});
        } catch (...) {
          object->~T();
          offset_ = previousOffset;
          throw;
        }
      }
      return object;
    } catch (...) {
      offset_ = previousOffset;
      throw;
    }
  }

  void Reset() noexcept {
    for (auto it = destructors_.rbegin(); it != destructors_.rend(); ++it) {
      it->destroy(it->object);
    }
    destructors_.clear();
    offset_ = 0;
  }

  std::size_t Capacity() const {
    return capacity_;
  }

  std::size_t Used() const {
    return offset_;
  }

  std::size_t HighWater() const {
    return highWater_;
  }

  std::size_t TrackedObjectCount() const {
    return destructors_.size();
  }

private:
  struct DestructorRecord {
    void* object = nullptr;
    void (*destroy)(void*) noexcept = nullptr;
  };

  template <typename T>
  static void Destroy(void* object) noexcept {
    static_cast<T*>(object)->~T();
  }

  static std::size_t AlignUp(std::size_t offset, std::size_t alignment) {
    if (alignment == 0 || (alignment & (alignment - 1)) != 0) {
      throw std::invalid_argument("ScratchBuffer alignment must be a nonzero power of two");
    }
    std::size_t mask = alignment - 1;
    if (offset > static_cast<std::size_t>(-1) - mask) {
      throw std::bad_alloc();
    }
    return (offset + mask) & ~mask;
  }

  static std::byte* AlignPointer(std::byte* pointer, std::size_t alignment) {
    std::uintptr_t raw = reinterpret_cast<std::uintptr_t>(pointer);
    std::uintptr_t mask = alignment - 1;
    std::uintptr_t aligned = (raw + mask) & ~mask;
    return reinterpret_cast<std::byte*>(aligned);
  }

  std::vector<std::byte> storage_;
  std::vector<DestructorRecord> destructors_;
  std::size_t capacity_ = 0;
  std::byte* base_ = nullptr;
  std::size_t offset_ = 0;
  std::size_t highWater_ = 0;
};

} // namespace base
} // namespace rayrender

#endif
