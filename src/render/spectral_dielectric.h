#ifndef RAYRENDER_RENDER_SPECTRAL_DIELECTRIC_H
#define RAYRENDER_RENDER_SPECTRAL_DIELECTRIC_H

#include "../base/base.h"
#include "../base/spectrum.h"
#include "../math/vectypes.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace rayrender {
namespace render {

class Scene;

using RegionId = std::uint32_t;

constexpr RegionId ExteriorRegionId = 0;

using EtaSpectrumHandle = std::shared_ptr<const base::Spectrum>;

EtaSpectrumHandle ConstantEtaSpectrum(Float eta);
EtaSpectrumHandle SampledEtaSpectrum(
  std::vector<Float> wavelengthNm,
  std::vector<Float> eta
);
EtaSpectrumHandle CauchyEtaSpectrum(Float a, Float b = 0, Float c = 0);
EtaSpectrumHandle SellmeierEtaSpectrum(std::vector<Float> b, std::vector<Float> c);
EtaSpectrumHandle NamedEtaSpectrum(
  const base::NamedSpectrumRegistry& registry,
  const std::string& name
);
Float EvaluateEtaSpectrum(const EtaSpectrumHandle& spectrum, Float lambdaNm);
bool EtaSpectrumIsConstant(const EtaSpectrumHandle& spectrum);

enum class RegionSide {
  NegativeNormal,
  PositiveNormal,
  Inside
};

enum class DielectricBoundaryKind {
  PrioritySkipped,
  IndexMatchedNull,
  ScatteringInterface
};

struct DielectricRegion {
  RegionId id = ExteriorRegionId;
  int priority = 0;
  Float eta = 1;
  EtaSpectrumHandle etaSpectrum = ConstantEtaSpectrum(1);
  base::MediumHandle medium = base::MediumHandle::Invalid();
  std::string debugName;
};

struct RegionBoundaryAttachment {
  RegionId region = ExteriorRegionId;
  RegionSide side = RegionSide::NegativeNormal;
};

struct DielectricTransitionChange {
  RegionId region = ExteriorRegionId;
  bool presentAfter = false;
};

struct DielectricTransitionToken {
  std::vector<DielectricTransitionChange> changes;
  std::uint64_t expectedGeneration = 0;
  RegionId expectedActiveAfter = ExteriorRegionId;
};

struct ResolvedDielectricInterface {
  int outsideRegionId = static_cast<int>(ExteriorRegionId);
  int insideRegionId = static_cast<int>(ExteriorRegionId);
  Float etaOutside = 1;
  Float etaInside = 1;
  EtaSpectrumHandle etaOutsideSpectrum = ConstantEtaSpectrum(1);
  EtaSpectrumHandle etaInsideSpectrum = ConstantEtaSpectrum(1);
  bool ratioIsConstant = true;
  bool ratioIsUnity = true;

  Float Eta(Float lambdaNm) const;
  Float Eta() const {
    return etaInside / etaOutside;
  }
  bool IsUnity() const;
};

struct ResolvedDielectricTransition {
  DielectricBoundaryKind kind = DielectricBoundaryKind::PrioritySkipped;
  bool crossesToNegativeSide = false;
  RegionId activePositiveSide = ExteriorRegionId;
  RegionId activeNegativeSide = ExteriorRegionId;
  RegionId activeBefore = ExteriorRegionId;
  RegionId activeAfter = ExteriorRegionId;
  ResolvedDielectricInterface interface;
  DielectricTransitionToken token;

  bool IsNullTraversal() const;
};

class DielectricRegionTable {
public:
  DielectricRegionTable();

  RegionId AddRegion(
    Float eta,
    int priority,
    std::string debugName = std::string(),
    base::MediumHandle medium = base::MediumHandle::Invalid()
  );
  RegionId AddRegion(
    EtaSpectrumHandle eta,
    int priority,
    std::string debugName = std::string(),
    base::MediumHandle medium = base::MediumHandle::Invalid()
  );
  RegionId AddRegion(
    base::Spectrum eta,
    int priority,
    std::string debugName = std::string(),
    base::MediumHandle medium = base::MediumHandle::Invalid()
  );
  RegionId AddRegion(DielectricRegion region);

  const DielectricRegion& Get(RegionId id) const;
  const DielectricRegion& Exterior() const;
  bool IsValid(RegionId id) const;
  std::size_t Size() const;
  bool VacuumOnly() const;

private:
  std::vector<DielectricRegion> regions_;
};

class DielectricPathState {
public:
  DielectricPathState() = default;
  explicit DielectricPathState(const DielectricRegionTable* table);

  static DielectricPathState FromInitialRegions(
    const DielectricRegionTable* table,
    const std::vector<RegionId>& regions
  );
  static DielectricPathState FromPointContainment(
    const DielectricRegionTable* table,
    const Scene& scene,
    const point3f& p
  );

  bool HasTable() const;
  bool Contains(RegionId id) const;
  RegionId ActiveRegionId() const;
  const DielectricRegion& ActiveRegion() const;
  base::MediumHandle ActiveMedium() const;
  std::uint64_t Generation() const;
  const std::vector<RegionId>& Members() const;

  ResolvedDielectricTransition Analyze(
    const std::vector<RegionBoundaryAttachment>& boundaries,
    normal3f geometricNormal,
    const vec3f& rayDirection
  ) const;
  void Commit(const DielectricTransitionToken& token);

private:
  DielectricPathState WithBoundarySide(
    const std::vector<RegionBoundaryAttachment>& boundaries,
    bool negativeSide
  ) const;
  void SetMembership(RegionId id, bool present);

  const DielectricRegionTable* table_ = nullptr;
  std::vector<RegionId> members_;
  std::uint64_t generation_ = 0;
};

} // namespace render
} // namespace rayrender

#endif
