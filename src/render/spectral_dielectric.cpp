#include "spectral_dielectric.h"

#include "spectral_scene.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace rayrender {
namespace render {
namespace {

bool RegionPresentOnSide(RegionSide side, bool negativeSide) {
  switch (side) {
  case RegionSide::NegativeNormal:
  case RegionSide::Inside:
    return negativeSide;
  case RegionSide::PositiveNormal:
    return !negativeSide;
  }
  return false;
}

bool AllPositiveFiniteEtaSamples(const EtaSpectrumHandle& spectrum) {
  if (!spectrum || !spectrum->IsValid()) {
    return false;
  }
  const Float samples[] = {
    base::LambdaMin,
    static_cast<Float>(550),
    base::LambdaMax
  };
  for (Float lambda : samples) {
    try {
      Float eta = (*spectrum)(lambda);
      if (!(eta > 0) || !std::isfinite(eta)) {
        return false;
      }
    } catch (const std::exception&) {
      return false;
    }
  }
  return true;
}

bool ProveSpectrumConstant(const base::Spectrum& spectrum, Float* value) {
  if (const base::ConstantSpectrum* constant = spectrum.GetIf<base::ConstantSpectrum>()) {
    *value = constant->Value();
    return true;
  }
  if (const base::CauchyIORSpectrum* cauchy = spectrum.GetIf<base::CauchyIORSpectrum>()) {
    if (cauchy->IsConstant()) {
      *value = cauchy->A();
      return true;
    }
    return false;
  }
  if (const base::SellmeierIORSpectrum* sellmeier =
        spectrum.GetIf<base::SellmeierIORSpectrum>()) {
    if (sellmeier->IsConstantOne()) {
      *value = 1;
      return true;
    }
    return false;
  }
  if (const base::PiecewiseLinearSpectrum* sampled =
        spectrum.GetIf<base::PiecewiseLinearSpectrum>()) {
    const std::vector<Float>& values = sampled->Values();
    if (values.empty()) {
      return false;
    }
    for (Float v : values) {
      if (v != values.front()) {
        return false;
      }
    }
    *value = values.front();
    return true;
  }
  if (const base::DenselySampledSpectrum* dense =
        spectrum.GetIf<base::DenselySampledSpectrum>()) {
    const std::vector<Float>& values = dense->Values();
    if (values.empty()) {
      return false;
    }
    for (Float v : values) {
      if (v != values.front()) {
        return false;
      }
    }
    *value = values.front();
    return true;
  }
  return false;
}

bool ValuesHaveConstantRatio(
  const std::vector<Float>& outside,
  const std::vector<Float>& inside,
  Float* ratio
) {
  if (outside.empty() || outside.size() != inside.size() || outside.front() == 0) {
    return false;
  }
  Float candidate = inside.front() / outside.front();
  for (std::size_t i = 0; i < outside.size(); ++i) {
    if (outside[i] == 0 || inside[i] != candidate * outside[i]) {
      return false;
    }
  }
  *ratio = candidate;
  return true;
}

bool LambdasMatch(
  const std::vector<Float>& lhs,
  const std::vector<Float>& rhs
) {
  return lhs == rhs;
}

bool ProveCauchyConstantRatio(
  const base::CauchyIORSpectrum& outside,
  const base::CauchyIORSpectrum& inside,
  Float* ratio
) {
  if (outside.A() == 0) {
    return false;
  }
  Float candidate = inside.A() / outside.A();
  if (inside.B() == candidate * outside.B() && inside.C() == candidate * outside.C()) {
    *ratio = candidate;
    return true;
  }
  return false;
}

bool VectorsEqual(const std::vector<Float>& lhs, const std::vector<Float>& rhs) {
  return lhs == rhs;
}

bool ProveConstantEtaRatio(
  const EtaSpectrumHandle& outside,
  const EtaSpectrumHandle& inside,
  Float* ratio
) {
  if (!outside || !inside) {
    return false;
  }
  if (outside.get() == inside.get()) {
    *ratio = 1;
    return true;
  }

  Float outsideConstant = 0;
  Float insideConstant = 0;
  if (ProveSpectrumConstant(*outside, &outsideConstant) &&
      ProveSpectrumConstant(*inside, &insideConstant) &&
      outsideConstant != 0) {
    *ratio = insideConstant / outsideConstant;
    return true;
  }

  const base::CauchyIORSpectrum* outsideCauchy =
    outside->GetIf<base::CauchyIORSpectrum>();
  const base::CauchyIORSpectrum* insideCauchy =
    inside->GetIf<base::CauchyIORSpectrum>();
  if (outsideCauchy != nullptr && insideCauchy != nullptr) {
    return ProveCauchyConstantRatio(*outsideCauchy, *insideCauchy, ratio);
  }

  const base::SellmeierIORSpectrum* outsideSellmeier =
    outside->GetIf<base::SellmeierIORSpectrum>();
  const base::SellmeierIORSpectrum* insideSellmeier =
    inside->GetIf<base::SellmeierIORSpectrum>();
  if (outsideSellmeier != nullptr && insideSellmeier != nullptr &&
      VectorsEqual(outsideSellmeier->B(), insideSellmeier->B()) &&
      VectorsEqual(outsideSellmeier->C(), insideSellmeier->C())) {
    *ratio = 1;
    return true;
  }

  const base::PiecewiseLinearSpectrum* outsideSampled =
    outside->GetIf<base::PiecewiseLinearSpectrum>();
  const base::PiecewiseLinearSpectrum* insideSampled =
    inside->GetIf<base::PiecewiseLinearSpectrum>();
  if (outsideSampled != nullptr && insideSampled != nullptr &&
      LambdasMatch(outsideSampled->Lambdas(), insideSampled->Lambdas())) {
    return ValuesHaveConstantRatio(outsideSampled->Values(), insideSampled->Values(), ratio);
  }

  const base::DenselySampledSpectrum* outsideDense =
    outside->GetIf<base::DenselySampledSpectrum>();
  const base::DenselySampledSpectrum* insideDense =
    inside->GetIf<base::DenselySampledSpectrum>();
  if (outsideDense != nullptr && insideDense != nullptr &&
      outsideDense->LambdaMinValue() == insideDense->LambdaMinValue() &&
      outsideDense->LambdaMaxValue() == insideDense->LambdaMaxValue()) {
    return ValuesHaveConstantRatio(outsideDense->Values(), insideDense->Values(), ratio);
  }

  return false;
}

void AssignInterfaceEta(
  ResolvedDielectricInterface& interface,
  const DielectricRegion& outside,
  const DielectricRegion& inside
) {
  interface.etaOutsideSpectrum = outside.etaSpectrum;
  interface.etaInsideSpectrum = inside.etaSpectrum;
  interface.etaOutside = EvaluateEtaSpectrum(interface.etaOutsideSpectrum, static_cast<Float>(550));
  interface.etaInside = EvaluateEtaSpectrum(interface.etaInsideSpectrum, static_cast<Float>(550));
  Float ratio = 0;
  interface.ratioIsConstant =
    ProveConstantEtaRatio(interface.etaOutsideSpectrum, interface.etaInsideSpectrum, &ratio);
  interface.ratioIsUnity = interface.ratioIsConstant && ratio == 1;
}

} // namespace

EtaSpectrumHandle ConstantEtaSpectrum(Float eta) {
  return std::make_shared<const base::Spectrum>(base::ConstantSpectrum(eta));
}

EtaSpectrumHandle SampledEtaSpectrum(
  std::vector<Float> wavelengthNm,
  std::vector<Float> eta
) {
  for (Float value : eta) {
    if (!(value > 0) || !std::isfinite(value)) {
      throw std::invalid_argument("sampled eta values must be positive finite");
    }
  }
  base::SpectrumDataPolicy policy;
  policy.order = base::SpectrumOrderPolicy::RequireSorted;
  policy.extrapolation = base::SpectrumExtrapolationPolicy::Constant;
  policy.validation = base::SpectrumValueValidation::Finite;
  base::SpectrumMetadata metadata;
  metadata.semantic = base::SpectrumSemantic::GlassEta;
  metadata.valueUnit = "dimensionless";
  return std::make_shared<const base::Spectrum>(
    base::PiecewiseLinearSpectrum(std::move(wavelengthNm), std::move(eta), policy, metadata)
  );
}

EtaSpectrumHandle CauchyEtaSpectrum(Float a, Float b, Float c) {
  return std::make_shared<const base::Spectrum>(base::CauchyIORSpectrum(a, b, c));
}

EtaSpectrumHandle SellmeierEtaSpectrum(std::vector<Float> b, std::vector<Float> c) {
  return std::make_shared<const base::Spectrum>(
    base::SellmeierIORSpectrum(std::move(b), std::move(c))
  );
}

EtaSpectrumHandle NamedEtaSpectrum(
  const base::NamedSpectrumRegistry& registry,
  const std::string& name
) {
  const base::SpectrumMetadata* metadata = registry.Metadata(name);
  if (metadata == nullptr) {
    throw std::out_of_range("Unknown named eta spectrum: " + name);
  }
  if (metadata->semantic != base::SpectrumSemantic::GlassEta &&
      metadata->semantic != base::SpectrumSemantic::OpticalEta) {
    throw std::invalid_argument("Named eta spectrum must have glass_eta or optical_eta semantic");
  }
  return std::make_shared<const base::Spectrum>(registry.GetOrThrow(name));
}

Float EvaluateEtaSpectrum(const EtaSpectrumHandle& spectrum, Float lambdaNm) {
  if (!spectrum || !spectrum->IsValid()) {
    throw std::runtime_error("eta spectrum is not valid");
  }
  Float eta = (*spectrum)(lambdaNm);
  if (!(eta > 0) || !std::isfinite(eta)) {
    throw std::runtime_error("eta spectrum evaluated to a non-positive or non-finite value");
  }
  return eta;
}

bool EtaSpectrumIsConstant(const EtaSpectrumHandle& spectrum) {
  if (!spectrum) {
    return false;
  }
  Float value = 0;
  return ProveSpectrumConstant(*spectrum, &value);
}

Float ResolvedDielectricInterface::Eta(Float lambdaNm) const {
  try {
    Float cachedOutside =
      EvaluateEtaSpectrum(etaOutsideSpectrum, static_cast<Float>(550));
    Float cachedInside =
      EvaluateEtaSpectrum(etaInsideSpectrum, static_cast<Float>(550));
    if (cachedOutside != etaOutside || cachedInside != etaInside) {
      return etaInside / etaOutside;
    }
  } catch (const std::exception&) {
    return etaInside / etaOutside;
  }
  return EvaluateEtaSpectrum(etaInsideSpectrum, lambdaNm) /
         EvaluateEtaSpectrum(etaOutsideSpectrum, lambdaNm);
}

bool ResolvedDielectricInterface::IsUnity() const {
  try {
    Float cachedOutside =
      EvaluateEtaSpectrum(etaOutsideSpectrum, static_cast<Float>(550));
    Float cachedInside =
      EvaluateEtaSpectrum(etaInsideSpectrum, static_cast<Float>(550));
    if (cachedOutside != etaOutside || cachedInside != etaInside) {
      return etaInside == etaOutside;
    }
  } catch (const std::exception&) {
    return etaInside == etaOutside;
  }
  return ratioIsUnity;
}

bool ResolvedDielectricTransition::IsNullTraversal() const {
  return kind == DielectricBoundaryKind::PrioritySkipped ||
         kind == DielectricBoundaryKind::IndexMatchedNull;
}

DielectricRegionTable::DielectricRegionTable() {
  DielectricRegion exterior;
  exterior.id = ExteriorRegionId;
  exterior.priority = std::numeric_limits<int>::max();
  exterior.eta = 1;
  exterior.etaSpectrum = ConstantEtaSpectrum(1);
  exterior.debugName = "exterior";
  regions_.push_back(std::move(exterior));
}

RegionId DielectricRegionTable::AddRegion(
  Float eta,
  int priority,
  std::string debugName,
  base::MediumHandle medium
) {
  DielectricRegion region;
  region.id = static_cast<RegionId>(regions_.size());
  region.priority = priority;
  region.eta = eta;
  region.etaSpectrum = ConstantEtaSpectrum(eta);
  region.medium = medium;
  region.debugName = std::move(debugName);
  return AddRegion(std::move(region));
}

RegionId DielectricRegionTable::AddRegion(
  EtaSpectrumHandle eta,
  int priority,
  std::string debugName,
  base::MediumHandle medium
) {
  DielectricRegion region;
  region.id = static_cast<RegionId>(regions_.size());
  region.priority = priority;
  region.etaSpectrum = std::move(eta);
  region.eta = EvaluateEtaSpectrum(region.etaSpectrum, static_cast<Float>(550));
  region.medium = medium;
  region.debugName = std::move(debugName);
  return AddRegion(std::move(region));
}

RegionId DielectricRegionTable::AddRegion(
  base::Spectrum eta,
  int priority,
  std::string debugName,
  base::MediumHandle medium
) {
  return AddRegion(
    std::make_shared<const base::Spectrum>(std::move(eta)),
    priority,
    std::move(debugName),
    medium
  );
}

RegionId DielectricRegionTable::AddRegion(DielectricRegion region) {
  if (region.id == ExteriorRegionId) {
    region.id = static_cast<RegionId>(regions_.size());
  }
  if (region.id != regions_.size()) {
    throw std::invalid_argument("dielectric regions must be added in RegionId order");
  }
  if (!region.etaSpectrum) {
    region.etaSpectrum = ConstantEtaSpectrum(region.eta);
  }
  if (!AllPositiveFiniteEtaSamples(region.etaSpectrum)) {
    throw std::invalid_argument("dielectric region eta spectrum must be positive finite");
  }
  region.eta = EvaluateEtaSpectrum(region.etaSpectrum, static_cast<Float>(550));
  if (!(region.eta > 0) || !std::isfinite(region.eta)) {
    throw std::invalid_argument("dielectric region eta must be positive finite");
  }
  regions_.push_back(std::move(region));
  return regions_.back().id;
}

const DielectricRegion& DielectricRegionTable::Get(RegionId id) const {
  if (!IsValid(id)) {
    throw std::out_of_range("invalid dielectric region id");
  }
  return regions_[id];
}

const DielectricRegion& DielectricRegionTable::Exterior() const {
  return regions_[ExteriorRegionId];
}

bool DielectricRegionTable::IsValid(RegionId id) const {
  return id < regions_.size();
}

std::size_t DielectricRegionTable::Size() const {
  return regions_.size();
}

bool DielectricRegionTable::VacuumOnly() const {
  for (const DielectricRegion& region : regions_) {
    if (region.medium.IsValid()) {
      return false;
    }
  }
  return true;
}

DielectricPathState::DielectricPathState(const DielectricRegionTable* table)
  : table_(table) {}

DielectricPathState DielectricPathState::FromInitialRegions(
  const DielectricRegionTable* table,
  const std::vector<RegionId>& regions
) {
  DielectricPathState state(table);
  for (RegionId region : regions) {
    state.SetMembership(region, true);
  }
  return state;
}

DielectricPathState DielectricPathState::FromPointContainment(
  const DielectricRegionTable* table,
  const Scene& scene,
  const point3f& p
) {
  DielectricPathState state(table);
  if (table == nullptr) {
    return state;
  }
  for (std::size_t i = 0; i < scene.PrimitiveCount(); ++i) {
    PrimitiveHandle handle = PrimitiveHandle::FromIndex(
      static_cast<PrimitiveHandle::IndexType>(i),
      1
    );
    const PrimitiveBinding& binding = scene.GetPrimitiveBinding(handle);
    std::vector<RegionBoundaryAttachment> boundaries = binding.DielectricBoundaries();
    if (boundaries.empty()) {
      continue;
    }
    const Shape& shape = scene.GetShape(scene.GetPrimitiveShape(handle));
    if (!shape.Capabilities().supportsContainment) {
      throw std::runtime_error("dielectric region boundary does not support point containment");
    }
    bool contains = shape.Contains(p);
    for (const RegionBoundaryAttachment& boundary : boundaries) {
      if (!table->IsValid(boundary.region)) {
        throw std::runtime_error("dielectric region boundary references an invalid region");
      }
      bool present = boundary.side == RegionSide::PositiveNormal ? !contains : contains;
      state.SetMembership(boundary.region, present);
    }
  }
  return state;
}

bool DielectricPathState::HasTable() const {
  return table_ != nullptr;
}

bool DielectricPathState::Contains(RegionId id) const {
  return std::find(members_.begin(), members_.end(), id) != members_.end();
}

RegionId DielectricPathState::ActiveRegionId() const {
  if (table_ == nullptr || members_.empty()) {
    return ExteriorRegionId;
  }
  RegionId active = ExteriorRegionId;
  int activePriority = table_->Exterior().priority;
  bool haveActive = false;
  for (RegionId id : members_) {
    const DielectricRegion& region = table_->Get(id);
    if (!haveActive || region.priority < activePriority) {
      active = id;
      activePriority = region.priority;
      haveActive = true;
    } else if (region.priority == activePriority) {
      throw std::runtime_error("equal-priority active dielectric regions are invalid");
    }
  }
  return active;
}

const DielectricRegion& DielectricPathState::ActiveRegion() const {
  if (table_ == nullptr) {
    throw std::runtime_error("dielectric path state has no region table");
  }
  return table_->Get(ActiveRegionId());
}

base::MediumHandle DielectricPathState::ActiveMedium() const {
  if (table_ == nullptr) {
    return base::MediumHandle::Invalid();
  }
  return ActiveRegion().medium;
}

std::uint64_t DielectricPathState::Generation() const {
  return generation_;
}

const std::vector<RegionId>& DielectricPathState::Members() const {
  return members_;
}

ResolvedDielectricTransition DielectricPathState::Analyze(
  const std::vector<RegionBoundaryAttachment>& boundaries,
  normal3f geometricNormal,
  const vec3f& rayDirection
) const {
  if (table_ == nullptr) {
    throw std::runtime_error("cannot analyze dielectric transition without a region table");
  }
  if (boundaries.empty()) {
    throw std::invalid_argument("dielectric transition requires at least one boundary");
  }

  DielectricPathState positive = WithBoundarySide(boundaries, false);
  DielectricPathState negative = WithBoundarySide(boundaries, true);
  RegionId activePositive = positive.ActiveRegionId();
  RegionId activeNegative = negative.ActiveRegionId();
  bool crossesToNegative = dot(rayDirection, geometricNormal) < 0;
  RegionId activeBefore = crossesToNegative ? activePositive : activeNegative;
  RegionId activeAfter = crossesToNegative ? activeNegative : activePositive;
  RegionId currentActive = ActiveRegionId();
  if (currentActive != activeBefore) {
    throw std::runtime_error("dielectric active-before state mismatch");
  }

  ResolvedDielectricTransition transition;
  transition.crossesToNegativeSide = crossesToNegative;
  transition.activePositiveSide = activePositive;
  transition.activeNegativeSide = activeNegative;
  transition.activeBefore = activeBefore;
  transition.activeAfter = activeAfter;
  transition.interface.outsideRegionId = static_cast<int>(activePositive);
  transition.interface.insideRegionId = static_cast<int>(activeNegative);
  AssignInterfaceEta(
    transition.interface,
    table_->Get(activePositive),
    table_->Get(activeNegative)
  );
  if (activePositive == activeNegative) {
    transition.kind = DielectricBoundaryKind::PrioritySkipped;
  } else if (transition.interface.IsUnity()) {
    transition.kind = DielectricBoundaryKind::IndexMatchedNull;
  } else {
    transition.kind = DielectricBoundaryKind::ScatteringInterface;
  }
  transition.token.expectedGeneration = generation_;
  transition.token.expectedActiveAfter = activeAfter;

  const DielectricPathState& afterState = crossesToNegative ? negative : positive;
  for (const RegionBoundaryAttachment& boundary : boundaries) {
    if (!table_->IsValid(boundary.region) || boundary.region == ExteriorRegionId) {
      throw std::runtime_error("dielectric boundary references an invalid region");
    }
    bool presentAfter = afterState.Contains(boundary.region);
    if (Contains(boundary.region) != presentAfter) {
      transition.token.changes.push_back({boundary.region, presentAfter});
    }
  }
  return transition;
}

void DielectricPathState::Commit(const DielectricTransitionToken& token) {
  if (token.expectedGeneration != generation_) {
    throw std::runtime_error("dielectric transition generation mismatch");
  }
  for (const DielectricTransitionChange& change : token.changes) {
    if (change.region == ExteriorRegionId || table_ == nullptr || !table_->IsValid(change.region)) {
      throw std::runtime_error("invalid dielectric transition token");
    }
  }
  for (const DielectricTransitionChange& change : token.changes) {
    SetMembership(change.region, change.presentAfter);
  }
  ++generation_;
  if (ActiveRegionId() != token.expectedActiveAfter) {
    throw std::runtime_error("dielectric transition committed to unexpected active region");
  }
}

DielectricPathState DielectricPathState::WithBoundarySide(
  const std::vector<RegionBoundaryAttachment>& boundaries,
  bool negativeSide
) const {
  DielectricPathState result(*this);
  for (const RegionBoundaryAttachment& boundary : boundaries) {
    if (table_ == nullptr || !table_->IsValid(boundary.region) ||
        boundary.region == ExteriorRegionId) {
      throw std::runtime_error("dielectric boundary references an invalid region");
    }
    result.SetMembership(boundary.region, RegionPresentOnSide(boundary.side, negativeSide));
  }
  return result;
}

void DielectricPathState::SetMembership(RegionId id, bool present) {
  if (id == ExteriorRegionId || table_ == nullptr || !table_->IsValid(id)) {
    throw std::runtime_error("invalid dielectric region membership");
  }
  auto found = std::find(members_.begin(), members_.end(), id);
  if (present) {
    if (found == members_.end()) {
      members_.push_back(id);
      std::sort(members_.begin(), members_.end());
    }
  } else if (found != members_.end()) {
    members_.erase(found);
  }
}

} // namespace render
} // namespace rayrender
