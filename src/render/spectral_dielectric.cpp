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

} // namespace

bool ResolvedDielectricInterface::IsUnity() const {
  return etaInside == etaOutside;
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
  region.medium = medium;
  region.debugName = std::move(debugName);
  return AddRegion(std::move(region));
}

RegionId DielectricRegionTable::AddRegion(DielectricRegion region) {
  if (region.id == ExteriorRegionId) {
    region.id = static_cast<RegionId>(regions_.size());
  }
  if (region.id != regions_.size()) {
    throw std::invalid_argument("dielectric regions must be added in RegionId order");
  }
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
  transition.interface.etaOutside = table_->Get(activePositive).eta;
  transition.interface.etaInside = table_->Get(activeNegative).eta;
  transition.interface.ratioIsConstant = true;
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
